/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "streamlinebundlecreator.h"

#include "../../datastructures/streamlinelist.h"

namespace voreen {

StreamlineBundleCreator::StreamlineBundleCreator()
    : AsyncComputeProcessor()
    , streamlineInport_(Port::INPORT, "streamlineInport", "Streamlines Input")
    , streamlineBundleOutport_(Port::OUTPORT, "streamlineBundleOutport", "Streamline Bundle Output")
    , streamlineNoiseOutport_(Port::OUTPORT, "streamlineNoiseOutport", "Streamline Noise Output")
    , enabled_("enabled", "Enabled", true)
    , maxAverageDistanceThresholdProp_("maxAverageDistanceThreshold", "Max. Average Distance Threshold (mm)", 1.0f, 0.0f, 100.0f)
    , minNumStreamlinesPerBundleProp_("minNumStreamlinesPerBundle", "Minimal number of Streamlines per Bundle (%)", 1.0f, 0.0f, 100.0f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_ADVANCED)
    , resampleSizeProp_("resampleSize", "Streamline Resample Size", 20, 2, 100, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_ADVANCED)
{
    addPort(streamlineInport_);
    addPort(streamlineBundleOutport_);
    addPort(streamlineNoiseOutport_);

    addProperty(enabled_);
        enabled_.setGroupID("output");
    setPropertyGroupGuiName("output", "Output");

        //streamline bundles
    addProperty(maxAverageDistanceThresholdProp_);
        maxAverageDistanceThresholdProp_.setTracking(false);
        maxAverageDistanceThresholdProp_.setGroupID("streamlinebundles");
    addProperty(minNumStreamlinesPerBundleProp_);
        minNumStreamlinesPerBundleProp_.setTracking(false);
        minNumStreamlinesPerBundleProp_.setGroupID("streamlinebundles");
    addProperty(resampleSizeProp_);
        resampleSizeProp_.setTracking(false);
        resampleSizeProp_.setGroupID("streamlinebundles");
    setPropertyGroupGuiName("streamlinebundles", "Streamline Bundle Settings");
}

bool StreamlineBundleCreator::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized");
        return false;
    }

    if(!streamlineInport_.isReady()) {
        setNotReadyErrorMessage("No streamlines");
        return false;
    }

    if(!streamlineBundleOutport_.isReady() && !streamlineNoiseOutport_.isReady()) {
        setNotReadyErrorMessage("No output connected");
        return false;
    }

    return true;
}

void StreamlineBundleCreator::adjustPropertiesToInput() {
    const StreamlineListBase* streamlines = streamlineInport_.getData();
    if(!streamlines) {
        return;
    }

    tgt::vec3 length = streamlines->getOriginalSpacing() * tgt::vec3(streamlines->getOriginalDimensions());
    maxAverageDistanceThresholdProp_.setMaxValue(tgt::length(length));
    maxAverageDistanceThresholdProp_.adaptDecimalsToRange(3);
}

StreamlineBundleCreatorInput StreamlineBundleCreator::prepareComputeInput() {

    if(!enabled_.get()) {
        streamlineBundleOutport_.setData(streamlineInport_.getData(), false);
        streamlineNoiseOutport_.setData(nullptr);
        throw InvalidInputException("", InvalidInputException::S_IGNORE);
    }

    auto streamlines = streamlineInport_.getThreadSafeData();
    if(!streamlines) {
        throw InvalidInputException("No streamlines", InvalidInputException::S_ERROR);
    }

    return StreamlineBundleCreatorInput {
            std::move(streamlines),
            resampleSizeProp_.get(),
            maxAverageDistanceThresholdProp_.get(),
            minNumStreamlinesPerBundleProp_.get(),
    };
}

StreamlineBundleCreatorOutput StreamlineBundleCreator::compute(StreamlineBundleCreatorInput input, ProgressReporter& progressReporter) const {

    PortDataPointer<StreamlineListBase> streamlines = std::move(input.streamlines);

    // Create output. This is done by cloning the input list to copy it's properties like spacing
    // and removing the contained streamlines.
    StreamlineListBase* emptyList = streamlines->clone();
    emptyList->clearStreamlines();

    std::unique_ptr<StreamlineListBase> streamlineBundleOutput(emptyList);
    std::unique_ptr<StreamlineListBase> streamlineNoiseOutput(emptyList->clone());

    // MDF metric used by QuickBundles.
    // Returns MDS distance between to streamlines in physical space.
    auto calculateMDF = [&] (const Streamline& s, const Streamline& t) -> float {
        // We actually calculate the average distance here, since our streamline calculation stores
        // the elements in sorted order. Thus, we can save some time and ignore the flipped distance calculation.

        float sum = 0.0f;

        for(int i = 0; i < input.resampleSize; i++) {
            sum += tgt::distance(s.getElementAt(i).position_, t.getElementAt(i).position_);
        }

        return sum / input.resampleSize;
    };

    // Step 1: Execute quickbundle algorithm.
    std::vector<StreamlineBundle> bundles;

    size_t numStreamlines = streamlines->getStreamlines().size();
    for(size_t i = 0; i < numStreamlines; i++) {

        // Resample the streamline.
        Streamline t = streamlines->getStreamlines()[i].resample(input.resampleSize);

        // Look for the bundle the streamlines fits best in.
        float m = std::numeric_limits<float>::max();
        size_t l = 0;
        for(size_t k = 0; k < bundles.size(); k++) {

            Streamline centroid = bundles[k].getCentroid();
            float distance = calculateMDF(centroid, t);

            if(distance < m) {
                m = distance;
                l = k;
            }
        }

        if(m < input.distanceThreshold)
            bundles[l].addStreamline(std::move(t));
        else
            bundles.emplace_back(StreamlineBundle(std::move(t)));

        progressReporter.setProgress(1.0f * i / numStreamlines);
    }

    // Step 2: Mark all Streamlines not being absorbed as noise.
    for(const StreamlineBundle& bundle : bundles) {
        if(bundle.getStreamlines().size() < input.minNumStreamlinesThreshold) {
            for(const Streamline& streamline : bundle.getStreamlines()) {
                streamlineNoiseOutput->addStreamline(streamline);
            }
        }
        else {
            streamlineBundleOutput->addStreamline(bundle.getCentroid());
        }
    }

    return StreamlineBundleCreatorOutput {
            std::move(streamlineBundleOutput),
            std::move(streamlineNoiseOutput)
    };
}

void StreamlineBundleCreator::processComputeOutput(StreamlineBundleCreatorOutput output) {
    streamlineBundleOutport_.setData(output.streamlineBundles.release());
    streamlineNoiseOutport_.setData(output.streamlineNoise.release());
}

/*
This Function could be used for merging bundles later on. In case it's not being considered, feel free to remove the code.

float StreamlineBundleCreator::calculateMAM(const Streamline& s, const Streamline& t) const {

    float sumS = 0.0f;
    float sumT = 0.0f;

    float minDistanceS, minDistanceT;

    for(size_t i = 0; i < resampleSize_; i++) {
        minDistanceS = minDistanceT = std::numeric_limits<float>::max();
        for(size_t j = 0; j < resampleSize_; j++) {
            minDistanceS = std::min(minDistanceS, tgt::length((s.getElementAt(i).position_ - t.getElementAt(j).position_) * output_->getOriginalSpacing()));
            minDistanceT = std::min(minDistanceT, tgt::length((s.getElementAt(i).position_ - t.getElementAt(i).position_) * output_->getOriginalSpacing()));
        }
        sumS += minDistanceS;
        sumT += minDistanceT;
    }

    switch(distanceFunction_) {
    case MAM_MEAN:
        return (sumS + sumT) * 0.5f / resampleSize_;
    case MAM_MIN:
        return std::min(sumS, sumT) / resampleSize_;
    case MAM_MAX:
        return std::max(sumS, sumT) / resampleSize_;
    default:
        tgtAssert(false, "unsupported distance function");
        return std::numeric_limits<float>::max();
    }
}
*/

}   // namespace
