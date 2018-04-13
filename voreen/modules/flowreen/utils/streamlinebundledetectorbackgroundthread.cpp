/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "streamlinebundledetectorbackgroundthread.h"

#include "modules/flowreen/datastructures/streamlinebundle.h"

namespace voreen {

StreamlineBundleDetectorBackgroundThread::StreamlineBundleDetectorBackgroundThread(
    StreamlineCreator* processor,
    StreamlineList* output,
    float distanceThreshold,
    size_t noiseThreshold,
    size_t stepSize)
    : ProcessorBackgroundThread<StreamlineCreator>(processor)
    , output_(output)
    , distanceThreshold_(distanceThreshold)
    , noiseThreshold_(noiseThreshold)
    , resampleSize_(stepSize)
{
}

StreamlineBundleDetectorBackgroundThread::~StreamlineBundleDetectorBackgroundThread() {
}

void StreamlineBundleDetectorBackgroundThread::threadMain() {

    size_t updateProcess = std::min(output_->getStreamlines().size() / 7 + 1, static_cast<size_t>(100));

    // Step 1: Execute quickbundle algorithm.
    std::vector<StreamlineBundle> bundles;

    for(size_t i = 0; i < output_->getStreamlines().size(); i++) {

        // Resample the streamline.
        const Streamline& t = output_->getStreamlines()[i].resample(resampleSize_);

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

        if(m < distanceThreshold_)
            bundles[l].addStreamline(i, t);
        else
            bundles.push_back(StreamlineBundle(i, t));

        interruptionPoint();
        if (i % updateProcess == 0) {
            processor_->bundleDetectionProgressProp_.setProgress(static_cast<float>(i) / output_->getStreamlines().size());
        }
    }

    // Step 2: Mark all Streamlines not being absorbed as noise.
    for(StreamlineBundle& bundle : bundles) {
        if(bundle.getStreamlines().size() < noiseThreshold_) {
            for(size_t i : bundle.getStreamlines()) {
                output_->setStreamlineNoiseFlag(i);
            }
        }
        else
            // Add the bundle to the output.
            output_->addStreamlineBundle(bundle);
    }

    // set progress to 100
    processor_->bundleDetectionProgressProp_.setProgress(1.0f);
}

//---------------------------------------------------------------------------
//          Helpers
//---------------------------------------------------------------------------
float StreamlineBundleDetectorBackgroundThread::calculateMDF(const Streamline& s, const Streamline& t) const {

    // We actually calculate the average distance here, since our streamline calculation stores
    // the elements in sorted order. Thus, we can save some time and ignore the flipped distance calculation.

    float sum = 0.0f;

    for(size_t i = 0; i < resampleSize_; i++) {
        sum += tgt::length((s.getElementAt(i).position_ - t.getElementAt(i).position_) * output_->getOriginalSpacing());
    }

    return sum / resampleSize_;
}
/*
This Function could be used for merging bundles later on. In case it's not being considered, feel free to remove the code.

float StreamlineBundleDetectorBackgroundThread::calculateMAM(const Streamline& s, const Streamline& t) const {

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
