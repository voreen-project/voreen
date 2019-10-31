/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "streamlinecreator.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "../../datastructures/streamlinelist.h"

#include <random>

namespace voreen {

StreamlineCreator::StreamlineCreator()
    : AsyncComputeProcessor()
    , volInport_(Port::INPORT, "volInport", "Flow Volume Input (vec3)")
    , seedMask_(Port::INPORT, "seedMaskPort", "Seed Mask (optional)")
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamlines Output")
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 5000, 1, 100000)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , streamlineLengthThresholdProp_("streamlineLengthThresholdProp", "Threshold of Streamline Length: ", tgt::ivec2(10, 100), 2, 1000)
    , absoluteMagnitudeThresholdProp_("absoluteMagnitudeThreshold", "Threshold of Magnitude (absolute)", tgt::vec2(0.0f, 1000.0f), 0.0f, 9999.99f)
    , fitAbsoluteMagnitudeProp_("fitAbsoluteMagnitude", "Fit absolute Threshold to Input", false)
    , relativeMagnitudeThresholdProp_("relativeMagnitudeThreshold", "Threshold of Magnitude (relative)", tgt::vec2(0.0f, 100.0f), 0.0f, 100.0f, Processor::VALID)
    , stopIntegrationAngleThresholdProp_("stopIntegrationAngleThreshold", "Stop Integration on Angle", 180, 0, 180, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_ADVANCED)
    , filterModeProp_("filterModeProp","Filtering:",Processor::INVALID_RESULT,false,Property::LOD_DEVELOPMENT)
    /*
    , detectStreamlineBundlesProp_("generateTubesProp", "Detect Streamline Bundles", false)
    , bundleDetectionProgressProp_("detectBundlesProgressProp", "Progress:")
    , maxAverageDistanceThresholdProp_("maxAverageDistanceThreshold", "Max. Average Distance Threshold (mm)", 1.0f, 0.0f, 100.0f)
    , minNumStreamlinesPerBundleProp_("minNumStreamlinesPerBundle", "Minimal number of Streamlines per Bundle (%)", 1.0f, 0.0f, 100.0f, Processor::INVALID_RESULT, NumericProperty<float>::STATIC, Property::LOD_ADVANCED)
    , resampleSizeProp_("resampleSize", "Streamline Resample Size", 20, 2, 100, Processor::INVALID_RESULT, NumericProperty<int>::STATIC, Property::LOD_ADVANCED)
     */
{
    volInport_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(volInport_);
    addPort(seedMask_);
    addPort(streamlineOutport_);

    addProperty(numSeedPoints_);
        numSeedPoints_.setGroupID("streamline");
    addProperty(seedTime_);
        seedTime_.setGroupID("streamline");
    addProperty(streamlineLengthThresholdProp_);
        streamlineLengthThresholdProp_.setGroupID("streamline");
    addProperty(absoluteMagnitudeThresholdProp_);
        absoluteMagnitudeThresholdProp_.onChange(MemberFunctionCallback<StreamlineCreator>(this,&StreamlineCreator::adjustRelativeThreshold));
        absoluteMagnitudeThresholdProp_.setGroupID("streamline");
    addProperty(fitAbsoluteMagnitudeProp_);
        fitAbsoluteMagnitudeProp_.setGroupID("streamline");
    addProperty(relativeMagnitudeThresholdProp_);
        relativeMagnitudeThresholdProp_.setReadOnlyFlag(true);
        relativeMagnitudeThresholdProp_.setGroupID("streamline");
    addProperty(stopIntegrationAngleThresholdProp_);
        stopIntegrationAngleThresholdProp_.setGroupID("streamline");
    addProperty(filterModeProp_);
        filterModeProp_.addOption("linear","Linear", VolumeRAM::LINEAR);
        filterModeProp_.addOption("nearest","Nearest", VolumeRAM::NEAREST);
        filterModeProp_.setGroupID("streamline");
    setPropertyGroupGuiName("streamline", "Streamline Settings");
    /*
        //streamline bundles
    addProperty(detectStreamlineBundlesProp_);
        detectStreamlineBundlesProp_.setGroupID("streamlinebundles");
    addProperty(bundleDetectionProgressProp_);
        bundleDetectionProgressProp_.setGroupID("streamlinebundles");
    addProperty(maxAverageDistanceThresholdProp_);
        maxAverageDistanceThresholdProp_.setGroupID("streamlinebundles");
    addProperty(minNumStreamlinesPerBundleProp_);
        minNumStreamlinesPerBundleProp_.setGroupID("streamlinebundles");
    addProperty(resampleSizeProp_);
        resampleSizeProp_.setGroupID("streamlinebundles");
    setPropertyGroupGuiName("streamlinebundles", "Streamline Bundle Settings");
    */
    // adjust properties
    adjustRelativeThreshold();
}


bool StreamlineCreator::isReady() const {
    if (!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if (!volInport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }

    // Note: Seed Mask is optional!

    return true;
}

std::vector<std::reference_wrapper<Port>> StreamlineCreator::getCriticalPorts() {
    auto criticalPorts = AsyncComputeProcessor<ComputeInput, ComputeOutput>::getCriticalPorts();
    criticalPorts.erase(std::remove_if(criticalPorts.begin(), criticalPorts.end(), [this] (const std::reference_wrapper<Port>& port){
        return port.get().getID() == seedMask_.getID();
    }), criticalPorts.end());
    return criticalPorts;
}

void StreamlineCreator::adjustPropertiesToInput() {

    if(!volInport_.hasData()) {
        return;
    }

    VolumeMinMaxMagnitude* data = volInport_.getData()->getDerivedData<VolumeMinMaxMagnitude>();

    absoluteMagnitudeThresholdProp_.setMinValue(data->getMinMagnitude());
    absoluteMagnitudeThresholdProp_.setMaxValue(data->getMaxMagnitude());
    if(fitAbsoluteMagnitudeProp_.get()) {
        absoluteMagnitudeThresholdProp_.set(tgt::vec2(data->getMinMagnitude(), data->getMaxMagnitude()));
    }

    /*
    tgt::vec3 length = volInport_.getData()->getSpacing() * tgt::vec3(volInport_.getData()->getDimensions());
    maxAverageDistanceThresholdProp_.setMaxValue(tgt::length(length));
    */
}

StreamlineCreatorInput StreamlineCreator::prepareComputeInput() {

    // Set up random generator.
    std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));

    auto flowVolume = volInport_.getThreadSafeData();
    if(!flowVolume) {
        throw InvalidInputException("No volume", InvalidInputException::S_ERROR);
    }
    tgt::Bounds roi = flowVolume->getBoundingBox(false).getBoundingBox(false);

    auto seedMask = seedMask_.getThreadSafeData();
    std::vector<tgt::vec3> seedPoints;
    seedPoints.reserve(numSeedPoints_.get());
    if (seedMask) {
        tgt::Bounds roiBounds = roi;
        tgt::Bounds seedMaskBounds = seedMask->getBoundingBox(false).getBoundingBox(false);

        roiBounds.intersectVolume(seedMaskBounds);
        if(!roiBounds.isDefined()) {
            throw InvalidInputException("Seed Mask does not overlap with ensemble ROI", InvalidInputException::S_ERROR);
        }

        VolumeRAMRepresentationLock seedMaskLock(seedMask);

        VolumeMinMax* vmm = seedMask->getDerivedData<VolumeMinMax>();
        if(vmm->getMinNormalized() == 0.0f && vmm->getMaxNormalized() == 0.0f) {
            throw InvalidInputException("Seed Mask is empty", InvalidInputException::S_ERROR);
        }

        tgt::mat4 seedMaskPhysicalToVoxelMatrix = seedMask->getPhysicalToVoxelMatrix();

        tgt::svec3 llf = tgt::round(seedMaskPhysicalToVoxelMatrix * roiBounds.getLLF());
        tgt::svec3 urb = tgt::round(seedMaskPhysicalToVoxelMatrix * roiBounds.getURB());

        std::vector<tgt::vec3> maskVoxels;
        for(size_t z=llf.z; z < urb.z; z++) {
            for(size_t y=llf.y; y < urb.y; y++) {
                for(size_t x=llf.x; x < urb.x; x++) {
                    if(seedMaskLock->getVoxelNormalized(x, y, z) != 0.0f) {
                        maskVoxels.push_back(tgt::vec3(x, y, z));
                    }
                }
            }
        }

        if (maskVoxels.empty()) {
            throw InvalidInputException("No seed points found in ROI", InvalidInputException::S_ERROR);
        }

        // If we have more seed mask voxel than we want to have seed points, reduce the list size.
        float probability = static_cast<float>(numSeedPoints_.get()) / maskVoxels.size();
        tgt::mat4 seedMaskVoxelToPhysicalMatrix = seedMask->getVoxelToPhysicalMatrix();
        for(const tgt::vec3& seedPoint : maskVoxels) {
            // Determine for each seed point, if we will keep it.
            if(probability >= 1.0f || rnd() < probability) {
                seedPoints.push_back(seedMaskVoxelToPhysicalMatrix * seedPoint);
            }
        }

        LINFO("Restricting seed points to volume mask using " << seedPoints.size() << " seeds");
    }
    else {
        // Without a seed mask, we uniformly sample the whole space enclosed by the roi.
        for (int k = 0; k<numSeedPoints_.get(); k++) {
            tgt::vec3 seedPoint;
            seedPoint = tgt::vec3(rnd(), rnd(), rnd());
            seedPoint = tgt::vec3(roi.getLLF()) + seedPoint * tgt::vec3(roi.diagonal());
            seedPoints.push_back(seedPoint);
        }
    }

    tgtAssert(!seedPoints.empty(), "no seed points found");
    if (seedPoints.empty()) {
        throw InvalidInputException("No seed points found", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<StreamlineListBase> output(new StreamlineList(flowVolume));

    return StreamlineCreatorInput {
        streamlineLengthThresholdProp_.get(),
        absoluteMagnitudeThresholdProp_.get(),
        stopIntegrationAngleThresholdProp_.get() * tgt::PIf / 180.0f,
        filterModeProp_.getValue(),
        volInport_.getThreadSafeData(),
        seedMask_.getThreadSafeData(),
        std::move(seedPoints),
        std::move(output)
    };
}

StreamlineCreatorOutput StreamlineCreator::compute(StreamlineCreatorInput input, ProgressReporter& progressReporter) const {

    const VolumeBase* flowVolume = input.flowVolume;
    VolumeRAMRepresentationLock representation(flowVolume);
    //const VolumeBase* seedMask_ = input.seedMask; // Currently not used.
    std::vector<tgt::vec3> seedPoints = std::move(input.seedPoints);
    std::unique_ptr<StreamlineListBase> output = std::move(input.output);

    IntegrationInput integrationInput {
        *representation,
        flowVolume->getRealWorldMapping(),
        input.filterMode,
        input.streamlineLengthThreshold,
        input.absoluteMagnitudeThreshold,
        input.stopIntegrationAngleThreshold
    };

    // TODO: parallelize, as soon as everything is working
//#ifdef VRN_MODULE_OPENMP
//#pragma omp parallel for
//#endif
    for (size_t i=0; i<seedPoints.size(); i++) {
        tgt::vec3& start = seedPoints[i];

        Streamline line = computeStreamlineRungeKutta(start, integrationInput);

        if (line.getNumElements() >= input.streamlineLengthThreshold.x &&
            line.getNumElements() <= input.streamlineLengthThreshold.y) {
            output->addStreamline(line);
        }

        progressReporter.setProgress(1.0f * i / seedPoints.size());
    }

    return StreamlineCreatorOutput {
        std::move(output)
    };
}

void StreamlineCreator::processComputeOutput(StreamlineCreatorOutput output) {
    streamlineOutport_.setData(output.streamlines.release());
}

void StreamlineCreator::adjustRelativeThreshold() {
    //adjust read only property
    tgt::vec2 range(absoluteMagnitudeThresholdProp_.getMinValue(), absoluteMagnitudeThresholdProp_.getMaxValue());
    tgt::vec2 value = absoluteMagnitudeThresholdProp_.get();

    relativeMagnitudeThresholdProp_.set(  (value - tgt::vec2(range.x)) / tgt::vec2(range.y - range.x) * 100.f);
}


Streamline StreamlineCreator::computeStreamlineRungeKutta(const tgt::vec3& start, const IntegrationInput& input) const {

    const tgt::vec3 dimAsVec3 = tgt::vec3(input.representation->getDimensions() - tgt::svec3::one);
    const float h = 0.5f; //stepwidth;
    const size_t maxNumElements = input.streamlineLengthThreshold.y;

    tgt::vec3 r(start), r_(start);
    tgt::vec3 k1(0.0f), k2(0.0f), k3(0.0f), k4(0.0f);
    tgt::vec3 k1_(0.0f), k2_(0.0f), k3_(0.0f), k4_(0.0f);
    tgt::vec3 velR = getVelocity(r, input), velR_ = getVelocity(r_, input);

    Streamline line;            // the streamline, which will be returned
    line.addElementAtEnd(Streamline::StreamlineElement(start, velR));

    bool lookupPos = true;  // integrate along the streamline in positive direction?
    bool lookupNeg = true;  // integrate along the streamline in negative direction?

    //we look in both directions at the same time
    while (line.getNumElements() < maxNumElements && (lookupPos || lookupNeg)) {

        //look in direction
        if (lookupPos) {
            //no magnitude
            if (velR == tgt::vec3::zero) {
                lookupPos = false;
            }
            else {
                // get next point
                k1 = tgt::normalize(velR) * h; //v != zero
                k2 = getVelocity(r + (k1 / 2.0f), input);
                if (k2 != tgt::vec3::zero) k2 = tgt::normalize(k2) * h;
                k3 = getVelocity(r + (k2 / 2.0f), input);
                if (k3 != tgt::vec3::zero) k3 = tgt::normalize(k3) * h;
                k4 = getVelocity(r + k3, input);
                if (k4 != tgt::vec3::zero) k4 = tgt::normalize(k4) * h;
                r += ((k1 / 6.0f) + (k2 / 3.0f) + (k3 / 3.0f) + (k4 / 6.0f));

                //is new r valid?
                if ((r != tgt::clamp(r, tgt::vec3::zero, dimAsVec3)) || r == line.getLastElement().position_) { // in case of no progress on streamline in this direction...
                    lookupPos = false;
                }
                else {//check length
                    tgt::vec3 oldVelR = velR;
                    velR = getVelocity(r, input);
                    float magnitudeR = tgt::length(velR);
                    if ((magnitudeR < input.absoluteMagnitudeThreshold.x) ||
                        (magnitudeR > input.absoluteMagnitudeThreshold.y)) {
                        lookupPos = false;
                    }
                    else if(std::acos(tgt::dot(oldVelR, velR) / (tgt::length(oldVelR) * magnitudeR)) > input.stopIntegrationAngleThreshold) {
                        lookupPos = false;
                    }
                    else {
                        line.addElementAtEnd(Streamline::StreamlineElement(r, velR));
                        if (line.getNumElements() == maxNumElements)
                            break;
                    }
                }
            }
        }

        // look previous direction
        if (lookupNeg) {
            //no magnitude
            if (velR_ == tgt::vec3::zero) {
                lookupNeg = false;
            }
            else {
                // get next point
                k1_ = tgt::normalize(velR_) * h; // velR_ != zero
                k2_ = getVelocity(r_ - (k1_ / 2.0f), input);
                if (k2_ != tgt::vec3::zero) k2_ = tgt::normalize(k2_) * h;
                k3_ = getVelocity(r_ - (k2_ / 2.0f), input);
                if (k3_ != tgt::vec3::zero) k3_ = tgt::normalize(k3_) * h;
                k4_ = getVelocity(r_ - k3_, input);
                if (k4_ != tgt::vec3::zero) k4_ = tgt::normalize(k4_) * h;
                r_ -= ((k1_ / 6.0f) + (k2_ / 3.0f) + (k3_ / 3.0f) + (k4_ / 6.0f));

                //is new r valid?
                if ((r_ != tgt::clamp(r_, tgt::vec3::zero, dimAsVec3)) || (r_ == line.getFirstElement().position_)) { // in case of no progress on streamline in this direction...
                    lookupNeg = false;
                }
                else { //check length
                    tgt::vec3 oldVelR_ = velR;
                    velR_ = getVelocity(r_, input);
                    float magnitudeR_ = tgt::length(velR_);
                    if ((magnitudeR_ < input.absoluteMagnitudeThreshold.x) ||
                        (magnitudeR_ > input.absoluteMagnitudeThreshold.y)) {
                        lookupNeg = false;
                    }
                    else if(std::acos(tgt::dot(oldVelR_, velR_) / (tgt::length(oldVelR_) * magnitudeR_)) > input.stopIntegrationAngleThreshold) {
                        lookupNeg = false;
                    }
                    else {
                        line.addElementAtFront(Streamline::StreamlineElement(r_, velR_));
                        //if (line.getNumElements() == maxNumElements)
                        //    break;
                    }
                }
            }
        }
    }

    return line;
}

tgt::vec3 StreamlineCreator::getVelocity(const tgt::vec3& pos, const IntegrationInput& input) const {

        tgt::vec3 voxel = tgt::vec3::zero;
        if(input.filterMode == VolumeRAM::NEAREST) {
            for (size_t channel = 0; channel < input.representation->getNumChannels(); channel++) {
                voxel[channel] = input.rwm.normalizedToRealWorld(
                        input.representation->getVoxelNormalized(pos, channel));
            }
        }
        else if(input.filterMode == VolumeRAM::LINEAR) {
            for (size_t channel = 0; channel < input.representation->getNumChannels(); channel++) {
                voxel[channel] = input.rwm.normalizedToRealWorld(
                        input.representation->getVoxelNormalizedLinear(pos, channel));
            }
        }
        else {
            tgtAssert(false, "unhandled filter mode")
        }
        return voxel;
    }

}   // namespace
