/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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
#include "../../utils/flowutils.h"

#include <random>

namespace voreen {

const std::string StreamlineCreator::loggerCat_("flowanalysis.StreamlineCreator");

StreamlineCreator::StreamlineCreator()
    : AsyncComputeProcessor()
    , volumeInport_(Port::INPORT, "volInport", "Flow Volume Input (vec3)")
    , seedMask_(Port::INPORT, "seedMaskPort", "Seed Mask (optional)")
    , streamlineOutport_(Port::OUTPORT, "streamlineOutport", "Streamlines Output")
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 5000, 1, 200000)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , streamlineLengthThreshold_("streamlineLengthThreshold", "Restrict streamline length", tgt::ivec2(10, 1000), 2, 10000)
    , absoluteMagnitudeThreshold_("absoluteMagnitudeThreshold", "Threshold of Magnitude (absolute)", tgt::vec2(0.0f, 1000.0f), 0.0f, 9999.99f)
    , stopOutsideMask_("stopOutsideMask", "Stop outside Mask", false)
    , fitAbsoluteMagnitudeThreshold_("fitAbsoluteMagnitude", "Fit absolute Threshold to Input", false)
    , stopIntegrationAngleThreshold_("stopIntegrationAngleThreshold", "Stop Integration on Angle", 180, 0, 180, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_ADVANCED)
    , filterMode_("filterModeProp", "Filtering:", Processor::INVALID_RESULT, false, Property::LOD_DEVELOPMENT)
    , transformVelocities_("transformVelocities", "Transform Velocities", false, InvalidationLevel::INVALID_RESULT, Property::LOD_ADVANCED)
    , velocityUnitConversion_("velocityUnitConversion", "Input Velocity Unit")
    , integrationSteps_("integrationSteps", "Integration Steps", 100, 1, 10000, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
{
    volumeInport_.addCondition(new PortConditionVolumeChannelCount(3));
    addPort(volumeInport_);
    addPort(seedMask_);
    addPort(streamlineOutport_);

    addProperty(numSeedPoints_);
        numSeedPoints_.setTracking(false);
        numSeedPoints_.setGroupID("streamline");
    addProperty(seedTime_);
        seedTime_.setTracking(false);
        seedTime_.setGroupID("streamline");
    addProperty(streamlineLengthThreshold_);
        streamlineLengthThreshold_.setTracking(false);
        streamlineLengthThreshold_.setGroupID("streamline");
    addProperty(absoluteMagnitudeThreshold_);
        absoluteMagnitudeThreshold_.setTracking(false);
        absoluteMagnitudeThreshold_.setNumDecimals(2);
        absoluteMagnitudeThreshold_.setGroupID("streamline");
    addProperty(stopOutsideMask_);
        stopOutsideMask_.setGroupID("streamline");
    addProperty(fitAbsoluteMagnitudeThreshold_);
        ON_CHANGE(fitAbsoluteMagnitudeThreshold_, StreamlineCreator, adjustPropertiesToInput);
        fitAbsoluteMagnitudeThreshold_.setGroupID("streamline");
    addProperty(stopIntegrationAngleThreshold_);
        stopIntegrationAngleThreshold_.setTracking(false);
        stopIntegrationAngleThreshold_.setGroupID("streamline");
    addProperty(filterMode_);
        filterMode_.addOption("linear", "Linear", VolumeRAM::LINEAR);
        filterMode_.addOption("nearest", "Nearest", VolumeRAM::NEAREST);
        filterMode_.setGroupID("streamline");
    addProperty(transformVelocities_);
        transformVelocities_.setGroupID("streamline");
    addProperty(velocityUnitConversion_);
        // Chose the values such that multiplying with real world values we get mm(/s)
        // which (for some reason) is the default voreen length unit.
        velocityUnitConversion_.addOption("km/s", "km/s", 1000000.0f);
        velocityUnitConversion_.addOption("m/s", "m/s", 1000.0f);
        //velocityUnitConversion_.addOption("dm/s", "dm/s", 100.0f); // Quite unusual.
        velocityUnitConversion_.addOption("cm/s", "cm/s", 10.0f);
        velocityUnitConversion_.addOption("mm/s", "mm/s", 1.0f);
        velocityUnitConversion_.set("m/s");
        velocityUnitConversion_.setGroupID("streamline");
    addProperty(integrationSteps_);
        integrationSteps_.setTracking(false);
        integrationSteps_.setGroupID("streamline");
    setPropertyGroupGuiName("streamline", "Streamline Settings");
}


bool StreamlineCreator::isReady() const {
    if (!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if (!volumeInport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }

    // Note: Seed Mask is optional!

    return true;
}

void StreamlineCreator::adjustPropertiesToInput() {

    const VolumeBase* volume = volumeInport_.getData();
    if(!volume) {
        return;
    }

    if(fitAbsoluteMagnitudeThreshold_.get()) {
        VolumeMinMaxMagnitude* data = volume->getDerivedData<VolumeMinMaxMagnitude>();

        absoluteMagnitudeThreshold_.setMinValue(data->getMinMagnitude());
        absoluteMagnitudeThreshold_.setMaxValue(data->getMaxMagnitude());
        absoluteMagnitudeThreshold_.set(tgt::vec2(data->getMinMagnitude(), data->getMaxMagnitude()));
    }
    else {
        absoluteMagnitudeThreshold_.setMinValue(0.0f);
        absoluteMagnitudeThreshold_.setMaxValue(5000.0f);
    }
}

StreamlineCreatorInput StreamlineCreator::prepareComputeInput() {

    auto flowVolume = volumeInport_.getThreadSafeData();
    if(!flowVolume) {
        throw InvalidInputException("No volume", InvalidInputException::S_ERROR);
    }

    // Set up random generator.
    std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));

    tgt::Bounds roi = flowVolume->getBoundingBox().getBoundingBox();
    auto numSeedPoints = static_cast<size_t>(numSeedPoints_.get());

    auto seedMask = seedMask_.getData();
    std::vector<tgt::vec3> seedPoints;
    seedPoints.reserve(numSeedPoints_.get());
    if (seedMask) {
        roi.intersectVolume(seedMask->getBoundingBox().getBoundingBox());
        if(!roi.isDefined()) {
            throw InvalidInputException("Seed Mask does not overlap with ensemble ROI", InvalidInputException::S_ERROR);
        }

        VolumeRAMRepresentationLock seedMaskLock(seedMask);
        tgt::mat4 seedMaskVoxelToWorldMatrix = seedMask->getVoxelToWorldMatrix();
        tgt::svec3 dim = seedMaskLock->getDimensions();
        for(size_t z=0; z < dim.z; z++) {
            for(size_t y=0; y < dim.y; y++) {
                for(size_t x=0; x < dim.x; x++) {
                    if(seedMaskLock->getVoxelNormalized(x, y, z) != 0.0f) {
                        tgt::vec3 pos = seedMaskVoxelToWorldMatrix * tgt::vec3(x, y, z);
                        if(roi.containsPoint(pos)) {
                            seedPoints.emplace_back(pos);
                        }
                    }
                }
            }
        }

        if (seedPoints.empty()) {
            throw InvalidInputException("No seed points found in ROI", InvalidInputException::S_ERROR);
        }

        tgt::shuffle(seedPoints.begin(), seedPoints.end(), std::mt19937(seedTime_.get()));
        seedPoints.resize(std::min(seedPoints.size(), numSeedPoints));

        LINFO("Restricting seed points to volume mask using " << seedPoints.size() << " seeds");
    }
    else {
        // Without a seed mask, we uniformly sample the whole space enclosed by the roi.
        for (size_t k = 0; k<numSeedPoints; k++) {
            // Since argument evaluation order is unspecified in c++, we need to ensure the order manually.
            float x = rnd();
            float y = rnd();
            float z = rnd();

            tgt::vec3 seedPoint(x, y, z);
            seedPoint = roi.getLLF() + seedPoint * roi.diagonal();
            seedPoints.push_back(seedPoint);
        }
    }

    tgtAssert(!seedPoints.empty(), "no seed points found");
    if (seedPoints.empty()) {
        throw InvalidInputException("No seed points found", InvalidInputException::S_ERROR);
    }

    auto mask = stopOutsideMask_.get() ? seedMask_.getData() : nullptr;
    std::unique_ptr<StreamlineListBase> output(new StreamlineList(flowVolume));

    return StreamlineCreatorInput {
            streamlineLengthThreshold_.get(),
            absoluteMagnitudeThreshold_.get(),
            velocityUnitConversion_.getValue(),
            integrationSteps_.get(),
            stopIntegrationAngleThreshold_.get() * tgt::PIf / 180.0f,
            filterMode_.getValue(),
            transformVelocities_.get(),
            volumeInport_.getThreadSafeData(),
            mask,
            std::move(seedPoints),
            std::move(output)
    };
}

StreamlineCreatorOutput StreamlineCreator::compute(StreamlineCreatorInput input, ProgressReporter& progressReporter) const {

    PortDataPointer<VolumeBase> flowVolume = std::move(input.flowVolume);
    VolumeRAMRepresentationLock representation(flowVolume);
    const VolumeBase* seedMask = input.seedMask;
    std::vector<tgt::vec3> seedPoints = std::move(input.seedPoints);

    size_t lowerLengthThreshold = input.streamlineLengthThreshold.x;
    size_t upperLengthThreshold = input.streamlineLengthThreshold.y;

    // Temp. requirements.
    const tgt::mat4 worldToVoxelMatrix = flowVolume->getWorldToVoxelMatrix();
    const tgt::mat4 velocityTransformationMatrix = input.transformVelocities ? flowVolume->getPhysicalToWorldMatrix().getRotationalPart() : tgt::mat4::identity;
    const tgt::Bounds roi = flowVolume->getBoundingBox().getBoundingBox();
    const RealWorldMapping rwm = flowVolume->getRealWorldMapping();

    const float stepSize = 1.0f / input.integrationSteps;

    std::function<bool(const tgt::vec3&)> bounds = [&] (const tgt::vec3& position) {
        return roi.containsPoint(position);
    };

    if(seedMask) {
        tgt::mat4 worldToVoxel = seedMask->getWorldToVoxelMatrix();
        VolumeRAMRepresentationLock lock(seedMask);
        bounds = [=] (const tgt::vec3& position) {
            if(!roi.containsPoint(position)) {
                return false;
            }
            return lock->getVoxelNormalized(worldToVoxel * position) != 0.0f;
        };
    }

    const IntegrationInput integrationInput{
            stepSize,
            input.velocityUnitConversion,
            upperLengthThreshold * input.integrationSteps,
            input.absoluteMagnitudeThreshold,
            input.stopIntegrationAngleThreshold,
            bounds
    };

    const SpatialSampler sampler(*representation, flowVolume->getRealWorldMapping(), input.filterMode, worldToVoxelMatrix, velocityTransformationMatrix);

    ThreadedTaskProgressReporter progress(progressReporter, seedPoints.size());
    bool aborted = false;

    std::vector<Streamline> streamlines(seedPoints.size());

#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
    for (long i=0; i<static_cast<long>(seedPoints.size()); i++) {
        if (aborted) {
            continue;
        }
#else
    for(size_t i=0; i<seedPoints.size(); i++) {
#endif
        streamlines[i] = integrateStreamline(seedPoints[i], sampler, integrationInput);

        if (progress.reportStepDone()) {
#ifdef VRN_MODULE_OPENMP
            #pragma omp critical
            aborted = true;
#else
            aborted = true;
            break;
#endif
        }
    }

    if (aborted) {
        throw boost::thread_interrupted();
    }

    // Add sequentially to output to guarantee deterministic results.
    std::unique_ptr<StreamlineListBase> output = std::move(input.output);
    for(const Streamline& streamline : streamlines) {
        if (streamline.getNumElements() >= lowerLengthThreshold) {
            output->addStreamline(streamline);
        }
    }

    return StreamlineCreatorOutput {
        std::move(output)
    };
}

void StreamlineCreator::processComputeOutput(StreamlineCreatorOutput output) {
    streamlineOutport_.setData(output.streamlines.release());
}

Streamline StreamlineCreator::integrateStreamline(const tgt::vec3& start, const SpatialSampler& sampler, const IntegrationInput& input) const {

    const float epsilon = 1e-5f; // std::numeric_limits<float>::epsilon() is not enough.

    // Position.
    tgt::vec3 r(start);
    tgt::vec3 r_(start);

    // Velocity.
    tgt::vec3 v = sampler.sample(r);
    tgt::vec3 v_ = v;

    // Return an empty line in case the initial velocity was zero already.
    if(v == tgt::vec3::zero) {
        return Streamline();
    }

    // Resulting streamline.
    Streamline line;
    line.addElementAtEnd(Streamline::StreamlineElement(r, v));

    bool lookupPositiveDirection = true;
    bool lookupNegativeDirection = true;

    // Look up positive and negative direction in alternating fashion.
    while (lookupPositiveDirection || lookupNegativeDirection) {

        if (lookupPositiveDirection) {

            // Execute 4th order Runge-Kutta step.
            tgt::vec3 k1 = v * input.stepSize * input.velocityUnitConversion;
            tgt::vec3 k2 = sampler.sample(r + (k1 / 2.0f)) * input.stepSize * input.velocityUnitConversion;;
            tgt::vec3 k3 = sampler.sample(r + (k2 / 2.0f)) * input.stepSize * input.velocityUnitConversion;
            tgt::vec3 k4 = sampler.sample(r + k3) * input.stepSize * input.velocityUnitConversion;
            tgt::vec3 dr = ((k1 / 6.0f) + (k2 / 3.0f) + (k3 / 3.0f) + (k4 / 6.0f));

            // Update position.
            r += dr;

            // Check constrains.
            lookupPositiveDirection &= input.bounds(r);
            lookupPositiveDirection &= (r != line.getLastElement().position_); // Progress in current direction?

            v = sampler.sample(r);
            float magnitude = tgt::length(v);
            lookupPositiveDirection &= (v != tgt::vec3::zero)
                    && (magnitude > input.absoluteMagnitudeThreshold.x - epsilon)
                    && (magnitude < input.absoluteMagnitudeThreshold.y + epsilon)
                    && std::acos(std::abs(tgt::dot(line.getLastElement().velocity_, v)) /
                                 (tgt::length(line.getLastElement().velocity_) * magnitude)) <= input.stopIntegrationAngleThreshold;

            if (lookupPositiveDirection) {
                line.addElementAtEnd(Streamline::StreamlineElement(r, v));
                if (line.getNumElements() >= input.upperLengthThreshold) {
                    break;
                }
            }
        }

        if (lookupNegativeDirection) {

            // Execute 4th order Runge-Kutta step.
            tgt::vec3 k1 = v_ * input.stepSize * input.velocityUnitConversion;
            tgt::vec3 k2 = sampler.sample(r_ + (k1 / 2.0f)) * input.stepSize * input.velocityUnitConversion;;
            tgt::vec3 k3 = sampler.sample(r_ + (k2 / 2.0f)) * input.stepSize * input.velocityUnitConversion;
            tgt::vec3 k4 = sampler.sample(r_ + k3) * input.stepSize * input.velocityUnitConversion;
            tgt::vec3 dr_ = ((k1 / 6.0f) + (k2 / 3.0f) + (k3 / 3.0f) + (k4 / 6.0f));

            // Update position.
            r_ -= dr_;

            // Check constrains.
            lookupNegativeDirection &= input.bounds(r_);
            lookupNegativeDirection &= (r_ != line.getFirstElement().position_); // Progress in current direction?

            v_ = sampler.sample(r_);
            float magnitude = tgt::length(v_);
            lookupNegativeDirection &= (v_ != tgt::vec3::zero)
                    && (magnitude > input.absoluteMagnitudeThreshold.x - epsilon)
                    && (magnitude < input.absoluteMagnitudeThreshold.y + epsilon)
                    && std::acos(std::abs(tgt::dot(line.getFirstElement().velocity_, v_)) /
                                 (tgt::length(line.getFirstElement().velocity_) * magnitude)) <= input.stopIntegrationAngleThreshold;

            if (lookupNegativeDirection) {
                line.addElementAtFront(Streamline::StreamlineElement(r_, v_));
                if (line.getNumElements() >= input.upperLengthThreshold) {
                    break;
                }
            }
        }
    }

    return line;
}

}   // namespace
