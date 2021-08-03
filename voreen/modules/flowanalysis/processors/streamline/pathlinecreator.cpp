/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "pathlinecreator.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "../../datastructures/streamlinelist.h"
#include "../../utils/flowutils.h"

#include <random>

namespace voreen {

const std::string PathlineCreator::loggerCat_("flowanalysis.PathlineCreator");

PathlineCreator::PathlineCreator()
    : AsyncComputeProcessor()
    , volumeListInport_(Port::INPORT, "volInport", "Flow Volume Input (vec3)")
    , seedMask_(Port::INPORT, "seedMaskPort", "Seed Mask (optional)")
    , pathlineOutport_(Port::OUTPORT, "pathlineOutport", "Pathlines Output")
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 5000, 1, 200000)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , enableReseeding_("enableReseeding", "Enable Reseeding", false)
    , reseedingInterval_("reseedingInterval", "Reseeding interval (steps)", 1, 1, 100)
    , reseedingIntervalUnitDisplay_("reseedingIntervalUnitDisplay", "Reseeding interval (ms)", 0.0f, 0.0f, 1.0f)
    , absoluteMagnitudeThreshold_("absoluteMagnitudeThreshold", "Threshold of Magnitude (absolute)", tgt::vec2(0.0f, 1000.0f), 0.0f, 9999.99f)
    , stopOutsideMask_("stopOutsideMask", "Stop outside Mask", false)
    , fitAbsoluteMagnitudeThreshold_("fitAbsoluteMagnitude", "Fit absolute Threshold to Input", false)
    , temporalResolution_("temporalResolution", "Temporal Resolution (ms)", 10.0f, 0.1f, 1000.0f)
    , filterMode_("filterModeProp", "Filtering:", Processor::INVALID_RESULT, false, Property::LOD_DEVELOPMENT)
    , velocityUnitConversion_("velocityUnitConversion", "Input Velocity Unit")
    , temporalIntegrationSteps_("temporalIntegrationSteps", "Temporal Integration Steps", 5, 1, 40, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
{
    volumeListInport_.addCondition(new PortConditionVolumeListEnsemble());
    volumeListInport_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeChannelCount(3)));
    addPort(volumeListInport_);
    addPort(seedMask_);
    addPort(pathlineOutport_);

    addProperty(numSeedPoints_);
        numSeedPoints_.setTracking(false);
        numSeedPoints_.setGroupID("pathline");
    addProperty(seedTime_);
        seedTime_.setTracking(false);
        seedTime_.setGroupID("pathline");
    addProperty(enableReseeding_);
        ON_CHANGE_LAMBDA(enableReseeding_, [this] {
            reseedingInterval_.setReadOnlyFlag(!enableReseeding_.get());
            reseedingIntervalUnitDisplay_.setReadOnlyFlag(!enableReseeding_.get());
        });
        enableReseeding_.setGroupID("pathline");
    addProperty(reseedingInterval_);
        reseedingInterval_.setReadOnlyFlag(!enableReseeding_.get());
        reseedingInterval_.setNumDecimals(2);
        reseedingInterval_.setTracking(false);
        reseedingInterval_.setGroupID("pathline");
        ON_CHANGE_LAMBDA(reseedingInterval_, [this] {
            float value = reseedingInterval_.get() * temporalResolution_.get() / temporalIntegrationSteps_.get();
            reseedingIntervalUnitDisplay_.setMinValue(value);
            reseedingIntervalUnitDisplay_.setMaxValue(value);
            reseedingIntervalUnitDisplay_.set(value);
            reseedingIntervalUnitDisplay_.adaptDecimalsToRange(3);
        });
    addProperty(reseedingIntervalUnitDisplay_);
        reseedingIntervalUnitDisplay_.setReadOnlyFlag(!enableReseeding_.get());
        reseedingIntervalUnitDisplay_.setReadOnlyFlag(true);
        reseedingIntervalUnitDisplay_.setGroupID("pathline");
    addProperty(absoluteMagnitudeThreshold_);
        absoluteMagnitudeThreshold_.setTracking(false);
        absoluteMagnitudeThreshold_.setNumDecimals(2);
        absoluteMagnitudeThreshold_.setGroupID("pathline");
    addProperty(stopOutsideMask_);
        stopOutsideMask_.setGroupID("pathline");
    addProperty(fitAbsoluteMagnitudeThreshold_);
        ON_CHANGE(fitAbsoluteMagnitudeThreshold_, PathlineCreator, adjustPropertiesToInput);
        fitAbsoluteMagnitudeThreshold_.setGroupID("pathline");
    addProperty(temporalResolution_);
        temporalResolution_.setTracking(false);
        temporalResolution_.setGroupID("pathline");
    addProperty(filterMode_);
        filterMode_.addOption("linear", "Linear", VolumeRAM::LINEAR);
        filterMode_.addOption("nearest", "Nearest", VolumeRAM::NEAREST);
        filterMode_.setGroupID("pathline");
    addProperty(velocityUnitConversion_);
        // Chose the values such that multiplying with real world values we get mm(/s)
        // which (for some reason) is the default voreen length unit.
        velocityUnitConversion_.addOption("km/s", "km/s", 1000000.0f);
        velocityUnitConversion_.addOption("m/s", "m/s", 1000.0f);
        //velocityUnitConversion_.addOption("dm/s", "dm/s", 100.0f); // Quite unusual.
        velocityUnitConversion_.addOption("cm/s", "cm/s", 10.0f);
        velocityUnitConversion_.addOption("mm/s", "mm/s", 1.0f);
        velocityUnitConversion_.set("m/s");
        velocityUnitConversion_.setGroupID("pathline");
    addProperty(temporalIntegrationSteps_);
        temporalIntegrationSteps_.setTracking(false);
        temporalIntegrationSteps_.setGroupID("pathline");
        ON_CHANGE_LAMBDA(temporalIntegrationSteps_, [this] {
            reseedingInterval_.setMaxValue(temporalIntegrationSteps_.get());
        });
    setPropertyGroupGuiName("pathline", "Pathline Settings");
}


bool PathlineCreator::isReady() const {
    if (!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if (!volumeListInport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }

    // Note: Seed Mask is optional!

    return true;
}

std::vector<std::reference_wrapper<Port>> PathlineCreator::getCriticalPorts() {
    auto criticalPorts = AsyncComputeProcessor::getCriticalPorts();
    criticalPorts.erase(std::remove_if(criticalPorts.begin(), criticalPorts.end(), [this] (const std::reference_wrapper<Port>& port){
        return port.get().getID() == seedMask_.getID();
    }), criticalPorts.end());
    return criticalPorts;
}

void PathlineCreator::adjustPropertiesToInput() {

    const VolumeList* volumeList = volumeListInport_.getData();
    if(!volumeList || volumeList->empty()) {
        return;
    }

    stopOutsideMask_.setReadOnlyFlag(!seedMask_.hasData());

    if(fitAbsoluteMagnitudeThreshold_.get()) {
        LWARNING("Calculating Min/Max Magnitudes. This may take a while...");

        float minMagnitude = std::numeric_limits<float>::max();
        float maxMagnitude = 0.0f;

        for (size_t i = 0; i < volumeList->size(); i++) {
            VolumeMinMaxMagnitude* vmmm = volumeList->at(i)->getDerivedData<VolumeMinMaxMagnitude>();
            minMagnitude = std::min(minMagnitude, vmmm->getMinMagnitude());
            maxMagnitude = std::max(maxMagnitude, vmmm->getMaxMagnitude());
        }

        absoluteMagnitudeThreshold_.setMinValue(minMagnitude);
        absoluteMagnitudeThreshold_.setMaxValue(maxMagnitude);
        absoluteMagnitudeThreshold_.set(tgt::vec2(minMagnitude, maxMagnitude));
    }
    else {
        absoluteMagnitudeThreshold_.setMinValue(0.0f);
        absoluteMagnitudeThreshold_.setMaxValue(5000.0f);
    }
}

PathlineCreatorInput PathlineCreator::prepareComputeInput() {

    auto flowVolumes = volumeListInport_.getThreadSafeData();
    if(!flowVolumes) {
        throw InvalidInputException("No input", InvalidInputException::S_ERROR);
    }

    if(flowVolumes->size() < 2) {
        throw InvalidInputException("Need at least two time steps", InvalidInputException::S_ERROR);
    }

    // Set up random generator.
    std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));

    const VolumeBase* referenceVolume = flowVolumes->first();
    VolumeRAMRepresentationLock reference(referenceVolume);

    tgt::Bounds roi = referenceVolume->getBoundingBox().getBoundingBox();
    auto numSeedPoints = static_cast<size_t>(numSeedPoints_.get());

    auto seedMask = seedMask_.getData();
    std::vector<tgt::vec3> seedPoints;
    seedPoints.reserve(numSeedPoints);
    if (seedMask) {
        tgt::Bounds seedMaskBounds = seedMask->getBoundingBox().getBoundingBox();

        roi.intersectVolume(seedMaskBounds);
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
    else  {
        // Without a seed mask, we uniformly sample the whole space enclosed by the roi.
        for (size_t k = 0; k < numSeedPoints; k++) {
            // Since argument evaluation order is unspecified in c++, we need to ensure the order manually.
            float x = rnd();
            float y = rnd();
            float z = rnd();

            tgt::vec3 seedPoint(x, y, z);
            seedPoint = roi.getLLF() + seedPoint * roi.diagonal();
            seedPoints.push_back(seedPoint);
        }
    }

    if (seedPoints.empty()) {
        throw InvalidInputException("No seed points found", InvalidInputException::S_ERROR);
    }

    auto mask = stopOutsideMask_.get() ? seedMask_.getData() : nullptr;
    std::unique_ptr<StreamlineListBase> output(new StreamlineList(referenceVolume));

    return PathlineCreatorInput {
            absoluteMagnitudeThreshold_.get(),
            velocityUnitConversion_.getValue(),
            temporalResolution_.get() / 1000.0f, // Convert to s.
            temporalIntegrationSteps_.get(),
            enableReseeding_.get(),
            reseedingInterval_.get(),
            filterMode_.getValue(),
            std::move(flowVolumes),
            mask,
            std::move(seedPoints),
            std::move(output)
    };
}

PathlineCreatorOutput PathlineCreator::compute(PathlineCreatorInput input, ProgressReporter& progressReporter) const {

    // Input.
    auto flowVolumes = std::move(input.flowVolumes);
    const VolumeBase* referenceVolume = flowVolumes->first();
    const VolumeBase* seedMask = input.seedMask;
    const auto seedPoints = std::move(input.seedPoints);

    // Output.
    std::unique_ptr<StreamlineListBase> output = std::move(input.output);

    // Temp. requirements.
    const tgt::mat4 worldToVoxelMatrix = referenceVolume->getWorldToVoxelMatrix();
    const tgt::Bounds roi = referenceVolume->getBoundingBox().getBoundingBox();
    const RealWorldMapping rwm = referenceVolume->getRealWorldMapping();

    const float totalTime = input.temporalResolution * (flowVolumes->size() - 1);
    const float dt = input.temporalResolution / input.temporalIntegrationSteps;

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

    const IntegrationInput integrationInput {
            dt,
            input.velocityUnitConversion,
            input.absoluteMagnitudeThreshold,
            bounds
    };

    VolumeRAMRepresentationLock currVol(flowVolumes->first());
    if(!*currVol) {
        LERROR("RAM representation not available");
        return PathlineCreatorOutput { nullptr };
    }

    std::list<Streamline> pathlines;
    for(size_t i=0; i<flowVolumes->size() - 1; i++) {

        // We ensure the current and next time frame has a RAM representation.
        VolumeRAMRepresentationLock nextVol(flowVolumes->at(i+1));
        if(!*nextVol) {
            LERROR("RAM representation not available");
            return PathlineCreatorOutput { nullptr };
        }

        // Temporal integration loop.
        for(int t=0; t<input.temporalIntegrationSteps; t++) {
            float alpha = t * dt / input.temporalResolution;
            SpatioTemporalSampler sampler(*currVol, *nextVol, alpha, rwm, input.filterMode, worldToVoxelMatrix);

            // Seeding. (If first step or reseeding interval is covered.
            if((i==0 && t==0) || (input.enableReseeding && t % input.reseedingSteps == 0)) {
                for (const tgt::vec3& seedPoint : seedPoints) {
                    tgt::vec3 velocity = sampler.sample(seedPoint);

                    Streamline pathline;
                    pathline.addElementAtEnd(Streamline::StreamlineElement(seedPoint, velocity, 0.0f, i*input.temporalResolution+t*dt));
                    pathlines.emplace_back(pathline);
                }
            }

            // Tracing.
            for (auto iter = pathlines.begin(); iter != pathlines.end();) {
                Streamline& pathline = *iter;
                bool continueIntegration = integrationStep(pathline, sampler, integrationInput);
                if (!continueIntegration) {
                    output->addStreamline(pathline);
                    iter = pathlines.erase(iter);
                }
                else {
                    iter++;
                }
            }

            std::swap(currVol, nextVol);
            progressReporter.setProgress((i * input.temporalIntegrationSteps + t) * dt / totalTime);
        }
    }

    // Add remaining path lines (Those have the full length).
    for(const Streamline& pathline : pathlines) {
        output->addStreamline(pathline);
    }

    progressReporter.setProgress(1.0f);

    return PathlineCreatorOutput {
            std::move(output)
    };
}

void PathlineCreator::processComputeOutput(PathlineCreatorOutput output) {
    pathlineOutport_.setData(output.pathlines.release());
}

bool PathlineCreator::integrationStep(Streamline& pathline, const SpatioTemporalSampler& sampler, const IntegrationInput& input) const {

    const float epsilon = 1e-5f; // std::numeric_limits<float>::epsilon() is not enough.

    // Position.
    tgt::vec3 r = pathline.getLastElement().position_;

    // Velocity.
    tgt::vec3 v = pathline.getLastElement().velocity_;

    // Execute 4th order Runge-Kutta step.
    tgt::vec3 k1 = v * input.stepSize * input.velocityUnitConversion;
    tgt::vec3 k2 = sampler.sample(r + (k1 / 2.0f)) * input.stepSize * input.velocityUnitConversion;
    tgt::vec3 k3 = sampler.sample(r + (k2 / 2.0f)) * input.stepSize * input.velocityUnitConversion;
    tgt::vec3 k4 = sampler.sample(r + k3) * input.stepSize * input.velocityUnitConversion;
    tgt::vec3 dr = (k1 / 6.0f) + (k2 / 3.0f) + (k3 / 3.0f) + (k4 / 6.0f);

    // Note: dr can be zero in one frame, but might change in a following frame.

    // Update position.
    r += dr;

    // Ran out of bounds?
    if(!input.bounds(r)) {
        return false;
    }

    // Sample local velocity.
    v = sampler.sample(r);
    float magnitude = tgt::length(v);

    // Magnitude within range?
    if(magnitude < input.absoluteMagnitudeThreshold.x - epsilon ||
        magnitude > input.absoluteMagnitudeThreshold.y + epsilon) {
        return false;
    }

    // Add element to pathline.
    float time = pathline.getFirstElement().time_ + pathline.getNumElements() * input.stepSize;
    pathline.addElementAtEnd(Streamline::StreamlineElement(r, v, 0.0f, time));

    return true;
}

}   // namespace
