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

#include "pathlinecreator.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "../../datastructures/streamlinelist.h"
#include "../../utils/flowutils.h"

#include <random>

namespace voreen {

const std::string PathlineCreator::loggerCat_("flowreen.PathlineCreator");

PathlineCreator::PathlineCreator()
    : AsyncComputeProcessor()
    , volumeListInport_(Port::INPORT, "volInport", "Flow Volume Input (vec3)")
    , seedMask_(Port::INPORT, "seedMaskPort", "Seed Mask (optional)")
    , pathlineOutport_(Port::OUTPORT, "pathlineOutport", "Pathlines Output")
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 5000, 1, 200000)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , absoluteMagnitudeThreshold_("absoluteMagnitudeThreshold", "Threshold of Magnitude (absolute)", tgt::vec2(0.0f, 1000.0f), 0.0f, 9999.99f)
    , fitAbsoluteMagnitudeThreshold_("fitAbsoluteMagnitude", "Fit absolute Threshold to Input", false)
    , temporalResolution_("temporalResolution", "Temporal Resolution (ms)", 3.1f, 0.1f, 1000.0f)
    , filterMode_("filterModeProp", "Filtering:", Processor::INVALID_RESULT, false, Property::LOD_DEVELOPMENT)
    , velocityUnitConversion_("velocityUnitConversion", "Input Velocity Unit")
    , temporalIntegrationSteps_("temporalIntegrationSteps", "Temporal Integration Steps",5, 1, 20, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
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
    addProperty(absoluteMagnitudeThreshold_);
        absoluteMagnitudeThreshold_.setTracking(false);
        absoluteMagnitudeThreshold_.setNumDecimals(2);
        absoluteMagnitudeThreshold_.setGroupID("pathline");
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
        //velocityUnitConversion_.addOption("dm/s", "dm/s", 100.0f); // Really unusual.
        velocityUnitConversion_.addOption("cm/s", "cm/s", 10.0f);
        velocityUnitConversion_.addOption("mm/s", "mm/s", 1.0f);
        velocityUnitConversion_.set("m/s");
        velocityUnitConversion_.setGroupID("pathline");
    addProperty(temporalIntegrationSteps_);
        temporalIntegrationSteps_.setTracking(false);
        temporalIntegrationSteps_.setGroupID("pathline");
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

    const VolumeBase* referenceVolume = flowVolumes->first();
    VolumeRAMRepresentationLock reference(referenceVolume);

    tgt::mat4 worldToVoxelMatrix = referenceVolume->getPhysicalToVoxelMatrix();
    tgt::Bounds roi = referenceVolume->getBoundingBox().getBoundingBox();
    RealWorldMapping rwm = referenceVolume->getRealWorldMapping();
    rwm.setScale(rwm.getScale() * velocityUnitConversion_.getValue()); // Now we have mm/s.
    SpatialSampler sampler(*reference, rwm, filterMode_.getValue(), worldToVoxelMatrix);

    // Set up random generator.
    std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));

    std::list<Streamline> pathlines;
    auto seedMask = seedMask_.getData();
    if (seedMask) {
        tgt::Bounds roiBounds = roi;
        tgt::Bounds seedMaskBounds = seedMask->getBoundingBox().getBoundingBox();

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
        tgt::mat4 seedMaskVoxelToWorldMatrix = seedMask->getVoxelToWorldMatrix();
        for(const tgt::vec3& seedPoint : maskVoxels) {
            // Determine for each seed point, if we will keep it.
            if(probability >= 1.0f || rnd() < probability) {
                tgt::vec3 position = worldToVoxelMatrix * seedMaskVoxelToWorldMatrix * seedPoint;
                tgt::vec3 velocity = sampler.sample(position);

                Streamline pathline;
                pathline.addElementAtEnd(Streamline::StreamlineElement(position, velocity));
                pathlines.push_back(pathline);
            }
        }

        LINFO("Restricting seed points to volume mask using " << pathlines.size() << " seeds");
    }
    else  {
        // Without a seed mask, we uniformly sample the whole space enclosed by the roi.
        for (int k = 0; k < numSeedPoints_.get(); k++) {
            tgt::vec3 position(rnd(), rnd(), rnd());
            position = roi.getLLF() + position * roi.diagonal();
            tgt::vec3 velocity = sampler.sample(position);

            Streamline pathline;
            pathline.addElementAtEnd(Streamline::StreamlineElement(position, velocity));
            pathlines.push_back(pathline);
        }
    }

    if (pathlines.empty()) {
        throw InvalidInputException("No seed points found", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<StreamlineListBase> output(new StreamlineList(referenceVolume));

    return PathlineCreatorInput {
            absoluteMagnitudeThreshold_.get(),
            velocityUnitConversion_.getValue(),
            temporalResolution_.get() / 1000.0f, // Convert to s.
            temporalIntegrationSteps_.get(),
            filterMode_.getValue(),
            std::move(flowVolumes),
            seedMask_.getThreadSafeData(),
            std::move(pathlines),
            std::move(output)
    };
}

PathlineCreatorOutput PathlineCreator::compute(PathlineCreatorInput input, ProgressReporter& progressReporter) const {

    // Input.
    auto flowVolumes = std::move(input.flowVolumes);
    const VolumeBase* referenceVolume = flowVolumes->first();
    //const VolumeBase* seedMask_ = input.seedMask; // Currently not used.

    // Output.
    std::list<Streamline> pathlines = std::move(input.pathlines);
    std::unique_ptr<StreamlineListBase> output = std::move(input.output);

    // Temp. requirements.
    RealWorldMapping rwm = referenceVolume->getRealWorldMapping();     // This maps to some unknown unit per second.
    rwm.setScale(rwm.getScale() * input.velocityUnitConversion); // Now we have mm/s.
    tgt::mat4 worldToVoxelMatrix = referenceVolume->getWorldToVoxelMatrix();
    tgt::Bounds roi = referenceVolume->getBoundingBox().getBoundingBox();

    const float totalTime = input.temporalResolution * (flowVolumes->size() - 1);
    const float dt = input.temporalResolution / input.temporalIntegrationSteps;

    const IntegrationInput integrationInput {
            dt,
            roi,
            input.absoluteMagnitudeThreshold * input.velocityUnitConversion,
    };

    for(size_t i=0; i<flowVolumes->size() - 1; i++) {

        // We ensure the current and next time frame has as RAM representation.
        VolumeRAMRepresentationLock volume0(flowVolumes->at(i+0));
        VolumeRAMRepresentationLock volume1(flowVolumes->at(i+1));

        // Temporal integration loop. Here we actually change the field.
        for(int t=0; t<input.temporalIntegrationSteps; t++) {
            float alpha = t * dt / input.temporalResolution;
            SpatioTemporalSampler sampler(*volume0, *volume1, alpha, rwm, input.filterMode, worldToVoxelMatrix);

            // Iterate pathlines in reverse, such that we can remove at the end.
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
    tgt::vec3 velR = pathline.getLastElement().velocity_;

    // Execute 4th order Runge-Kutta step.
    tgt::vec3 k1 = velR * input.stepSize;
    tgt::vec3 k2 = sampler.sample(r + (k1 / 2.0f)) * input.stepSize;
    tgt::vec3 k3 = sampler.sample(r + (k2 / 2.0f)) * input.stepSize;
    tgt::vec3 k4 = sampler.sample(r + k3) * input.stepSize;
    tgt::vec3 dr = (k1 / 6.0f) + (k2 / 3.0f) + (k3 / 3.0f) + (k4 / 6.0f);

    // Note: dr can be zero in one frame, but might change in a following frame.

    // Update position.
    r += dr;

    // Ran out of bounds?
    if(!input.bounds.containsPoint(r)) {
        return false;
    }

    // Sample local velocity.
    velR = sampler.sample(r);
    float magnitude = tgt::length(velR);

    // Magnitude within range?
    if(magnitude < input.absoluteMagnitudeThreshold.x - epsilon ||
        magnitude > input.absoluteMagnitudeThreshold.y + epsilon) {
        return false;
    }

    // Add element to pathline.
    pathline.addElementAtEnd(Streamline::StreamlineElement(r, velR, 0.0f, pathline.getNumElements() * input.stepSize));

    return true;
}

}   // namespace
