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

#include "referencevolumecreator.h"

#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/utils/statistics.h"
#include "../utils/utils.h"

namespace voreen {

const std::string ReferenceVolumeCreator::loggerCat_("voreen.ensembleanalysis.ReferenceVolumeCreator");

ReferenceVolumeCreator::ReferenceVolumeCreator()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , selectedField_("selectedField", "Selected Field")
    , time_("time", "Time", 0.0f, 0.0f, 1000000.0f)
    , referenceMethod_("referenceMethode", "Reference Method")
    , referenceRun_("referenceRun", "Reference Run")
    , sampleRegion_("sampleRegion", "Sample Region")
    , outputDimensions_("outputDimensions", "Output Dimensions", tgt::ivec3(200), tgt::ivec3(1), tgt::ivec3(1000))
{
    // Ports
    addPort(inport_);
    addPort(outport_);

    addProperty(selectedField_);
    addProperty(time_);
    time_.setTracking(false);
    addProperty(referenceMethod_);
    referenceMethod_.addOption("run", "Select Run");
    referenceMethod_.addOption("zero", "Zero volume");
    referenceMethod_.addOption("mean", "Global Mean");
    //referenceMethod_.addOption("median", "Global Median");
    ON_CHANGE_LAMBDA(referenceMethod_, [this] {
        referenceRun_.setReadOnlyFlag(referenceMethod_.get() != "run");
    });
    addProperty(referenceRun_);
    addProperty(sampleRegion_);
    sampleRegion_.addOption("bounds", "Ensemble Bounds");
    sampleRegion_.addOption("common", "Common Bounds");
    sampleRegion_.addOption("roi", "Region of Interest");
    addProperty(outputDimensions_);
}

ReferenceVolumeCreator::~ReferenceVolumeCreator() {
}

ReferenceVolumeCreatorInput ReferenceVolumeCreator::prepareComputeInput() {
    PortDataPointer<EnsembleDataset> ensemble = inport_.getThreadSafeData();
    if (!ensemble) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    // Get required information about reference volume format.
    tgt::ivec3 newDims = outputDimensions_.get();
    size_t numChannels = ensemble->getNumChannels(selectedField_.get());
    const std::string& baseType = ensemble->getBaseType(selectedField_.get());

    // Create output volume.
    VolumeFactory factory;
    std::string format = factory.getFormat(baseType, numChannels);
    std::unique_ptr<VolumeRAM> outputVolume(factory.create(format, newDims));
    outputVolume->clear(); // Sets all to zero.

    tgt::Bounds bounds;
    if (sampleRegion_.get() == "bounds") {
        bounds = ensemble->getBounds();
    }
    else if (sampleRegion_.get() == "common") {
        bounds = ensemble->getCommonBounds();
    }
    else if (sampleRegion_.get() == "roi") {
        bounds = ensemble->getRoi();
    }
    else {
        throw InvalidInputException("Unknown sample region", InvalidInputException::S_ERROR);
    }

    return ReferenceVolumeCreatorInput{
            std::move(ensemble),
            std::move(outputVolume),
            bounds,
            selectedField_.get(),
            time_.get(),
            referenceMethod_.get(),
            static_cast<size_t>(referenceRun_.getSelectedIndex())
    };
}

ReferenceVolumeCreatorOutput ReferenceVolumeCreator::compute(ReferenceVolumeCreatorInput input, ProgressReporter& progress) const {

    auto ensemble = std::move(input.ensemble);
    const tgt::Bounds& bounds = input.bounds;
    const std::string& field = input.field;
    const size_t numRuns = ensemble->getRuns().size();
    const size_t numChannels = ensemble->getNumChannels(field);
    std::unique_ptr<VolumeRAM> output = std::move(input.outputVolume);
    const tgt::ivec3 newDims = output->getDimensions();

    if(input.referenceMethod == "run") {

        size_t referenceTimeStep = ensemble->pickTimeStep(input.referenceRun, input.time);
        const VolumeBase* refVolume = ensemble->getRuns()[input.referenceRun].timeSteps_[referenceTimeStep].fieldNames_.at(input.field);
        VolumeRAMRepresentationLock referenceVolume(refVolume);
        tgt::mat4 refPhysicalToVoxel = refVolume->getPhysicalToVoxelMatrix();

        tgt::ivec3 pos = tgt::ivec3::zero;
        for (pos.z = 0; pos.z < newDims.z; ++pos.z) {
            for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
                for (pos.x = 0; pos.x < newDims.x; ++pos.x) {

                    // Map sample position to physical space.
                    tgt::vec3 sample = mapRange(tgt::vec3(pos), tgt::vec3::zero, tgt::vec3(newDims), bounds.getLLF(), bounds.getURB());

                    // Map to voxel space.
                    tgt::svec3 sampleInRefVoxelSpace = refPhysicalToVoxel * sample;

                    // Ignore, if out of bounds.
                    if(tgt::clamp(sampleInRefVoxelSpace, tgt::svec3::zero, referenceVolume->getDimensions() - tgt::svec3::one) != sampleInRefVoxelSpace) {
                        continue;
                    }

                    // Sample the volume.
                    for(size_t channel=0; channel<numChannels; channel++) {
                        float value = referenceVolume->getVoxelNormalized(sampleInRefVoxelSpace, channel);
                        output->setVoxelNormalized(value, pos, channel);
                    }
                }
            }
            progress.setProgress((pos.z+1.0f) / newDims.z);
        }
    }
    else if(input.referenceMethod == "zero") {
        // Do nothing, the volume was already initialized all zero.
    }
    else if(input.referenceMethod == "mean") {

        for (size_t r = 0; r < numRuns; r++) {
            size_t t = ensemble->pickTimeStep(r, input.time);
            const VolumeBase* vol = ensemble->getRuns()[r].timeSteps_[t].fieldNames_.at(field);
            VolumeRAMRepresentationLock lock(vol);
            tgt::mat4 physicalToVoxel = vol->getPhysicalToVoxelMatrix();

            tgt::ivec3 pos = tgt::ivec3::zero;
            for (pos.z = 0; pos.z < newDims.z; ++pos.z) {
                for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
                    for (pos.x = 0; pos.x < newDims.x; ++pos.x) {

                        // Map sample position to physical space.
                        tgt::vec3 sample = mapRange(tgt::vec3(pos), tgt::vec3::zero, tgt::vec3(newDims), bounds.getLLF(), bounds.getURB());

                        // Map to voxel space.
                        tgt::svec3 sampleInVoxelSpace = physicalToVoxel * sample;

                        // Ignore, if out of source bounds.
                        if(tgt::clamp(sampleInVoxelSpace, tgt::svec3::zero, lock->getDimensions() - tgt::svec3::one) != sampleInVoxelSpace) {
                            continue;
                        }

                        for (size_t channel = 0; channel < numChannels; channel++) {
                            float oldValue = output->getVoxelNormalized(pos, channel);
                            float newValue = lock->getVoxelNormalized(sampleInVoxelSpace, channel) / numRuns;
                            output->setVoxelNormalized(oldValue + newValue, pos, channel);
                        }
                    }
                }
            }
            progress.setProgress((r+1.0f) / numRuns);
        }
    }
    else if(input.referenceMethod == "median") {
        LERROR("Not yet implemented");
        // TODO: implement efficiently
        /*
        tgt::ivec3 pos = tgt::ivec3::zero;
        for (pos.z = 0; pos.z < newDims.z; ++pos.z) {
            for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
                for (pos.x = 0; pos.x < newDims.x; ++pos.x) {

                    // Map sample position to physical space.
                    tgt::vec3 sample = mapRange(tgt::vec3(pos), tgt::vec3::zero, tgt::vec3(newDims), roi.getLLF(),
                                                roi.getURB());

                    std::vector<tgt::vec3> samples;

                    Statistics samples(false);
                    for (size_t r = 0; r < numRuns; r++) {
                        if (r != input.referenceRun) {

                            size_t t = ensemble->pickTimeStep(r, input.time);
                            const VolumeBase* vol = ensemble->getRuns()[r].timeSteps_[t].fieldNames_.at(field);
                            VolumeRAMRepresentationLock lock(vol);
                            tgt::vec3 sampleInVoxelSpace = vol->getPhysicalToVoxelMatrix() * sample;

                            tgt::vec3 voxelDiff = tgt::vec3::zero;
                            for (size_t channel = 0; channel < numChannels; channel++) {
                                voxelDiff[channel] =
                                        lock->getVoxelNormalized(sampleInVoxelSpace, channel) - referenceVoxel[channel];
                            }

                            samples.addSample(tgt::length(voxelDiff));
                        }
                    }

                    output->voxel(pos) = samples.getStdDev();
                }
            }
            progress.setProgress(1.0f * pos.z / newDims.z);
        }
         */
    }

    tgt::vec3 spacing = bounds.diagonal() / tgt::vec3(newDims);
    std::unique_ptr<Volume> volume(new Volume(output.release(), spacing, bounds.getLLF()));
    volume->getMetaDataContainer().addMetaData("time", new FloatMetaData(input.time));
    volume->getMetaDataContainer().addMetaData("field", new StringMetaData(field));

    progress.setProgress(1.0f);

    return ReferenceVolumeCreatorOutput{
            std::move(volume)
    };
}

void ReferenceVolumeCreator::processComputeOutput(ReferenceVolumeCreatorOutput output) {
    outport_.setData(output.volume.release(), true);
}


void ReferenceVolumeCreator::adjustPropertiesToInput() {
    if(!inport_.hasData()) return;

    const EnsembleDataset* ensemble = inport_.getData();

    referenceRun_.reset();
    referenceRun_.setOptions(std::deque<Option<std::string>>());
    for(const EnsembleDataset::Run& run : ensemble->getRuns()) {
        referenceRun_.addOption(run.name_, run.name_);
    }

    selectedField_.reset();
    selectedField_.setOptions(std::deque<Option<std::string>>());
    for(const std::string& field : ensemble->getCommonFieldNames()) {
        selectedField_.addOption(field, field);
    }

    time_.setMinValue(ensemble->getStartTime());
    time_.setMaxValue(ensemble->getEndTime());
    //time_.set(ensemble->getStartTime());
}

Processor* ReferenceVolumeCreator::create() const {
    return new ReferenceVolumeCreator();
}

} // namespace voreen
