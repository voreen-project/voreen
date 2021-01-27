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

#include "ensemblemeancreator.h"

#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/utils/statistics.h"
#include "../utils/utils.h"

namespace voreen {

const std::string EnsembleMeanCreator::loggerCat_("voreen.ensembleanalysis.EnsembleMeanCreator");

EnsembleMeanCreator::EnsembleMeanCreator()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , selectedField_("selectedField", "Selected Field")
    , time_("time", "Time", 0.0f, 0.0f, 1000000.0f)
    , sampleRegion_("sampleRegion", "Sample Region")
    , outputDimensions_("outputDimensions", "Output Dimensions", tgt::ivec3(200), tgt::ivec3(2), tgt::ivec3(1000))
{
    // Ports
    addPort(inport_);
    ON_CHANGE(inport_, EnsembleMeanCreator, adjustToEnsemble);
    addPort(outport_);

    addProperty(selectedField_);
    addProperty(time_);
    time_.setTracking(false);
    addProperty(sampleRegion_);
    sampleRegion_.addOption("bounds", "Ensemble Bounds");
    sampleRegion_.addOption("common", "Common Bounds");
    addProperty(outputDimensions_);
}

EnsembleMeanCreator::~EnsembleMeanCreator() {
}

EnsembleMeanCreatorInput EnsembleMeanCreator::prepareComputeInput() {
    PortDataPointer<EnsembleDataset> ensemble = inport_.getThreadSafeData();
    if (!ensemble) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    // Get required information about mean volume format.
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
    else {
        throw InvalidInputException("Unknown sample region", InvalidInputException::S_ERROR);
    }

    // Clear old data.
    outport_.clear();

    return EnsembleMeanCreatorInput{
            std::move(ensemble),
            std::move(outputVolume),
            bounds,
            selectedField_.get(),
            time_.get()
    };
}

EnsembleMeanCreatorOutput EnsembleMeanCreator::compute(EnsembleMeanCreatorInput input, ProgressReporter& progress) const {

    auto ensemble = std::move(input.ensemble);
    const tgt::Bounds& bounds = input.bounds;
    std::string field = std::move(input.field);
    const size_t numMembers = ensemble->getMembers().size();
    const size_t numChannels = ensemble->getNumChannels(field);
    std::unique_ptr<VolumeRAM> output = std::move(input.outputVolume);
    const tgt::ivec3 newDims = output->getDimensions();

    for (size_t r = 0; r < numMembers; r++) {
        size_t t = ensemble->getMembers()[r].getTimeStep(input.time);
        const VolumeBase* vol = ensemble->getMembers()[r].getTimeSteps()[t].getVolume(field);
        VolumeRAMRepresentationLock lock(vol);
        tgt::mat4 worldToVoxel = vol->getWorldToVoxelMatrix();

        tgt::ivec3 pos = tgt::ivec3::zero;
        for (pos.z = 0; pos.z < newDims.z; ++pos.z) {
            for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
                for (pos.x = 0; pos.x < newDims.x; ++pos.x) {

                    // Map sample position to world space.
                    tgt::vec3 sample = mapRange(tgt::vec3(pos), tgt::vec3::zero, tgt::vec3(newDims), bounds.getLLF(), bounds.getURB());

                    // Map to voxel space.
                    tgt::svec3 sampleInVoxelSpace = worldToVoxel * sample;

                    // Ignore, if out of source bounds.
                    if(tgt::clamp(sampleInVoxelSpace, tgt::svec3::zero, lock->getDimensions() - tgt::svec3::one) != sampleInVoxelSpace) {
                        continue;
                    }

                    for (size_t channel = 0; channel < numChannels; channel++) {
                        float oldValue = output->getVoxelNormalized(pos, channel);
                        float newValue = lock->getVoxelNormalized(sampleInVoxelSpace, channel) / numMembers;
                        output->setVoxelNormalized(oldValue + newValue, pos, channel);
                    }
                }
            }
        }
        progress.setProgress((r+1.0f) / numMembers);
    }

    tgt::vec3 spacing = bounds.diagonal() / tgt::vec3(newDims);
    std::unique_ptr<Volume> volume(new Volume(output.release(), spacing, bounds.getLLF()));
    volume->setTimestep(input.time);
    volume->setModality(Modality(field));

    progress.setProgress(1.0f);

    return EnsembleMeanCreatorOutput{
            std::move(volume)
    };
}

void EnsembleMeanCreator::processComputeOutput(EnsembleMeanCreatorOutput output) {
    outport_.setData(output.volume.release(), true);
}


void EnsembleMeanCreator::adjustToEnsemble() {
    if(!inport_.hasData()) return;

    const EnsembleDataset* ensemble = inport_.getData();

    selectedField_.reset();
    selectedField_.setOptions(std::deque<Option<std::string>>());
    for(const std::string& field : ensemble->getCommonFieldNames()) {
        selectedField_.addOption(field, field);
    }

    time_.setMinValue(ensemble->getStartTime());
    time_.setMaxValue(ensemble->getEndTime());
    //time_.set(ensemble->getStartTime());
}

Processor* EnsembleMeanCreator::create() const {
    return new EnsembleMeanCreator();
}

} // namespace voreen
