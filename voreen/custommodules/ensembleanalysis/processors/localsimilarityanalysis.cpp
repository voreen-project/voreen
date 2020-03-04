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

#include "localsimilarityanalysis.h"

#include "voreen/core/utils/statistics.h"
#include "../utils/utils.h"

namespace voreen {

const std::string LocalSimilarityAnalysis::loggerCat_("voreen.ensembleanalysis.LocalSimilarityAnalysis");

LocalSimilarityAnalysis::LocalSimilarityAnalysis()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , ensembleInport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , referencePort_(Port::INPORT, "referenceport", "Reference Volume Port")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , selectedField_("selectedField", "Selected Field")
    , time_("time", "Time", 0.0f, 0.0f, 1000000.0f)
{
    // Ports
    addPort(ensembleInport_);
    addPort(referencePort_);
    addPort(outport_);

    addProperty(selectedField_);
    addProperty(time_);
    time_.setTracking(false);
}

LocalSimilarityAnalysis::~LocalSimilarityAnalysis() {
}

LocalSimilarityAnalysisInput LocalSimilarityAnalysis::prepareComputeInput() {
    PortDataPointer<EnsembleDataset> ensemble = ensembleInport_.getThreadSafeData();
    if (!ensemble) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    if(ensemble->getRuns().empty()) {
        throw InvalidInputException("Need at least a single run", InvalidInputException::S_ERROR);
    }

    const VolumeBase* reference = referencePort_.getData();
    if(!reference) {
        throw InvalidInputException("No reference volume", InvalidInputException::S_ERROR);
    }

    VolumeRAMRepresentationLock referenceVolume(referencePort_.getData());
    if(ensemble->getNumChannels(selectedField_.get()) != referenceVolume->getNumChannels()) {
        throw InvalidInputException("Reference Volume channel count is different from selected field", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<VolumeRAM_Float> outputVolume(new VolumeRAM_Float(referenceVolume->getDimensions()));
    outputVolume->clear();

    return LocalSimilarityAnalysisInput{
            std::move(ensemble),
            std::move(referenceVolume),
            std::move(outputVolume),
            reference->getPhysicalToVoxelMatrix(),
            selectedField_.get(),
            time_.get()
    };
}

LocalSimilarityAnalysisOutput LocalSimilarityAnalysis::compute(LocalSimilarityAnalysisInput input, ProgressReporter& progress) const {

    auto ensemble = std::move(input.ensemble);
    const std::string& field = input.field;
    const size_t numRuns = ensemble->getRuns().size();
    const size_t numChannels = ensemble->getNumChannels(field);
    std::unique_ptr<VolumeRAM_Float> output = std::move(input.outputVolume);
    const tgt::ivec3 newDims = output->getDimensions();
    const tgt::Bounds& roi = ensemble->getRoi();

    const VolumeRAMRepresentationLock& referenceVolume = input.referenceVolume;
    const tgt::mat4& refPhysicalToVoxel = input.physicalToVoxel;

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
                    tgt::vec3 sample = mapRange(tgt::vec3(pos), tgt::vec3::zero, tgt::vec3(newDims), roi.getLLF(), roi.getURB());

                    // Map to voxel space.
                    tgt::ivec3 sampleInVoxelSpace = physicalToVoxel * sample;

                    // Ignore, if out of bounds.
                    if (tgt::clamp(sampleInVoxelSpace, tgt::ivec3::zero, newDims - tgt::ivec3::one) != sampleInVoxelSpace) {
                        continue;
                    }

                    float length = 0.0f;
                    for(size_t channel=0; channel<numChannels; channel++) {
                        float value = lock->getVoxelNormalized(sampleInVoxelSpace, channel) - referenceVolume->getVoxelNormalized(pos, channel);
                        length += value * value;
                    }

                    output->voxel(pos) += length / (numRuns - 1.0f);
                }
            }
        }

        // Update progress.
        progress.setProgress(1.0f * r / numRuns);
    }

    // Convert variance to standard deviation.
    for (size_t i = 0; i < output->getNumVoxels(); i++) {
        float variance = output->voxel(i);
        output->voxel(i) = std::sqrt(variance);
    }

    tgt::vec3 spacing = roi.diagonal() / tgt::vec3(newDims);
    std::unique_ptr<Volume> volume(new Volume(output.release(), spacing, roi.getLLF()));
    volume->getMetaDataContainer().addMetaData("time", new FloatMetaData(input.time));
    volume->getMetaDataContainer().addMetaData("field", new StringMetaData(field));

    progress.setProgress(1.0f);

    return LocalSimilarityAnalysisOutput{
        std::move(volume)
    };
}

void LocalSimilarityAnalysis::processComputeOutput(LocalSimilarityAnalysisOutput output) {
    outport_.setData(output.volume.release(), true);
}


void LocalSimilarityAnalysis::adjustPropertiesToInput() {
    if(!ensembleInport_.hasData()) return;

    const EnsembleDataset* ensemble = ensembleInport_.getData();

    selectedField_.reset();
    selectedField_.setOptions(std::deque<Option<std::string>>());
    for(const std::string& field : ensemble->getCommonFieldNames()) {
        selectedField_.addOption(field, field);
    }

    time_.setMinValue(ensemble->getStartTime());
    time_.setMaxValue(ensemble->getEndTime());
    //time_.set(ensemble->getEndTime());
}

Processor* LocalSimilarityAnalysis::create() const {
    return new LocalSimilarityAnalysis();
}

} // namespace voreen
