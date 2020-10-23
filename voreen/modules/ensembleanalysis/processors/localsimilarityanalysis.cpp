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
    , time_("time", "Time", 0.0f, 0.0f, std::numeric_limits<float>::max())
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

    if(ensemble->getMembers().empty()) {
        throw InvalidInputException("Need at least a single run", InvalidInputException::S_ERROR);
    }

    const VolumeBase* referenceVolume = referencePort_.getData();
    if(!referenceVolume) {
        throw InvalidInputException("No reference volume", InvalidInputException::S_ERROR);
    }

    if(ensemble->getNumChannels(selectedField_.get()) != referenceVolume->getNumChannels()) {
        throw InvalidInputException("Reference Volume channel count is different from selected field", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<VolumeRAM_Float> outputVolume(new VolumeRAM_Float(referenceVolume->getDimensions()));
    outputVolume->clear();

    return LocalSimilarityAnalysisInput{
            std::move(ensemble),
            referenceVolume,
            std::move(outputVolume),
            selectedField_.get(),
            time_.get()
    };
}

LocalSimilarityAnalysisOutput LocalSimilarityAnalysis::compute(LocalSimilarityAnalysisInput input, ProgressReporter& progress) const {

    auto ensemble = std::move(input.ensemble);
    const std::string& field = input.field;
    const size_t numMembers = ensemble->getMembers().size();
    const size_t numChannels = ensemble->getNumChannels(field);
    std::unique_ptr<VolumeRAM_Float> output = std::move(input.outputVolume);
    const tgt::svec3 dims = output->getDimensions();

    VolumeRAMRepresentationLock referenceVolume(input.referenceVolume);
    tgt::mat4 refVoxelToWorld = input.referenceVolume->getVoxelToWorldMatrix();

    for (size_t r = 0; r < numMembers; r++) {
        size_t t = ensemble->getMembers()[r].getTimeStep(input.time);
        const VolumeBase* vol = ensemble->getMembers()[r].getTimeSteps()[t].getVolume(field);
        VolumeRAMRepresentationLock lock(vol);
        tgt::Bounds bounds = vol->getBoundingBox().getBoundingBox();
        tgt::mat4 worldToVoxel = vol->getWorldToVoxelMatrix();

        tgt::svec3 pos = tgt::svec3::zero;
        for (pos.z = 0; pos.z < dims.z; ++pos.z) {
            for (pos.y = 0; pos.y < dims.y; ++pos.y) {
                for (pos.x = 0; pos.x < dims.x; ++pos.x) {

                    // Transform sample into world space.
                    tgt::vec3 sample = refVoxelToWorld * tgt::vec3(pos);

                    // Ignore, if out of bounds.
                    if(!bounds.containsPoint(sample)) {
                        continue;
                    }

                    // Transform to local voxel space.
                    sample = worldToVoxel * sample;

                    float length = 0.0f;
                    for(size_t channel=0; channel<numChannels; channel++) {
                        float value = lock->getVoxelNormalized(sample, channel) - referenceVolume->getVoxelNormalized(pos, channel);
                        length += value * value;
                    }

                    output->voxel(pos) += length / numMembers;
                }
            }
        }

        // Update progress.
        progress.setProgress(1.0f * r / numMembers);
    }

    // Convert variance to standard deviation.
    for (size_t i = 0; i < output->getNumVoxels(); i++) {
        float variance = output->voxel(i);
        output->voxel(i) = std::sqrt(variance);
    }

    std::unique_ptr<Volume> volume(new Volume(output.release(), input.referenceVolume->getSpacing(), input.referenceVolume->getOffset()));
    volume->setMetaDataValue<FloatMetaData>("time", input.time);
    volume->setMetaDataValue<StringMetaData>("field", field);

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
