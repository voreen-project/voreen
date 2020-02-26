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
    , inport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , selectedField_("selectedField", "Selected Field")
    , referenceRun_("referenceRun", "Reference Run")
    , time_("time", "Time", 0.0f, 0.0f, 1000000.0f)
    , outputDimensions_("outputDimensions", "Output Dimensions", tgt::ivec3(200), tgt::ivec3(10), tgt::ivec3(1000))
{
    // Ports
    addPort(inport_);
    addPort(outport_);

    addProperty(selectedField_);
    addProperty(referenceRun_);
    addProperty(time_);
    addProperty(outputDimensions_);
}

LocalSimilarityAnalysis::~LocalSimilarityAnalysis() {
}

LocalSimilarityAnalysisInput LocalSimilarityAnalysis::prepareComputeInput() {
    PortDataPointer<EnsembleDataset> ensemble = inport_.getThreadSafeData();
    if (!ensemble) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    if(ensemble->getRuns().size() < 2) {
        throw InvalidInputException("Need at least 2 runs", InvalidInputException::S_ERROR);
    }

    tgt::ivec3 newDims = outputDimensions_.get();
    std::unique_ptr<VolumeRAM_Float> volumeData(new VolumeRAM_Float(newDims));
    volumeData->clear();

    size_t referenceRun = -1;
    for(size_t i = 0; i < ensemble->getRuns().size(); i++) {
        if(ensemble->getRuns()[i].name_ == referenceRun_.get()) {
            referenceRun = i;
            break;
        }
    }

    tgtAssert(referenceRun != -1, "Selected run not available");

    return LocalSimilarityAnalysisInput{
            std::move(ensemble),
            std::move(volumeData),
            time_.get(),
            referenceRun,
            selectedField_.get()
    };
}

LocalSimilarityAnalysisOutput LocalSimilarityAnalysis::compute(LocalSimilarityAnalysisInput input, ProgressReporter& progress) const {

    auto ensemble = std::move(input.ensemble);
    const std::string& field = input.field;
    const size_t numRuns = ensemble->getRuns().size();
    const size_t numChannels = ensemble->getNumChannels(field);
    std::unique_ptr<VolumeRAM_Float> output = std::move(input.volumeData);
    const tgt::ivec3 newDims = output->getDimensions();

    const tgt::Bounds& roi = ensemble->getRoi();
    size_t referenceTimeStep = ensemble->pickTimeStep(input.referenceRun, input.time);

    const VolumeBase* refVol = ensemble->getRuns()[input.referenceRun].timeSteps_[referenceTimeStep].fieldNames_.at(input.field);
    VolumeRAMRepresentationLock referenceVolume(refVol);
    tgt::mat4 refPhysicalToVoxel = refVol->getPhysicalToVoxelMatrix();
    //RealWorldMapping refRwm = refVol->getRealWorldMapping(); // TODO: apply rwm?

    tgt::ivec3 pos = tgt::ivec3::zero;
    for (pos.z = 0; pos.z < newDims.z; ++pos.z) {
        for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
            for (pos.x = 0; pos.x < newDims.x; ++pos.x) {

                // Map sample position to physical space.
                tgt::vec3 sample = mapRange(tgt::vec3(pos), tgt::vec3::zero, tgt::vec3(newDims), roi.getLLF(), roi.getURB());

                // We first sample our reference volume.
                tgt::vec3 sampleInRefVoxelSpace = refPhysicalToVoxel * sample;
                tgt::vec3 referenceVoxel = tgt::vec3::zero;
                for(size_t channel=0; channel<numChannels; channel++) {
                    referenceVoxel[channel] = referenceVolume->getVoxelNormalized(sampleInRefVoxelSpace, channel);
                }

                // Now we calculate the difference
                Statistics samples(false);
                for(size_t r = 0; r<numRuns; r++) {
                    if(r != input.referenceRun) {

                        size_t t = ensemble->pickTimeStep(r, input.time);
                        const VolumeBase* vol = ensemble->getRuns()[r].timeSteps_[t].fieldNames_.at(field);
                        VolumeRAMRepresentationLock lock(vol);
                        tgt::vec3 sampleInVoxelSpace = vol->getPhysicalToVoxelMatrix() * sample;

                        tgt::vec3 voxelDiff = tgt::vec3::zero;
                        for(size_t channel=0; channel<numChannels; channel++) {
                            voxelDiff[channel] = lock->getVoxelNormalized(sampleInVoxelSpace, channel) - referenceVoxel[channel];
                        }

                        samples.addSample(tgt::length(voxelDiff));
                    }
                }

                output->voxel(pos) = samples.getStdDev();
            }
        }

        // Update progress.
        progress.setProgress(1.0f * pos.z / newDims.z);
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
    time_.set(ensemble->getStartTime());
}

Processor* LocalSimilarityAnalysis::create() const {
    return new LocalSimilarityAnalysis();
}

float LocalSimilarityAnalysis::calculateVariance(const std::vector<float>& voxelData) const {
    Statistics statistics(false);

    for(float voxelRunValue : voxelData) {
        statistics.addSample(voxelRunValue);
    }

    return statistics.getVariance();
}

} // namespace voreen
