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

#include "ensemblevarianceanalysis.h"

#include "voreen/core/utils/statistics.h"
#include "../utils/utils.h"

namespace voreen {

const std::string EnsembleVarianceAnalysis::loggerCat_("voreen.ensembleanalysis.LocalSimilarityAnalysis");

EnsembleVarianceAnalysis::EnsembleVarianceAnalysis()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , ensembleInport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , ensembleMeanPort_(Port::INPORT, "meanport", "Ensemble Mean Volume Port")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , selectedField_("selectedField", "Selected Field")
    , vectorMagnitudeThreshold_("vectorMagnitudeThreshold", "Vector Magnitude Threshold", 0.0f, 0.0f, 10000.0f)
    , vectorComponent_("vectorStrategy", "Strategy")
    , time_("time", "Time", 0.0f, 0.0f, std::numeric_limits<float>::max())
{
    // Ports
    addPort(ensembleInport_);
    ON_CHANGE(ensembleInport_, EnsembleVarianceAnalysis, adjustToEnsemble);
    addPort(ensembleMeanPort_);
    addPort(outport_);

    addProperty(selectedField_);
    ON_CHANGE_LAMBDA(selectedField_, [this] {
        const EnsembleDataset* ensemble = ensembleInport_.getData();
        if(!ensemble) {
            return;
        }
        bool scalarField = ensemble->getNumChannels(selectedField_.get()) == 1;
        vectorMagnitudeThreshold_.setReadOnlyFlag(scalarField);
        vectorComponent_.setReadOnlyFlag(scalarField);
        if(scalarField) {
            tgt::vec2 range = ensemble->getMagnitudeRange(selectedField_.get());
            vectorMagnitudeThreshold_.setMinValue(range.x);
            vectorMagnitudeThreshold_.setMaxValue(range.y);
            vectorMagnitudeThreshold_.set(range.x);
        }
    });
    addProperty(vectorMagnitudeThreshold_);
    addProperty(vectorComponent_);
    vectorComponent_.addOption("both", "Both", BOTH);
    vectorComponent_.addOption("magnitude", "Magnitude", MAGNITUDE);
    vectorComponent_.addOption("direction", "Direction", DIRECTION);
    addProperty(time_);
    time_.setTracking(false);
}

EnsembleVarianceAnalysis::~EnsembleVarianceAnalysis() {
}

EnsembleVarianceAnalysisInput EnsembleVarianceAnalysis::prepareComputeInput() {
    PortDataPointer<EnsembleDataset> ensemble = ensembleInport_.getThreadSafeData();
    if (!ensemble) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    if(ensemble->getMembers().empty()) {
        throw InvalidInputException("Empty ensemble", InvalidInputException::S_ERROR);
    }

    const VolumeBase* meanVolume = ensembleMeanPort_.getData();
    if(!meanVolume) {
        throw InvalidInputException("No mean volume", InvalidInputException::S_ERROR);
    }

    if(ensemble->getNumChannels(selectedField_.get()) != meanVolume->getNumChannels()) {
        throw InvalidInputException("Mean Volume channel count is different from selected field", InvalidInputException::S_ERROR);
    }

    if(ensemble->getNumChannels(selectedField_.get()) > 4) {
        throw InvalidInputException("Only up to 4 channels supported", InvalidInputException::S_ERROR);
    }

    // We output float so there is no need to find a proper real world mapping.
    std::unique_ptr<VolumeRAM_Float> outputVolume(new VolumeRAM_Float(meanVolume->getDimensions()));
    outputVolume->clear();

    return EnsembleVarianceAnalysisInput{
            std::move(ensemble),
            meanVolume,
            std::move(outputVolume),
            selectedField_.get(),
            vectorMagnitudeThreshold_.get(),
            vectorComponent_.getValue(),
            time_.get()
    };
}

EnsembleVarianceAnalysisOutput EnsembleVarianceAnalysis::compute(EnsembleVarianceAnalysisInput input, ProgressReporter& progress) const {

    auto ensemble = std::move(input.ensemble);
    std::string field = std::move(input.field);
    const size_t numMembers = ensemble->getMembers().size();
    const size_t numChannels = ensemble->getNumChannels(field);
    std::unique_ptr<VolumeRAM_Float> output = std::move(input.outputVolume);
    const tgt::svec3 dims = output->getDimensions();

    VolumeRAMRepresentationLock meanVolume(input.meanVolume);
    RealWorldMapping rwmMean = input.meanVolume->getRealWorldMapping();
    tgt::mat4 meanVoxelToWorld = input.meanVolume->getVoxelToWorldMatrix();

    for (size_t r = 0; r < numMembers; r++) {
        size_t t = ensemble->getMembers()[r].getTimeStep(input.time);
        const VolumeBase* vol = ensemble->getMembers()[r].getTimeSteps()[t].getVolume(field);
        if(!vol) {
            continue;
        }

        RealWorldMapping rwmCurr = vol->getRealWorldMapping();
        VolumeRAMRepresentationLock lock(vol);
        tgt::Bounds bounds = vol->getBoundingBox().getBoundingBox();
        tgt::mat4 worldToVoxel = vol->getWorldToVoxelMatrix();

        tgt::svec3 pos = tgt::svec3::zero;
        for (pos.z = 0; pos.z < dims.z; ++pos.z) {
            for (pos.y = 0; pos.y < dims.y; ++pos.y) {
                for (pos.x = 0; pos.x < dims.x; ++pos.x) {

                    // Transform sample into world space.
                    tgt::vec3 sample = meanVoxelToWorld * tgt::vec3(pos);

                    // Ignore, if out of bounds.
                    if(!bounds.containsPoint(sample)) {
                        continue;
                    }

                    // Transform to local voxel space.
                    sample = worldToVoxel * sample;

                    if(numChannels == 1 || input.vectorComponent == BOTH) {
                        float length = 0.0f;
                        for (size_t channel = 0; channel < numChannels; channel++) {
                            float value = rwmCurr.normalizedToRealWorld(lock->getVoxelNormalized(sample, channel))
                                        - rwmMean.normalizedToRealWorld(meanVolume->getVoxelNormalized(pos, channel));
                            length += value * value;
                        }

                        output->voxel(pos) += length / numMembers;
                    }
                    else if(input.vectorComponent == MAGNITUDE) {
                        float lengthSqCurrent = 0.0f;
                        float lengthSqMean = 0.0f;
                        for (size_t channel = 0; channel < numChannels; channel++) {
                            float value = rwmCurr.normalizedToRealWorld(lock->getVoxelNormalized(sample, channel));
                            lengthSqCurrent += value * value;

                            value = rwmMean.normalizedToRealWorld(meanVolume->getVoxelNormalized(sample, channel));
                            lengthSqMean += value * value;
                        }
                        float magnitude = std::abs(std::sqrt(lengthSqCurrent) - std::sqrt(lengthSqMean));
                        output->voxel(pos) += magnitude;
                    }
                    else if(input.vectorComponent == DIRECTION) {
                        tgt::vec4 v1 = tgt::vec4::zero;
                        tgt::vec4 v2 = tgt::vec4::zero;
                        for (size_t channel = 0; channel < numChannels; channel++) {
                            v1[channel] = rwmCurr.normalizedToRealWorld(lock->getVoxelNormalized(sample, channel));
                            v2[channel] = rwmMean.normalizedToRealWorld(meanVolume->getVoxelNormalized(sample, channel));
                        }

                        // Test if any magnitude of both vectors is too small to be considered for direction calculation.
                        if(tgt::lengthSq(v1) >= input.vectorMagnitudeThreshold * input.vectorMagnitudeThreshold &&
                            tgt::lengthSq(v2) >= input.vectorMagnitudeThreshold * input.vectorMagnitudeThreshold) {

                            v1 = tgt::normalize(v1);
                            v2 = tgt::normalize(v2);

                            float dot = tgt::dot(v1, v2);
                            float angle = std::acos(tgt::clamp(dot, -1.0f, 1.0f)) / tgt::PIf;
                            output->voxel(pos) += angle / numMembers;
                        }
                    }
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

    std::unique_ptr<Volume> volume(new Volume(output.release(), input.meanVolume->getSpacing(), input.meanVolume->getOffset()));
    volume->setTimestep(input.time);
    volume->setModality(Modality(field));

    progress.setProgress(1.0f);

    return EnsembleVarianceAnalysisOutput{
        std::move(volume)
    };
}

void EnsembleVarianceAnalysis::processComputeOutput(EnsembleVarianceAnalysisOutput output) {
    outport_.setData(output.volume.release(), true);
}


void EnsembleVarianceAnalysis::adjustToEnsemble() {
    if(!ensembleInport_.hasData()) return;

    const EnsembleDataset* ensemble = ensembleInport_.getData();

    selectedField_.setOptions(std::deque<Option<std::string>>());
    for(const std::string& field : ensemble->getUniqueFieldNames()) {
        selectedField_.addOption(field, field);
    }
    selectedField_.selectByIndex(0); // Revert to first field.

    time_.setMinValue(ensemble->getStartTime());
    time_.setMaxValue(ensemble->getEndTime());
    //time_.set(ensemble->getEndTime());
}

Processor* EnsembleVarianceAnalysis::create() const {
    return new EnsembleVarianceAnalysis();
}

} // namespace voreen
