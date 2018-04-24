/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "similaritydatavolume.h"

#include "voreen/core/utils/statistics.h"

namespace voreen {
    
const std::string SimilartyDataVolume::loggerCat_("voreen.viscontest2018.SimilartyDataVolume");

SimilartyDataVolume::SimilartyDataVolume()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , resampleFactor_("resampleFactor", "Resample Factor", 1.0f, 0.1f, 1.0f)
    , time_("time", "Time", 0.0f, 0.0f, 1000000.0f)
    , similarityMethod_("similarityMethod", "Similarity Method")
    , selectedChannel_("selectedChannel", "Selected Channel")
    , group1_("similartyDataVolumeGroup1", "Group 1")
    , group2_("similartyDataVolumeGroup2", "Group 2")
{
    // Ports
    addPort(inport_);
    ON_CHANGE(inport_, SimilartyDataVolume, adjustToEnsemble);
    addPort(outport_);

    addProperty(resampleFactor_);
    addProperty(time_);
    addProperty(similarityMethod_);
    similarityMethod_.addOption("variance", "Variance");
    similarityMethod_.addOption("minmax", "Min/Max Comparison");
    similarityMethod_.select("variance");
    addProperty(selectedChannel_);

    addProperty(group1_);
    addProperty(group2_);
}

SimilartyDataVolume::~SimilartyDataVolume() {
}

SimilarityDataVolumeCreatorInput SimilartyDataVolume::prepareComputeInput() {
    const EnsembleDataset* inputPtr = inport_.getData();
    if (!inputPtr)
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);

    const EnsembleDataset& input = *inputPtr;

    // TODO: check if selected time is simulated in each selected run.

    tgt::ivec3 newDims = tgt::vec3(input.getDimensions()) * resampleFactor_.get();
    VolumeRAM_Float* volumeData = new VolumeRAM_Float(newDims, true);

    std::string runGroup1;
    for(int index : group1_.getSelectedRowIndices())
        runGroup1 += input.getRuns()[index].name_ + " ";

    std::string runGroup2;
    for(int index : group2_.getSelectedRowIndices())
        runGroup2 += input.getRuns()[index].name_ + " ";

    return SimilarityDataVolumeCreatorInput{
            input,
            volumeData,
            resampleFactor_.get(),
            time_.get(),
            runGroup1,
            runGroup2,
            selectedChannel_.get()
    };
}

SimilarityDataVolumeCreatorOutput SimilartyDataVolume::compute(SimilarityDataVolumeCreatorInput input, ProgressReporter& progress) const {

    progress.setProgress(0.0f);

    const std::string& channel = input.channel;

    tgt::vec3 ratio(1.0f / input.resampleFactor);
    tgt::ivec3 newDims = input.volumeData->getDimensions();

    float progressIncrement = 0.95f / newDims.z;

    tgt::vec3 d_a = tgt::vec3(newDims - tgt::ivec3::one) / 2.0f;
    tgt::vec3 d_b = tgt::vec3(input.dataset.getDimensions() - tgt::svec3::one) / 2.0f;

    tgt::ivec3 pos = tgt::ivec3::zero; // iteration variable
    tgt::vec3 nearest; // stores the new position of the target volume

    for (pos.z = 0; pos.z < newDims.z; ++pos.z) {
        nearest.z = (static_cast<float>(pos.z) - d_a.z) * ratio.z + d_b.z;

        for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
            nearest.y = (static_cast<float>(pos.y) - d_a.y) * ratio.y + d_b.y;

            for (pos.x = 0; pos.x < newDims.x; ++pos.x) {
                nearest.x = (static_cast<float>(pos.x) - d_a.x) * ratio.x + d_b.x;

                std::vector<float> samples(input.dataset.getRuns().size());
                for(size_t r = 0; r<input.dataset.getRuns().size(); r++) {

                    // Filter unused runs to save time.
                    if(std::find(group1_.getSelectedRowIndices().begin(), group1_.getSelectedRowIndices().end(), r)
                       == group1_.getSelectedRowIndices().end() && group2_.getSelectedRowIndices().end() ==
                       std::find(group2_.getSelectedRowIndices().begin(), group2_.getSelectedRowIndices().end(), r)
                       )
                        continue;

                    const EnsembleDataset::Run& run = input.dataset.getRuns()[r];

                    size_t t = input.dataset.pickTimeStep(r, input.time);
                    const VolumeBase* volume = run.timeSteps_[t].channels_.at(channel);
                    const VolumeRAM_Float* volumeData = dynamic_cast<const VolumeRAM_Float*>(volume->getRepresentation<VolumeRAM>());

                    samples[r] = input.dataset.pickSample(volumeData, volume->getSpacing(), nearest);
                }

                // Apply group logic.
                samples = applyGroupLogic(samples);

                // Apply similarity method.
                if(similarityMethod_.get() == "variance") {
                    input.volumeData->voxel(pos) = calculateVariance(samples);
                }
                else if(similarityMethod_.get() == "minmax") {
                    input.volumeData->voxel(pos) = calculateMinMaxDiff(samples);
                }
            }
        }

        // Update progress.
        progress.setProgress(std::min(progress.getProgress() + progressIncrement, 1.0f));
    }

    progress.setProgress(1.0f);

    Volume* volume = new Volume(input.volumeData, input.dataset.getSpacing(), tgt::vec3::zero);
    volume->getMetaDataContainer().addMetaData("time", new FloatMetaData(input.time));
    volume->getMetaDataContainer().addMetaData("channel", new StringMetaData(channel));
    volume->getMetaDataContainer().addMetaData("group1", new StringMetaData(input.runGroup1));
    volume->getMetaDataContainer().addMetaData("group2", new StringMetaData(input.runGroup2));
    return SimilarityDataVolumeCreatorOutput{volume};
}

void SimilartyDataVolume::processComputeOutput(SimilarityDataVolumeCreatorOutput output) {
    outport_.setData(output.volume, true);
}

void SimilartyDataVolume::adjustToEnsemble() {
    if(!inport_.hasData()) return;

    const EnsembleDataset* ensemble = inport_.getData();

    group1_.reset();
    group2_.reset();
    for(const EnsembleDataset::Run& run : ensemble->getRuns()) {
        group1_.addRow(run.name_);
        group2_.addRow(run.name_);
    }

    selectedChannel_.setOptions(std::deque<Option<std::string>>());
    for(const std::string& channel : ensemble->getCommonChannels()) {
        selectedChannel_.addOption(channel, channel);
    }

    time_.setMinValue(ensemble->getStartTime());
    time_.setMaxValue(ensemble->getEndTime());
    time_.set(ensemble->getStartTime());
}

Processor* SimilartyDataVolume::create() const {
    return new SimilartyDataVolume();
}

const std::vector<float> SimilartyDataVolume::applyGroupLogic(const std::vector<float>& rawVoxelData) const {
    std::vector<float> modifiedVoxelData;
    // compare two groups
    if(!group1_.getSelectedRowIndices().empty() && !group2_.getSelectedRowIndices().empty()) {
        Statistics statistics(false);
        for(int selectedIndex : group1_.getSelectedRowIndices()) {
            statistics.addSample(rawVoxelData.at(selectedIndex));
        }
        modifiedVoxelData.push_back(statistics.getMean());

        statistics.reset();
        for(int selectedIndex : group2_.getSelectedRowIndices()) {
            statistics.addSample(rawVoxelData.at(selectedIndex));
        }
        modifiedVoxelData.push_back(statistics.getMean());
    }
    // normal run filtering
    else if(!group1_.getSelectedRowIndices().empty()) {
        for(int selectedIndex : group1_.getSelectedRowIndices()) {
            modifiedVoxelData.push_back(rawVoxelData.at(selectedIndex));
        }
    }
    return rawVoxelData;
}

float SimilartyDataVolume::calculateVariance(const std::vector<float>& voxelData) const {
    Statistics statistics(true);

    for(float voxelRunValue : voxelData) {
        statistics.addSample(voxelRunValue);
    }

    return statistics.getVariance();
}

float SimilartyDataVolume::calculateMinMaxDiff(const std::vector<float>& voxelData) const{
    float min = std::numeric_limits<float>::max();
    float max = std::numeric_limits<float>::lowest();

    for(float voxelRunValue : voxelData) {
        min = std::min(min, voxelRunValue);
        max = std::max(max, voxelRunValue);
    }

    return std::abs(min - max);
}

} // namespace voreen
