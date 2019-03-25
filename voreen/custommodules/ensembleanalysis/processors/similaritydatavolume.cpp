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
    , outputDimensions_("outputDimensions", "Output Dimensions", tgt::ivec3(300), tgt::ivec3(10), tgt::ivec3(1000))
    , time_("time", "Time", 0.0f, 0.0f, 1000000.0f)
    , similarityMethod_("similarityMethod", "Similarity Method")
    , selectedChannel_("selectedChannel", "Selected Channel")
    , comparisonMethod_("comparisonMethod", "Comparison Method")
    , groupBehaviour_("groupBehaviour", "Group behaviour")
    , singleRunSelection_("singleRunSelection", "Compare run ...")
    , group1_("similartyDataVolumeGroup1", "Group 1")
    , group2_("similartyDataVolumeGroup2", "Group 2")
{
    // Ports
    addPort(inport_);
    ON_CHANGE(inport_, SimilartyDataVolume, adjustToEnsemble);
    addPort(outport_);

    addProperty(outputDimensions_);
    addProperty(time_);
    addProperty(similarityMethod_);
    similarityMethod_.addOption("variance", "Variance");
    similarityMethod_.addOption("minmax", "Min/Max Comparison");
    similarityMethod_.select("variance");
    addProperty(selectedChannel_);

    addProperty(comparisonMethod_);
    comparisonMethod_.addOption("selection", "Compare selection");
    comparisonMethod_.addOption("oneToGroup", "Compare one to group");
    comparisonMethod_.addOption("groupToGroup", "Compare group to group");
    comparisonMethod_.addOption("oneToItself", "Compare one to itself");
    ON_CHANGE(comparisonMethod_, SimilartyDataVolume, onComparisonMethodChange);
    onComparisonMethodChange();

    addProperty(singleRunSelection_);
    addProperty(group1_);
    addProperty(group2_);

//    ON_CHANGE(singleRunSelection_, SimilartyDataVolume, updateProperties);
//    ON_CHANGE(group1_, SimilartyDataVolume, updateProperties);
//    ON_CHANGE(group2_, SimilartyDataVolume, updateProperties);

    addProperty(groupBehaviour_);
    groupBehaviour_.addOption("groupAvg", "Compare group average");
    groupBehaviour_.addOption("groupMin", "Compare group min value");
    groupBehaviour_.addOption("groupMax", "Compare group max value");

}

SimilartyDataVolume::~SimilartyDataVolume() {
}

void SimilartyDataVolume::onComparisonMethodChange() {
    if(comparisonMethod_.getKey() == "selection") {
        group1_.setVisibleFlag(true);
        group1_.setGuiName("Run selection");
        group2_.setVisibleFlag(false);
        group2_.setGuiName("Selection");
        groupBehaviour_.setVisibleFlag(false);
        singleRunSelection_.setVisibleFlag(false);
    }
    else if(comparisonMethod_.getKey() == "oneToGroup") {
        group1_.setVisibleFlag(false);
        group1_.setGuiName("Compare one run ...");
        group2_.setVisibleFlag(true);
        group2_.setGuiName("... to group of runs");
        groupBehaviour_.setVisibleFlag(true);
        singleRunSelection_.setVisibleFlag(true);
    }
    else if(comparisonMethod_.getKey() == "groupToGroup") {
        group1_.setVisibleFlag(true);
        group1_.setGuiName("Select group of runs");
        group2_.setVisibleFlag(true);
        group2_.setGuiName("Select other group of runs");
        groupBehaviour_.setVisibleFlag(true);
        singleRunSelection_.setVisibleFlag(false);
    }
    else if(comparisonMethod_.getKey() == "oneToItself") {
        group1_.setVisibleFlag(false);
        group1_.setGuiName("Compare one run ...");
        group2_.setVisibleFlag(false);
        group2_.setGuiName("... to group of runs");
        groupBehaviour_.setVisibleFlag(false);
        singleRunSelection_.setVisibleFlag(true);
    }
}

SimilarityDataVolumeCreatorInput SimilartyDataVolume::prepareComputeInput() {
    const EnsembleDataset* inputPtr = inport_.getData();
    if (!inputPtr)
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);

    const EnsembleDataset& input = *inputPtr;

    // TODO: check if selected time is simulated in each selected run.

    isReadyToCompute();

    tgt::ivec3 newDims = outputDimensions_.get();
    std::unique_ptr<VolumeRAM_Float> volumeData(new VolumeRAM_Float(newDims, true));
    memset(volumeData->getData(), 0, volumeData->getNumBytes());

    std::string runGroup1;
    for(int index : group1_.getSelectedRowIndices())
        runGroup1 += input.getRuns()[index].name_ + " ";

    std::string runGroup2;
    for(int index : group2_.getSelectedRowIndices())
        runGroup2 += input.getRuns()[index].name_ + " ";

    return SimilarityDataVolumeCreatorInput{
            input,
            std::move(volumeData),
            time_.get(),
            runGroup1,
            runGroup2,
            selectedChannel_.get()
    };
}

SimilarityDataVolumeCreatorOutput SimilartyDataVolume::compute(SimilarityDataVolumeCreatorInput input, ProgressReporter& progress) const {

    progress.setProgress(0.0f);

    const std::string& channel = input.channel;
    const tgt::Bounds& roi = input.dataset.getRoi();

    tgt::ivec3 newDims = input.volumeData->getDimensions();
    tgt::vec3 spacing = roi.diagonal() / tgt::vec3(newDims);

    std::unique_ptr<Volume> volume(new Volume(input.volumeData.release(), spacing, roi.getLLF()));

    float progressIncrement = 0.95f / newDims.z;
    tgt::ivec3 pos = tgt::ivec3::zero; // iteration variable
    for (pos.z = 0; pos.z < newDims.z; ++pos.z) {
        for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
            for (pos.x = 0; pos.x < newDims.x; ++pos.x) {

                tgt::vec3 sample = volume->getVoxelToPhysicalMatrix() * tgt::vec3(pos);

                std::vector<float> samples(input.dataset.getRuns().size());
                for(size_t r = 0; r<input.dataset.getRuns().size(); r++) {

                    // Filter unused runs to save time.
//                    if(std::find(group1_.getSelectedRowIndices().begin(), group1_.getSelectedRowIndices().end(), r)
//                       == group1_.getSelectedRowIndices().end() && group2_.getSelectedRowIndices().end() ==
//                       std::find(group2_.getSelectedRowIndices().begin(), group2_.getSelectedRowIndices().end(), r)
//                       )
//                        continue;
                    const EnsembleDataset::Run& run = input.dataset.getRuns()[r];

                    if(comparisonMethod_.get() == "oneToItself") {
                        if(r != singleRunSelection_.getSelectedIndex()) continue;

                        samples.resize(2);
                        size_t t = input.dataset.pickTimeStep(r, input.time);
                        const VolumeBase* volume = run.timeSteps_[t].channels_.at(channel);
                        const VolumeRAM_Float* volumeData = dynamic_cast<const VolumeRAM_Float*>(volume->getRepresentation<VolumeRAM>());

                        const VolumeBase* volumeStart = run.timeSteps_[0].channels_.at(channel);
                        const VolumeRAM_Float* volumeDataStart = dynamic_cast<const VolumeRAM_Float*>(volumeStart->getRepresentation<VolumeRAM>());



                        samples[r] = volumeData->getVoxelNormalizedLinear(volume->getPhysicalToVoxelMatrix() * sample);
                        samples[r+1] = volumeDataStart->getVoxelNormalizedLinear(volumeStart->getPhysicalToVoxelMatrix() * sample);
                        break;
                    }
                    else {
                        size_t t = input.dataset.pickTimeStep(r, input.time);
                        const VolumeBase* volume = run.timeSteps_[t].channels_.at(channel);
                        const VolumeRAM_Float* volumeData = dynamic_cast<const VolumeRAM_Float*>(volume->getRepresentation<VolumeRAM>());

                        samples[r] = volumeData->getVoxelNormalizedLinear(volume->getPhysicalToVoxelMatrix() * sample);
                    }
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

    volume->getMetaDataContainer().addMetaData("time", new FloatMetaData(input.time));
    volume->getMetaDataContainer().addMetaData("channel", new StringMetaData(channel));
    std::string group1string = comparisonMethod_.getKey() == "oneToGroup" ? singleRunSelection_.get() : input.runGroup1;
    std::string group2string = comparisonMethod_.getKey() == "selection" ? "" : input.runGroup2;
    volume->getMetaDataContainer().addMetaData("group1", new StringMetaData(group1string));
    volume->getMetaDataContainer().addMetaData("group2", new StringMetaData(group2string));
    return SimilarityDataVolumeCreatorOutput{std::move(volume)};
}

void SimilartyDataVolume::processComputeOutput(SimilarityDataVolumeCreatorOutput output) {
    outport_.setData(output.volume.release(), true);
}

void SimilartyDataVolume::updateProperties() {
    tgt::vec2 timeSpan(0, std::numeric_limits<float>::max());

    std::vector<int> runIndices;
    if(comparisonMethod_.getKey() == "selection") {
        runIndices = group1_.getSelectedRowIndices();
    }
    else if(comparisonMethod_.getKey() == "oneToGroup") {
        runIndices = group2_.getSelectedRowIndices();
        runIndices.push_back(singleRunSelection_.getSelectedIndex());
    }
    else if(comparisonMethod_.getKey() == "groupToGroup") {
        runIndices = group1_.getSelectedRowIndices();
        for(int index : group2_.getSelectedRowIndices()) {
            runIndices.push_back(index);
        }
    }
    for(int index : runIndices) {
        EnsembleDataset::Run run = inport_.getData()->getRuns().at(index);
        EnsembleDataset::TimeStep firstTimeStep = run.timeSteps_.front();
        EnsembleDataset::TimeStep lastTimeStep = run.timeSteps_.back();
        timeSpan.x = std::max(timeSpan.x, firstTimeStep.time_);
        timeSpan.y = std::min(timeSpan.y, lastTimeStep.time_);
    }

    time_.setMinValue(timeSpan.x);
    time_.setMaxValue(timeSpan.y);
}

tgt::vec3 SimilartyDataVolume::getSpacing() const {

    tgt::vec3 commonSpacing(std::numeric_limits<float>::max());

    std::vector<int> runIndices;
    if(comparisonMethod_.getKey() == "selection") {
        runIndices = group1_.getSelectedRowIndices();
    }
    else if(comparisonMethod_.getKey() == "oneToGroup") {
        runIndices = group2_.getSelectedRowIndices();
        runIndices.push_back(singleRunSelection_.getSelectedIndex());
    }
    else if(comparisonMethod_.getKey() == "groupToGroup") {
        runIndices = group1_.getSelectedRowIndices();
        for(int index : group2_.getSelectedRowIndices()) {
            runIndices.push_back(index);
        }
    }
    for(int index : runIndices) {
        EnsembleDataset::Run run = inport_.getData()->getRuns().at(index);
        const VolumeBase* firstVolume = run.timeSteps_.front().channels_.begin()->second;
        tgt::vec3 spacing = firstVolume->getSpacing();
        commonSpacing.x = std::min(commonSpacing.x, spacing.x);
        commonSpacing.y = std::min(commonSpacing.y, spacing.y);
        commonSpacing.z = std::min(commonSpacing.z, spacing.z);
    }

    return commonSpacing;
}

void SimilartyDataVolume::adjustToEnsemble() {
    if(!inport_.hasData()) return;

    const EnsembleDataset* ensemble = inport_.getData();

    group1_.reset();
    group2_.reset();
    singleRunSelection_.reset();
    singleRunSelection_.setOptions(std::deque<Option<std::string>>());
    std::vector<int> selection;
    for(const EnsembleDataset::Run& run : ensemble->getRuns()) {
        group1_.addRow(run.name_);
        selection.push_back(static_cast<int>(selection.size()));
        group2_.addRow(run.name_);
        singleRunSelection_.addOption(run.name_, run.name_);
    }
    group1_.setSelectedRowIndices(selection);
    group2_.setSelectedRowIndices(selection);

    selectedChannel_.reset();
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
    if(comparisonMethod_.getKey() == "selection") {
        for(int selectedIndex : group1_.getSelectedRowIndices()) {
            modifiedVoxelData.push_back(rawVoxelData.at(selectedIndex));
        }
        return modifiedVoxelData;
    }
    else if(comparisonMethod_.getKey() == "oneToGroup") {
        modifiedVoxelData.push_back(rawVoxelData.at(singleRunSelection_.getSelectedIndex()));
        if(groupBehaviour_.getKey() == "groupAvg") {
            Statistics statistics(false);
            for(int selectedIndex : group2_.getSelectedRowIndices()) {
                statistics.addSample(rawVoxelData.at(selectedIndex));
            }
            modifiedVoxelData.push_back(statistics.getMean());
        }
        else if(groupBehaviour_.getKey() == "groupMin") {
            float minValue = std::numeric_limits<float>::max();
            for(int selectedIndex : group2_.getSelectedRowIndices()) {
                minValue = std::min(minValue, rawVoxelData.at(selectedIndex));
            }
            modifiedVoxelData.push_back(minValue);
        }
        else if(groupBehaviour_.getKey() == "groupMax") {
            float maxValue = std::numeric_limits<float>::lowest();
            for(int selectedIndex : group2_.getSelectedRowIndices()) {
                maxValue = std::max(maxValue, rawVoxelData.at(selectedIndex));
            }
            modifiedVoxelData.push_back(maxValue);
        }
        return modifiedVoxelData;
    }
    else if(comparisonMethod_.getKey() == "groupToGroup") {
        if(groupBehaviour_.getKey() == "groupAvg") {
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
        else if(groupBehaviour_.getKey() == "groupMin") {
            float minValue = std::numeric_limits<float>::max();
            for(int selectedIndex : group1_.getSelectedRowIndices()) {
                minValue = std::min(minValue, rawVoxelData.at(selectedIndex));
            }
            modifiedVoxelData.push_back(minValue);

            minValue = std::numeric_limits<float>::max();
            for(int selectedIndex : group2_.getSelectedRowIndices()) {
                minValue = std::min(minValue, rawVoxelData.at(selectedIndex));
            }
            modifiedVoxelData.push_back(minValue);
        }
        else if(groupBehaviour_.getKey() == "groupMax") {
            float maxValue = std::numeric_limits<float>::lowest();
            for(int selectedIndex : group1_.getSelectedRowIndices()) {
                maxValue = std::max(maxValue, rawVoxelData.at(selectedIndex));
            }
            modifiedVoxelData.push_back(maxValue);

            maxValue = std::numeric_limits<float>::lowest();
            for(int selectedIndex : group2_.getSelectedRowIndices()) {
                maxValue = std::max(maxValue, rawVoxelData.at(selectedIndex));
            }
            modifiedVoxelData.push_back(maxValue);
        }
        return modifiedVoxelData;
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

bool SimilartyDataVolume::isReadyToCompute() const {
    if(comparisonMethod_.getKey() == "selection") {
        if(group1_.getSelectedRowIndices().size() == 0) {
            throw InvalidInputException("No selection", InvalidInputException::S_WARNING);
        }
        if(group1_.getSelectedRowIndices().size() == 1) {
            throw InvalidInputException("Select at least two runs", InvalidInputException::S_WARNING);
        }
    }
    else if(comparisonMethod_.getKey() == "oneToGroup") {
        if(singleRunSelection_.getSelectedIndex() < 0) {
            throw InvalidInputException("No run selected", InvalidInputException::S_WARNING);
        }
        if(group2_.getSelectedRowIndices().size() == 0) {
            throw InvalidInputException("No group selected", InvalidInputException::S_WARNING);
        }
    }
    else if(comparisonMethod_.getKey() == "groupToGroup") {
        if(group1_.getSelectedRowIndices().size() == 0) {
            throw InvalidInputException("Group 1 is empty", InvalidInputException::S_WARNING);
        }
        if(group2_.getSelectedRowIndices().size() == 0) {
            throw InvalidInputException("Group 2 is empty", InvalidInputException::S_WARNING);
        }
    }

    return true;
}

} // namespace voreen
