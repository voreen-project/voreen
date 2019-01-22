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

#include "ensemblesimilarityplot.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/textureunit.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/glsl.h"

#include "../datastructures/ensembledataset.h"
#include "../utils/utils.h"

namespace voreen {

const std::string EnsembleSimilarityPlot::loggerCat_("voreen.viscontest2018.EnsembleSimilarityPlot");

EnsembleSimilarityPlot::EnsembleSimilarityPlot()
        : RenderProcessor()
        , inport_(Port::INPORT, "ensemble.inport", "Ensemble Input", false)
        , outport_(Port::OUTPORT, "", "Similarity Plot Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
        , sampleRate_(32)
        , similarityRange_(std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest())
        , calculateButton_("applyButton", "Calculate Similarities")
        , isoValueProperty_("isovalue.property", "Iso-Value", 0.0f, 0.0f, 0.0f)
        , channelProperty_("channel", "Selected Channel")
        , progressProperty_("progress", "Progress")
        , selectedTimeStep_("selectedTimeStep", "Selected Time Step", 0, 0, std::numeric_limits<int>::max(), Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
        , autoCalculateProperty_("autocalculatesimilarityplot", "Auto-calculate plot")
        , timeStepPosition_(-1.0f)
        , maxNumTimeSteps_(0)
{
    addPort(inport_);
    addPort(outport_);
    addProperty(calculateButton_);
    addProperty(progressProperty_);
    addProperty(channelProperty_);
    addProperty(isoValueProperty_);
    addProperty(autoCalculateProperty_);
    addProperty(selectedTimeStep_);

    addProgressBar(&progressProperty_);

    ON_CHANGE(calculateButton_, EnsembleSimilarityPlot, createSimilarityPlot);
    ON_CHANGE(inport_, EnsembleSimilarityPlot, onInportChanged);
    ON_CHANGE(channelProperty_, EnsembleSimilarityPlot, initIsoValueProperty);
    ON_CHANGE(isoValueProperty_, EnsembleSimilarityPlot, onIsoValueChanged);
}

EnsembleSimilarityPlot::~EnsembleSimilarityPlot() {
}

Processor* EnsembleSimilarityPlot::create() const {
    return new EnsembleSimilarityPlot();
}

void EnsembleSimilarityPlot::onInportChanged() {
    similarityPlot_.clear();

    // Disable properties
    calculateButton_.setReadOnlyFlag(true);

    if(inport_.hasData() && !inport_.getData()->getRuns().empty()) {
        initProperties();
        calculateButton_.setReadOnlyFlag(false);
    }
}

void EnsembleSimilarityPlot::initProperties() {
    if(!inport_.hasData()) return;

    channelProperty_.setOptions(std::deque<Option<std::string>>());
    for(const std::string& channel : inport_.getData()->getCommonChannels()) {
        channelProperty_.addOption(channel, channel, channel);
    }

    selectedTimeStep_.setMinValue(0);
    selectedTimeStep_.setMaxValue(inport_.getData()->getMinNumTimeSteps());

    initIsoValueProperty();
}

void EnsembleSimilarityPlot::initIsoValueProperty() {
    float minValue = inport_.getData()->getValueRange(channelProperty_.get()).x;
    float maxValue = inport_.getData()->getValueRange(channelProperty_.get()).y;
    isoValueProperty_.setMinValue(minValue);
    isoValueProperty_.setMaxValue(maxValue);
    isoValueProperty_.set((maxValue - minValue) / 2.0f);
}

void EnsembleSimilarityPlot::onIsoValueChanged() {
    if(autoCalculateProperty_.get() == true) createSimilarityPlot();
}

void EnsembleSimilarityPlot::process() {

    if(!inport_.hasData()) return;

    outport_.activateTarget();
    outport_.clearTarget();

    glDepthFunc(GL_ALWAYS);

    drawSimilarityTimeSeries();
    drawTimeStepSelection();

    glLineWidth(1.0f);
    glDepthFunc(GL_LESS);
    IMode.color(tgt::vec4::one);

    outport_.deactivateTarget();
}

void EnsembleSimilarityPlot::drawSimilarityTimeSeries() {

    int numberOfTimesteps = maxNumTimeSteps_;
    float timestepWidth = 2.0f / (numberOfTimesteps - 1.0f);

    glLineWidth(4.0f);
    IMode.begin(tgt::ImmediateMode::LINES);

    const tgt::vec3 colors [] = {{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 1.0f}, {1.0f, 1.0f, 0.0f}};

    size_t i = 0;
    for(std::pair<std::string, std::vector<std::map<std::string, float>>> runSimilarity : similarityPlot_) {

        IMode.color(colors[i]);
        for(int timestepIndex = 1; timestepIndex < numberOfTimesteps; timestepIndex++) {

            if(runSimilarity.second.size() < timestepIndex + 1) continue;

            float minSimilarity = similarityRange_.first;
            float maxSimilarity = similarityRange_.second;

            float currSimilarityValue = (runSimilarity.second[timestepIndex][channelProperty_.get()] - minSimilarity) / (maxSimilarity - minSimilarity);
            float prevSimilarityValue = (runSimilarity.second[timestepIndex - 1][channelProperty_.get()] - minSimilarity) / (maxSimilarity - minSimilarity);


            float x1 = (timestepIndex - 1.0f) * timestepWidth - 1.0f;
            float y1 = (prevSimilarityValue * 2.0f) - 1.0f;
            float x2 = timestepIndex * timestepWidth - 1.0f;
            float y2 = (currSimilarityValue * 2.0f) - 1.0f;

            IMode.vertex(x1, y1);
            IMode.vertex(x2, y2);
        }

        i++;
    }

    IMode.end();
}

void EnsembleSimilarityPlot::drawTimeStepSelection() {
    glDepthFunc(GL_ALWAYS);
    IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
    glLineWidth(5.0f);
    IMode.begin(tgt::ImmediateMode::LINES);
        IMode.vertex(timeStepPosition_, 1.0f);
        IMode.vertex(timeStepPosition_, -1.0f);
    IMode.end();
    glLineWidth(1.0f);
    glDepthFunc(GL_LESS);
}

void EnsembleSimilarityPlot::createSimilarityPlot() {
    if(!inport_.hasData()) return;

    similarityPlot_.clear();
    calculateButton_.setReadOnlyFlag(true);
    setProgress(0.0f);

    similarityRange_.first = std::numeric_limits<float>::max();
    similarityRange_.second = std::numeric_limits<float>::lowest();

    const EnsembleDataset* dataset = inport_.getData();
    TimeStepMapper timeStepMapper(dataset);

    float progressStep = 1.0f / (timeStepMapper.getMaxNumTimeSteps() * inport_.getData()->getRuns().size());
    maxNumTimeSteps_ = timeStepMapper.getMaxNumTimeSteps();

    timeStepMapper.iterate([&](std::string runName, const EnsembleDataset::TimeStep timestep, float time) {
        std::map<std::string, float> channelSimilarity;

        const VolumeBase* channelVolume = static_cast<const VolumeBase*>(timestep.channels_.at(channelProperty_.get()));
        float normalizationFactor = 1.0f; //numberOfTimesteps;
        float similarity = calculateSimilarity(channelVolume, isoValueProperty_.get(), normalizationFactor);
        similarityRange_.first = std::min(similarityRange_.first, similarity);
        similarityRange_.second = std::max(similarityRange_.second, similarity);
        channelSimilarity[channelProperty_.get()] = similarity;

        similarityPlot_[runName].push_back(channelSimilarity);
        setProgress(std::min(getProgress() + progressStep, 1.0f));

    });

    setProgress(1.0f);
    calculateButton_.setReadOnlyFlag(false);

    invalidate();
}

float EnsembleSimilarityPlot::calculateSimilarity(const VolumeBase* timestepVolume, float isoValue, float normalizationFactor = 1.0f) {
    const VolumeRAM_Float* volumeRepresentation = dynamic_cast<const VolumeRAM_Float*>(timestepVolume->getRepresentation<VolumeRAM>());
    tgt::ivec3 volumeDimensions = tgt::vec3(volumeRepresentation->getDimensions());
    float similarity = 0.0f;

    tgt::ivec3 llf = inport_.getData()->getRoi().getLLF();
    tgt::ivec3 urb = inport_.getData()->getRoi().getURB();

    size_t numberOfSamplePoints = ((volumeDimensions.x) * (volumeDimensions.y) * (volumeDimensions.z)) / sampleRate_;
    for(size_t z = llf.z; z < urb.z; z+=sampleRate_) {
        for(size_t y = llf.y; y < urb.y; y+=sampleRate_) {
            for(size_t x = llf.x; x < urb.x; x+=sampleRate_) {
                float voxelValue = volumeRepresentation->voxel(x,y,z);

                if(voxelValue > isoValue) {
                    similarity += 16.0f / static_cast<float>(numberOfSamplePoints);
                }
            }
        }
    }

    return similarity;
}

void EnsembleSimilarityPlot::mouseEvent(tgt::MouseEvent* e) {

    if(e->action() != tgt::MouseEvent::PRESSED)
        return;

    viewPortWidth_ = static_cast<float>(e->viewport().x);

    timeStepPosition_ = e->x() / viewPortWidth_;
    selectedTimeStep_.set(static_cast<int>(timeStepPosition_ * inport_.getData()->getMinNumTimeSteps()));

    // Convert to Normalized Device Coordinates.
    timeStepPosition_ = timeStepPosition_*2.0f - 1.0f;

    invalidate();
}

void EnsembleSimilarityPlot::onEvent(tgt::Event* e) {
    tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);

    if (!event)
        RenderProcessor::onEvent(e);
    else
        mouseEvent(event);
}

} // namespace voreen
