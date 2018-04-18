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

#include "fieldparallelplotviewer.h"
#include "fieldparallelplotcreator.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/textureunit.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/glsl.h"

#include "../ports/conditions/portconditionensemble.h"
#include "../utils/ensemblehash.h"
#include "../utils/utils.h"

namespace voreen {

static int PIXEL_BRUSH_ENABLE_DISTANCE = 3;
static const tgt::vec2 NO_SELECTION(-2.0f); // needs to be any value not inside range [-1, 1]

const std::string FieldParallelPlotViewer::loggerCat_("voreen.viscontest2018.FieldParallelPlotViewer");

FieldParallelPlotViewer::FieldParallelPlotViewer()
    : RenderProcessor()
    , plotDataInport_(Port::INPORT, "plotdatainport", "Field Plot Data Input")
    , ensembleInport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , outport_(Port::OUTPORT, "", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , volumeOutport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , transferFunc_("transferFunction", "Transfer Function for Plot")
    , renderedChannel_("channel", "Rendered Channel")
    , renderedRuns_("renderedRuns", "Rendered Runs")
    , volumeTransferFunc_("volumeTransferFunction", "Transfer Function for Volume")
    , valueRange_("valueRange", "Value Range", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
    , timeInterval_("timeInterval", "Selected Time Interval", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
    , selectedRuns_("selectedRuns", "Selected Runs")
    , hasLogarithmicDensity_("logarithmicDensity", "Logarithmic Density", true, Processor::VALID)
    , plotShader_("shaderprop", "Shader:", "fieldplot.frag", "passthrough.vert", "", Processor::INVALID_PROGRAM, Property::LOD_DEBUG)
    , timeStepPosition_(0.0f)
    , viewPortWidth_(0.0f)
    , requirePlotDataUpdate_(false)
    , requireConsistencyCheck_(false)
    , consistent_(false)
    , plotData_(nullptr)
    // UI
    , selectionStart_(NO_SELECTION)
    , selectionEnd_(NO_SELECTION)
    , isSelectionMode_(false)
{
    // Ports
    addPort(plotDataInport_);
        ON_CHANGE(plotDataInport_, FieldParallelPlotViewer, loadPlotData);
    addPort(ensembleInport_);
        ON_CHANGE_LAMBDA(ensembleInport_, [this] { requireConsistencyCheck_ = true; });
    addPort(outport_);
    addPort(volumeOutport_);

    // Properties
    addProperty(transferFunc_);
        transferFunc_.setGroupID("rendering");
    addProperty(hasLogarithmicDensity_);
        hasLogarithmicDensity_.setGroupID("rendering");
        ON_CHANGE_LAMBDA(hasLogarithmicDensity_, [this] { requirePlotDataUpdate_ = true; });
    addProperty(renderedChannel_);
        renderedChannel_.setGroupID("rendering");
        renderedChannel_.setReadOnlyFlag(true);
        ON_CHANGE(renderedChannel_, FieldParallelPlotViewer, switchChannel);
    addProperty(renderedRuns_);
        renderedRuns_.setGroupID("rendering");
    setPropertyGroupGuiName("rendering", "Rendering");

    addProperty(volumeTransferFunc_);
        volumeTransferFunc_.setGroupID("linking");
    addProperty(valueRange_);
        valueRange_.setGroupID("linking");
        valueRange_.setReadOnlyFlag(true);
        ON_CHANGE_LAMBDA(valueRange_, [this] { requirePlotDataUpdate_ = true; });
    addProperty(timeInterval_);
        timeInterval_.setGroupID("linking");
        timeInterval_.setReadOnlyFlag(true);
    addProperty(selectedRuns_);
        selectedRuns_.setGroupID("linking");
        selectedRuns_.setReadOnlyFlag(true);
    setPropertyGroupGuiName("linking", "Linking");

    addProperty(plotShader_);
        transferFunc_.setDomainFittingStrategy(TransFunc1DKeysProperty::FIT_DOMAIN_ALWAYS);
}

FieldParallelPlotViewer::~FieldParallelPlotViewer() {
}

Processor* FieldParallelPlotViewer::create() const {
    return new FieldParallelPlotViewer();
}

void FieldParallelPlotViewer::initialize() {
    RenderProcessor::initialize();
}

void FieldParallelPlotViewer::deinitialize() {
    currentPlot_.reset();
    plotTexture_.reset();

    RenderProcessor::deinitialize();
}

std::string FieldParallelPlotViewer::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = RenderProcessor::generateHeader();

    header += "#define NUM_RUNS " + (ensembleInport_.hasData() ? std::to_string(ensembleInport_.getData()->getRuns().size()) : "0") + "\n";
    header += transferFunc_.get()->getShaderDefines();

    return header;
}

bool FieldParallelPlotViewer::rebuildShader() {
    plotShader_.setHeader(generateHeader());
    return plotShader_.rebuild();
}

void FieldParallelPlotViewer::switchChannel() {

    size_t numRuns = ensembleInport_.getData()->getRuns().size();
    size_t firstSlice = renderedChannel_.getSelectedIndex() * numRuns;
    tgt::svec3 offset(0, 0, firstSlice);
    tgt::svec3 dimensions(plotDataInport_.getData()->getWidth(), plotDataInport_.getData()->getHeight(), numRuns);
    VolumeRAM_Float* slices = dynamic_cast<const VolumeRAM_Float*>(plotData_->getRepresentation<VolumeRAM>())->getSubVolume(dimensions, offset);
    tgtAssert(slices, "slices could not be read");

    // Set unfiltered slices to outport.
    channelSlices_ = new Volume(slices, tgt::svec3::one, tgt::svec3::zero);
    volumeOutport_.setData(channelSlices_, true);

    requirePlotDataUpdate_ = true;
}

void FieldParallelPlotViewer::updatePlotData() {
    if(!ensembleInport_.hasData()) return;

    // Apply range mask to selected plot data.
    VolumeRAM_Float* slices = dynamic_cast<const VolumeRAM_Float*>(channelSlices_->getRepresentation<VolumeRAM>())->clone();

    tgt::vec2 valueRange(slices->min(), slices->max());

    valueRange_.setMinValue(valueRange.x);
    valueRange_.setMaxValue(valueRange.y);
    //valueRange_.set(valueRange); // do NOT set, because it's gonna be linked
    valueRange_.blockCallbacks(false);

    // Apply threshold
    applyThreshold(slices);

    currentPlot_.reset(new Volume(slices, tgt::vec3::one, tgt::vec3::zero));
    plotTexture_.reset(new tgt::Texture(tgt::ivec3(slices->getDimensions()), GL_RED, GL_R32F, GL_FLOAT, tgt::Texture::NEAREST, tgt::Texture::CLAMP_TO_EDGE, static_cast<GLubyte*>(slices->getData()), false));
    plotTexture_->uploadTexture();
    plotTexture_->setCpuTextureData(0, false);
    
    // Fit transfer function to volume.
    transferFunc_.setVolume(currentPlot_.get());

    requirePlotDataUpdate_ = false;
}

void FieldParallelPlotViewer::loadPlotData() {

    // Reset old data.
    plotData_ = nullptr;

    if(!plotDataInport_.hasData())
        return;

    try {
        VolumeBase* plotData = plotDataInport_.getData()->getVolume();

        //if(!plotData->hasRepresentation<VolumeDisk>())
        //    throw VoreenException("Plot data has no Disk representation.");
        
        if(!plotData->hasMetaData(FieldParallelPlotCreator::META_DATA_HASH))
            throw VoreenException("Plot Dataset could not be loaded - Attribute " + FieldParallelPlotCreator::META_DATA_HASH + " missing.");

        plotData_ = plotData;

    } catch(VoreenException e) {
        VoreenApplication::app()->showMessageBox("Failed loading", e.what(), true);
        return;
    }

    // Trigger consistency check.
    requireConsistencyCheck_ = true;
}

bool FieldParallelPlotViewer::checkConsistency() {
    requireConsistencyCheck_ = false;

    selectedRuns_.reset();
    renderedRuns_.reset();
    renderedChannel_.setOptions(std::deque<Option<std::string>>());
    renderedChannel_.setReadOnlyFlag(true);
    valueRange_.setReadOnlyFlag(true);

    if (!ensembleInport_.isReady() || !plotData_)
        return false;

    const EnsembleDataset* dataset = ensembleInport_.getData();
    tgtAssert(dataset, "dataset must not be null");

    const MetaDataBase* ensembleHash = plotData_->getMetaData(FieldParallelPlotCreator::META_DATA_HASH);
    tgtAssert(ensembleHash, "Plot does not contain Hash");

    // Compare hashes.
    if (EnsembleHash(*dataset).getHash() != ensembleHash->toString())
        return false;

    renderedChannel_.blockCallbacks(true);
    for (const std::string& channel : dataset->getCommonChannels())
        renderedChannel_.addOption(channel, channel, channel);
    renderedChannel_.blockCallbacks(false);

    timeInterval_.setMinValue(dataset->getStartTime());
    timeInterval_.setMaxValue(dataset->getEndTime());
    timeInterval_.set(tgt::vec2(dataset->getStartTime(), dataset->getEndTime()));

    std::vector<int> renderedRunsIndices;
    for(const EnsembleDataset::Run& run : dataset->getRuns()) {
        renderedRunsIndices.push_back(static_cast<int>(renderedRunsIndices.size()));
        renderedRuns_.addRow(run.name_, run.color_);
        selectedRuns_.addRow(run.name_, run.color_);
    }
    renderedRuns_.setSelectedRowIndices(renderedRunsIndices);

    // Switch to default channel (first).
    switchChannel();

    // Update UI.
    renderedChannel_.setReadOnlyFlag(false);
    valueRange_.setReadOnlyFlag(false);

    // Build shader for ensemble.
    return rebuildShader();
}

void FieldParallelPlotViewer::applyThreshold(VolumeRAM_Float* volume) {
    for (size_t z = 0; z < volume->getDimensions().z; z++) {
        for(size_t x = 0; x < volume->getDimensions().x; x++) {
            for(size_t y = 0; y < volume->getDimensions().y; y++) {
                float& voxel = volume->voxel(x, y, z);
                if(voxel < valueRange_.get().x || voxel > valueRange_.get().y)
                    voxel = 0.0f;
                if(hasLogarithmicDensity_.get() && voxel != 0.0f)
                    voxel = voxel > 0 ? std::log(voxel) : -std::log(std::abs(voxel));
            }
        }
    }
    volume->invalidate();
}

void FieldParallelPlotViewer::mouseEvent(tgt::MouseEvent* e) {

    float positionX = (e->x() / static_cast<float>(e->viewport().x)) * 2.0f - 1.0f;
    float positionY = ((e->y() / static_cast<float>(e->viewport().y)) * 2.0f - 1.0f) * -1.0f;

    switch (e->action()) {
    case tgt::MouseEvent::PRESSED:
        if (selectionStart_ == NO_SELECTION || selectionEnd_ == NO_SELECTION) {
            isSelectionMode_ = true;
            selectionStart_ = tgt::vec2(positionX, positionY);
        }
        else {
            selectionStart_ = NO_SELECTION;
            selectionEnd_ = NO_SELECTION;
            updateSelection();
        }
        break;
    case tgt::MouseEvent::MOTION:
        if (isSelectionMode_) {
            selectionEnd_ = tgt::vec2(positionX, positionY);
        }
        break;
    case tgt::MouseEvent::RELEASED:
        isSelectionMode_ = false;
        updateSelection();
        break;
    default:
        e->ignore();
        return;
    }

    e->accept();

    // Trigger redraw.
    invalidate();
}

void FieldParallelPlotViewer::onEvent(tgt::Event* e) {

    if(!isReady())
        return;

    tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);

    if (!event)
        RenderProcessor::onEvent(e);
    else
        mouseEvent(event);
}

bool FieldParallelPlotViewer::isReady() const {
    return plotDataInport_.isReady() && ensembleInport_.isReady();
}

void FieldParallelPlotViewer::process() {

    if (requireConsistencyCheck_) {
        consistent_ = checkConsistency();
        if(!consistent_)
            VoreenApplication::app()->showMessageBox("Mismatching data", "The plot was probably not generated by the ensemble dataset passed to the FieldParallelPlotViewer", true);
    }

    if (!consistent_)
        return;

    if(requirePlotDataUpdate_)
        updatePlotData();

    outport_.activateTarget();
    outport_.clearTarget();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    /////////////////////////////////////////
    // draw plot

    tgt::Shader* shader = plotShader_.getShader();
    shader->activate();
    setGlobalShaderParameters(shader);
    LGL_ERROR;

    // Setup plot texture
    tgt::TextureUnit plotUnit;
    plotUnit.activate();
    plotTexture_->bind();
    shader->setUniform("plotData_", plotUnit.getUnitNumber());

    // Setup selection.
    std::vector<GLint> runs(ensembleInport_.getData()->getRuns().size());
    for (size_t i = 0; i < runs.size(); i++)
        runs[i] = 0;

    for (int i : renderedRuns_.getSelectedRowIndices())
        runs[i] = 1;

    shader->setUniform("renderedRuns_", runs.data(), runs.size());

    // Setup transfer function.
    tgt::TextureUnit transferUnit;
    transferUnit.activate();
    transferFunc_.get()->getTexture()->bind();
    transferFunc_.get()->setUniform(shader, "transFuncParams_", "transFuncTex_", transferUnit.getUnitNumber());
    LGL_ERROR;

    renderQuad();
    LGL_ERROR;

    shader->deactivate();

    /////////////////////////////////////////
    // draw selection

    if(selectionStart_ != NO_SELECTION && selectionEnd_ != NO_SELECTION) {
        // draw run quad
        glDepthFunc(GL_ALWAYS);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        IMode.color(1.0f, 1.0f, 1.0f, 0.2f);
        IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(selectionStart_.x, selectionStart_.y);
        IMode.vertex(selectionEnd_.x, selectionStart_.y);
        IMode.vertex(selectionEnd_.x, selectionEnd_.y);
        IMode.vertex(selectionStart_.x, selectionEnd_.y);
        IMode.end();
        glLineWidth(1.0f);
        glDisable(GL_BLEND);
        glDepthFunc(GL_LESS);
        glBlendFunc(GL_ONE, GL_ZERO);

    }

    float startTime = ensembleInport_.getData()->getStartTime();
    float endTime = ensembleInport_.getData()->getEndTime();
    float duration = endTime - startTime;
    if(timeInterval_.get().x != startTime || timeInterval_.get().y != endTime) {
        tgt::vec2 timeBounds(-1, 1);
        timeBounds.x = timeInterval_.get().x / duration * 2.0f - 1.0f;
        timeBounds.y = timeInterval_.get().y / duration * 2.0f - 1.0f;

        // draw time min and max lines
        glDepthFunc(GL_ALWAYS);
        IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
        IMode.begin(tgt::ImmediateMode::LINES);
        IMode.vertex(timeBounds.x, -1);
        IMode.vertex(timeBounds.x, 1);
        IMode.vertex(timeBounds.y, 1);
        IMode.vertex(timeBounds.y, -1);
        IMode.end();
        glLineWidth(1.0f);
        glDepthFunc(GL_LESS);
    }

    /////////////////////////////////////////

    IMode.color(tgt::vec4::one);

    tgt::TextureUnit::setZeroUnit();
    outport_.deactivateTarget();
    LGL_ERROR;

}

void FieldParallelPlotViewer::updateSelection() {

    bool selecting = selectionStart_ != NO_SELECTION && selectionEnd_ != NO_SELECTION;

    float minX = std::min(selectionStart_.x, selectionEnd_.x) + 1;
    float maxX = std::max(selectionStart_.x, selectionEnd_.x) + 1;
    float minY = std::min(selectionStart_.y, selectionEnd_.y) + 1;
    float maxY = std::max(selectionStart_.y, selectionEnd_.y) + 1;

    /**
     * update run selection
     */
    size_t width = plotData_->getDimensions().x;
    size_t height = plotData_->getDimensions().y;

    size_t selectedMinPixelX;
    size_t selectedMaxPixelX;
    size_t selectedMinPixelY;
    size_t selectedMaxPixelY;

    if(selecting) {
        selectedMinPixelX = static_cast<size_t>(width * minX / 2.0f);
        selectedMaxPixelX = static_cast<size_t>(width * maxX / 2.0f);
        selectedMinPixelY = static_cast<size_t>(height * minY / 2.0f);
        selectedMaxPixelY = static_cast<size_t>(height * maxY / 2.0f);
    }
    else {
        selectedMinPixelX = 0;
        selectedMaxPixelX = plotData_->getDimensions().x;
        selectedMinPixelY = 0;
        selectedMaxPixelY = plotData_->getDimensions().y;
    }

    const VolumeRAM_Float* slices = dynamic_cast<const VolumeRAM_Float*>(channelSlices_->getRepresentation<VolumeRAM>());
    tgtAssert(slices, "slices null");

    std::vector<int> relevantRuns;
    for(int run : renderedRuns_.getSelectedRowIndices()) {
        bool isRelevant = false;
        for(size_t y = selectedMinPixelY; y < selectedMaxPixelY; y++) {
            for(size_t x = selectedMinPixelX; x < selectedMaxPixelX; x++) {
                // We assume that a value from 0.0f means no data has been set
                // which might include some rare cases with ranges including negative numbers
                // were we generate a false negative. For now, we leave it that way.
                if(slices->voxel(x,y,run) != 0.0f) {
                    relevantRuns.push_back(run);
                    isRelevant = true;
                    break;
                }
            }
            if(isRelevant) break;
        }
    }

    selectedRuns_.setSelectedRowIndices(relevantRuns);

    /**
     * update time selection
     */

    float startTime = ensembleInport_.getData()->getStartTime();
    float endTime = ensembleInport_.getData()->getEndTime();
    float duration = endTime - startTime;

    float selectedMinTime;
    float selectedMaxTime;

    if(selecting) {
        selectedMinTime = startTime + duration * minX / 2.0f;
        selectedMaxTime = startTime + duration * maxX / 2.0f;
    }
    else {
        selectedMinTime = startTime;
        selectedMaxTime = endTime;
    }

    timeInterval_.set(tgt::vec2(selectedMinTime, selectedMaxTime));

    /**
     * update volume transfer function
     */
    if (selecting)
        volumeTransferFunc_.get()->setThreshold(mapRange(selectedMinPixelY, size_t(0), slices->getDimensions().y, 0.0f, 1.0f), mapRange(selectedMaxPixelY, size_t(0), slices->getDimensions().y, 0.0f, 1.0f));
    else
        volumeTransferFunc_.get()->setThreshold(0.0f, 1.0f);
    volumeTransferFunc_.invalidate();
}


} // namespace voreen
