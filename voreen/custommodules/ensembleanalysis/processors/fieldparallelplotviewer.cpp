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

#include "modules/plotting/utils/plotlibrary/plotlibrary.h"
#include "modules/plotting/utils/plotlibrary/plotlibraryopengl.h"

#include "../ports/conditions/portconditionensemble.h"
#include "../utils/ensemblehash.h"
#include "../utils/utils.h"

namespace voreen {

static const int PIXEL_BRUSH_ENABLE_DISTANCE = 3;
static const tgt::ivec2 MARGINS(75, 50);
static const tgt::vec2 NO_SELECTION(-2.0f); // needs to be any value not inside range [-1, 1]

const std::string FieldParallelPlotViewer::loggerCat_("voreen.ensembleanalysis.FieldParallelPlotViewer");

FieldParallelPlotViewer::FieldParallelPlotViewer()
    : RenderProcessor()
    , plotDataInport_(Port::INPORT, "plotdatainport", "Field Plot Data Input")
    , ensembleInport_(Port::INPORT, "ensembleinport", "Ensemble Data Input")
    , outport_(Port::OUTPORT, "imageOutport", "Image Output", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , volumeOutport_(Port::OUTPORT, "volumehandle.volumehandle", "Volume Output")
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    , transferFunc_("transferFunction", "Transfer Function for Plot")
    , renderedField_("channel", "Rendered Field")
    , renderedRuns_("renderedRuns", "Rendered Runs")
    , volumeTransferFunc_("volumeTransferFunction", "Transfer Function for Volume")
    , valueRange_("valueRange", "Value Range", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
    , timeInterval_("timeInterval", "Selected Time Interval", tgt::vec2(0.0f, 0.0f), 0.0f, 0.0f)
    , selectedRuns_("selectedRuns", "Selected Runs")
    , hasLogarithmicDensity_("logarithmicDensity", "Logarithmic Density", false)
    , zoomX_("zoomX", "Zoom X", tgt::vec2(0.0, 1.0), 0.0f, 1.0f)
    , zoomY_("zoomY", "Zoom Y", tgt::vec2(0.0, 1.0), 0.0f, 1.0f)
    , xUnit_("xUnit", "x-Unit", "time")
    , yUnit_("yUnit", "y-Unit", "value")
    , fontSize_("fontSize", "Font Size", 10, 1, 30)
    , plotShader_("shaderprop", "Plot Shader", "fieldplot.frag", "passthrough.vert", "", Processor::INVALID_RESULT, Property::LOD_DEBUG)
    , timeStepPosition_(0.0f)
    , viewPortWidth_(0.0f)
    , plotData_(nullptr)
    // UI
    , plotLib_(new PlotLibraryOpenGl())
    , selectionStart_(NO_SELECTION)
    , selectionEnd_(NO_SELECTION)
    , isSelectionMode_(false)
{
    // Ports
    addPort(ensembleInport_);
    addPort(plotDataInport_);
    addPort(outport_);
    addPort(volumeOutport_);
    addPrivateRenderPort(privatePort_);

    // Properties
    addProperty(transferFunc_);
        transferFunc_.setGroupID("rendering");
    addProperty(hasLogarithmicDensity_);
        ON_CHANGE(hasLogarithmicDensity_, FieldParallelPlotViewer, adjustPropertiesToInput);
        hasLogarithmicDensity_.setGroupID("rendering");
    addProperty(zoomX_);
        zoomX_.setGroupID("rendering");
    addProperty(zoomY_);
        zoomY_.setGroupID("rendering");
    addProperty(xUnit_);
        xUnit_.setGroupID("rendering");
    addProperty(yUnit_);
        yUnit_.setGroupID("rendering");
    addProperty(fontSize_);
        fontSize_.setGroupID("rendering");
    addProperty(renderedField_);
        renderedField_.setGroupID("rendering");
        renderedField_.setReadOnlyFlag(true);
        ON_CHANGE(renderedField_, FieldParallelPlotViewer, switchField);
    addProperty(renderedRuns_);
        renderedRuns_.setGroupID("rendering");
    setPropertyGroupGuiName("rendering", "Rendering");

    addProperty(volumeTransferFunc_);
        volumeTransferFunc_.setDomainFittingStrategy(TransFunc1DKeysProperty::FIT_DOMAIN_ALWAYS);
        volumeTransferFunc_.setGroupID("linking");
    addProperty(valueRange_);
        //ON_CHANGE(valueRange_, FieldParallelPlotViewer, adjustPropertiesToInput);
        valueRange_.setGroupID("linking");
        valueRange_.setReadOnlyFlag(true);
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

    RenderProcessor::deinitialize();
}

std::string FieldParallelPlotViewer::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = RenderProcessor::generateHeader(version);

    header += "#define NUM_RUNS " + (ensembleInport_.hasData() ? std::to_string(ensembleInport_.getData()->getRuns().size()) : "0") + "\n";
    header += transferFunc_.get()->getShaderDefines();

    return header;
}

bool FieldParallelPlotViewer::rebuildShader() {
    plotShader_.setHeader(generateHeader());
    return plotShader_.rebuild();
}

void FieldParallelPlotViewer::adjustPropertiesToInput() {

    plotData_ = nullptr;
    selectedRuns_.reset();
    renderedRuns_.reset();
    renderedField_.setOptions(std::deque<Option<std::string>>());
    renderedField_.setReadOnlyFlag(true);
    valueRange_.setReadOnlyFlag(true);

    if (!ensembleInport_.isReady())
        return;

    const EnsembleDataset* dataset = ensembleInport_.getData();

    if(!plotDataInport_.isReady())
        return;

    const VolumeBase* plotData = plotDataInport_.getData()->getVolume();
    const MetaDataBase* ensembleHash = plotData->getMetaData(FieldParallelPlotCreator::META_DATA_HASH);
    tgtAssert(ensembleHash, "Plot does not contain Hash");

    // Compare hashes.
    if (EnsembleHash(*dataset).getHash() != ensembleHash->toString())
        //LERROR("HASH does not match - currently ignored, application likely to crash");
        return;

    plotData_ = plotData;

    renderedField_.blockCallbacks(true);
    for (const std::string& fieldName : dataset->getCommonFieldNames())
        renderedField_.addOption(fieldName, fieldName, fieldName);
    renderedField_.blockCallbacks(false);

    timeInterval_.setMinValue(dataset->getStartTime());
    timeInterval_.setMaxValue(dataset->getEndTime());
    timeInterval_.set(tgt::vec2(dataset->getStartTime(), dataset->getEndTime()));

    std::vector<int> renderedRunsIndices;
    for(const EnsembleDataset::Run& run : dataset->getRuns()) {
        renderedRuns_.addRow(run.getName(), run.getColor());
        selectedRuns_.addRow(run.getName(), run.getColor());
        renderedRunsIndices.push_back(static_cast<int>(renderedRunsIndices.size()));
    }
    renderedRuns_.setSelectedRowIndices(renderedRunsIndices);

    // Switch to default field (first).
    switchField();

    // Update UI.
    renderedField_.setReadOnlyFlag(false);
    valueRange_.setReadOnlyFlag(false);

    // Build shader for ensemble.
    rebuildShader();

    // Apply range mask to selected plot data.
    VolumeRAM_Float* slices = dynamic_cast<const VolumeRAM_Float*>(fieldSlices_->getRepresentation<VolumeRAM>())->clone();

    tgt::vec2 valueRange(slices->min(), slices->max());

    valueRange_.setMinValue(valueRange.x);
    valueRange_.setMaxValue(valueRange.y);
    //valueRange_.set(valueRange); // do NOT set, because it's gonna be linked
    valueRange_.blockCallbacks(false);

    // Apply threshold
    applyThreshold(slices);

    currentPlot_.reset(new Volume(slices, tgt::vec3::one, tgt::vec3::zero));

    // Fit transfer function to volume.
    transferFunc_.setVolume(currentPlot_.get());
    volumeTransferFunc_.get()->setDomain(valueRange);
}

void FieldParallelPlotViewer::switchField() {

    size_t numRuns = ensembleInport_.getData()->getRuns().size();
    size_t firstSlice = renderedField_.getSelectedIndex() * numRuns;
    tgt::svec3 offset(0, 0, firstSlice);
    tgt::svec3 dimensions(plotDataInport_.getData()->getWidth(), plotDataInport_.getData()->getHeight(), numRuns);
    VolumeRAM_Float* slices = dynamic_cast<const VolumeRAM_Float*>(plotData_->getRepresentation<VolumeRAM>())->getSubVolume(dimensions, offset);
    tgtAssert(slices, "slices could not be read");

    // Set unfiltered slices to outport.
    fieldSlices_ = new Volume(slices, tgt::svec3::one, tgt::svec3::zero);
    volumeOutport_.setData(fieldSlices_, true);
}

void FieldParallelPlotViewer::applyThreshold(VolumeRAM_Float* volume) {
    for (size_t z = 0; z < volume->getDimensions().z; z++) {
        for(size_t x = 0; x < volume->getDimensions().x; x++) {
            for(size_t y = 0; y < volume->getDimensions().y; y++) {
                float& voxel = volume->voxel(x, y, z);
                if(voxel < valueRange_.get().x || voxel > valueRange_.get().y)
                    voxel = 0;
                else if(hasLogarithmicDensity_.get())
                    voxel = std::log(voxel+1); // shift the function such that we obtain only positive values
            }
        }
    }
    volume->invalidate();
}

void FieldParallelPlotViewer::mouseEvent(tgt::MouseEvent* e) {

    int ex = tgt::clamp(e->x(), MARGINS.x, e->viewport().x - MARGINS.x);
    int ey = tgt::clamp(e->y(), MARGINS.y, e->viewport().y - MARGINS.y);

    // TODO: optimize
    if(ex != e->x() || ey != e->y())
        return;

    float mx = mapRange(ex, 0, e->viewport().x, -1.0f,  1.0f);
    float my = mapRange(ey, 0, e->viewport().y,  1.0f, -1.0f);

    switch (e->action()) {
    case tgt::MouseEvent::PRESSED:
        if (selectionStart_ == NO_SELECTION || selectionEnd_ == NO_SELECTION) {
            isSelectionMode_ = true;
            selectionStart_ = tgt::vec2(mx, my);
        }
        else {
            selectionStart_ = NO_SELECTION;
            selectionEnd_ = NO_SELECTION;
            updateSelection();
        }
        break;
    case tgt::MouseEvent::MOTION:
        if (isSelectionMode_) {
            selectionEnd_ = tgt::vec2(mx, my);
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

    if(!ensembleInport_.isReady()) {
        setNotReadyErrorMessage("No Ensemble Input");
        return false;
    }

    if(!plotDataInport_.isReady()) {
        setNotReadyErrorMessage("No Plot Data Input");
        return false;
    }

    if(ensembleInport_.getData()->getCommonFieldNames().empty()) {
        setNotReadyErrorMessage("No common fields available");
        return false;
    }

    if(!plotData_) {
        setNotReadyErrorMessage("No Plot data available");
        return false;
    }

    return true;
}

void FieldParallelPlotViewer::beforeProcess() {
    RenderProcessor::beforeProcess();


}

void FieldParallelPlotViewer::process() {
    // Perform the actual rendering by embedding the rendered
    // plot inside a coordinate system generated by the plot library
    renderPlot();
    renderAxes();
}

void FieldParallelPlotViewer::renderPlot() {

    if(privatePort_.getSize() != outport_.getSize())
        privatePort_.resize(outport_.getSize());

    privatePort_.activateTarget();
    privatePort_.clearTarget();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    tgt::Shader* shader = plotShader_.getShader();
    shader->activate();
    setGlobalShaderParameters(shader);
    LGL_ERROR;

    // Setup plot texture.
    tgt::TextureUnit plotUnit;
    plotUnit.activate();
    currentPlot_->getRepresentation<VolumeGL>()->getTexture()->bind();
    shader->setUniform("plotData_", plotUnit.getUnitNumber());

    // Set zoom region.
    shader->setUniform("rangeX_", zoomX_.get());
    shader->setUniform("rangeY_", zoomY_.get());

    // Setup selection.
    std::vector<GLint> runs(ensembleInport_.getData()->getRuns().size(), 0);
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

    IMode.color(tgt::vec4::one);
    tgt::TextureUnit::setZeroUnit();
    privatePort_.deactivateTarget();
}

void FieldParallelPlotViewer::renderAxes() {

    outport_.activateTarget();
    outport_.clearTarget();

    // Set Plot status.
    plotLib_->setWindowSize(outport_.getSize());
    plotLib_->setAxesWidth(1.0f);
    plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
    plotLib_->setLineWidth(1.0f);
    plotLib_->setMaxGlyphSize(1.0f);
    plotLib_->setMarginBottom(MARGINS.y);
    plotLib_->setMarginTop(MARGINS.y);
    plotLib_->setMarginLeft(MARGINS.x);
    plotLib_->setMarginRight(MARGINS.x);
    plotLib_->setMinimumScaleStep(32, PlotLibrary::X_AXIS);
    plotLib_->setMinimumScaleStep(32, PlotLibrary::Y_AXIS);
    tgt::vec2 valueRange = ensembleInport_.getData()->getValueRange(renderedField_.get());
    tgt::vec2 zoomedTimeRange = applyZoomToRange(tgt::vec2(timeInterval_.getMinValue(), timeInterval_.getMaxValue()), zoomX_.get());
    tgt::vec2 zoomedValueRange = applyZoomToRange(tgt::vec2(valueRange.x, valueRange.y), zoomY_.get());

    plotLib_->setDomain(Interval<plot_t>(zoomedTimeRange.x, zoomedTimeRange.y), PlotLibrary::X_AXIS);
    plotLib_->setDomain(Interval<plot_t>(zoomedValueRange.x, zoomedValueRange.y), PlotLibrary::Y_AXIS);

    if (plotLib_->setRenderStatus()) {
        plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->renderAxes();
        plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
        plotLib_->setFontSize(fontSize_.get());
        plotLib_->setFontColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->renderAxisScales(PlotLibrary::X_AXIS, false);
        plotLib_->renderAxisScales(PlotLibrary::Y_AXIS, false);
        plotLib_->setFontSize(fontSize_.get() + 2);
        plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, xUnit_.get());
        plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, yUnit_.get());
    }
    plotLib_->resetRenderStatus();

    // Plot data
    float xMinMarginNDC = mapRange(MARGINS.x, 0, outport_.getSize().x, -1.0f, 1.0f);
    float xMaxMarginNDC = mapRange(outport_.getSize().x-MARGINS.x, 0, outport_.getSize().x, -1.0f, 1.0f);
    float yMinMarginNDC = mapRange(MARGINS.y, 0, outport_.getSize().y, -1.0f, 1.0f);
    float yMaxMarginNDC = mapRange(outport_.getSize().y-MARGINS.y, 0, outport_.getSize().y, -1.0f, 1.0f);

    tgt::TextureUnit::setZeroUnit();
    tgt::Texture* texture = privatePort_.getColorTexture();
    if(texture) {
        texture->enable();
        texture->bind();
    }

    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.texcoord(0.0f, 0.0f); IMode.vertex(xMinMarginNDC, yMinMarginNDC);
        IMode.texcoord(1.0f, 0.0f); IMode.vertex(xMaxMarginNDC, yMinMarginNDC);
        IMode.texcoord(1.0f, 1.0f); IMode.vertex(xMaxMarginNDC, yMaxMarginNDC);
        IMode.texcoord(0.0f, 1.0f); IMode.vertex(xMinMarginNDC, yMaxMarginNDC);
    IMode.end();

    if(texture)
        texture->disable();

    ////////////////////////////////
    /// Draw selection.

    // draw run quad
    if(selectionStart_ != NO_SELECTION && selectionEnd_ != NO_SELECTION) {

        glDepthFunc(GL_ALWAYS);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        IMode.color(0.0f, 0.0f, 0.0f, 0.3f);
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(selectionStart_.x, selectionStart_.y);
            IMode.vertex(selectionEnd_.x, selectionStart_.y);
            IMode.vertex(selectionEnd_.x, selectionEnd_.y);
            IMode.vertex(selectionStart_.x, selectionEnd_.y);
        IMode.end();
        glDisable(GL_BLEND);
        glDepthFunc(GL_LESS);
        glBlendFunc(GL_ONE, GL_ZERO);

        if(!isSelectionMode_) {

            float startTime = zoomedTimeRange.x;
            float endTime = zoomedTimeRange.y;
            tgt::vec2 timeBounds = mapRange(timeInterval_.get(), tgt::vec2(startTime), tgt::vec2(endTime), -tgt::vec2::one, tgt::vec2::one);

            // Apply margins
            timeBounds *= ((outport_.getSize().x - 2.0f*MARGINS.x) / outport_.getSize().x);

            // draw time min and max lines
            glDepthFunc(GL_ALWAYS);
            IMode.color(0.0f, 0.0f, 0.0f, 1.0f);
            glLineWidth(3.0f);
            IMode.begin(tgt::ImmediateMode::LINES);
                IMode.vertex(timeBounds.x, yMinMarginNDC);
                IMode.vertex(timeBounds.x, yMaxMarginNDC);
                IMode.vertex(timeBounds.y, yMaxMarginNDC);
                IMode.vertex(timeBounds.y, yMinMarginNDC);
            IMode.end();
            glLineWidth(1.0f);
            glDepthFunc(GL_LESS);
        }

    }

    // Reset state.
    IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
    outport_.deactivateTarget();
    LGL_ERROR;
}

tgt::vec2 FieldParallelPlotViewer::applyZoomToRange(const tgt::vec2& range, const tgt::vec2& zoom) {
    float lowerLimit = range.x + (range.y - range.x) * zoom.x;
    float upperLimit = range.y * zoom.y;
    return tgt::vec2(lowerLimit, upperLimit);
}

void FieldParallelPlotViewer::updateSelection() {

    const EnsembleDataset* ensemble = ensembleInport_.getData();
    if(!ensemble)
        return;

    bool selecting = selectionStart_ != NO_SELECTION && selectionEnd_ != NO_SELECTION;

    tgt::vec2 selectionX = tgt::vec2(std::min(selectionStart_.x, selectionEnd_.x), std::max(selectionStart_.x, selectionEnd_.x));
    tgt::vec2 selectionY = tgt::vec2(std::min(selectionStart_.y, selectionEnd_.y), std::max(selectionStart_.y, selectionEnd_.y));

    // Map to whole canvas size.
    selectionX = selectionX * (outport_.getSize().x / (outport_.getSize().x - 2.0f*MARGINS.x));
    selectionY = selectionY * (outport_.getSize().y / (outport_.getSize().y - 2.0f*MARGINS.y));

    /**
     * update run selection
     */
    size_t width  = plotData_->getDimensions().x;
    size_t height = plotData_->getDimensions().y;

    tgt::ivec2 selectedPixelX(0, width);
    tgt::ivec2 selectedPixelY(0, height);

    if(selecting) {
        selectedPixelX = tgt::ivec2(mapRange(selectionX, -tgt::vec2::one, tgt::vec2::one, tgt::vec2::zero, tgt::vec2(width)));
        selectedPixelY = tgt::ivec2(mapRange(selectionY, -tgt::vec2::one, tgt::vec2::one, tgt::vec2::zero, tgt::vec2(height)));

        const VolumeRAM_Float* slices = dynamic_cast<const VolumeRAM_Float*>(fieldSlices_->getRepresentation<VolumeRAM>());
        tgtAssert(slices, "slices null");

        std::vector<int> relevantRuns;
        for (int run : renderedRuns_.getSelectedRowIndices()) {
            bool isRelevant = false;
            for (int y = selectedPixelY.x; y < selectedPixelY.y; y++) {
                for (int x = selectedPixelX.x; x < selectedPixelX.y; x++) {
                    // We assume that a value from 0.0f means no data has been set
                    // which might include some rare cases with ranges including negative numbers
                    // were we generate a false negative. For now, we leave it that way.
                    if (slices->voxel(x, y, run) != 0.0f) {
                        relevantRuns.push_back(run);
                        isRelevant = true;
                        break;
                    }
                }
                if (isRelevant) break;
            }
        }

        selectedRuns_.setSelectedRowIndices(relevantRuns);
    }
    else {
        selectedRuns_.setSelectedRowIndices(renderedRuns_.get());
    }

    /**
     * update time selection
     */
    float startTime = ensemble->getStartTime();
    float endTime = ensemble->getEndTime();

    tgt::vec2 zoomedTimeRange = applyZoomToRange(tgt::vec2(startTime, endTime), zoomX_.get());
    startTime = zoomedTimeRange.x;
    endTime = zoomedTimeRange.y;

    if(selecting) {
        tgt::vec2 selectedTime = mapRange(selectionX, -tgt::vec2::one, tgt::vec2::one, tgt::vec2(startTime), tgt::vec2(endTime));
        timeInterval_.set(selectedTime);
    }

    /**
     * update volume transfer function
     */
    if (selecting) {
        volumeTransferFunc_.get()->setThreshold(mapRange(selectionY, -tgt::vec2::one, tgt::vec2::one, tgt::vec2::zero, tgt::vec2::one));
    }
    else {
        volumeTransferFunc_.get()->setThreshold(0.0f, 1.0f);
    }
    volumeTransferFunc_.invalidate();
}

} // namespace voreen
