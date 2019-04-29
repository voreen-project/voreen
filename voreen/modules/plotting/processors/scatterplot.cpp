/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "scatterplot.h"

#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/utils/multisampler.h"
#include "../interaction/plotcamerainteractionhandler.h"

#include "tgt/logmanager.h"
#include "tgt/texturemanager.h"

namespace voreen {

const std::string ScatterPlot::loggerCat_("voreen.ScatterPlot");

ScatterPlot::ScatterPlot()
    : PlotProcessor(PlotEntitySettings::SCATTER, false)
    , lighting_("lighting", "Enable Lighting", false)
{
    plotEntitiesProp_.setGuiName("Scatter Data");
    addZLabelProperties();
    addProperty(selectionPlaneColor_);
    addProperty(renderMousePosition_);
    addProperty(renderLegend_);
    addProperty(orthographicCamera_);
    addProperty(lighting_);
    addProperty(renderXHelperLines_);
    addProperty(renderYHelperLines_);
    addProperty(renderZHelperLines_);
    addProperty(xScaleStep_);
    addProperty(yScaleStep_);
    addProperty(zScaleStep_);
    addProperty(marginLeft_);
    addProperty(marginRight_);
    addProperty(marginBottom_);
    addProperty(marginTop_);
    addProperty(camera_);

    // group properties
    lighting_.setGroupID("general");

    addEventProperty(eventHighlight_);
    addEventProperty(eventLabel_);
    addEventProperty(eventZoom_);
    addEventProperty(eventHighlightAdditive_);
    addEventProperty(eventLabelAdditive_);
    addEventProperty(eventZoomAdditive_);
    addEventProperty(mousePositionUpdateEvent_);
    addEventProperty(mouseEventEnterExit_);

    addInteractionHandler(cameraHandler_);
    addInteractionHandler(plotCameraHandler_);
    orthographicCamera_.onChange(MemberFunctionCallback<ScatterPlot>(this, &ScatterPlot::toggleProperties));
}

Processor* ScatterPlot::create() const {
    return new ScatterPlot();
}

void ScatterPlot::render() {
    {
        Multisampler m(outport_);
        plotLib_->setUsePlotPickingManager(false);
        setPlotStatus();
        if (plotLib_->setRenderStatus()) {
            renderAxes();
            renderData();
            renderPlotLabel();
            createPlotLabels();
            plotLib_->renderPlotLabels();
            renderMousePosition();
            renderLegends();
            if (threeDimensional_)
                renderSelectionPlanes();
        }
        plotLib_->resetRenderStatus();
        if (!threeDimensional_)
            renderSelectedRegion();
    }
    plotPickingManager_.activateTarget();
    plotPickingManager_.clearTarget();
    if (enablePicking_.get()) {
        plotLib_->setUsePlotPickingManager(true);
        if (plotLib_->setRenderStatus())
            renderData();
        plotLib_->resetRenderStatus();
    }
    plotPickingManager_.deactivateTarget();
}

void ScatterPlot::setPlotStatus(){
    plotLib_->setWindowSize(outport_.getSize());
    plotLib_->setAxesWidth(axesWidth_.get());
    plotLib_->setCamera(camera_.get());
    plotLib_->setMinimumScaleStep(xScaleStep_.get(), PlotLibrary::X_AXIS);
    plotLib_->setMinimumScaleStep(yScaleStep_.get(), PlotLibrary::Y_AXIS);
    plotLib_->setMinimumScaleStep(zScaleStep_.get(), PlotLibrary::Z_AXIS);
    plotLib_->setMarginBottom(marginBottom_.get());
    plotLib_->setMarginTop(marginTop_.get());
    plotLib_->setMarginLeft(marginLeft_.get());
    plotLib_->setMarginRight(marginRight_.get());
    plotLib_->setDimension(threeDimensional_ ? PlotLibrary::THREE : PlotLibrary::TWO);
    plotLib_->setLightingFlag(lighting_.get());
    plotLib_->setOrthographicCameraFlag(orthographicCamera_.get());
    plotLib_->setDomain(selectionProp_.getZoom().xZoom_, PlotLibrary::X_AXIS);
    plotLib_->setDomain(selectionProp_.getZoom().yZoom_, PlotLibrary::Y_AXIS);
    plotLib_->setDomain(selectionProp_.getZoom().zZoom_, PlotLibrary::Z_AXIS);
}

void ScatterPlot::renderData() {
    plotLib_->setHighlightColor(highlightColor_.get());
    std::vector<PlotEntitySettings>::const_iterator it;
    for (it = plotEntitiesProp_.get().begin(); it < plotEntitiesProp_.get().end(); ++it) {
        plotLib_->setDrawingColor(it->getFirstColor());
        plotLib_->setColorMap(it->getColorMap());
        plotLib_->setMinGlyphSize(it->getMinGlyphSize());
        plotLib_->setMaxGlyphSize(it->getMaxGlyphSize());
        plotLib_->setGlyphStyle(it->getGlyphStyle());
        if (it->getUseTextureFlag())
            loadTextures(it->getTexturePath());
        else
            loadTextures("");
        threeDimensional_ ?
            plotLib_->renderScatter(data_, plotEntitiesProp_.getXColumnIndex(),
                                        plotEntitiesProp_.getYColumnIndex(),
                                        it->getMainColumnIndex(),
                                        it->getOptionalColumnIndex(),
                                        it->getSecondOptionalColumnIndex())
          : plotLib_->renderScatter(data_, plotEntitiesProp_.getXColumnIndex(),
                                        it->getMainColumnIndex(),
                                        -1,
                                        it->getOptionalColumnIndex(),
                                        it->getSecondOptionalColumnIndex());
    }
    LGL_ERROR;
}

void ScatterPlot::renderAxes() {
    //render axes
    if (renderAxes_.get()) {
        plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->setFontColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->renderAxes();
        plotLib_->setFontSize(10);
        if (renderScales_.get()) {
            if (threeDimensional_) {
                plotLib_->renderAxisScales(PlotLibrary::X_AXIS, renderXHelperLines_.get(), getXLabel());
                plotLib_->renderAxisScales(PlotLibrary::Y_AXIS, renderYHelperLines_.get(), getYLabel());
                plotLib_->renderAxisScales(PlotLibrary::Z_AXIS, renderZHelperLines_.get(), getZLabel());
            }
            else {
                plotLib_->renderAxisScales(PlotLibrary::X_AXIS, renderXHelperLines_.get());
                plotLib_->renderAxisScales(PlotLibrary::Y_AXIS, renderYHelperLines_.get());
                plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, getXLabel());
                plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, getYLabel());
            }
        }
    }
}

void ScatterPlot::readFromInport() {
    if (dynamic_cast<const PlotData*>(inport_.getData())) {
        data_ = *dynamic_cast<const PlotData*>(inport_.getData());
        dataProp_.set(&data_);
        plotEntitiesProp_.setPlotData(&data_);
        selectionProp_.setPlotData(&data_);
        plotPickingManager_.setColumnCount(data_.getColumnCount());
    }
    else {
        LWARNINGC("ScatterPlot", "ScatterPlot can only handle PlotData objects");
        data_ = PlotData(0,0);
        dataProp_.set(&data_);
        plotEntitiesProp_.setPlotData(&data_);
        selectionProp_.setPlotData(&data_);
        plotPickingManager_.setColumnCount(data_.getColumnCount());
    }
}

void ScatterPlot::calcDomains() {
    if (plotEntitiesProp_.dataValid()) {
        Interval<plot_t> xDomain = data_.getInterval(plotEntitiesProp_.getXColumnIndex());
        Interval<plot_t> yDomain;
        Interval<plot_t> zDomain;
        if (threeDimensional_) {
            yDomain = data_.getInterval(plotEntitiesProp_.getYColumnIndex());
            //get z domain
            std::vector<PlotEntitySettings>::const_iterator it = plotEntitiesProp_.get().begin();
            for (; it < plotEntitiesProp_.get().end(); ++it)
                zDomain.unionWith(data_.getInterval(it->getMainColumnIndex()));
        }
        else {
            std::vector<PlotEntitySettings>::const_iterator it = plotEntitiesProp_.get().begin();
            for (; it < plotEntitiesProp_.get().end(); ++it)
                yDomain.unionWith(data_.getInterval(it->getMainColumnIndex()));
        }
        selectionProp_.setBaseZoomState(PlotZoomState(xDomain, yDomain, zDomain));
    }
}

void ScatterPlot::toggleProperties() {
    threeDimensional_ = plotEntitiesProp_.getYColumnIndex() != -1;
    marginLeft_.setVisibleFlag(!threeDimensional_);
    marginRight_.setVisibleFlag(!threeDimensional_);
    marginBottom_.setVisibleFlag(!threeDimensional_);
    marginTop_.setVisibleFlag(!threeDimensional_);
    setPropertyGroupVisible("margins", !threeDimensional_);
    cameraHandler_->setVisible(threeDimensional_ && !orthographicCamera_.get());
    cameraHandler_->setEnabled(threeDimensional_ && !orthographicCamera_.get());
    plotCameraHandler_->setVisible(threeDimensional_ && orthographicCamera_.get());
    plotCameraHandler_->setEnabled(threeDimensional_ && orthographicCamera_.get());
    //zLabel_.setVisibleFlag(threeDimensional_);
    renderZHelperLines_.setVisibleFlag(threeDimensional_);
    zScaleStep_.setVisibleFlag(threeDimensional_);
    selectionProp_.setDimension(threeDimensional_);
    orthographicCamera_.setVisibleFlag(threeDimensional_);
    lighting_.setVisibleFlag(threeDimensional_);
    lighting_.set(lighting_.get() && threeDimensional_);
    if (orthographicCamera_.get())
        camera_.set(tgt::Camera(tgt::normalize(camera_.get().getPosition()), tgt::vec3(0,0,0), tgt::vec3(0,0,1)));
}

#ifdef VRN_MODULE_DEVIL
void ScatterPlot::loadTextures(std::string path) {
    // try loading textures
    if (!path.empty())
        plotLib_->setTexture(TexMgr.load(path),path);
    else
        plotLib_->setTexture(0,"");
    LGL_ERROR;
}
#else
void ScatterPlot::loadTextures(std::string /*path*/) {
}

#endif //VRN_MODULE_DEVIL

} // namespace voreen
