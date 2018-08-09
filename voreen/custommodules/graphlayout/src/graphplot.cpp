/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "../include/graphplot.h"
#include "modules/plotting/datastructures/plotrow.h"

namespace voreen {

const std::string GraphPlot::loggerCat_("voreen.GraphPlot");

GraphPlot::GraphPlot()
    : PlotProcessor(PlotEntitySettings::NONE, false)
    , connectionPort_(Port::INPORT, "connectionPort")
    , connectionData_(2, 0)
{
    plotEntitiesProp_.setGuiName("No Data");
    addPort(connectionPort_);

    addProperty(marginLeft_);
    addProperty(marginRight_);
    addProperty(marginTop_);
    addProperty(marginBottom_);
}

void GraphPlot::process() {
    if (inport_.isReady() && connectionPort_.isReady() && (inport_.hasChanged() || connectionPort_.hasChanged()))
        readFromInport();
    render();
}

void GraphPlot::render() {
    outport_.activateTarget();
    plotLib_->setUsePlotPickingManager(false);
    calcDomains();
    setPlotStatus();
    if (plotLib_->setRenderStatus()) {
        renderAxes();
        renderData();
    }
    plotLib_->resetRenderStatus();
    outport_.deactivateTarget();
/*    plotPickingManager_.activateTarget();
    plotPickingManager_.clearTarget();
    if (enablePicking_.get()) {
        plotLib_.setUsePlotPickingManager(true);
        if (plotLib_.setOpenGLStatus())
            renderData();
        plotLib_.resetOpenGLStatus();
    }
    plotPickingManager_.deactivateTarget();*/
}

void GraphPlot::setPlotStatus() {
    plotLib_->setDimension(PlotLibrary::TWO);
    plotLib_->setWindowSize(outport_.getSize());
    plotLib_->setDrawingColor(tgt::Color(0,0,0,1));
    plotLib_->setMarginBottom(marginBottom_.get());
    plotLib_->setMarginTop(marginTop_.get());
    plotLib_->setMarginLeft(marginLeft_.get());
    plotLib_->setMarginRight(marginRight_.get());
    plotLib_->setAxesWidth(axesWidth_.get());
    calcDomains();
}

void GraphPlot::renderData() {
    plotLib_->setDrawingColor(tgt::Color(1,0,0,1));
    plotLib_->setFillColor(tgt::Color(0,0,1,1));
    plotLib_->setMinGlyphSize(10);
    plotLib_->setMaxGlyphSize(15);

    plotLib_->renderNodeGraph(data_, connectionData_, 1, 2, 3, 4);
    LGL_ERROR;
}

void GraphPlot::renderAxes() {
    if (renderAxes_.get()) {
        plotLib_->setDrawingColor(tgt::Color(0,0,0,1));
        plotLib_->renderAxes();
        plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
        if (renderScales_.get()) {
            plotLib_->setFontSize(10);
            plotLib_->setFontColor(tgt::Color(0,0,0,1));
            plotLib_->renderAxisScales(PlotLibrary::X_AXIS, renderXHelperLines_.get());
            plotLib_->renderAxisScales(PlotLibrary::Y_AXIS, renderYHelperLines_.get());
            plotLib_->setFontSize(12);
            plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, getXLabel());
            plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, getYLabel());
        }
    }}

void GraphPlot::readFromInport() {
    // create local copy of data and assign it to property
    const PlotData* nodeData = dynamic_cast<const PlotData*>(inport_.getData());
    const PlotData* connectionData = dynamic_cast<const PlotData*>(connectionPort_.getData());

    if (nodeData && connectionData) {
        data_ = *nodeData;
        connectionData_ = *connectionData;
        dataProp_.set(&data_);
        plotEntitiesProp_.setPlotData(&data_);
    }
    else {
        LWARNINGC("GraphPlot", "GraphPlot can only handle PlotData objects");
        data_ = PlotData(0,0);
        dataProp_.set(&data_);
        plotEntitiesProp_.setPlotData(&data_);
        selectionProp_.setPlotData(&data_);
        plotPickingManager_.setColumnCount(data_.getColumnCount());
    }
}

void GraphPlot::calcDomains() {
    Interval<plot_t> xInterval = data_.getInterval(1);
    Interval<plot_t> yInterval = data_.getInterval(2);

    if (xInterval.size() > yInterval.size()) {
        plot_t diff = xInterval.size() / yInterval.size();
        yInterval.enlarge(diff);
    }
    else {
        plot_t diff = yInterval.size() / xInterval.size();
        xInterval.enlarge(diff);
    }

    plotLib_->setDomain(xInterval, PlotLibrary::X_AXIS);
    plotLib_->setDomain(yInterval, PlotLibrary::Y_AXIS);
}

void GraphPlot::toggleProperties() {
}

void GraphPlot::createPlotLabels() {
}

} // namespace voreen
