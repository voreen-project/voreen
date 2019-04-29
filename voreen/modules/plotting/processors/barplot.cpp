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

#include "barplot.h"

#include "../datastructures/plotrow.h"
#include "voreen/core/utils/multisampler.h"

#include <iomanip>

namespace voreen {

const std::string BarPlot::loggerCat_("voreen.BarPlot");

BarPlot::BarPlot()
    : PlotProcessor(PlotEntitySettings::BAR, false)
    , barWidth_("barWidth", "Bar Width", 70, 1, 100)
    , shear_("shear", "Shear", tgt::vec2(35,70), tgt::vec2(0,0), tgt::vec2(100,100))
    , squeezeFactor_("squeezeFactor", "Squeeze Factor", 0.9,0.1,0.95)
    , barMode_("barMode", "Bar Mode")
{
    plotEntitiesProp_.setGuiName("Bar Data");
    barMode_.addOption("groupedBar", "Grouped Bar", PlotLibrary::GROUPED);
    barMode_.addOption("stackedBar", "Stacked Bar", PlotLibrary::STACKED);
    barMode_.addOption("mergedBar", "Merged Bar", PlotLibrary::MERGED);
    addProperty(barMode_);
    selectionProp_.setEnableZoomToFlags(false);
    addProperty(barWidth_);
    addProperty(shear_);
    addProperty(squeezeFactor_);
    squeezeFactor_.setVisibleFlag(false);
    addProperty(renderYHelperLines_);
    addProperty(yScaleStep_);
    addProperty(marginLeft_);
    addProperty(marginRight_);
    addProperty(marginBottom_);
    addProperty(marginTop_);

    // group properties
    barMode_.setGroupID("bar");
    barWidth_.setGroupID("bar");
    shear_.setGroupID("bar");
    squeezeFactor_.setGroupID("bar");
    setPropertyGroupGuiName("bar", "Bar Settings");

    addEventProperty(eventHighlight_);
    addEventProperty(eventLabel_);
    addEventProperty(eventHighlightAdditive_);
    addEventProperty(eventLabelAdditive_);

    barMode_.onChange(MemberFunctionCallback<BarPlot>(this, &BarPlot::toggleProperties));
    barMode_.onChange(MemberFunctionCallback<BarPlot>(this, &BarPlot::calcDomains));
}

Processor* BarPlot::create() const {
    return new BarPlot();
}

void BarPlot::select(tgt::MouseEvent* e, bool highlight, bool label, bool zoom, bool additive) {
    mousePosition_ = tgt::ivec2(e->x(),e->viewport().y-e->y());
    if (e->action() == tgt::MouseEvent::PRESSED) {
        mousePressed_ = true;
        mousePressedPosition_ = mousePosition_;
        e->accept();
    }
    else if (e->action() == tgt::MouseEvent::RELEASED) {
        mousePressed_ = false;
        if (tgt::length(mousePosition_ - mousePressedPosition_) < 4 && pickingBuffer_.hasRenderTarget()) {
            tgt::ivec2 cell = plotPickingManager_.getCellAtPos(tgt::ivec2(e->x(), e->viewport().y-e->y()));
            if (cell != tgt::ivec2(-1,-1)) {
                PlotSelectionEntry entry(PlotSelection(cell), highlight, label, zoom);
                selectionProp_.add(entry, additive);
            }
            else {
                selectionProp_.clearFlags(highlight, label, zoom);
            }
        }
        e->accept();
    }
    invalidate();
}

void BarPlot::render() {
    {
        Multisampler m(outport_);
        plotLib_->setUsePlotPickingManager(false);
        setPlotStatus();
        if (plotLib_->setRenderStatus()) {
            renderAxes();
            renderData();
            createPlotLabels();
            plotLib_->renderPlotLabels();
            renderPlotLabel();
        }
        plotLib_->resetRenderStatus();
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

void BarPlot::setPlotStatus() {
    plotLib_->setDimension(PlotLibrary::FAKETHREE);
    plotLib_->setWindowSize(outport_.getSize());
    plotLib_->setBarGroupingMode(barMode_.getValue());
    plotLib_->setBarWidth(static_cast<double>(barWidth_.get())/100.0);
    plotLib_->setAxesWidth(axesWidth_.get());
    plotLib_->setLineWidth(0.4f);
    plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
    plotLib_->setHighlightColor(highlightColor_.get());
    plotLib_->setShear(tgt::vec2(shear_.get().x/1000.f, shear_.get().y/1000.f));
    plotLib_->setSqueezeFactor(static_cast<double>(squeezeFactor_.get()));
    plotLib_->setMarginBottom(marginBottom_.get());
    plotLib_->setMarginTop(marginTop_.get());
    plotLib_->setMarginLeft(marginLeft_.get());
    plotLib_->setMarginRight(marginRight_.get());
    plotLib_->setMinimumScaleStep(yScaleStep_.get(), PlotLibrary::Y_AXIS);
    calcDomains();
}

void BarPlot::renderData() {
    std::vector<int> columns;
    std::vector<PlotEntitySettings>::const_iterator it;
    std::vector<tgt::Color> colors;
    for (it = plotEntitiesProp_.get().begin(); it < plotEntitiesProp_.get().end(); ++it) {
        columns.push_back(it->getMainColumnIndex());
        colors.push_back(it->getFirstColor());
    }
    plotLib_->setColorMap(ColorMap::createFromVector(colors));
    plotLib_->renderBars(data_, columns);
    LGL_ERROR;
}

void BarPlot::renderAxes() {
    if (renderAxes_.get()) {
        plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
        plotLib_->renderAxes();
        plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
        if (renderScales_.get()) {
            plotLib_->setFontSize(10);
            plotLib_->setFontColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
            plotLib_->renderAxisLabelScales(data_,plotEntitiesProp_.getXColumnIndex(), false);
            plotLib_->renderAxisScales(PlotLibrary::Y_AXIS, renderYHelperLines_.get());
            plotLib_->setFontSize(12);
            plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, getXLabel());
            plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, getYLabel());
        }
    }
}

void BarPlot::readFromInport() {
    // create local copy of data and assign it to property
    if (dynamic_cast<const PlotData*>(inport_.getData())) {
        data_ = *dynamic_cast<const PlotData*>(inport_.getData());
        if (!data_.sorted())
            data_.sortRows();
        dataProp_.set(&data_);
        plotEntitiesProp_.setPlotData(&data_);
        selectionProp_.setPlotData(&data_);
        plotPickingManager_.setColumnCount(data_.getColumnCount());
    }
    else {
        LWARNINGC("BarPlot", "BarPlot can only handle PlotData objects");
        data_ = PlotData(0,0);
        dataProp_.set(&data_);
        plotEntitiesProp_.setPlotData(&data_);
        selectionProp_.setPlotData(&data_);
        plotPickingManager_.setColumnCount(data_.getColumnCount());
    }
}

void BarPlot::calcDomains() {
    if (plotEntitiesProp_.dataValid()) {
        //the x domain is simply [0,number of rows-1]
        Interval<plot_t> xDomain = Interval<plot_t>(0,data_.getRowsCount()-1,false,false);
        //get y domain
        Interval<plot_t> yDomain;
        std::vector<PlotEntitySettings>::const_iterator it = plotEntitiesProp_.get().begin();
        if (barMode_.getValue() == PlotLibrary::GROUPED || barMode_.getValue() == PlotLibrary::MERGED) {
            for (; it < plotEntitiesProp_.get().end(); ++it) {
                Interval<plot_t> lineDomain = data_.getInterval(it->getMainColumnIndex());
                yDomain.unionWith(lineDomain);
            }
        }
        else {
            std::vector<int> columns;
            for (; it < plotEntitiesProp_.get().end(); ++it)
                columns.push_back(it->getMainColumnIndex());
            yDomain = data_.getSumInterval(columns);
        }
        xDomain = Interval<plot_t>(xDomain.getLeft()-0.5, xDomain.getRight()+0.5, false, false);
        yDomain = Interval<plot_t>(std::min(0.0, yDomain.getLeft()), std::max(0.0, yDomain.getRight()), true, true);
        plotLib_->setDomain(xDomain, PlotLibrary::X_AXIS);
        plotLib_->setDomain(yDomain, PlotLibrary::Y_AXIS);
    }
}

void BarPlot::toggleProperties() {
    squeezeFactor_.setVisibleFlag((barMode_.getValue() == PlotLibrary::MERGED));
}

void BarPlot::createPlotLabels() {
    plotLib_->resetPlotLabels();
    PlotSelectionProperty::LabelSelectionIterator lit = selectionProp_.getLabelsBegin();
    if (lit == selectionProp_.getLabelsEnd() || plotEntitiesProp_.get().empty()) // no labels
        return;
    std::stringstream ss;
    double singleBarWidth = barWidth_.get()/(100.0*static_cast<double>(plotEntitiesProp_.get().size()));
    for (; lit != selectionProp_.getLabelsEnd(); ++lit) {
        if ((*lit).isTablePosition()) {
            tgt::ivec2 cell = (*lit).getTablePosition();
            size_t currentEntityIndex = 0; // the index of the column in PlotEntitiesProperty::value_
            while (currentEntityIndex < plotEntitiesProp_.get().size() && cell.y != plotEntitiesProp_.get().at(currentEntityIndex).getMainColumnIndex())
                ++currentEntityIndex;
            if (currentEntityIndex < plotEntitiesProp_.get().size()) { // the cell is rendered
                if (cell.x >= 0 && cell.x < data_.getRowsCount() && cell.y >= 0 && cell.y < data_.getColumnCount()) {
                    plot_t value = data_.getRow(cell.x).getValueAt(cell.y);
                    tgt::dvec3 viewportCoords;
                    if (barMode_.getValue() == PlotLibrary::GROUPED) {
                        viewportCoords = plotLib_->convertPlotCoordinatesToViewport3(tgt::dvec3(
                                cell.x-barWidth_.get()/200.0+static_cast<double>(currentEntityIndex+0.5)*singleBarWidth, value/2, 1.0));
                    }
                    else if (barMode_.getValue() == PlotLibrary::MERGED) {
                        viewportCoords = plotLib_->convertPlotCoordinatesToViewport3(tgt::dvec3(cell.x, value, 1.0));
                    }
                    else { //STACKED
                        plot_t yPos = value / 2;
                        for (size_t i = 0; i < currentEntityIndex; ++i) {
                            yPos += data_.getRow(cell.x).getValueAt(plotEntitiesProp_.get().at(i).getMainColumnIndex());
                        }
                        viewportCoords = plotLib_->convertPlotCoordinatesToViewport3(tgt::dvec3(cell.x, yPos, 1.0));
                    }
                    ss.str("");
                    ss.clear();
                    ss << std::fixed << std::setprecision(2) << value;
                    plotLib_->addPlotLabel(ss.str(), viewportCoords, tgt::Color(0.f, 0.f, 0.f, 1.f), 10, SmartLabel::CENTERED);
                }
            }
        }
    }
}

} // namespace voreen
