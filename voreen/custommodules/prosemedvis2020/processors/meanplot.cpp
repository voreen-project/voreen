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

#include "meanplot.h"

#include "modules/plotting/datastructures/plotrow.h"
#include "voreen/core/utils/multisampler.h"
#include "tgt/tgt_math.h"

#include <algorithm>
#include <math.h>
#include <iomanip>
#include <string>

namespace voreen {

	MeanPlot::MeanPlot()
		: PlotProcessor(PlotEntitySettings::LINE, false)

		, colorArea_("colorArea", "Color area", true)
		, heatMap_("hearMap", "Timestamp (Heatmap)", 0, 0, 200)
		, lineWidth_("lineWidth", "Line Width", 2.f, 1.f, 5.f)
		, pointSize_("pointSize", "Point Size", 2.f, 1.f, 9.f)
		, renderLineLabel_("renderLineLabel", "Line Labels", false)
		, showStdLabels_("showStdLabel", "Standard Deviation (Line) Labels", false)

	{
		removePort(inport_);
		inport_ = PlotPort(Port::INPORT, "inport", "inport port", true, Processor::INVALID_RESULT);

		addPort(inport_);

		plotEntitiesProp_.setGuiName("Line Data");
		addProperty(selectionPlaneColor_);
		addProperty(renderMousePosition_);
		addProperty(discreteStep_);
		addProperty(renderXHelperLines_);
		addProperty(renderYHelperLines_);
		addProperty(renderLineLabel_);
		addProperty(showStdLabels_);
		addProperty(colorArea_);
		addProperty(heatMap_);
		addProperty(lineWidth_);
		addProperty(pointSize_);
		addProperty(xScaleStep_);
		addProperty(yScaleStep_);
		addProperty(marginLeft_);
		addProperty(marginRight_);
		addProperty(marginBottom_);
		addProperty(marginTop_);

		// group properties
		renderLineLabel_.setGroupID("line");
		showStdLabels_.setGroupID("line");
		colorArea_.setGroupID("line");
		lineWidth_.setGroupID("line");
		heatMap_.setGroupID("line");
		pointSize_.setGroupID("line");
		setPropertyGroupGuiName("line", "Line Settings");

		addEventProperty(eventHighlight_);
		addEventProperty(eventLabel_);
		addEventProperty(eventZoom_);
		addEventProperty(eventHighlightAdditive_);
		addEventProperty(eventLabelAdditive_);
		addEventProperty(eventZoomAdditive_);
		addEventProperty(mousePositionUpdateEvent_);
		addEventProperty(mouseEventEnterExit_);

		//if one of the following properties is changed we handle it like plot entities property is changed
		renderLineLabel_.onChange(MemberFunctionCallback<MeanPlot>(this, &MeanPlot::regenDisplayLists));
		lineWidth_.onChange(MemberFunctionCallback<MeanPlot>(this, &MeanPlot::regenDisplayLists));
		heatMap_.onChange(MemberFunctionCallback<MeanPlot>(this, &MeanPlot::regenDisplayLists));
		pointSize_.onChange(MemberFunctionCallback<MeanPlot>(this, &MeanPlot::regenDisplayLists));
	}

	Processor* MeanPlot::create() const {
		return new MeanPlot();
	}

	void MeanPlot::render() {
		{
			Multisampler m(outport_);

			plotLib_->setUsePlotPickingManager(false);
			setPlotStatus();
			if (plotLib_->setRenderStatus()) {
				renderAxes();
				renderData();
				createLineLabels();
				plotLib_->renderLineLabels();
				createPlotLabels();
				plotLib_->renderPlotLabels();
				renderPlotLabel();
				renderMousePosition();
			}
			plotLib_->resetRenderStatus();

			renderSelectedRegion();
		}

		plotPickingManager_.activateTarget();
		plotPickingManager_.clearTarget();
		if (enablePicking_.get()) {
			plotLib_->setUsePlotPickingManager(true);
			renderData();
			plotLib_->resetRenderStatus();
		}
		plotPickingManager_.deactivateTarget();
	}

	void MeanPlot::setPlotStatus() {
		plotLib_->setWindowSize(outport_.getSize());
		plotLib_->setAxesWidth(axesWidth_.get());
		plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
		plotLib_->setLineWidth(lineWidth_.get());
		plotLib_->setMaxGlyphSize(pointSize_.get());
		plotLib_->setMarginBottom(marginBottom_.get());
		plotLib_->setMarginTop(marginTop_.get());
		plotLib_->setMarginLeft(marginLeft_.get());
		plotLib_->setMarginRight(marginRight_.get());
		plotLib_->setMinimumScaleStep(xScaleStep_.get(), PlotLibrary::X_AXIS);
		plotLib_->setMinimumScaleStep(yScaleStep_.get(), PlotLibrary::Y_AXIS);
		plotLib_->setDomain(selectionProp_.getZoom().xZoom_, PlotLibrary::X_AXIS);
		plotLib_->setDomain(selectionProp_.getZoom().yZoom_, PlotLibrary::Y_AXIS);
	}


	MeanPlot::MeanData MeanPlot::readFromInport(const PlotData* p) {

		const int rows = p->getRowsCount();
		const int columns = p->getColumnCount();
		const int data_size = columns - 1;

		MeanData meandata;

		// Set name:
		meandata.name_ = p->getColumnLabel(1);//getTracerName(p);

		meandata.mean_ = new plot_t[rows];
		meandata.std_ = new plot_t[rows];
		std::fill(meandata.mean_, meandata.mean_ + rows, 0);
		std::fill(meandata.std_, meandata.std_ + rows, 0);


		// Calculate Mean
		for (int j = 1; j < columns; j++) {
			for (int i = 0; i < rows; i++) {
				meandata.mean_[i] += p->getRow(i).getValueAt(j);
			}
		}

		for (int i = 0; i < rows; i++) {
			meandata.mean_[i] = meandata.mean_[i] / data_size;
		}

		// Calculate Std (standard deviation)
		for (int j = 1; j < columns; j++) {
			for (int i = 0; i < rows; i++) {

				PlotRowValue value = p->getRow(i);
				meandata.std_[i] += std::pow(value.getValueAt(j) - meandata.mean_[i], 2);
			}
		}

		for (int i = 0; i < rows; i++) {
			meandata.std_[i] = std::sqrt(meandata.std_[i] / data_size);
		}

		meandata.std_sum_ = 0;
		for (int i = 0; i < rows; i++) {
			meandata.std_sum_ += meandata.std_[i];
		}

		return meandata;
	}

	void MeanPlot::readFromInport() {

		// Cast to vector of PlotData*
		std::vector<const PlotData*> data_vec;

		for (const PlotBase* pb : inport_.getAllData()) {
			if (dynamic_cast<const PlotData*>(pb)) {
				if (pb->getColumnCount() != 0
					&& pb->getDataColumnCount() != 0) {
					// Only add plot data if it contains any Data-Columns
					const PlotData* pd = dynamic_cast<const PlotData*>(pb);
					data_vec.push_back(pd);
				}
			}
			else {
				return;
			}
		}

		// create and configure PlotData
		const PlotData* p_first = data_vec.front();
		int amount_columns = static_cast<int>(COLUMNS_PER_TRACER_ * data_vec.size());
		int rowCount = p_first->getRowsCount();

		// adjust heatMap range
		if (rowCount > 0) {
			heatMap_.setMinValue(0);
			heatMap_.setMaxValue(rowCount - 1);
		}

		// create new PlotData
		data_ = PlotData(1, amount_columns);
		data_.setColumnLabel(0, "Time");

		std::vector<MeanData> data_mean;

		for (int i = 0; i < data_vec.size(); i++) {
			data_mean.push_back(readFromInport(data_vec.at(i)));
		}

		// sort values by range between +Std and -Std
		std::sort(data_mean.begin(), data_mean.end(), meandata_sort());

		// Set column names
		for (int i = 0; i < data_mean.size(); ++i) {
			const std::string& tracerName = data_mean.at(i).name_;
			data_.setColumnLabel(i * COLUMNS_PER_TRACER_ + 1, "+Std (" + tracerName + ")");
			data_.setColumnLabel(i * COLUMNS_PER_TRACER_ + 2, "Mean + Std (" + tracerName + ")");
			data_.setColumnLabel(i * COLUMNS_PER_TRACER_ + 3, "Mean (" + tracerName + ")");
			data_.setColumnLabel(i * COLUMNS_PER_TRACER_ + 4, "Mean - Std (" + tracerName + ")");
			data_.setColumnLabel(i * COLUMNS_PER_TRACER_ + 5, "-Std (" + tracerName + ")");
		}

		// insert values (mean and std)
		for (int i = 0; i < rowCount; i++) {

			std::vector<PlotCellValue> vec;
			vec.push_back(PlotCellValue(p_first->getRow(i).getValueAt(0)));

			for (MeanData m : data_mean) {
				vec.push_back(PlotCellValue(m.std_[i]));
				vec.push_back(PlotCellValue(m.mean_[i] + m.std_[i]));
				vec.push_back(PlotCellValue(m.mean_[i]));
				vec.push_back(PlotCellValue(m.mean_[i] - m.std_[i]));
				vec.push_back(PlotCellValue(-1 * m.std_[i]));
			}

			data_.insert(vec);
		}

		for (MeanData m : data_mean) {
			delete[]m.mean_; m.mean_ = nullptr;
			delete[]m.std_; m.std_ = nullptr;
		}

		// adjust label colors to colormap
		ColorMap colorMap_ = plotEntitiesProp_.getColorMap();

		for (int i = 1; i <= amount_columns; i++) {
			int index = (i - 1) / COLUMNS_PER_TRACER_;
			data_.setColumnColorHint(i, colorMap_.getColorAtIndex(index >= colorMap_.getColorCount() ? 0 : index));
		}

		plotEntitiesProp_.setPlotData(&data_);
		dataProp_.set(&data_);
		selectionProp_.setPlotData(&data_);
		plotPickingManager_.setColumnCount(data_.getColumnCount());
		inportHasPlotFunction_ = false;
		discreteStep_.setVisibleFlag(false);
	}

	std::string MeanPlot::getTracerName(const PlotData* data) {

		if (data->getColumnCount() <= 1) {
			return "";
		}
		// Name format: Region::pointId::Tracer_IDmemberid
		std::string fullName = data->getColumnLabel(1);
		// Cut of member-id
		fullName = fullName.substr(0, fullName.find("_ID"));
		const std::string region = fullName.substr(0, fullName.find("::"));
		fullName = fullName.substr(region.length() + 2);
		const std::string tracer = fullName.substr(fullName.find("::") + 2);

		
		return region + "::" + tracer;
	}

	void MeanPlot::renderData() {

		plotLib_->setHighlightColor(highlightColor_.get());
		std::vector<PlotEntitySettings> vec = plotEntitiesProp_.get();

		// iterate through plot entity settings
		for (int i = 0; i < vec.size(); i++) {

			PlotEntitySettings it = vec[i];
			tgt::Color color = getColor(i);

			plotLib_->setDrawingColor(color);
			plotLib_->setHighlightColor(color);
			plotLib_->setFontColor(color);

			color.a = 0.5;
			plotLib_->setFillColor(color);
			plotLib_->setLineStyle(it.getLineStyle());

			int main_column_index = it.getMainColumnIndex();

			// draw line
			if (isMeanColumn(i, it, true))
				plotLib_->renderLine(data_, 0, main_column_index);
			// color areas between mean's line and std's lines
			else if (colorArea_.get())
				plotLib_->renderErrorline(data_, 0, getTracerInLoop(i, true) + 3, main_column_index);
		}

		// draw vertical line to show heatMap's timestamp
		const int heatMapStep = heatMap_.get();
		float heatMapTime = 0;

		if (! inport_.getAllData().empty() 
			&& dynamic_cast<const PlotData*>(inport_.getAllData().front())) {
			const PlotData* p_first = dynamic_cast<const PlotData*>(inport_.getAllData().front());
			heatMapTime = static_cast<float>(p_first->getRow(heatMapStep).getValueAt(0));
		}

		if (heatMapTime != 0) {

			Interval<plot_t> i = selectionProp_.getZoom().yZoom_;
			IMode.color(tgt::Color(0.f, 0.f, 0.f, 1.f));
			IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
			IMode.vertex(heatMapTime, static_cast<float>(i.getLeft()));
			IMode.vertex(heatMapTime, static_cast<float>(i.getRight()));
			IMode.end();
		}

		LGL_ERROR;
	}

	void MeanPlot::renderAxes() {
		if (renderAxes_.get()) {
			plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
			plotLib_->renderAxes();
			plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
			if (renderScales_.get()) {
				plotLib_->setFontSize(10);
				plotLib_->setFontColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
				if (data_.getColumnType(plotEntitiesProp_.getXColumnIndex()) == PlotBase::STRING)
					plotLib_->renderAxisLabelScales(data_, plotEntitiesProp_.getXColumnIndex(), renderXHelperLines_.get());
				else
					plotLib_->renderAxisScales(PlotLibrary::X_AXIS, renderXHelperLines_.get());
				plotLib_->renderAxisScales(PlotLibrary::Y_AXIS, renderYHelperLines_.get());
				plotLib_->setFontSize(12);
				plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, getXLabel());
				plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, getYLabel());
			}
		}
	}


	void MeanPlot::calcDomains() {
		if (plotEntitiesProp_.dataValid() && !inportHasPlotFunction_) {
			Interval<plot_t> xDomain = data_.getInterval(plotEntitiesProp_.getXColumnIndex());
			Interval<plot_t> yDomain = Interval<plot_t>();
			std::vector<PlotEntitySettings>::const_iterator it = plotEntitiesProp_.get().begin();
			for (; it < plotEntitiesProp_.get().end(); ++it) {
				if ((it->getMainColumnIndex() - 1) % 5 == 1
					|| (it->getMainColumnIndex() - 1) % 5 == 2
					|| (it->getMainColumnIndex() - 1) % 5 == 3) {
					// Skip the +Std and -Std columns as they are not plotted
					yDomain.unionWith(data_.getInterval(it->getMainColumnIndex()));
				}
			}
			yDomain.enlarge(1.1);
			selectionProp_.setBaseZoomState(PlotZoomState(xDomain, yDomain));
		}
	}

	tgt::Color MeanPlot::getColor(int index) {

		ColorMap colorMap_ = plotEntitiesProp_.getColorMap();

		int index_color = getTracerInLoop(index, false);
		return colorMap_.getColorAtIndex(index_color >= colorMap_.getColorCount() ? 0 : index_color);
	}

	void MeanPlot::createLineLabels() {
		plotLib_->resetLineLabels();

		if (!renderLineLabel_.get())
			return;

		plot_t xr = selectionProp_.getZoom().xZoom_.getRight();
		// we have to find last rendered x value
		int row = data_.getRowsCount() - 1;
		bool xAxisIsString = (data_.getColumnType(plotEntitiesProp_.getXColumnIndex()) == PlotBase::STRING);
		if (xAxisIsString)
			row = static_cast<int>(floor(xr));
		else {
			while (row > 0 && data_.getRow(row - 1).getValueAt(plotEntitiesProp_.getXColumnIndex()) > xr)
				--row;
		}
		if (row < 1)
			return;

		std::vector<PlotEntitySettings> vec = plotEntitiesProp_.get();

		for (int i = 0; i < vec.size(); i++) {

			PlotEntitySettings it = vec[i];

			if (!isMeanColumn(i, it, showStdLabels_.get()))
				continue;

			tgt::Color color = getColor(i);

			tgt::dvec2 last(xAxisIsString ? static_cast<double>(row) :
				data_.getRow(row).getValueAt(plotEntitiesProp_.getXColumnIndex()),
				data_.getRow(row).getValueAt(it.getMainColumnIndex()));
			tgt::dvec2 lastButOne(xAxisIsString ? static_cast<double>(row) - 1 :
				data_.getRow(row - 1).getValueAt(plotEntitiesProp_.getXColumnIndex()),
				data_.getRow(row - 1).getValueAt(it.getMainColumnIndex()));
			double dydx = (last.y - lastButOne.y) / (last.x - lastButOne.x);
			last.y = dydx * (xr - last.x) + last.y;
			if (tgt::isNaN(last.y) && !tgt::isNaN(lastButOne.y))
				last.y = lastButOne.y;
			if (selectionProp_.getZoom().yZoom_.contains(last.y)) {
				plotLib_->addLineLabel(data_.getColumnLabel(it.getMainColumnIndex()),
					tgt::dvec3(10, 0, 0) + plotLib_->convertPlotCoordinatesToViewport3(
						tgt::dvec3(xr, last.y, 0)),
					color, 12, SmartLabel::MIDDLELEFT);
			}
		}
	}

	void MeanPlot::createPlotLabels() {
		plotLib_->resetPlotLabels();
		PlotSelectionProperty::LabelSelectionIterator lit = selectionProp_.getLabelsBegin();
		if (lit == selectionProp_.getLabelsEnd()) // no labels
			return;
		std::stringstream ss;
		// iterate label selection
		for (; lit != selectionProp_.getLabelsEnd(); ++lit) {
			if (lit->isTablePosition()) {
				tgt::ivec2 cell = lit->getTablePosition();
				if (cell.x >= -1 && cell.x < data_.getRowsCount() && cell.y >= 0 && cell.y < data_.getColumnCount()) {
					int start = 0;
					int end = data_.getRowsCount();
					if (cell.x != -1) {
						start = cell.x;
						end = cell.x + 1;
					}
					for (int i = start; i < end; ++i) {
						// check if line or error is labeled
						std::vector<PlotEntitySettings>::const_iterator eit;
						for (eit = plotEntitiesProp_.get().begin(); eit < plotEntitiesProp_.get().end(); ++eit) {
							if (eit->getMainColumnIndex() == lit->getTablePosition().y
								|| eit->getStickTopColumnIndex() == lit->getTablePosition().y
								|| eit->getCandleBottomColumnIndex() == lit->getTablePosition().y
								|| eit->getCandleTopColumnIndex() == lit->getTablePosition().y
								|| eit->getStickBottomColumnIndex() == lit->getTablePosition().y) {
								const PlotRowValue& row = data_.getRow(i);
								plot_t x = (data_.getColumnType(plotEntitiesProp_.getXColumnIndex()) == PlotBase::STRING ?
									i : row.getValueAt(plotEntitiesProp_.getXColumnIndex()));
								plot_t y = row.getValueAt(lit->getTablePosition().y);
								if (selectionProp_.getZoom().xZoom_.contains(x) && selectionProp_.getZoom().yZoom_.contains(y)) {
									tgt::dvec3 viewportCoords = plotLib_->convertPlotCoordinatesToViewport3(tgt::dvec3(x, y, 0));
									ss.str("");
									ss.clear();
									ss << std::fixed << std::setprecision(4) << "x: " << x << std::endl << "y: " << y;
									plotLib_->addPlotLabel(ss.str(), viewportCoords, tgt::Color(0.f, 0.f, 0.f, 1.f), 10, SmartLabel::CENTERED);
								}
							}
							else if (eit->getOptionalColumnIndex() == lit->getTablePosition().y && eit->getErrorbarFlag()) {
								const PlotRowValue& row = data_.getRow(i);
								plot_t x = row.getValueAt(plotEntitiesProp_.getXColumnIndex());
								plot_t y = row.getValueAt(eit->getMainColumnIndex());
								plot_t error = row.getValueAt(lit->getTablePosition().y);
								if (selectionProp_.getZoom().xZoom_.contains(x) && selectionProp_.getZoom().yZoom_.contains(y)) {
									tgt::dvec3 viewportCoords = plotLib_->convertPlotCoordinatesToViewport3(tgt::dvec3(x, y, 0));
									ss.str("");
									ss.clear();
									ss << std::fixed << std::setprecision(4) << "Error: " << error;
									plotLib_->addPlotLabel(ss.str(), viewportCoords, tgt::Color(0.f, 0.f, 0.f, 1.f), 10, SmartLabel::CENTERED);
								}
							}
						}
					}
				}
			}
		}
	}

} // namespace voreen
