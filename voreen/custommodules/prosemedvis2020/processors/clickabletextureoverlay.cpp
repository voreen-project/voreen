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

//header file
#include "clickabletextureoverlay.h"
//needed headers (used in process())
#include "tgt/textureunit.h"
#include "tgt/filesystem.h"
#include "tgt/texturereaderdevil.h"
#include "../modules/ensembleanalysis/utils/utils.h"
#include "../modules/plotting/datastructures/plotdata.h"
#include "../modules/plotting/datastructures/plotbase.h"
#include "../modules/plotting/datastructures/plotrow.h"
#include "../modules/plotting/datastructures/colormap.h"
#include "tgt/immediatemode/immediatemode.h"

#include <set>
#include <float.h>
#include <math.h> 
#include <stdexcept>
#include <exception>
#include <iomanip>
#include <sstream>


namespace voreen {

	ClickableTextureOverlay::ClickableTextureOverlay()
		: RenderProcessor()
		, inport_(Port::INPORT, "inport", "Unmodified Image")
		, outport_(Port::OUTPORT, "outport", "Modified Image")
		, concentrationOutport_(Port::OUTPORT, "concentrationOutport", "TracerConcentration")
		, plotDataInport_(Port::INPORT, "plotDataInport", "PlotData")
		, pickingBuffer_(Port::OUTPORT, "picking", "Picking", false)
		, fileProp_("filesoure", "RegionTextures", "Select File", VoreenApplication::app()->getUserDataPath(),
			"CSV (*.csv);;Text (*.txt)", FileDialogProperty::OPEN_FILE)
		, selRegionProp_("region_prop", "Selected region")
		, selTimestepProp_("selTimestep_prop", "Timestep")
		, selTracerProp_("selTracerProp", "Tracer")
		, colorMapProp_("colorMapProp", "Color Map")
		, heatMapFontProp_("heatMapFontProp", "Font")
		, selectedTimeProp_("selTime_prop", "Time")
		, timeUnitProp_("timeUnit_prop", "Time-Unit")
		, measurementUnitProp_("measurementUnit_prop", "Measurement-Unit")
		, regions_()
		, manualUpperBoundProp_("manualUpperBound_prop", "Limit max. value", false)
		, upperBoundProp_("upperBound_prop", "Max. value")
{	
	
	//register ports
	addPort(inport_);
	addPort(outport_);
	addPort(concentrationOutport_);
	addPrivateRenderPort(pickingBuffer_);
	addPort(plotDataInport_);
	//register properties
	addProperty(fileProp_);
	ON_CHANGE(fileProp_, ClickableTextureOverlay, directoryChanged);
	addProperty(selRegionProp_);
	addProperty(selTimestepProp_);
	addProperty(selTracerProp_);
	addProperty(colorMapProp_);
	addProperty(heatMapFontProp_);
	mouseDirectory = "";
	selRegionProp_.setEditable(false);
	addProperty(selectedTimeProp_);
	selectedTimeProp_.setEditable(false);
	addProperty(timeUnitProp_);
	addProperty(measurementUnitProp_);
	addProperty(manualUpperBoundProp_);
	addProperty(upperBoundProp_);
}


void ClickableTextureOverlay::initialize() {
	RenderProcessor::initialize();

	// load fragment shader 'clickabletextureoverlay.frag'
	shader_ = ShdrMgr.loadSeparate("passthrough.vert", "clickabletextureoverlay.frag",
		generateHeader(), // see RenderProcessor
		false);           // do not set default flags
}

void ClickableTextureOverlay::deinitialize() {
	// free shader
	ShdrMgr.dispose(shader_);
	shader_ = 0;

	RenderProcessor::deinitialize();
}

void ClickableTextureOverlay::process() {

	// Check if plotdata changed:
	if (plotDataInport_.hasChanged()) {
		plotDataChanged();
	}

	// Resize frame buffer accordingly.
	if (pickingBuffer_.getSize() != outport_.getSize()) {
		pickingBuffer_.resize(outport_.getSize());
	}

	// Update the displayed time-value
	if (plotDataInport_.hasData()
		&& dynamic_cast<const PlotData*>(plotDataInport_.getData()) != nullptr) {

		const PlotData* plotData = dynamic_cast<const PlotData*>(plotDataInport_.getData());	
		const int selStep = selTimestepProp_.get();
		std::string newValue = "";
		if (selStep >= 0 && selStep < plotData->getRowsCount()) {
			std::stringstream stream;
			stream.precision(2);
			stream << std::noshowpoint << plotData->getRow(selStep).getValueAt(0);
			stream << " " << timeUnitProp_.get();
			newValue = stream.str();
		}
		selectedTimeProp_.set(newValue);
	}


	const ColorMap& cMap = colorMapProp_.get();
	
		// activate and clear output render target
		outport_.activateTarget();
		outport_.clearTarget();
		shader_->activate();
		setGlobalShaderParameters(shader_);
		shader_->setUniform("color_", tgt::Vector3f(1.0, 1.0, 1.0));
		inport_.getRenderTarget()->getColorTexture()->enable();
		inport_.getRenderTarget()->getColorTexture()->bind();
		renderQuad();
		inport_.getRenderTarget()->getColorTexture()->disable();
		if (!regions_.empty() && mouseDirectory != "") {
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			std::string currentTexturePath = "";
			for (SamplePointConfig& region : regions_) {
				// Get color for region:
				try {
					const plot_t maxVal = manualUpperBoundProp_.get() ? upperBoundProp_.get() : maxMeanBuffer_;
					const plot_t meanVal = meanBuffer_.at(selTracerProp_.get()).at(region.name_).at(selTimestepProp_.get());
					const plot_t meanValNorm = meanVal / maxVal;
					const tgt::Color regionColor = cMap.getColorAtPosition(static_cast<float>(meanValNorm));
					shader_->setUniform("color_", regionColor.xyz());
				} catch (std::out_of_range& e){
					// Unable to get color from buffer --> use default-color
					shader_->setUniform("color_", tgt::Color(0, 0, 0, 1).xyz());
				}

				currentTexturePath = mouseDirectory + "/" + region.textureFilename_;
                tgt::Texture* texture = TexMgr.load(currentTexturePath);
				texture->setFilter(tgt::Texture::Filter::NEAREST);
				texture->enable();
				texture->bind();
				renderQuad();
				texture->disable();
			}
			glDisable(GL_BLEND);
			glBlendFunc(GL_ONE, GL_ZERO);
		}
		// cleanup
		shader_->deactivate();
		outport_.deactivateTarget();

		if (concentrationOutport_.hasRenderTarget()) {
			drawConcentration();
		}
}

void ClickableTextureOverlay::directoryChanged() {
	
	selRegionProp_.set("");
	
	regions_ = SamplePointConfigLoader().loadSamplepointConfigFile(fileProp_.get());

	// Cut of filename from path
	size_t pos = fileProp_.get().find_last_of("\\/");
	mouseDirectory = (std::string::npos == pos) ? "" : fileProp_.get().substr(0, pos);
} 

void ClickableTextureOverlay::mouseEvent(tgt::MouseEvent* e) {
	e->accept();
	
	//Left Click, selects a region using a pickingBuffer 
	if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT && e->action() == tgt::MouseEvent::PRESSED) {
		//get coordinates and map them 
		int x = e->x();
		int y = e->y();
		RenderTarget* targetOutport = outport_.getRenderTarget();
		targetOutport->getColorTexture()->setFilter(tgt::Texture::Filter::NEAREST);
		tgt::ivec2 pixelOutport = mapRange(tgt::ivec2(x,y), tgt::ivec2::zero, e->viewport(), tgt::ivec2::zero, targetOutport->getSize());
		tgt::vec4 texelOutport = targetOutport->getColorAtPos(tgt::ivec2(pixelOutport.x, targetOutport->getSize().y - pixelOutport.y - 1));
		std::string currentTexturePath = "";
		tgt::Texture* texture = NULL;
		//renders every region onto the otherwise black pickingBuffer
		for (SamplePointConfig& region : regions_) {
			pickingBuffer_.activateTarget();
			pickingBuffer_.clearTarget();
			currentTexturePath = mouseDirectory + "/" + region.textureFilename_;
			texture = TexMgr.load(currentTexturePath);
			texture->setFilter(tgt::Texture::Filter::NEAREST);
			texture->enable();
			texture->bind();
			renderQuad();
            texture->disable();
			RenderTarget* targetPicking = pickingBuffer_.getRenderTarget();
			targetPicking->getColorTexture()->setFilter(tgt::Texture::Filter::NEAREST);
			//get coordinates and map them
			tgt::ivec2 pixelPicking = mapRange(tgt::ivec2(x, y), tgt::ivec2::zero, e->viewport(), tgt::ivec2::zero, targetPicking->getSize());
			tgt::vec4 texelPicking = targetPicking->getColorAtPos(tgt::ivec2(pixelPicking.x, targetPicking->getSize().y - pixelPicking.y - 1));
			//if clicked region is not black, its a hit
			if (texelPicking.r != 0 || texelPicking.g != 0 || texelPicking.b != 0) {
				selRegionProp_.set(region.name_);
			}
			pickingBuffer_.deactivateTarget();
		}	
	}

	//change timesteps when mousewheel is used
	if (e->button() == tgt::MouseEvent::MOUSE_WHEEL_UP && e->action() == tgt::MouseEvent::WHEEL) {
		selTimestepProp_.set((selTimestepProp_.get() + 1) % (getNumTimesteps() - 1));
	}

	if (e->button() == tgt::MouseEvent::MOUSE_WHEEL_DOWN && e->action() == tgt::MouseEvent::WHEEL) {
		selTimestepProp_.set(selTimestepProp_.get() == 0 ? getNumTimesteps() - 1 : selTimestepProp_.get() - 1);
	}
	
}


void ClickableTextureOverlay::onEvent(tgt::Event* e) {
	tgt::MouseEvent* event = dynamic_cast<tgt::MouseEvent*>(e);
	if (event) {
		mouseEvent(event);
	}

	RenderProcessor::onEvent(e);
}

std::vector<std::string> voreen::ClickableTextureOverlay::getTracerNames() const {
	if (!plotDataInport_.hasData()) {
		return std::vector<std::string>();
	}
	
	const PlotBase* plotData = plotDataInport_.getData();
	std::set<std::string> names;

	for (int i = plotData->getKeyColumnCount(); i < plotData->getColumnCount(); ++i) {
		// Get tracername from columname (REGION::ID::TRACERNAME_MOUSEID)
		std::string colLabel = plotData->getColumnLabel(i);

		// Get begin pos (Pos after second occurence of ::)
		size_t beginPos = colLabel.find("::");
		if (beginPos == std::string::npos || beginPos + 2 >= colLabel.size()) {
			continue;
		}
		beginPos = colLabel.find("::", beginPos + 2);
		if (beginPos == std::string::npos || beginPos + 2 >= colLabel.size()) {
			continue;
		}
		beginPos += 2;

		// Get end-Pos (Position before last occurence of '_')
		size_t endPos = colLabel.find_last_of("_");
		if (endPos == std::string::npos || endPos <= beginPos) {
			continue;
		}

		std::string tracerName = colLabel.substr(beginPos, endPos - beginPos);
		names.insert(tracerName);
	}

	return std::vector<std::string>(names.begin(), names.end());
}

std::pair<plot_t, plot_t> voreen::ClickableTextureOverlay::getMeanValue(const int timestepNr, const std::string& tracerName, const std::string& regionName) const {
	if (!plotDataInport_.hasData() 
		|| dynamic_cast<const PlotData*>(plotDataInport_.getData()) == nullptr) {
		return { DBL_MAX, -1 };
	}

	const PlotData* plotData = dynamic_cast<const PlotData*>(plotDataInport_.getData());

	// Find row
	if (timestepNr < 0 || timestepNr >= plotData->getRowsCount()) {
		return { DBL_MAX, -1 };
	}
	const PlotRowValue& row = plotData->getRow(timestepNr);

	// Find columns and calc mean:
	plot_t mean = 0;
	size_t count = 0;
	std::vector<int> columnIndices;
	for (int c = plotData->getKeyColumnCount(); c < plotData->getColumnCount(); ++c) {
		const std::string& colName = plotData->getColumnLabel(c);
		// Check if column-name strart with Region::
		if (!colName.rfind(regionName + "::", 0) == 0) {
			continue;
		}
		// Check if column contains Tracername
		if (colName.find("::" + tracerName + "_", regionName.size() + 2) == std::string::npos) {
			continue;
		}

		// Check that cell is valid
		if (row.getCellAt(c).isNull() || row.getCellAt(c).isTag()) {
			// TODO: Ignore cell or exit with error?
			LWARNING("Skipping cell " + colName + "::" + std::to_string(c) + " during mean calculation!");
			continue;
		}
		mean += row.getCellAt(c).getValue();
		// Save column-index for std-calculation
		columnIndices.push_back(c);
		++count;
	}

	if (count == 0) {
		return { DBL_MAX, -1 };
	}

	mean = mean / static_cast<plot_t>(count);

	// Calculate std:
	plot_t std = 0;
	for (int c : columnIndices) {
		const plot_t diff = row.getCellAt(c).getValue() - mean;
		std += diff * diff;
	}
	std = std / static_cast<plot_t>(columnIndices.size());
	std = sqrt(std);

	return { mean, std };
}

bool voreen::ClickableTextureOverlay::isReady() const {
	return inport_.isReady();
}

void voreen::ClickableTextureOverlay::plotDataChanged() {
	if (!plotDataInport_.hasData()) {
		// No PlotData --> Clear properties
		selTracerProp_.reset();
		selTracerProp_.setOptions(std::deque<voreen::Option<std::string>>());
		return;
	}

	// Update the Tracer-selection-property
	std::vector<std::string> tracerNames = getTracerNames();
	std::string oldSelection = selTracerProp_.get();
	selTracerProp_.reset();
	selTracerProp_.setOptions(std::deque<voreen::Option<std::string>>());
	for (std::string& tracer : tracerNames) {
		selTracerProp_.addOption(tracer, tracer, tracer);
	}
	// Try to reconstruct the old selection
	if (selTracerProp_.hasKey(oldSelection)) {
		selTracerProp_.select(oldSelection);
	}

	// Set range of selTimestepProp_:
	selTimestepProp_.setMinValue(0);
	selTimestepProp_.setMaxValue(getNumTimesteps() - 1);

	// Calculate and buffer mean-values and biggest mean-value
	maxMeanBuffer_ = -DBL_MAX;
	minMeanBuffer_ = DBL_MAX;
	maxStdBuffer_ = -DBL_MAX;
	minStdBuffer_ = DBL_MAX;
	const int numTimesetps = getNumTimesteps();
	meanBuffer_ = std::map<std::string, std::map<std::string, std::vector<plot_t>>>();
	stdBuffer_ = std::map<std::string, std::map<std::string, std::vector<plot_t>>>();
	for (std::string tracer : tracerNames) {
		meanBuffer_.insert({ tracer, std::map<std::string, std::vector<plot_t>>() });
		stdBuffer_.insert({ tracer, std::map<std::string, std::vector<plot_t>>() });
		for (SamplePointConfig region : regions_) {
			meanBuffer_.at(tracer).insert({region.name_, std::vector<plot_t>()});
			stdBuffer_.at(tracer).insert({ region.name_, std::vector<plot_t>() });
			for (int t = 0; t < numTimesetps; ++t) {
				const std::pair<plot_t, plot_t> res = getMeanValue(t, tracer, region.name_);
				const plot_t mean = std::get<0>(res);
				const plot_t std = std::get<1>(res);
				if (mean == DBL_MAX && std == -1) {
					// No values available --> skip
					continue;
				}
				meanBuffer_.at(tracer).at(region.name_).push_back(mean);
				stdBuffer_.at(tracer).at(region.name_).push_back(std);
				if (mean > maxMeanBuffer_) {
					maxMeanBuffer_ = mean;
				}
				if (mean < minMeanBuffer_) {
					minMeanBuffer_ = mean;
				}
				if (std > maxStdBuffer_) {
					maxStdBuffer_ = std;
				}
				if (std <minStdBuffer_) {
					minStdBuffer_ = std;
				}
			}
		}
	}

	upperBoundProp_.setMinValue(minMeanBuffer_);
	upperBoundProp_.setMaxValue(maxMeanBuffer_);
}

int voreen::ClickableTextureOverlay::getNumTimesteps() const {
	if (!plotDataInport_.hasData()
		|| dynamic_cast<const PlotData*>(plotDataInport_.getData()) == nullptr) {
		return -1;
	}

	const PlotData* plotData = dynamic_cast<const PlotData*>(plotDataInport_.getData());
	return plotData->getRowsCount();
}

void ClickableTextureOverlay::drawConcentration() {

	concentrationOutport_.activateTarget();
	concentrationOutport_.clearTarget();

	std::vector<std::string> tracers = getTracerNames();

	if (regions_.size() == 0 || tracers.size() == 0) {
		concentrationOutport_.deactivateTarget();
		return;
	}

	if (concentrationOutport_.getSize() != outport_.getSize()) {
		concentrationOutport_.resize(outport_.getSize());
	}

	const ColorMap& cMap = colorMapProp_.get();

	//calculate offsets and distances between squares based on windowsize
	int numberOfRegions = regions_.size();
	int numberOfTracers = tracers.size();
	int height = concentrationOutport_.getSize().x;
	int width = concentrationOutport_.getSize().y;
	int xInterval = height / (numberOfRegions + 2);
	int yInterval = width / (numberOfTracers + 4);
	int xOffset = 0;
	int yOffset = height / (4 * numberOfTracers);
	int stdInterval = yInterval / 5;


	// Render legend for mean-ColorMap
	const float legendWidth = 0.6 * concentrationOutport_.getSize().x;
	const float legendHeight = yInterval * 0.5;
	const float legendXOffset = 0.2 * concentrationOutport_.getSize().x;
	const float legendYOffset = 0.25 * yInterval;
	const plot_t maxVal = manualUpperBoundProp_.get() ? upperBoundProp_.get() : maxMeanBuffer_;
	renderLegend(cMap, legendWidth, legendHeight, legendXOffset, legendYOffset, minMeanBuffer_, maxVal, 
		"Mean (" + measurementUnitProp_.get() + ")");
	
	heatMapFontProp_.get()->setTextAlignment(tgt::Font::BottomLeft);
	heatMapFontProp_.get()->setLineWidth(xInterval);

	// Render tracer-names
	for (int r = 0; r < tracers.size(); ++r) {
		const std::string& tracerName = tracers.at(r);
		const tgt::vec3 pos(xOffset, (r + 1) * (yInterval + yOffset) + 0.5 * yInterval, 0);
		heatMapFontProp_.get()->render(pos, tracerName, concentrationOutport_.getSize());
	}

	for (int xDir = 1; xDir < (numberOfRegions + 1) ; xDir++) {
		// Render region-name
		const std::string& region = regions_.at(xDir - 1).name_;
		const tgt::vec3 pos(xDir * (xInterval + xOffset), (numberOfTracers + 1) * (yInterval + yOffset), 0);
		heatMapFontProp_.get()->render(pos, region, concentrationOutport_.getSize());
		const plot_t maxVal = manualUpperBoundProp_.get() ? upperBoundProp_.get() : maxMeanBuffer_;

		for (int yDir = 1; yDir < (numberOfTracers + 1) ; yDir++) {
			// Get color
			std::string& tracer = tracers.at(yDir - 1);
			tgt::Color meanSquareColor = tgt::vec4(1.0, 1.0, 1.0, 1.0);
			tgt::Color maxStdSquareColor = tgt::vec4(1.0, 1.0, 1.0, 1.0);
			tgt::Color minStdSquareColor = tgt::vec4(1.0, 1.0, 1.0, 1.0);
			try {
				//calculate color for the meanSquare
				const plot_t meanVal = meanBuffer_.at(tracer).at(region).at(selTimestepProp_.get());
				const plot_t meanValNorm = meanVal / maxVal;
				meanSquareColor = cMap.getColorAtPosition(static_cast<float>(meanValNorm));

				//calculate color for the maxStdSquare
				const plot_t maxStdVal = meanVal + stdBuffer_.at(tracer).at(region).at(selTimestepProp_.get());
				const plot_t maxStdValNorm = maxStdVal / maxVal;
				maxStdSquareColor = cMap.getColorAtPosition(static_cast<float>(maxStdValNorm));

				//calculate color for the minStdSquare
				const plot_t minStdVal = meanVal - stdBuffer_.at(tracer).at(region).at(selTimestepProp_.get());
				const plot_t minStdValNorm = minStdVal / maxVal;
				minStdSquareColor = cMap.getColorAtPosition(static_cast<float>(minStdValNorm));





			} catch (std::out_of_range& e) {
				// Unable to get mean from buffer --> skip (or use default-color??)
				continue;
			}

			//calculate coordinates for corners of squares
			tgt::vec2 A = tgt::vec2(xDir * xInterval + xDir * xOffset, yDir * yInterval + yDir * yOffset);
			tgt::vec2 B = tgt::vec2((xDir + 1) * xInterval + xDir * xOffset, yDir * yInterval + yDir * yOffset);
			tgt::vec2 C = tgt::vec2((xDir + 1) * xInterval + xDir * xOffset, yDir * yInterval + yDir * yOffset + stdInterval);
			tgt::vec2 D = tgt::vec2((xDir + 1)* xInterval + xDir * xOffset, yDir * yInterval + yDir * yOffset + 4 * stdInterval);
			tgt::vec2 E = tgt::vec2((xDir + 1)* xInterval + xDir * xOffset, (yDir + 1) * yInterval + yDir * yOffset);
			tgt::vec2 F = tgt::vec2(xDir * xInterval + xDir * xOffset, (yDir + 1) * yInterval + yDir * yOffset);
			tgt::vec2 G = tgt::vec2(xDir * xInterval + xDir * xOffset, yDir * yInterval + yDir * yOffset + 4 * stdInterval);
			tgt::vec2 H = tgt::vec2(xDir * xInterval + xDir * xOffset, yDir * yInterval + yDir * yOffset + stdInterval);

			tgt::vec2 outportSize = tgt::vec2(concentrationOutport_.getSize());

			//draw maxStdSquare
			IMode.color(maxStdSquareColor);
			IMode.begin(tgt::ImmediateMode::QUADS);
			IMode.vertex(mapRange(A, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.vertex(mapRange(B, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.vertex(mapRange(C, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.vertex(mapRange(H, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.end();

			//draw meanSquare
			IMode.color(meanSquareColor);
			IMode.begin(tgt::ImmediateMode::QUADS);
			IMode.vertex(mapRange(H, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.vertex(mapRange(C, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.vertex(mapRange(D, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.vertex(mapRange(G, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.end();
			
			//draw minStdSquare
			IMode.color(minStdSquareColor);
			IMode.begin(tgt::ImmediateMode::QUADS);
			IMode.vertex(mapRange(G, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.vertex(mapRange(D, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.vertex(mapRange(E, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.vertex(mapRange(F, tgt::vec2::zero, outportSize, -tgt::vec2::one, tgt::vec2::one));
			IMode.end();
		}
	}
	IMode.color(tgt::vec4::one);
	concentrationOutport_.deactivateTarget();
}

void voreen::ClickableTextureOverlay::renderLegend(const ColorMap& cMap, const float width, const float height, 
	const float xOffset, const float yOffset, const plot_t minValue, const plot_t maxValue, const std::string& label) {
	const float segmentWidth = width / static_cast<float>(cMap.getColorCount() - 1);
	// Draw gradient
	IMode.begin(tgt::ImmediateMode::QUADS);
	for (int i = 0; i < cMap.getColorCount() - 1; ++i) {
		IMode.color(cMap.getColorAtIndex(i));
		IMode.vertex(mapRange(tgt::vec2(xOffset + segmentWidth * i, yOffset + 0), tgt::vec2::zero, tgt::vec2(concentrationOutport_.getSize()), -tgt::vec2::one, tgt::vec2::one));
		IMode.vertex(mapRange(tgt::vec2(xOffset + segmentWidth * i, yOffset + height), tgt::vec2::zero, tgt::vec2(concentrationOutport_.getSize()), -tgt::vec2::one, tgt::vec2::one));
		IMode.color(cMap.getColorAtIndex(i + 1));
		IMode.vertex(mapRange(tgt::vec2(xOffset + segmentWidth * (i + 1), yOffset + height), tgt::vec2::zero, tgt::vec2(concentrationOutport_.getSize()), -tgt::vec2::one, tgt::vec2::one));
		IMode.vertex(mapRange(tgt::vec2(xOffset + segmentWidth * (i + 1), yOffset + 0), tgt::vec2::zero, tgt::vec2(concentrationOutport_.getSize()), -tgt::vec2::one, tgt::vec2::one));
	}
	IMode.end();
	// Draw bounding box for legend
	IMode.begin(tgt::ImmediateMode::LINE_LOOP);
	IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
	IMode.vertex(mapRange(tgt::vec2(xOffset, yOffset), tgt::vec2::zero, tgt::vec2(concentrationOutport_.getSize()), -tgt::vec2::one, tgt::vec2::one));
	IMode.vertex(mapRange(tgt::vec2(xOffset, yOffset + height), tgt::vec2::zero, tgt::vec2(concentrationOutport_.getSize()), -tgt::vec2::one, tgt::vec2::one));
	IMode.vertex(mapRange(tgt::vec2(xOffset + segmentWidth * (cMap.getColorCount() - 1), yOffset + height), tgt::vec2::zero, tgt::vec2(concentrationOutport_.getSize()), -tgt::vec2::one, tgt::vec2::one));
	IMode.vertex(mapRange(tgt::vec2(xOffset + segmentWidth * (cMap.getColorCount() - 1), yOffset), tgt::vec2::zero, tgt::vec2(concentrationOutport_.getSize()), -tgt::vec2::one, tgt::vec2::one));
	IMode.end();
	// Add text
	// Minimum value
	const tgt::vec3 posMin(0.70 * xOffset, yOffset + 0.5 * height, 0);
	std::stringstream streamMin;
	streamMin << std::fixed << std::setprecision(3) << minValue;
	heatMapFontProp_.get()->render(posMin, streamMin.str(), concentrationOutport_.getSize());
	// Maximum value
	const tgt::vec3 posMax(1.1f * xOffset + segmentWidth * (cMap.getColorCount() - 1), yOffset + 0.5 * height, 0);
	std::stringstream streamMax;
	if (manualUpperBoundProp_.get()) {
		streamMax << ">= ";
	}
	streamMax << std::fixed << std::setprecision(3) << maxValue;
	heatMapFontProp_.get()->render(posMax, streamMax.str(), concentrationOutport_.getSize());
	// Label
	const tgt::vec3 posLabel(0.1f * xOffset, yOffset + 0.5 * height, 0);
	heatMapFontProp_.get()->render(posLabel, label, concentrationOutport_.getSize());
}

} // namespace
