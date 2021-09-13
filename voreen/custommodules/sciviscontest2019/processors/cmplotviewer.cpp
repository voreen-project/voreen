#include "cmplotviewer.h"

#include "tgt/textureunit.h"

#include <iostream>

#include <stdlib.h> 

namespace voreen {

	static const tgt::ivec2 MARGINS(75, 50);

	CMPlotViewer::CMPlotViewer()
		: RenderProcessor()
		, inport_(Port::INPORT, "cmplotdata", "Particle Plot Input")
		, outport_(Port::OUTPORT, "outport", "Modified Image", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
		, rowSelector_("selectedrows", "Select Rows to Display")
		, lineThickness_("lineThickness", "Line Thickness", 1.0f, 0.5f, 50.0f)
		, rowSelection_("Selected Rows", "Selected Rows")
		, plotLib_(new PlotLibraryOpenGl())
		, rowSelected_(std::vector<bool>())
		, colorProp_("Select Color", "Select Color", tgt::vec4(0.0f, 0.0f, 0.0f, 1.0f))
		, deletebtn_("deletebtn", "Delete Row")
		, plotColors()
		, selectedRows_()
	{
		addPort(inport_);
		addPort(outport_);
		ON_CHANGE(inport_, CMPlotViewer, inDataChange);

		addProperty(lineThickness_);

		addProperty(rowSelector_);
		addProperty(rowSelection_);
		addProperty(colorProp_);
		addProperty(deletebtn_);
		ON_CHANGE(rowSelector_, CMPlotViewer, rowSelected);
		ON_CHANGE(deletebtn_, CMPlotViewer, deleteRow);
		ON_CHANGE(colorProp_, CMPlotViewer, colorSelection);

	}

	CMPlotViewer::~CMPlotViewer() {
		
	}

	void CMPlotViewer::inDataChange()
	{
		CMPlotData* data = (CMPlotData*)inport_.getData();
		if (data == NULL) {
			selectedRows_.clear();
			rowSelection_.reset();
			return;
		}
		rowSelected_.clear();
		rowSelected_.resize(data->getDataRowCount());

		rowSelection_.blockCallbacks(true);
		rowSelection_.reset();
		rowSelection_.addRow(data->getDataRow(0).name);
		rowSelection_.blockCallbacks(false);
		rowSelector_.blockCallbacks(true);
		rowSelector_.reset();
		std::vector<std::string> optionStrings = rowSelector_.getKeys();
		for (std::string option : optionStrings) {
			rowSelector_.removeOption(option);
		}
		for (int i = 0; i < data->getDataRowCount(); i++) {
			rowSelector_.addOption(data->getDataRow(i).name, data->getDataRow(i).name, data->getDataRow(i).name);
			rowSelected_[i] = false;
		}
		rowSelected_[0] = true;
		rowSelector_.blockCallbacks(false);


		plotColors.clear();
		plotColors.push_back(tgt::vec3(0.0f, 0.0f, 0.0f));
		selectedRows_.clear();
		selectedRows_.push_back(0);
	}

	void CMPlotViewer::rowSelected()
	{
		if (!rowSelected_[rowSelector_.getSelectedIndex()]) {
			rowSelection_.blockCallbacks(true);
			rowSelection_.addRow(((CMPlotData*)inport_.getData())->getDataRow(rowSelector_.getSelectedIndex()).name);
			rowSelection_.blockCallbacks(false);

			rowSelected_[rowSelector_.getSelectedIndex()] = true;
			plotColors.push_back(tgt::vec3(0.0f, 0.0f, 0.0f));
			selectedRows_.push_back(rowSelector_.getSelectedIndex());
		}
		
	}

	void CMPlotViewer::deleteRow()
	{
		
		std::vector<int> selected = rowSelection_.getSelectedRowIndices();
		if (selected.size() < 1) {
			return;
		}
		rowSelection_.blockCallbacks(true);

		for (int currentSelection = selected.size() - 1; currentSelection >= 0; currentSelection--) {

			rowSelection_.removeRow(selected[currentSelection]);

			plotColors.erase(plotColors.begin() + selected[currentSelection]);

			rowSelected_[selectedRows_[selected[currentSelection]]] = false;

			selectedRows_.erase(selectedRows_.begin() + selected[currentSelection]);

		}

		rowSelection_.blockCallbacks(false);

	}

	void CMPlotViewer::colorSelection() {

		std::vector<int> selected = rowSelection_.getSelectedRowIndices();
		if (selected.size() < 1) {
			return;
		}

		tgt::vec4 color = colorProp_.get();

		for (int selection : selected) {
			plotColors[selection] = tgt::vec3(color.r, color.g, color.b);
			//rowSelection_.setRowColor(selection, tgt::vec3(color.r, color.g, color.b)); //TODO: make color adjustable
		}

	}

	void CMPlotViewer::initialize() {
		RenderProcessor::initialize();
		
		glGenVertexArrays(1, &vao_);
		glGenBuffers(1, &ssbo_);

		glGenBuffers(1, &ebo_);

		GLuint elements[] = {
			0, 1, 2,
			2, 3, 0
		};

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);

		vertexShader_ = glCreateShader(GL_VERTEX_SHADER);
		fragmentShader_ = glCreateShader(GL_FRAGMENT_SHADER);

		shaderProgram_ = glCreateProgram();

	}

	void CMPlotViewer::deinitialize() {
		RenderProcessor::deinitialize();

		glDeleteProgram(shaderProgram_);
		glDeleteShader(fragmentShader_);
		glDeleteShader(vertexShader_);

		glDeleteBuffers(1, &ebo_);
		glDeleteBuffers(1, &ssbo_);

		glDeleteVertexArrays(1, &vao_);

	}

	struct lineSegment {
		float x1;
		float y1;
		float x2;
		float y2;
	};

	const GLchar* vertexSource = R"glsl(
#version 430 core

	#define featherpixels 1.0f
	
	layout(location = 0) uniform float lineThickness;
	layout(location = 1) uniform vec2 canvaspixel;
	layout(location = 2) uniform vec3 inColor;

	struct segment{
		vec2 start;
		vec2 end;
	};

	layout(binding = 0, std430) buffer posblock
	{
		segment segArray[];
	}segBlock;

	layout(location = 0) out float lineThicknessOut;
	layout(location = 1) out vec2 linecenter;
	layout(location = 3) out vec3 outColor;

	out vec4 gl_Position;
	
	vec4 transformLine(int lineid){
		segment currentSegment = segBlock.segArray[lineid];
		vec2 direction = normalize( vec2( currentSegment.end.x - currentSegment.start.x, currentSegment.end.y - currentSegment.start.y));
		vec2 normal = vec2(0.0f - direction.y, direction.x);

		//normal should always point upwards for convenience
		if (normal.y < 0){
			normal = normal * -1.0f;
		}
		vec2 linepoint = (gl_VertexID == 0 || gl_VertexID == 3 )? currentSegment.start : currentSegment.end;
		float linedown = (gl_VertexID == 3 || gl_VertexID == 2)? -1.0f : 1.0f;

		// aspect ratio normal correction
		normal.y = normal.y * (canvaspixel.x / canvaspixel.y);
		normal.x = normal.x * (canvaspixel.y / canvaspixel.x);
		normal = normalize(normal);
	
		// transform to screenspace for pixel accuracy
		vec2 normalTransformed = normal * 0.5f * canvaspixel;
		float scaleFactorNormal = ((lineThickness / 2.0f) + featherpixels) / length(normalTransformed) ;
		
		//in pixel coordinates
		linecenter = (linepoint + vec2(1,1)) * 0.5f * canvaspixel;	
		
		return vec4(linepoint + (normal * scaleFactorNormal * linedown) , 0.0f, 1.0f);
	}
	
	vec4 transformStart(int lineid){
		segment currentSegment = segBlock.segArray[lineid];
		vec2 direction = normalize( vec2( currentSegment.start.x - currentSegment.end.x, currentSegment.start.y - currentSegment.end.y));
		vec2 normal = vec2(0.0f - direction.y, direction.x);
		
		//normal should always point upwards for convenience
		if (normal.y < 0){
			normal = normal * -1.0f;
		}
		
		//in pixel coordinates
		linecenter = (currentSegment.start + vec2(1,1)) * 0.5f * canvaspixel;
		
		float dirTransform = (gl_VertexID == 0 || gl_VertexID == 3 )? 1.0f : 0.0f;
		float linedown = (gl_VertexID == 3 || gl_VertexID == 2)? -1.0f : 1.0f;

		// aspect ratio normal correction
		normal.y = normal.y * (canvaspixel.x / canvaspixel.y);
		normal.x = normal.x * (canvaspixel.y / canvaspixel.x);
		normal = normalize(normal);
		
		// transform to screenspace for pixel accuracy
		vec2 normalTransformed = normal * 0.5f * canvaspixel;
		float scaleFactorNormal = ((lineThickness / 2.0f) + featherpixels) / length(normalTransformed) ;
		
		// transform to screenspace for pixel accuracy
		vec2 directionTransformed = direction * 0.5f * canvaspixel;
		float scaleFactorDirection = ((lineThickness / 2.0f) + featherpixels) / length(directionTransformed) ;
		
		return vec4(currentSegment.start + (normal * scaleFactorNormal * linedown) + (direction * scaleFactorDirection * dirTransform), 0.0f, 1.0f);
	}
	
	vec4 transformEnd(int lineid){
		segment currentSegment = segBlock.segArray[lineid];
		vec2 direction = normalize( vec2( currentSegment.end.x - currentSegment.start.x, currentSegment.end.y - currentSegment.start.y));
		vec2 normal = vec2(0.0f - direction.y, direction.x);
		
		//normal should always point upwards for convenience
		if (normal.y < 0){
			normal = normal * -1.0f;
		}
		
		//in pixel coordinates
		linecenter = (currentSegment.end + vec2(1,1)) * 0.5f * canvaspixel;
		
		float dirTransform = (gl_VertexID == 0 || gl_VertexID == 3 )? 0.0f : 1.0f;
		float linedown = (gl_VertexID == 3 || gl_VertexID == 2)? -1.0f : 1.0f;
		
		// aspect ratio normal correction
		normal.y = normal.y * (canvaspixel.x / canvaspixel.y);
		normal.x = normal.x * (canvaspixel.y / canvaspixel.x);
		normal = normalize(normal);

		// transform to screenspace for pixel accuracy
		vec2 normalTransformed = normal * 0.5f * canvaspixel;
		float scaleFactorNormal = ((lineThickness / 2.0f) + featherpixels) / length(normalTransformed) ;
		
		// transform to screenspace for pixel accuracy
		vec2 directionTransformed = direction * 0.5f * canvaspixel;
		float scaleFactorDirection = ((lineThickness / 2.0f) + featherpixels) / length(directionTransformed) ;
		
		return vec4(currentSegment.end + (normal * scaleFactorNormal * linedown) + (direction * scaleFactorDirection * dirTransform), 0.0f, 1.0f);
	}
	
	void main()
	{
		int lineID = int(floor(gl_InstanceID / 3.0f));
		
		float drawID = mod(gl_InstanceID, 3);
		
		// line start butt
		if (drawID == 0){
			gl_Position = transformStart(lineID);
		} 
		// line end butt
		else if (drawID == 2){
			gl_Position = transformEnd(lineID);
		}
		// line
		else {
			gl_Position = transformLine(lineID);
		}
	
		//in pixels
		lineThicknessOut = lineThickness / 2.0f;

		outColor = inColor;
	}

	)glsl";
	const GLchar* fragmentSource = R"glsl(
	#version 430 core

	layout(location = 0) in float lineThickness;
	layout(location = 1) in vec2 linecenter;
	layout(location = 3) in vec3 color;

	out vec4 outColor;

	void main()
	{
		vec4 col = vec4(color, 1.0f);
		float d = length(linecenter - gl_FragCoord.xy);
		float w = lineThickness;
		if( w < d){
			col.w = 0.33f;
			if ( d-w > 1.5f){
				col.w = 0.0f;
			}
		} /*else {
		}*/
		outColor = col;
	}
	)glsl";

	void CMPlotViewer::process() {

		outport_.activateTarget();
		outport_.clearTarget();
		outport_.getRenderTarget()->activateTarget();

		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		renderAxes();
		tgt::ivec2 dataViewPort = tgt::ivec2(outport_.getSize().x - (MARGINS.x * 2), outport_.getSize().y - (MARGINS.y * 2));
		if (dataViewPort.x > 0 && dataViewPort.y > 0) {
			renderData(MARGINS, outport_.getSize());
		}

		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		outport_.deactivateTarget();
		tgt::TextureUnit::setZeroUnit();
		LGL_ERROR;
	}

	void CMPlotViewer::renderAxes() {
		// Set Plot status.

		CMPlotData* plotData = (CMPlotData*)inport_.getData();

		if (plotData == NULL) {
			LWARNING("Data in is NULL");
		}

		tgt::vec2 valueRange = tgt::vec2(0.0f, 1.0f);
		for (int i : selectedRows_) {
			CMPlotDataRow data = plotData->getDataRow(i);
			valueRange.x = std::min(valueRange.x, data.min);
			valueRange.y = std::max(valueRange.y, data.max);
		}

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
		tgt::vec2 xAxisRange = tgt::vec2(0, plotData->getDataRowLength()-1);
		tgt::vec2 yAxisRange = valueRange;

		plotLib_->setDomain(Interval<plot_t>(xAxisRange.x, xAxisRange.y), PlotLibrary::X_AXIS);
		plotLib_->setDomain(Interval<plot_t>(yAxisRange.x, yAxisRange.y), PlotLibrary::Y_AXIS);

		float fontSize = 10.0f;

		if (plotLib_->setRenderStatus()) {
			plotLib_->setDrawingColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
			plotLib_->renderAxes();
			plotLib_->setDrawingColor(tgt::Color(0, 0, 0, .5f));
			plotLib_->setFontSize(fontSize);
			plotLib_->setFontColor(tgt::Color(0.f, 0.f, 0.f, 1.f));
			plotLib_->renderAxisScales(PlotLibrary::X_AXIS, false);
			plotLib_->renderAxisScales(PlotLibrary::Y_AXIS, false);
			plotLib_->setFontSize(fontSize + 2);
			plotLib_->renderAxisLabel(PlotLibrary::X_AXIS, "Timestep");
			plotLib_->renderAxisLabel(PlotLibrary::Y_AXIS, "Count");
		}

		plotLib_->resetRenderStatus();
	}

	void CMPlotViewer::renderData(tgt::vec2 viewPortOffset, tgt::vec2 viewPortSize){

		float lineThickness = lineThickness_.get();

		struct lineSegment {
			tgt::vec2 start;
			tgt::vec2 end;
		};

		CMPlotData* plotData = (CMPlotData*)inport_.getData();

		if (plotData == NULL) {
			LWARNING("Data in is NULL");
		}

		// Create and compile the vertex shader
		glShaderSource(vertexShader_, 1, &vertexSource, NULL);
		glCompileShader(vertexShader_);

		// Create and compile the fragment shader
		glShaderSource(fragmentShader_, 1, &fragmentSource, NULL);
		glCompileShader(fragmentShader_);

		// Link the vertex and fragment shader into a shader program
		glAttachShader(shaderProgram_, vertexShader_);
		glAttachShader(shaderProgram_, fragmentShader_);
		glLinkProgram(shaderProgram_);
		glUseProgram(shaderProgram_);

		glBindVertexArray(vao_);

		glUniform1f(0, lineThickness);
		glUniform2f(1, outport_.getSize().x, outport_.getSize().y);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_);

		glDisable(GL_DEPTH_TEST);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		float min = 0.0f;
		float max = 1.0f;

		for (int i : selectedRows_) {
			CMPlotDataRow data = plotData->getDataRow(i);
			min = std::min(min, data.min);
			max = std::max(max, data.max);
		}
		float range = max - min;

		int ofSelected = 0;

		for (int currentRow : selectedRows_) {

			CMPlotDataRow data = plotData->getDataRow(currentRow);
			int datapoints = data.data.size() - 1;
			//int datapoints = 2;

			std::vector<lineSegment> segments = std::vector<lineSegment>(datapoints);

			//segments[0] = { tgt::vec2(-0.5f, 0.0f), tgt::vec2(0.0f, 0.5f) };
			//segments[0] = { tgt::vec2(0.0f, 0.5f), tgt::vec2(0.5f, 0.0f) };

			for (int i = 0; i < datapoints; i++) {
				tgt::vec2 start = tgt::vec2((((float)i) / (datapoints - 1)) * 2.0f - 1.0f, (data.data[i] / range) * 2.0f - 1.0f);
				tgt::vec2 end = tgt::vec2((((float)i + 1) / (datapoints - 1)) * 2.0f - 1.0f, (data.data[i + 1] / range) * 2.0f - 1.0f);

				// to viewport coordinates
				start = start * ((viewPortSize - (viewPortOffset * 2.0f)) / viewPortSize);
				end = end * ((viewPortSize - (viewPortOffset * 2.0f)) / viewPortSize);
				segments[i] = { start, end };
			}

			glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(lineSegment) * datapoints, segments.data(), GL_STATIC_DRAW);


			glBindBufferRange(GL_SHADER_STORAGE_BUFFER, 0, ssbo_, 0, sizeof(lineSegment) * datapoints);
			glUniform3fv(2, 1, &plotColors[ofSelected][0]);

			ofSelected++;

			glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, datapoints * 3);

		}

		glBlendFuncSeparate(GL_ONE, GL_ZERO, GL_ONE, GL_ZERO);

		glDisable(GL_BLEND);

		glEnable(GL_DEPTH_TEST);

		glUseProgram(0);

	}

}