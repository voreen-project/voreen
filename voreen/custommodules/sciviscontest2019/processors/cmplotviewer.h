#ifndef VRN_CMPLOTVIEWER_H
#define VRN_CMPLOTVIEWER_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/ports/renderport.h"
#include "../ports/cmplotport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include "voreen/core/interaction/camerainteractionhandler.h"

#include "modules/plotting/utils/plotlibrary/plotlibrary.h"
#include "modules/plotting/utils/plotlibrary/plotlibraryopengl.h"

namespace voreen {

class VRN_CORE_API CMPlotViewer : public RenderProcessor {

public:
	CMPlotViewer();
	~CMPlotViewer();

	virtual Processor*  create() const { return new CMPlotViewer(); }
	virtual std::string getClassName() const { return "CMPlotViewer"; }
	virtual std::string getCategory() const { return "Plotting"; }

protected:
	virtual void setDescriptions() { setDescription("Processor that draws particle plots."); }
	virtual void process();

	virtual void initialize();
	virtual void deinitialize();

private:

	void renderData(tgt::vec2 viewPortOffset, tgt::vec2 viewPortSize);
	void renderAxes();

	void inDataChange();
	void rowSelected();
	void deleteRow();
	void colorSelection();

	std::vector<bool> rowSelected_;

	CMPlotPort inport_;
	RenderPort outport_;

	FloatProperty lineThickness_;

	OptionProperty<std::string> rowSelector_;
	StringListProperty rowSelection_;
	ColorProperty colorProp_;
	ButtonProperty deletebtn_;

	std::unique_ptr<PlotLibrary> plotLib_;

	std::vector<tgt::vec3> plotColors;
	std::vector<int> selectedRows_;
	//std::vector<CMPlotDataRow&> selectedRows_;

	GLuint vao_;
	GLuint ssbo_;
	GLuint ebo_;
	GLuint vertexShader_;
	GLuint fragmentShader_;
	GLuint shaderProgram_;

	//static const std::string loggerCat_; ///< category used in logging

};

}

#endif