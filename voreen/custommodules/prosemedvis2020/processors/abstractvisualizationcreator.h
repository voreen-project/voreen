#ifndef VRN_ABSTRACTVISUALIZATIONCREATOR_H
#define VRN_ABSTRACTVISUALIZATIONCREATOR_H

#include "voreen/core/processors/processor.h"

#include "modules/plotting/ports/plotport.h"
#include "../ports/templatevesselgraphport.h"
#include "../ports/concretevesselgraphport.h"
#include "voreen/core/ports/genericport.h"
#include "voreen/core/ports/geometryport.h"
#include "custommodules/flowsimulation/ports/flowparametrizationport.h"

#include "../datastructures/templatevesselgraph.h"
#include "../datastructures/concretevesselgraph.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"
namespace voreen {

    class VRN_CORE_API AbstractVisualizationCreator : public Processor {
    public:
        AbstractVisualizationCreator();
        virtual Processor* create() const;

        virtual std::string getClassName() const { return "AbstractVisualizationCreator"; }
        virtual std::string getCategory() const { return "Volume Processing"; }

    private:
        TemplateVesselGraphPort templateVesselGraphInport;

        VolumeListPort controlVolumeInport;
        ConcreteVesselGraphPort controlConcreteVesselGraphInport;
        VolumeListPort diseasedVolumeInport;
        ConcreteVesselGraphPort diseasedConcreteVesselGraphInport;
        VolumeListPort treatedVolumeInport;
        ConcreteVesselGraphPort treatedConcreteVesselGraphInport;

        StringOptionProperty edgeLabel;
        FloatProperty edgeSlider;
		StringListProperty controlEdgeIds;
		StringListProperty diseasedEdgeIds;
		StringListProperty treatedEdgeIds;
        ButtonProperty calculateButton;

        PlotPort plotOutport;

        GeometryPort controlPointsOutport;
        GeometryPort diseasedPointsOutport;
        GeometryPort treatedPointsOutport;

        FlowParametrizationPort controlParameterOutport;
        FlowParametrizationPort diseasedParameterOutport;
        FlowParametrizationPort treatedParameterOutport;

		

        std::map<std::string, float> rememberValues;
        virtual void process();

        void setLabels();
        void remember();
        void calculate();
		void selectIdInProperty(std::string edgelabel);
    };
}

#endif

