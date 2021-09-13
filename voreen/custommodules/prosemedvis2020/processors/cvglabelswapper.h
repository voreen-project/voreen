#ifndef VRN_CVGLABELSWAPPER_H
#define VRN_CVGLABELSWAPPER_H

#include "voreen/core/processors/processor.h"
#include "../ports/concretevesselgraphport.h"
#include "../datastructures/concretevesselgraph.h"
#include "voreen/core/properties/buttonproperty.h"


namespace voreen {

	class VRN_CORE_API CvgLabelSwapper : public Processor {
	public:
		CvgLabelSwapper();
		virtual Processor* create() const;

		virtual std::string getClassName() const { return "CvgLabelSwapper"; }
		virtual std::string getCategory() const { return "Label Processing"; } 

	private:
		virtual void process();

		ConcreteVesselGraphPort concreteVesselGraphOutport;
		ConcreteVesselGraphPort concreteVesselGraphInport;
		StringOptionProperty edgeLabel1;
		StringOptionProperty edgeLabel2;

		
	};
}

#endif