#ifndef VRN_CMPLOTPORT_H
#define VRN_CMPLOTPORT_H

#include "voreen/core/ports/genericport.h"
#include "../datastructures/cmplotdata.h"

namespace voreen {

class CMPlotPort : public GenericPort<CMPlotData> {
	public:

	CMPlotPort(PortDirection direction, const std::string& name,
		bool allowMultipleConnections = false,
		Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);
	~CMPlotPort();

	virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const { return new CMPlotPort(direction, id); }
	virtual std::string getClassName() const { return "PlotPort"; }
	virtual tgt::col3 getColorHint() const;
};

}

#endif  // VRN_CMPLOTPORT_H