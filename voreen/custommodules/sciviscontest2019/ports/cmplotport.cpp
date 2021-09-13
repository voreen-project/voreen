#include "cmplotport.h"

namespace voreen {

	CMPlotPort::CMPlotPort(PortDirection direction, const std::string& name,
		bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
		: GenericPort<CMPlotData>(direction, name, name, allowMultipleConnections, invalidationLevel)
	{
		if (isOutport())
			setData(0, false);
	}

	CMPlotPort::~CMPlotPort() {
		if ((portData_) && (isOutport())) {
			delete portData_;
			portData_ = 0;
		}
	}

	tgt::col3 CMPlotPort::getColorHint() const {
		return tgt::col3(204, 0, 153);
	}

}