#include "parallelvectorsolutionpointsport.h"
#include <sstream>

namespace voreen
{
	ParallelVectorSolutionPointsPort::ParallelVectorSolutionPointsPort(PortDirection direction, const std::string &id, const std::string &guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
		: GenericPort<ParallelVectorSolutions>(direction, id, guiName, allowMultipleConnections, invalidationLevel) {}

	Port *ParallelVectorSolutionPointsPort::create(PortDirection direction, const std::string &id, const std::string &guiName) const
	{
		return new ParallelVectorSolutionPointsPort(direction, id, guiName);
	}
	std::string ParallelVectorSolutionPointsPort::getClassName() const
	{
		return "ParallelVectorSolutionPointsPort";
	}
	std::string ParallelVectorSolutionPointsPort::getContentDescriptionHTML() const
	{
		auto stream = std::stringstream();
		stream << Port::getContentDescriptionHTML();

		if (this->hasData())
		{
			stream << "<br/>Size: " << this->getData()->solutions.size();
		}

		return stream.str();
	}
} // namespace voreen
