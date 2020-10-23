#pragma once
#ifndef VRN_PARALLELVECTORSOLUTIONPOINTSPORT_H
#define VRN_PARALLELVECTORSOLUTIONPOINTSPORT_H

#include "voreen/core/ports/genericport.h"
#include <vector>

namespace voreen
{

	struct ParallelVectorSolutions
	{
		std::vector<tgt::vec3> solutions; // Points at which vectors are parallel
		std::vector<int32_t> triangleSolutionIndices; // Indices of solutions inside each triangle or -1 if there is no solution
		tgt::svec3 dimensions; // Dimensions of the original volumes

		ParallelVectorSolutions() = default;
		ParallelVectorSolutions(std::vector<tgt::vec3> &&solutions, std::vector<int32_t> &&solutionIndices, tgt::svec3 dimensions) : solutions(std::move(solutions)), triangleSolutionIndices(std::move(solutionIndices)), dimensions(dimensions) {}
	};

#ifdef DLL_TEMPLATE_INST
	template class VRN_CORE_API GenericPort<ParallelVectorSolutions>;
#endif

	class VRN_CORE_API ParallelVectorSolutionPointsPort : public GenericPort<ParallelVectorSolutions>
	{
	public:
		ParallelVectorSolutionPointsPort(PortDirection direction, const std::string &id, const std::string &guiName = {}, bool allowMultipleConnections = false, Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT);

		virtual Port *create(PortDirection direction, const std::string &id, const std::string &guiName = {}) const;
		virtual std::string getClassName() const;
		virtual std::string getContentDescriptionHTML() const;
	};
} // namespace voreen

#endif // VRN_PARALLELVECTORSOLUTIONPOINTSPORT_H
