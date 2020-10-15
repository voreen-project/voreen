#pragma once
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/intproperty.h"

#include "modules/flowanalysis/ports/parallelvectorsolutionpointsport.h"

namespace voreen
{
	class CorelineCreator : public Processor
	{
	public:
		CorelineCreator();
		virtual Processor *create() const { return new CorelineCreator(); }
		virtual std::string getClassName() const { return "CorelineCreator"; }
		virtual std::string getCategory() const { return "Geometry"; }

		static void Process( const ParallelVectorSolutions& solutions, int lengthThreshold, std::vector<std::vector<tgt::vec3>>& corelines );

	protected:
		virtual void process();

	private:
		ParallelVectorSolutionPointsPort _in;
		GeometryPort _out;
		IntProperty _lengthThreshold;
	};
} // namespace voreen
