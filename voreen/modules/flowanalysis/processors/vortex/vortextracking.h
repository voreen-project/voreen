#pragma once

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"

#include "modules/flowanalysis/ports/vortexport.h"

namespace voreen
{
	class VortexTracking : public Processor
	{
	public:
		VortexTracking();
		virtual Processor *create() const { return new VortexTracking(); }
		virtual std::string getClassName() const { return "VortexTracking"; }
		virtual std::string getCategory() const { return "Vortex Processing"; }

		static void Process(const Vortex &vortex, const std::vector<Vortex> &vortices, float maxDistanceSameCoreline, size_t& outTrackedVortexIndex);

	protected:
		virtual void process();

	private:
		VortexPort _inVortex;
		VortexListPort _inVortices;
		GeometryPort _outTrackedCoreline;
		FloatProperty _maxDistanceSameCoreline;
	};
} // namespace voreen
