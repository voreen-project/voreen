#pragma once
#ifndef VRN_APPROXIMATEPARALLELVECTORS_H
#define VRN_APPROXIMATEPARALLELVECTORS_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"
#include "modules/ensembleanalysis/ports/ensembledatasetport.h"

namespace voreen
{
	class ApproximateParallelVectors : public Processor
	{
	public:
		ApproximateParallelVectors();

		Processor* create() const override;
		std::string getClassName() const override;
		std::string getCategory() const override;

	private:
		void process() override;
		void updateButton();

		EnsembleDatasetPort _inportEnsemble;
		VolumePort _inportMask, _outportA, _outportB, _outportEps;

		IntProperty _propertyTimestep;
		FloatProperty _propertyEpsilon;
		BoolProperty _propertyUseAcceleration;
		ButtonProperty _propertyUpdateButton;
	};
}

#endif // VRN_APPROXIMATEPARALLELVECTORS_H