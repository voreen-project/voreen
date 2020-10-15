#pragma once
#ifndef VRN_UNCERTAINVECTORFIELDPROCESSOR_H
#define VRN_UNCERTAINVECTORFIELDPROCESSOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/string/stringlistproperty.h"
#include "voreen/core/ports/volumeport.h"

#include "modules/ensembleanalysis/ports/ensembledatasetport.h"

namespace voreen
{
	class UncertainVectorFieldProcessor : public Processor
	{
	public:
		UncertainVectorFieldProcessor();

		Processor* create() const override
		{
			return new UncertainVectorFieldProcessor();
		}
		std::string getClassName() const override
		{
			return "UncertainVectorFieldProcessor";
		}
		std::string getCategory() const override
		{
			return "Volume Processing";
		}

	private:
		void process() override
		{}
		void update();

		EnsembleDatasetPort _inportEnsemble;
		VolumePort _inportMask, _outportQ, _outportL;

		StringListProperty _propertySelectedMembers;
		IntProperty _propertyTimestep, _propertySampleCount;
		FloatProperty _propertyThresholdQ, _propertyThresholdL;
		ButtonProperty _propertyUpdate;
	};
}

#endif // VRN_UNCERTAINVECTORFIELDPROCESSOR_H