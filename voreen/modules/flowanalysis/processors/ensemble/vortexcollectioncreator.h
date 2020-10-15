#pragma once
#ifndef VRN_VORTEXCOLLECTIONCREATOR_H
#define VRN_VORTEXCOLLECTIONCREATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include "modules/ensembleanalysis/ports/ensembledatasetport.h"
#include "modules/flowanalysis/ports/vortexport.h"

#include <chrono>

namespace voreen
{
	class VortexCollectionCreator : public Processor
	{
	public:
		VortexCollectionCreator();

		Processor* create() const override
		{
			return new VortexCollectionCreator();
		}
		std::string getClassName() const override
		{
			return "VortexCollectionCreator";
		}
		std::string getCategory() const override
		{
			return "Vortex Processing";
		}
		bool isReady() const override
		{
			return _inportEnsemble.isReady();
		}

	private:
		void process() override {}
		void updateButton();

		EnsembleDatasetPort _inportEnsemble;
		VortexCollectionPort _outportVortexCollection;

		StringListProperty _propertySelectedMembers;
		IntIntervalProperty _propertyTimestepInterval;
		IntProperty _propertyCorelineLength;
		ButtonProperty _propertyUpdateButton;

		FileDialogProperty _propertyFileDialog;
		ButtonProperty _propertySaveButton;
	};
}

#endif // VRN_VORTEXCOLLECTIONCREATOR_H
