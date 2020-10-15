#pragma once
#ifndef VRN_VORTEXCOLLECTIONSOURCE_H
#define VRN_VORTEXCOLLECTIONSOURCE_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/filedialogproperty.h"

#include "modules/flowanalysis/ports/vortexport.h"

namespace voreen
{
	class VortexCollectionSource : public Processor
	{
	public:
		VortexCollectionSource();
		virtual Processor* create() const;

		virtual std::string getClassName() const;
		virtual std::string getCategory() const;

	private:
		virtual void process();

		VortexCollectionPort _outportVortexCollection;

		FileDialogProperty _propertyFileDialog;
		ButtonProperty _propertyLoadButton;
	};
}

#endif // VRN_VORTEXCOLLECTIONSOURCE_H
