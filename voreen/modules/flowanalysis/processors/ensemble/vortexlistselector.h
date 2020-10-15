#pragma once
#ifndef VRN_VORTEXLISTSELECTOR_H
#define VRN_VORTEXLISTSELECTOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include "modules/flowanalysis/ports/vortexport.h"

namespace voreen
{
	class VortexListSelector : public Processor
	{
	private:
		enum RotationOptions {
			OPTION_CW, //Clockwise
			OPTION_CCW, //Counterclockwise
			OPTION_B, //Both
		};

	public:
		VortexListSelector();

		Processor* create() const override
		{
			return new VortexListSelector();
		}
		std::string getClassName() const override
		{
			return "VortexListSelector";
		}
		std::string getCategory() const override
		{
			return "Vortex Processing";
		}

		bool isReady() const override
		{
			return _inportVortexCollection.isReady();
		}

		static void Process( const VortexCollection& vortices, const std::vector<int>& runs, int firstTimestep, int lastTimestep, int minLength,RotationOptions rot, std::vector<Vortex>& outVortexList );

	private:
		void process() override;
		void updatePropertyCorelineLength();		

		VortexCollectionPort _inportVortexCollection;
		VortexListPort _outportVortexList;
		GeometryPort _outportGeometry;

		StringListProperty _propertyRuns;
		IntIntervalProperty _propertyTimesteps;
		IntProperty _propertyCorelineLength;
		OptionProperty<RotationOptions> _Rotation;
	};
}

#endif // VRN_VORTEXLISTSELECTOR_H
