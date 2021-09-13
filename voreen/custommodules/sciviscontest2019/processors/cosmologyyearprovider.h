#ifndef VRN_COSMOLOGYYEARPROVIDER_H
#define VRN_COSMOLOGYYEARPROVIDER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/ports/textport.h"

namespace voreen {
	class CosmologyYearProvider : public Processor
	{
	public:
		CosmologyYearProvider();
		virtual ~CosmologyYearProvider();
		virtual Processor* create() const { return new CosmologyYearProvider(); };

		virtual std::string getClassName() const { return "CosmologyYearProvider"; }
		virtual std::string getCategory() const { return "Viscontest2019"; }
		virtual void setDescriptions() { setDescription("Provides the lookback time at timestep as Text"); }
		virtual CodeState   getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

	protected:
		virtual void initialize();
		virtual void deinitialize();

		virtual void process();

	private:

		//void timeStepChanged();

		TextPort outport_;
		FloatProperty timeStep_;

	};
}

#endif
