#ifndef VRN_CMPLOTCREATOR_H
#define VRN_CMPLOTCREATOR_H
#pragma once

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "../ports/cmparticleport.h"
#include "../ports/cmplotport.h"

namespace voreen {


	class VRN_CORE_API CMPlotCreator : public AsyncComputeProcessor<CMParticleData*, CMPlotData*> {
	public:
		CMPlotCreator();
		virtual ~CMPlotCreator();
		virtual Processor* create() const;

		virtual std::string getClassName() const { return "CMPlotCreator"; }
		virtual std::string getCategory() const { return "Plotting"; }
		virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

		virtual ComputeInput prepareComputeInput();
		virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
		virtual void processComputeOutput(ComputeOutput output);

	protected:

		//virtual bool isReady() const;

	private:

		CMParticlePort inport_;
		CMPlotPort outport_;

		CMPlotDataRow data;

		static const std::string loggerCat_;

	};

}

#endif // VRN_CMPLOTCREATOR_H 