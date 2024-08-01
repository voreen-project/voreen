/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

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