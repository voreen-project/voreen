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
