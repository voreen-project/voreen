/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include <algorithm>

#include "rotationaldirectionprocessor.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"


namespace voreen
{

	RotationalDirectionProcessor::RotationalDirectionProcessor()
		: Processor(), inport_(Port::INPORT, "inport", "Curl"), inport2_(Port::INPORT, "inport2", "vortex"), outport_(Port::OUTPORT, "outport", "vortex")

	{
		//inport_.addCondition(new PortConditionVolumeType("Matrix3(float)", "VolumeRAM_Mat3Float"));
		addPort(inport_);
		addPort(inport2_);
		addPort(outport_);
	}

	Processor *RotationalDirectionProcessor::create() const
	{
		return new RotationalDirectionProcessor();
	}

	void RotationalDirectionProcessor::Process( const VolumeRAM_3xFloat& curl, Vortex& vortex )
	{
		int cCounter = 0;
		int ccCounter = 0;

		for(const auto& position : vortex.coreline())
		{
			if (curl.voxel(position).z > 0.0) { //changex to ccounter to counter and vice versa
				cCounter = cCounter + 1; 
			}
			else { 
				ccCounter = ccCounter + 1; 
			}
		}

		vortex.setOrientation(ccCounter > cCounter? Vortex::Orientation::eCounterClockwise : Vortex::Orientation::eClockwise);
	}

	void RotationalDirectionProcessor::process()
	{
		auto input = inport_.getData();
		auto Vort = inport2_.getData();	
		auto CurlVolume = dynamic_cast<const VolumeRAM_3xFloat *>(input->getRepresentation<VolumeRAM>());		

		Vortex* Vort2 = new Vortex(Vortex::Orientation::eUnknown, Vort->coreline());
		RotationalDirectionProcessor::Process( *CurlVolume, *Vort2 );

		outport_.setData(Vort2);
	}

} // namespace voreen
