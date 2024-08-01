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

#ifndef VRN_MATERIALDERIVATIVE_H
#define VRN_MATERIALDERIVATIVE_H

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/properties/optionproperty.h"

namespace voreen {

class MaterialDerivative : public Processor {
public:
    MaterialDerivative();

	virtual Processor* create() const             { return new MaterialDerivative(); }
	virtual std::string getClassName() const      { return "MaterialDerivative"; }
	virtual std::string getCategory() const       { return "Volume Processing";     }
	virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL; }
protected:
	virtual void setDescriptions() {
		setDescription(
            "Computes the material derivative Dv/Dt = dv/dt + u * grad v<br>"
            "<ul>"
            "<li>V0: Volume <i>v</i> for timestep <i>t-1</i></li>"
            "<li>V1: Volume <i>v</i> for timestep <i>t</i></li>"
            "<li>Velocity: 3-component flow field</li>"
            "<li>Jacobian: Jacobian of the volume <i>v</i> for timestep <i>t</i></li>"
            "</ul>"
        );
	}

	virtual void process();    
private:
	// Ports
    VolumePort v0_;
	VolumePort v1_;
    VolumePort velocity_;
    VolumePort jacobian_;
    VolumePort outport_;

	static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMEINTERPOLATION_H