/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_PYROGENICVORTICITYMAPPER_H
#define VRN_PYROGENICVORTICITYMAPPER_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/properties/optionproperty.h"

namespace voreen {

class PyrogenicVorticityMapper : public Processor {
public:
    PyrogenicVorticityMapper();

	virtual Processor* create() const             { return new PyrogenicVorticityMapper(); }
	virtual std::string getClassName() const      { return "PyrogenicVorticityMapper"; }
	virtual std::string getCategory() const       { return "Volume Processing";     }
	virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL; }

	virtual bool isReady() const override;
protected:
	virtual void setDescriptions() {
		setDescription("Maps the pyrogenic vorticity to a float value");
	}

	virtual void process();    
private:
	// Ports
    VolumePort vorticityInport_;
	VolumePort materialDerivative_;
	VolumePort outputVolume_;
    VolumePort pyroVorticityVolume_;

    IntOptionProperty mappingAlgorithm_;

	BoolProperty enableFirstComponent_;
	BoolProperty enableSecondComponent_;
	BoolProperty enableThirdComponent_;

	FloatProperty epsilonProperty_;

	BoolProperty invertProperty_;

	static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMEINTERPOLATION_H