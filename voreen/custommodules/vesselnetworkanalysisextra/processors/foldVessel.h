/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_foldVessel_H
#define VRN_foldVessel_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "tgt/vector.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"


namespace voreen {

class foldVessel : public VolumeProcessor {
public:
    foldVessel();
    virtual ~foldVessel();

    virtual std::string getClassName() const         { return "foldVessel";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual std::string setDescriptions() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() { setDescription( "Creates test Volume with a given distribution"); }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual void process();
    virtual void mapToQuader(const VolumeRAM*, VolumeRAM_UInt8*);


private:

    VolumePort outport_;
    VolumePort inport_;

    BoolProperty enableProcessing_;
    BoolProperty bezierCurve_;
    FloatProperty threshold_;
    FloatProperty bezierX0_;
    FloatProperty bezierX1_;
    FloatProperty bezierY1_;
    FloatProperty bezierX2_;



    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMEFLOODFILL_H
