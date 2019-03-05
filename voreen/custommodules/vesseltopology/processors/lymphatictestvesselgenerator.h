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

#ifndef VRN_LYMPHATICTESTVESSELGENERATOR
#define VRN_LYMPHATICTESTVESSELGENERATOR

// This was very much an experiment and is not ready to use for anything at all, really.

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"

#include "tgt/vector.h"
#include <random>

namespace voreen {

class LymphaticTestVesselGenerator : public VolumeProcessor {
public:
    LymphaticTestVesselGenerator();
    virtual ~LymphaticTestVesselGenerator();

    virtual std::string getClassName() const         { return "LymphaticTestVesselGenerator";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual std::string setDescriptions() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() { setDescription( "Processor generates an artificial lymphatic vessel like binary segmentation and ground truth vessel graph"); }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual void process();

private:
    VolumePort outport_;

    BoolProperty enabledProp_;

    IntVec3Property dimensions_;
    FloatVec3Property spacing_;
    BoolProperty usePredeterminedSeed_;
    IntProperty predeterminedSeed_;

    std::random_device randomDevice_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_LYMPHATICTESTVESSELGENERATOR
