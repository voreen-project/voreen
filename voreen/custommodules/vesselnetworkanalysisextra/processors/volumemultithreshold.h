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

#ifndef VRN_VOLUMEMULTITHRESHOLD
#define VRN_VOLUMEMULTITHRESHOLD

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "tgt/vector.h"

namespace voreen {

/**
 * Used to create a list binary volumes from a list of
 * input volumes using global thresholding.
 *
 * !!!!!
 * Note: This processor currently includes a hack that stores the threshold used as the offset of the RWM!!
 * !!!!!
 */


class VolumeMultiThreshold : public VolumeProcessor {
public:
    VolumeMultiThreshold();
    virtual ~VolumeMultiThreshold();

    virtual std::string getClassName() const         { return "VolumeMultiThreshold";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual std::string setDescriptions() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() { setDescription("Volume processor that performs multiple binarization using a range of thresholds producing a VolumeList"); }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual void process();

private:
    VolumePort inport_;
    VolumeListPort outport_;

    FloatIntervalProperty binarizationThresholdRange_;
    IntProperty steps_;
    ButtonProperty start_;

    bool startPressed_;

    static const std::string loggerCat_;
};
} // namespace voreen

#endif
