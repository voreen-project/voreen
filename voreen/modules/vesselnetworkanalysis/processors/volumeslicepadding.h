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

#ifndef VRN_VOLUMESLICEPADDING
#define VRN_VOLUMESLICEPADDING

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "tgt/vector.h"

namespace voreen {

class VolumeSlicePadding : public VolumeProcessor {
public:
    VolumeSlicePadding();
    virtual ~VolumeSlicePadding();

    virtual std::string getClassName() const         { return "VolumeSlicePadding";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual std::string setDescriptions() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() {
        setDescription(
            "Processor that modifies an input volume by adding slices before/after the first/last slice. "
            "This is escpecially useful when trying to extract the vessel network topology from a 'nearly-2D' volume using <b>VesselGraphCreator</b> to avoid connections to the edge of the volume."
            );
    }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual void process();

private:
    VolumePort inport_;
    VolumePort outport_;

    BoolProperty enabledProp_;
    FloatProperty sliceVoxelValue_;
    IntProperty numSlicesBefore_;
    IntProperty numSlicesAfter_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMESLICEPADDING
