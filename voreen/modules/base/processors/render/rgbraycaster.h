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

#ifndef VRN_RGBARAYCASTER_H
#define VRN_RGBARAYCASTER_H

#include "voreen/core/processors/volumeraycaster.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * Performs raycasting of an RGB volume by applying a hue-based transfer function.
 */
class VRN_CORE_API RGBRaycaster : public VolumeRaycaster {
public:
    RGBRaycaster();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "RGBRaycaster";     }
    virtual std::string getCategory() const  { return "Raycasting";       }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING; }

protected:
    virtual void setDescriptions() {
        setDescription("Performs raycasting of an RGB volume data set by applying a hue-based transfer function.");
    }

    virtual void process();
    virtual void initialize();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void compile();

private:

    ShaderProperty shaderProp_;
    TransFunc1DKeysProperty transferFunc_;     ///< the property that controls the transfer-function
    BoolProperty applyColorModulation_;
    CameraProperty camera_;
};

} // namespace

#endif // VRN_RGBARAYCASTER_H
