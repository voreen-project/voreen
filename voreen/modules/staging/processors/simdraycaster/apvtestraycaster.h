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

#ifndef VRN_SIMPLERAYCASTER_H
#define VRN_SIMPLERAYCASTER_H

#include "voreen/core/processors/volumeraycaster.h"

#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/shaderproperty.h"

namespace voreen {

/**
 * Modified SimpleRaycaster for testing
 */
class VRN_CORE_API APVTestRaycaster : public VolumeRaycaster {
public:
    APVTestRaycaster();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Raycasting";      }
    virtual std::string getClassName() const  { return "APVTestRaycaster"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("Performs a simple single pass raycasting without lighting, which can be modified through a shader property.");
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void rebuildShader();
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

private:

    ShaderProperty shader_;           ///< the property that stores the used shader
    TransFunc1DKeysProperty transferFunc_;  ///< the property that controls the transfer function
    CameraProperty camera_;           ///< necessary for depth value calculation
};


} // namespace

#endif // VRN_SIMPLERAYCASTER_H
