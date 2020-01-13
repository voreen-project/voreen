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

#ifndef VRN_FLOWARROWRENDERER3D_H
#define VRN_FLOWARROWRENDERER3D_H
//super class
#include "voreen/core/processors/renderprocessor.h"
//ports
#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"
//properties
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/interaction/camerainteractionhandler.h" //< for camera interaction
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "tgt/shadermanager.h"

namespace voreen {

class FlowArrowRenderer3D : public RenderProcessor {

public:
    enum ArrowColorMode{
        ARROW_LENGTH = 0,
        ARROW_DIRECTION = 1
    };

    enum ArrowComponent{
        AC_ALL_DIM = 0,
        AC_ONLY_X_DIM = 1,
        AC_ONLY_Y_DIM = 2,
        AC_ONLY_Z_DIM = 3
    };

    FlowArrowRenderer3D();
    virtual ~FlowArrowRenderer3D();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "FlowArrowRenderer3D";  }
    virtual std::string getCategory() const  { return "Experimental Processor";     }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("Renders Arrows");
    }

    virtual void process();
    virtual bool isReady() const;
    virtual void initialize();
    virtual void deinitialize();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

    virtual void beforeProcess();
    virtual void compile();

    void changePropertyVisibility();
private:
    //ports
    VolumePort volumeInport_;
    RenderPort renderOutport_;
    //shader
    ShaderProperty shaderProp_;
    tgt::Shader* program_;
    //camera
    CameraProperty cameraProp_;
    CameraInteractionHandler* cameraHandler_;
    //other properties
    TransFunc1DKeysProperty tfProp_;
    FloatProperty arrowSize_;
    OptionProperty<FlowArrowRenderer3D::ArrowColorMode> arrowColorMode_;
    OptionProperty<FlowArrowRenderer3D::ArrowComponent> arrowComponentMode_;

    FloatProperty xVoxelOffset_;        ///< the shift inside a voxel
    FloatProperty yVoxelOffset_;        ///< the shift inside a voxel
    FloatProperty zVoxelOffset_;        ///< the shift inside a voxel

    bool gpuGsSupport_;
    std::string gpuErrorString_;

    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_FlowARROWRENDERER3D_H
