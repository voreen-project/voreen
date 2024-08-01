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

#ifndef VRN_COSMOLOGYARROWRENDERER3D_H
#define VRN_COSMOLOGYARROWRENDERER3D_H
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
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"

#include "tgt/shadermanager.h"
#include "../ports/cmparticleport.h"

namespace voreen {

class CosmologyArrowRenderer3D : public RenderProcessor {

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

    CosmologyArrowRenderer3D();
    virtual ~CosmologyArrowRenderer3D();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "CosmologyArrowRenderer3D"; }
    virtual std::string getCategory() const  { return "Viscontest2019";           }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;    }

protected:
    virtual void setDescriptions() {
        setDescription("Renders Arrows from velocity");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    virtual void compile();

    void changePropertyVisibility();

    void buildBuffers();

private:
    struct VertexLayout{
        tgt::vec3 pos;
        tgt::vec3 dir;
    };

    //ports
    CMParticlePort inport_;
    RenderPort     renderOutport_;
    // time
    IntProperty  timeStep_;

	BoolProperty useAlpha_;
	FloatProperty alphaFactor_;

    //shader
    ShaderProperty shaderProp_;
    ShaderProperty shaderPropFast_;
    //camera
    CameraProperty            cameraProp_;
    CameraInteractionHandler* cameraHandler_;
    //other properties
    TransFunc1DKeysProperty tfProp_;
    FloatProperty           arrowSize_;
    OptionProperty<CosmologyArrowRenderer3D::ArrowColorMode> arrowColorMode_;
    OptionProperty<CosmologyArrowRenderer3D::ArrowComponent> arrowComponentMode_;

    FloatProperty xVoxelOffset_;        ///< the shift inside a voxel
    FloatProperty yVoxelOffset_;        ///< the shift inside a voxel
    FloatProperty zVoxelOffset_;        ///< the shift inside a voxel

    BoolProperty reduceQuality_;
    BoolProperty adaptCamera_;
    FloatBoundingBoxProperty universeDimensions_;
    FloatMat4Property universeMatrix_;


    GLuint vbo_;
    int    vertexCount_;
    float  maxLength_;
    
    static const std::string loggerCat_;
};

}   //namespace

#endif // VRN_CosmologyARROWRENDERER3D_H
