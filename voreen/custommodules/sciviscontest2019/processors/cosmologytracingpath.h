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

#ifndef VRN_COSMOLOGYTRACINGPATH_H
#define VRN_COSMOLOGYTRACINGPATH_H


#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/renderport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/colorproperty.h"

#include "voreen/core/interaction/camerainteractionhandler.h"

#include "../ports/cmparticleport.h"

//use namespace voreen
namespace voreen {

/**
 * Sample render processor, which gray-scales an input image based on a user-defined parameter.
 * VRN_CORE_API is a macro needed for shared libs on windows (see voreencoreapi.h)
 */
class VRN_CORE_API CosmologyTracingPath : public RenderProcessor {

public:
    /**
     * Constructor
     */
    CosmologyTracingPath();
    ~CosmologyTracingPath();

    //------------------------------------------
    //  Pure virtual functions of base classes
    //------------------------------------------
    virtual Processor* create() const { return new CosmologyTracingPath();    }
    virtual std::string getClassName() const { return "CosmologyTracingPath"; }
    virtual std::string getCategory() const  { return "Viscontest2019";       }

protected:

    virtual void setDescriptions() { setDescription("Traces the path of particle in Input Port. Slows down with increasing particle count!"); }
    virtual void process();

    /**
     * Overwrites the base implementation of this function.
     * It is used to load the needed shader.
     * @see Processor
     */
    virtual void initialize(); 

    /**
     * Overwrites the base implementation of this function.
     * It is used to free the used shader.
     * @see Processor
     */
    virtual void deinitialize();

private:
    void setupBuffers();
    void setupBounds();
    void invalidateBuffers();

    //-------------
    //  members
    //-------------
    RenderPort outport_;            ///< output of the modified image
    CMParticlePort inport_;
    CameraProperty camera_;
    ShaderProperty shaderProp_;
    TransFunc1DKeysProperty transFunc_;
    FloatIntervalProperty timeStep_;
    BoolProperty useAlpha_;
    FloatProperty alphaFactor_;
    BoolProperty adaptCamera_;
	BoolProperty traceJumps_;

	ColorProperty baryonColor;
	ColorProperty darkMatterColor;
	ColorProperty starColor;
	ColorProperty windColor;
	ColorProperty gasColor;
	ColorProperty agnColor;


    FloatBoundingBoxProperty universeDimensions_;
    FloatMat4Property universeMatrix_;


    CameraInteractionHandler* cameraHandler_;


    std::vector<tgt::mat4> normMatricies_;

    GLuint vbo_;
    GLuint vao_;
    GLuint ibo_;

    int drawCount_;
    
    std::vector<tgt::Bounds> timeSliceBounds_;
    bool buffersInvalid_;

};

} // namespace

#endif // VRN_COSMOLOGYTRACINGPATH_H
