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

#ifndef VRN_PARTICLES_H
#define VRN_PARTICLES_H

//check for GL_COMPUTE_SHADER define flag (needed under windows without PCH)
#include "tgt/tgt_gl.h"

#ifdef GL_COMPUTE_SHADER //disable compilation for old gl headers

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/shaderproperty.h"

namespace voreen {

class CameraInteractionHandler;

class Particles : public RenderProcessor {

public:
    Particles();
    virtual ~Particles();

    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Compute Shader"; }
    virtual std::string getClassName() const { return "Particles"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING; }

    virtual bool isReady() const;
    virtual void initialize();
    virtual void deinitialize();
    virtual void compile();

protected:

    // struct needed for cpu/gpu compute shader interaction
    struct DataStruct{
        tgt::vec4 pos;
        tgt::vec4 vel;
        tgt::vec4 col;
    };

    virtual void setDescriptions() {
        setDescription("Image processor which performs screen space ambient occlusion.");
    }

    void process();

    float getRandomPositiveFloat() const;
    float getRandomFloat() const;
    virtual void timerEvent(tgt::TimeEvent* e);

    RenderPort outport_;

    CameraProperty camera_;
    IntProperty radius_;
    ShaderProperty computeShader_;
    ShaderProperty renderShader_;

    GLuint bufferObject_;
    GLuint vertexArrayID_;

    tgt::EventHandler eventHandler_;    // A local eventhandler which is added to the timer.
    tgt::Timer* timer_;
    CameraInteractionHandler* cameraHandler_;

    float time_;
};

} // namespace voreen

#endif //GL_COMPUTE_SHADER

#endif //VRN_PARTICLES_H
