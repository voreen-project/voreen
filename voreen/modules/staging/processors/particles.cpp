/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "particles.h"

#ifdef GL_COMPUTE_SHADER //disable compilation for old gl headers

#include "voreen/core/interaction/camerainteractionhandler.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/textureunit.h"
#include "tgt/timer.h"

#include <cmath>
#include <time.h>

//#define NUM_PARTICLES      1024*1024
#define NUM_PARTICLES      128*128
#define WORK_GROUP_SIZE    128
#define SPHERE_RADIUS      0.5f

using tgt::TextureUnit;

namespace voreen {

Particles::Particles()
    : RenderProcessor()
    , outport_(Port::OUTPORT, "renderOutport", "Outport", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , computeShader_("compute.prg", "Compute shader", ShaderFileList(tgt::ShaderObject::COMPUTE_SHADER, "particles.comp"))
    , renderShader_("render.prg", "Render shader",
            ShaderFileList(tgt::ShaderObject::FRAGMENT_SHADER, "particles.frag")
                          (tgt::ShaderObject::VERTEX_SHADER,   "particles.vert"))
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
    , radius_("radius", "Point Radius", 3, 1, 10)
    , bufferObject_(0)
    , vertexArrayID_(0)
{
    addPort(outport_);

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    addInteractionHandler(cameraHandler_);

    addProperty(camera_);
    addProperty(radius_);

    addProperty(computeShader_);
    addProperty(renderShader_);

    eventHandler_.addListenerToBack(this);
    std::srand((unsigned)time(NULL));
}

Particles::~Particles() {
    delete cameraHandler_;
}

Processor* Particles::create() const {
    return new Particles();
}

bool Particles::isReady() const {
    return outport_.isReady();
}

void Particles::initialize() {
    if(GpuCaps.getShaderModel() < tgt::GpuCapabilities::SHADER_MODEL_5) {
        LERROR("This system does not seem to support compute shaders, cannot initialize particles processor.");
        return;
    }

    RenderProcessor::initialize();

    compile();

    glGenVertexArrays(1, &vertexArrayID_);
    glBindVertexArray(vertexArrayID_);
    glGenBuffers( 1, &bufferObject_);
    glBindBuffer( GL_SHADER_STORAGE_BUFFER, bufferObject_);
    glBufferData( GL_SHADER_STORAGE_BUFFER, NUM_PARTICLES * sizeof(DataStruct), NULL, GL_STATIC_DRAW );

    DataStruct* points = (DataStruct*) glMapBuffer( GL_SHADER_STORAGE_BUFFER, GL_WRITE_ONLY );
    for( int i = 0; i < NUM_PARTICLES; i++ ) {
        points[i].pos.x = getRandomFloat() * 2.f * SPHERE_RADIUS;
        points[i].pos.y = getRandomFloat() * 2.f * SPHERE_RADIUS;
        points[i].pos.z = getRandomFloat() * 2.f * SPHERE_RADIUS;
        if(length(points[i].pos.xyz()) < SPHERE_RADIUS)
            points[i].pos.xyz() = normalize(points[i].pos.xyz()) * 2.5f * SPHERE_RADIUS;
        points[i].vel.x = getRandomFloat() * 0.01f * SPHERE_RADIUS;
        points[i].vel.y = getRandomFloat() * 0.01f * SPHERE_RADIUS;
        points[i].vel.z = getRandomFloat() * 0.01f * SPHERE_RADIUS;
    }

    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER );
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,  2,  bufferObject_);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0 );
    glBindVertexArray(0);

    timer_ = VoreenApplication::app()->createTimer(&eventHandler_);
    if (timer_)
        timer_->start(30);
}

void Particles::deinitialize() {
    if (timer_)
        timer_->stop();
    delete timer_;
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER,0);
    if(vertexArrayID_ != 0)
        glDeleteVertexArrays(1, &vertexArrayID_);
    if(bufferObject_ != 0)
        glDeleteBuffers(1, &bufferObject_);
    RenderProcessor::deinitialize();
}

float Particles::getRandomPositiveFloat() const {
    return (float)std::rand() / RAND_MAX;
}

float Particles::getRandomFloat() const {
    return 2.f * getRandomPositiveFloat()  - 1.f;
}

void Particles::compile() {
    std::string computeHeader = generateHeader();

    // there can be no direct output from a compute shader, so
    // we manually delete the FragData0 output variable from the header
    size_t outPos = computeHeader.find("FragData0");
    size_t secondNL = computeHeader.find('\n', outPos);
    size_t firstNL = computeHeader.rfind('\n', outPos);
    computeHeader.erase(firstNL, secondNL - firstNL);

    computeShader_.setHeader(computeHeader);
    computeShader_.rebuild();

    renderShader_.setHeader(generateHeader());
    renderShader_.rebuild();
}

void Particles::process() {

    // compile program if needed
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;

    if(!computeShader_.hasValidShader() || !renderShader_.hasValidShader())
        return;

    glBindVertexArray(vertexArrayID_);
    glBindBuffer(GL_ARRAY_BUFFER, bufferObject_);

    // activate compute shader
    tgt::Shader* computeProg = computeShader_.getShader();
    computeProg->activate();

    static float periodic = 0.f;
    periodic += 0.01f;
    periodic = std::fmod(periodic, 2.f * tgt::PIf);
    computeProg->setUniform("time_", periodic);

    glDispatchCompute( NUM_PARTICLES  / WORK_GROUP_SIZE, 1,  1 );
    glMemoryBarrier( GL_SHADER_STORAGE_BARRIER_BIT );
    computeProg->deactivate();

    outport_.activateTarget();
    outport_.clearTarget();

    // activate rendering shader
    tgt::Shader* renderProg = renderShader_.getShader();
    renderProg->activate();
    setGlobalShaderParameters(renderProg, &camera_.get());

    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(DataStruct), 0);
    // specify offset to jump over velocity vectors to get colors
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(DataStruct), (void*)(2 * sizeof(tgt::vec4)));

    // render points
    glPointSize(static_cast<float>(radius_.get()));
    glDrawArrays( GL_POINTS, 0, NUM_PARTICLES );
    glPointSize(1.f);

    renderProg->deactivate();
    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    LGL_ERROR;
}

void Particles:: timerEvent(tgt::TimeEvent* e) {
    e->accept();
    // FIXME calling the compute shader here does not work due to multiple
    // OpenGL context / VAO issues.
    invalidate();
}

} // voreen namespace

#endif // GL_COMPUTE_SHADER
