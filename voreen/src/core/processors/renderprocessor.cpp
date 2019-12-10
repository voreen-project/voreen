/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/processors/renderprocessor.h"

#include "tgt/glmath.h"
#include "tgt/camera.h"
#include "tgt/shadermanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/immediatemode/immediatemode.h"

#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/utils/glsl.h"

#include "voreen/core/properties/cameraproperty.h"

#include <sstream>

using tgt::vec2;
using tgt::vec3;
using tgt::vec4;
using tgt::Color;
using std::map;

namespace voreen {

const std::string RenderProcessor::loggerCat_("voreen.RenderProcessor");

GLuint RenderProcessor::fullscreenVBO__ = 0;
tgt::Shader* RenderProcessor::fullscreenShader_ = 0;
size_t RenderProcessor::fullScreenShaderInUse_ = 0;

RenderProcessor::RenderProcessor()
    : Processor()
{}

void RenderProcessor::initialize() {

    Processor::initialize();

    if (fullscreenVBO__ == 0) {
        // create geometry for fullscreen triangle
        glGenBuffers(1, &fullscreenVBO__);

        const float extent = 3.0f;
        vec2 fullscreenBuffer[] = {
            vec2(-1, -1), vec2(1, -1), vec2(1, 1), vec2(-1, 1)
        };

        glBindBuffer(GL_ARRAY_BUFFER, fullscreenVBO__);
        glBufferData(GL_ARRAY_BUFFER, sizeof(fullscreenBuffer), fullscreenBuffer, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    if (fullscreenShader_ == 0) {
        fullscreenShader_ = ShdrMgr.load("passthrough", generateHeader(), false);
    }

    // count instances that use the shader
    fullScreenShaderInUse_++;

    const std::vector<RenderPort*> pports = getPrivateRenderPorts();
    for (size_t i=0; i<pports.size(); ++i) {
        pports[i]->initialize();
    }
    LGL_ERROR;

    adjustRenderOutportSizes();
    LGL_ERROR;
}

void RenderProcessor::deinitialize() {
    const std::vector<RenderPort*> pports = getPrivateRenderPorts();
    for (size_t i=0; i<pports.size(); ++i) {
        pports[i]->deinitialize();
    }
    LGL_ERROR;

    // decrement counter of shader instances in use
    fullScreenShaderInUse_--;

    // dispose if no further instances use the shader
    if(fullscreenShader_ != 0 && fullScreenShaderInUse_ == 0) {
        ShdrMgr.dispose(fullscreenShader_);
        fullscreenShader_ = 0;
    }

    Processor::deinitialize();
}

void RenderProcessor::invalidate(int inv) {

    Processor::invalidate(inv);

    if (inv == Processor::VALID)
        return;

    if (!isInitialized())
        return;

    // invalidate result of render ports
    /*for (size_t i=0; i<getOutports().size(); ++i) {
        RenderPort* renderPort = dynamic_cast<RenderPort*>(getOutports()[i]);
        if (renderPort)
            renderPort->invalidateResult();
    }*/
}

void RenderProcessor::beforeProcess() {
    Processor::beforeProcess();

    tgtAssert(isInitialized(), "No initialized");
    manageRenderTargets();
    adjustRenderOutportSizes();
}

void RenderProcessor::manageRenderTargets() {
    const std::vector<Port*> outports = getOutports();
    for (size_t i=0; i<outports.size(); ++i) {
        RenderPort* rp = dynamic_cast<RenderPort*>(outports[i]);
        if (rp && !rp->getRenderTargetSharing()) {
            if (rp->isConnected()) {
                if (!rp->hasRenderTarget()) {
                    rp->initialize();
                }
            }
            else {
                if (rp->hasRenderTarget() && rp->getDeinitializeOnDisconnect())
                    rp->deinitialize();
            }
        }
    }
}

void RenderProcessor::adjustRenderOutportSizes() {
    // detect dimension of first connected render inport
    tgt::ivec2 inputDim(-1);
    const std::vector<Port*> inports = getInports();
    for (size_t i=0; i<inports.size() && inputDim == tgt::ivec2(-1); ++i) {
        RenderPort* rp = dynamic_cast<RenderPort*>(inports[i]);
        if (rp && rp->hasRenderTarget()){
            inputDim = rp->getSize();
            break;
        }
    }

    // - assign inport dimension to all connected render outports, which are not size receivers
    // - if a outport is a size receiver, assign its received rendering size to it
    tgt::ivec2 resizeDim = tgt::ivec2(-1);
    const std::vector<Port*> outports = getOutports();
    for (size_t i=0; i<outports.size(); ++i) {
        RenderPort* rp = dynamic_cast<RenderPort*>(outports[i]);
        if (!rp || !rp->hasRenderTarget())
            continue;

        if (rp->getRenderSizePropagation() == RenderPort::RENDERSIZE_DEFAULT) {
            if (inputDim != rp->getSize() && inputDim != tgt::ivec2(-1)) {
                rp->resize(inputDim);
                if (resizeDim == tgt::ivec2(-1))
                    resizeDim = inputDim;
            }
        }
        else if (rp->getRenderSizePropagation() == RenderPort::RENDERSIZE_RECEIVER) {
            RenderSizeReceiveProperty* receiveProp = rp->getSizeReceiveProperty();
            tgtAssert(receiveProp, "Render outport is size receiver, but has no SizeReceiveProperty");
            if (receiveProp->get() != rp->getSize()) {
                rp->resize(receiveProp->get());
                if (resizeDim == tgt::ivec2(-1))
                    resizeDim = receiveProp->get();
            }
        }
        else if (rp->getRenderSizePropagation() == RenderPort::RENDERSIZE_STATIC) {
            // nothing to do
        }
        else {
            LERROR("Render outport has invalid render size propgation mode: " << rp->getQualifiedName());
        }
    }

    // resize private render ports to resizeDim
    if (resizeDim != tgt::ivec2(-1)) {
        for (size_t i=0; i<privateRenderPorts_.size(); i++)
            privateRenderPorts_.at(i)->resize(resizeDim);
    }
}

void RenderProcessor::addPrivateRenderPort(RenderPort* port) {
    port->setProcessor(this);
    privateRenderPorts_.push_back(port);

    map<std::string, Port*>::const_iterator it = portMap_.find(port->getID());
    if (it == portMap_.end())
        portMap_.insert(std::make_pair(port->getID(), port));
    else {
        LERROR("Port with name " << port->getID() << " has already been inserted!");
        tgtAssert(false, std::string("Port with name " + port->getID() + " has already been inserted").c_str());
    }
}

void RenderProcessor::addPrivateRenderPort(RenderPort& port) {
    addPrivateRenderPort(&port);
}

void RenderProcessor::renderQuad() {
    GLuint vao;
    glGenVertexArrays(1, &vao);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, fullscreenVBO__);
    glVertexAttribPointer(0, 2, GL_FLOAT, false, 0, 0);
    glEnableVertexAttribArray(0);
    bool noShader = tgt::Shader::getCurrentProgram() == 0;
    if (noShader) {
        fullscreenShader_->activate();
        IMode.setMatstackUniforms(fullscreenShader_);
        fullscreenShader_->setUniform("colorTex_", 0);
    }
    IMode.setMatstackUniforms();

    glDepthFunc(GL_ALWAYS);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
    glDepthFunc(GL_LESS);

    if (noShader) {
        fullscreenShader_->deactivate();
    }
    glBindVertexArray(0);
    glDeleteVertexArrays(1, &vao);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// Parameters currently set:
// - screenDim_
// - screenDimRCP_
// - cameraPosition_ (camera position in world coordinates)
void RenderProcessor::setGlobalShaderParameters(tgt::Shader* shader, const tgt::Camera* camera, tgt::ivec2 screenDim) {
    shader->setIgnoreUniformLocationError(true);
    if (screenDim == tgt::ivec2(-1)) {
        RenderPort* rp = 0;
        for (size_t i=0; i<getOutports().size(); ++i) {
            rp = dynamic_cast<RenderPort*>(getOutports()[i]);
            if (rp && rp->hasRenderTarget())
                break;
        }
        if (rp) {
            screenDim = rp->getSize();
        }
    }
    if(screenDim != tgt::ivec2(-1)) {
        shader->setUniform("screenDim_", tgt::vec2(screenDim));
        shader->setUniform("screenDimRCP_",  tgt::vec2(1.f) / tgt::vec2(screenDim));
        LGL_ERROR;
    }

    // camera position in world coordinates, and corresponding transformation matrices
    if (camera) {
        shader->setUniform("cameraPosition_", camera->getPosition());
        shader->setUniform("viewMatrix_", camera->getViewMatrix());
        shader->setUniform("projectionMatrix_", camera->getProjectionMatrix(screenDim));
        tgt::mat4 viewInvert;
        if(camera->getViewMatrix().invert(viewInvert))
            shader->setUniform("viewMatrixInverse_", viewInvert);
        tgt::mat4 projInvert;
        if(camera->getProjectionMatrix(screenDim).invert(projInvert))
            shader->setUniform("projectionMatrixInverse_", projInvert);
    }

    IMode.setMatstackUniforms(shader);

    shader->setIgnoreUniformLocationError(false);

    LGL_ERROR;
}

std::string RenderProcessor::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    return GLSL::generateStandardShaderHeader(version);
}

const std::vector<RenderPort*>& RenderProcessor::getPrivateRenderPorts() const {
    return privateRenderPorts_;
}

} // namespace voreen
