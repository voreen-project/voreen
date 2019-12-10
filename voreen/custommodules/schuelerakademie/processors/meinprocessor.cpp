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

#include "meinprocessor.h"

#include "voreen/core/interaction/camerainteractionhandler.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "tgt/glmath.h"
#include "tgt/gpucapabilities.h"
#include "tgt/texturemanager.h"
#include "tgt/shadermanager.h"
#include "tgt/textureunit.h"

using tgt::vec3;
using tgt::mat4;
using tgt::TextureUnit;

//#define VRN_CULLING

namespace voreen {

//-----------------------------------------------------------------------
//          Hier Code einfuegen                                           
//-----------------------------------------------------------------------
void MeinProcessor::zeichneWuerfel(tgt::vec3 llf, tgt::vec3 urb) {
    //blaues Dreieck
    glBegin(GL_TRIANGLES);
        glColor4f(0.f,0.f,1.f,1.f); glVertex3f(-1.f,-1.f,0.f);
        glColor4f(0.f,0.f,1.f,1.f); glVertex3f(1.f,-1.f,0.f);
        glColor4f(0.f,0.f,1.f,1.f); glVertex3f(0.f,1.f,0.f);
    glEnd();
}

//----------------------------------------------------------------------
//                  Vorgegebener Code                                   
//----------------------------------------------------------------------

const std::string MeinProcessor::loggerCat_("voreen.MeinProcessor");

MeinProcessor::MeinProcessor()
    : RenderProcessor()
    , entryPort_(Port::OUTPORT, "image.entrypoints", "Entry-points Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , exitPort_(Port::OUTPORT, "image.exitpoints", "Exit-points Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport_(Port::INPORT, "volumeinport", "Volume Input")
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)))
{

    addProperty(camera_);

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    addInteractionHandler(cameraHandler_);

    entryPort_.setDeinitializeOnDisconnect(false);
    addPort(entryPort_);
    exitPort_.setDeinitializeOnDisconnect(false);
    addPort(exitPort_);
    addPort(inport_);
}

MeinProcessor::~MeinProcessor() {
    delete cameraHandler_;
}

Processor* MeinProcessor::create() const {
    return new MeinProcessor();
}

bool MeinProcessor::isReady() const {
    // We want to render if at least one of the outports is connected
    return ((entryPort_.isReady() || exitPort_.isReady()) && inport_.isReady());
}

void MeinProcessor::beforeProcess() {
    RenderProcessor::beforeProcess();

    RenderPort& refPort = (entryPort_.isReady() ? entryPort_ : exitPort_);

    camera_.adaptInteractionToScene(inport_.getData()->getBoundingBox().getBoundingBox());

    if (refPort.getRenderTarget()->getColorTexture()->getType() == GL_FLOAT) {
        entryPort_.changeFormat(GL_RGBA16);
        exitPort_.changeFormat(GL_RGBA16);
    }
}

void MeinProcessor::process() {
    const VolumeBase* input = inport_.getData();

    // set modelview and projection matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadMatrix(camera_.get().getProjectionMatrix(entryPort_.isReady() ? entryPort_.getSize() : exitPort_.getSize()));

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadMatrix(camera_.get().getViewMatrix());

    // render front texture, use temporary target if necessary
    if (entryPort_.isReady())
       renderGeometry(input, entryPort_, GL_LESS, 1.0f, GL_BACK);


    // render back texture
    if (exitPort_.isReady())
        renderGeometry(input, exitPort_, GL_GREATER, 0.0f, GL_FRONT);

    // restore OpenGL state
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;

    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

void MeinProcessor::renderGeometry(const VolumeBase* volume, RenderPort& outport, GLenum depthFunc, float clearDepth, GLenum cullFace) {

    // enable culling
#ifdef VRN_CULLING
    glEnable(GL_CULL_FACE);
#endif
    outport.activateTarget();
    glClearDepth(clearDepth);
    glDepthFunc(depthFunc);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glCullFace(cullFace);

    zeichneWuerfel(volume->getLLF(), volume->getURB());
    LGL_ERROR;

    outport.deactivateTarget();

    glDepthFunc(GL_LESS);
    glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glClearDepth(1.0f);
    LGL_ERROR;
}


void MeinProcessor::adjustRenderOutportSizes() {
    tgt::ivec2 size = tgt::ivec2(-1);

    RenderSizeReceiveProperty* entrySizeProp = entryPort_.getSizeReceiveProperty();
    RenderSizeReceiveProperty* exitSizeProp = exitPort_.getSizeReceiveProperty();

    if(entrySizeProp->get() == exitSizeProp->get()) {
        size = entrySizeProp->get();
    }
    else {
        if(entryPort_.isConnected() && exitPort_.isConnected()) {
            size = entrySizeProp->get();
            LWARNING("Requested size for entry- and exit-point ports differ! Using size of entry-port");
            //TODO: Check for inbound links
        }
        else if(entryPort_.isConnected())
            size = entrySizeProp->get();
        else if(exitPort_.isConnected())
            size = exitSizeProp->get();
        else {
            //size = entrySizeProp->get();
        }
    }

    if(size != tgt::ivec2(-1)) {
        entryPort_.resize(size);
        exitPort_.resize(size);
    }
}

} // namespace voreen
