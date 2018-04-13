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

#include "meshentryexitpoints.h"

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

namespace voreen {

const std::string MeshEntryExitPoints::loggerCat_("voreen.MeshEntryExitPoints");

MeshEntryExitPoints::MeshEntryExitPoints()
    : RenderProcessor()
    , entryPort_(Port::OUTPORT, "image.entrypoints", "Entry-points Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , exitPort_(Port::OUTPORT, "image.exitpoints", "Exit-points Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport_(Port::INPORT, "proxgeometry.geometry", "Proxy Geometry Input")
    , tmpPort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    , outputCoordinateSystem_("outputCoordinateSystem", "Output coordinate system")
    , cameraInsideVolumeTechnique_("cameraInsideVolumeTechnique", "Camera inside volume technique")
    , jitterEntryPoints_("jitterEntryPoints", "Jitter entry params", false)
    , useCulling_("useCulling", "Use culling", true)
    , jitterStepLength_("jitterStepLength", "Jitter step length", 0.005f, 0.0005f, 0.025f)
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)))
    , shaderProgram_(0)
    , shaderProgramJitter_(0)
    , shaderProgramInsideVolume_(0)
{
    outputCoordinateSystem_.addOption("texture", "Texture Coordinates");
    outputCoordinateSystem_.addOption("world", "World Coordinates");
    outputCoordinateSystem_.select("texture");
    addProperty(outputCoordinateSystem_);

    cameraInsideVolumeTechnique_.addOption("none", "None");
    cameraInsideVolumeTechnique_.addOption("cpu", "CPU (slow)");
    cameraInsideVolumeTechnique_.addOption("gpu", "GPU");
    cameraInsideVolumeTechnique_.select("gpu");
    addProperty(cameraInsideVolumeTechnique_);

    // jittering
    addProperty(jitterEntryPoints_);
    jitterEntryPoints_.onChange(MemberFunctionCallback<MeshEntryExitPoints>(this, &MeshEntryExitPoints::onJitterEntryPointsChanged));
    onJitterEntryPointsChanged(); // set initial state
    jitterStepLength_.setStepping(0.001f);
    jitterStepLength_.setNumDecimals(3);
    addProperty(jitterStepLength_);
    addProperty(useCulling_);

    addProperty(camera_);

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    addInteractionHandler(cameraHandler_);

    entryPort_.setDeinitializeOnDisconnect(false);
    addPort(entryPort_);
    exitPort_.setDeinitializeOnDisconnect(false);
    addPort(exitPort_);
    addPort(inport_);
    addPrivateRenderPort(&tmpPort_);
}

MeshEntryExitPoints::~MeshEntryExitPoints() {
    delete cameraHandler_;
}

Processor* MeshEntryExitPoints::create() const {
    return new MeshEntryExitPoints();
}

void MeshEntryExitPoints::initialize() {
    RenderProcessor::initialize();

    shaderProgram_ = ShdrMgr.load("eep_simple", generateHeader(), false);
    shaderProgramJitter_ = ShdrMgr.loadSeparate("passthrough.vert", "eep_jitter.frag", generateHeader(), false);
    shaderProgramInsideVolume_ = ShdrMgr.loadSeparate("passthrough.vert", "eep_insidevolume.frag", generateHeader(), false);
}

void MeshEntryExitPoints::deinitialize() {
    ShdrMgr.dispose(shaderProgram_);
    shaderProgram_ = 0;

    ShdrMgr.dispose(shaderProgramJitter_);
    shaderProgramJitter_ = 0;

    ShdrMgr.dispose(shaderProgramInsideVolume_);
    shaderProgramInsideVolume_ = 0;

    RenderProcessor::deinitialize();
}

bool MeshEntryExitPoints::isReady() const {
    // We want to render if at least one of the outports is connected
    return ((entryPort_.isReady() || exitPort_.isReady()) && inport_.isReady());
}

void MeshEntryExitPoints::beforeProcess() {
    RenderProcessor::beforeProcess();

    RenderPort& refPort = (entryPort_.isReady() ? entryPort_ : exitPort_);

    if(outputCoordinateSystem_.isSelected("world")) {
        if (refPort.getRenderTarget()->getColorTexture()->getGLDataType() != GL_FLOAT) {
            entryPort_.changeFormat(GL_RGBA16F);
            exitPort_.changeFormat(GL_RGBA16F);
            tmpPort_.changeFormat(GL_RGBA16F);
        }
    }
    else {
        if (refPort.getRenderTarget()->getColorTexture()->getGLDataType() == GL_FLOAT) {
            entryPort_.changeFormat(GL_RGBA16);
            exitPort_.changeFormat(GL_RGBA16);
            tmpPort_.changeFormat(GL_RGBA16);
        }
    }
}

void MeshEntryExitPoints::renderGeometry(const Geometry* geometry, RenderPort& outport, GLenum depthFunc, float clearDepth, GLenum cullFace) {
    // activate shader program
    shaderProgram_->activate();
    if(outputCoordinateSystem_.isSelected("texture"))
        shaderProgram_->setUniform("useTextureCoordinates_", true);
    else
        shaderProgram_->setUniform("useTextureCoordinates_", false);

    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(shaderProgram_, &cam);
    LGL_ERROR;

    // enable culling
    if(useCulling_.get())
        glEnable(GL_CULL_FACE);

    outport.activateTarget();
    glClearDepth(clearDepth);
    glDepthFunc(depthFunc);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glCullFace(cullFace);

    geometry->render();
    LGL_ERROR;

    outport.deactivateTarget();

    // deactivate shader program
    shaderProgram_->deactivate();

    glDepthFunc(GL_LESS);
    glDisable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glClearDepth(1.0f);
    LGL_ERROR;
}

void MeshEntryExitPoints::process() {
    const Geometry* input = inport_.getData();
    bool deleteInput = false;

    if (cameraInsideVolumeTechnique_.isSelected("cpu") ) {
        tgt::Bounds b = input->getBoundingBox();

        // clip proxy geometry against near-plane
        float nearPlaneDistToOrigin = tgt::dot(-(camera_.get().getPosition() + camera_.get().getLook() * camera_.get().getNearDist()), camera_.get().getLook())/* - camera_.get().getNearDist()*/ /* - 0.001f*/;
        double epsilon = static_cast<double>(tgt::length(input->getBoundingBox().diagonal())) * 1e-6;

        if(b.intersects(tgt::plane(-camera_.get().getLook(), nearPlaneDistToOrigin))) {
            Geometry* copy = input->clone().release();
            copy->clip(tgt::plane(-camera_.get().getLook(), nearPlaneDistToOrigin), epsilon);
            deleteInput = true;
            input = copy;
        }
    }

    // set modelview and projection matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadMatrix(camera_.get().getProjectionMatrix(entryPort_.isReady() ? entryPort_.getSize() : exitPort_.getSize()));

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadMatrix(camera_.get().getViewMatrix());

    if(const TriangleMeshGeometryBase* tmgb = dynamic_cast<const TriangleMeshGeometryBase*>(inport_.getData())) {
        MatStack.multMatrix(tmgb->getTransformationMatrix());
    }

    if (cameraInsideVolumeTechnique_.isSelected("gpu")) {
        if (jitterEntryPoints_.get()) {
            renderGeometry(input, exitPort_, GL_LESS, 1.0f, GL_BACK); // render first front face
            renderGeometry(input, entryPort_, GL_LESS, 1.0f, GL_FRONT); // render first back face
            fillEntryPoints(entryPort_, exitPort_, tmpPort_, input);
        }
        else {
            renderGeometry(input, tmpPort_, GL_LESS, 1.0f, GL_BACK); // render first front face
            renderGeometry(input, exitPort_, GL_LESS, 1.0f, GL_FRONT); // render first back face
            fillEntryPoints(exitPort_, tmpPort_, entryPort_, input);
        }
    }
    else {
        // render front texture, use temporary target if necessary
        if (entryPort_.isReady()) {
            if (cameraInsideVolumeTechnique_.isSelected("cpu"))
                glEnable(GL_DEPTH_CLAMP);

            if (!jitterEntryPoints_.get())
                renderGeometry(input, entryPort_, GL_LESS, 1.0f, GL_BACK);
            else
                renderGeometry(input, tmpPort_, GL_LESS, 1.0f, GL_BACK);

            glDisable(GL_DEPTH_CLAMP);
            LGL_ERROR;
        }
    }

    // render back texture
    if (exitPort_.isReady())
        renderGeometry(input, exitPort_, GL_GREATER, 0.0f, GL_FRONT);

    if (entryPort_.isReady()) {
        // jittering of entry points
        if (jitterEntryPoints_.get())
            jitterEntryPoints(tmpPort_, exitPort_, entryPort_);
    }

    // restore OpenGL state
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;

    TextureUnit::setZeroUnit();
    LGL_ERROR;

    if(deleteInput)
        delete input;
}

void MeshEntryExitPoints::jitterEntryPoints(RenderPort& entryPort, RenderPort& exitPort, RenderPort& outport) {
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    shaderProgramJitter_->activate();
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(shaderProgramJitter_, &cam);

    // bind texture
    TextureUnit entryUnit, entryDepthUnit, exitUnit;

    // bind entry points texture and depth texture (have been rendered to temporary port)
    entryPort.bindColorTexture(entryUnit.getEnum());
    shaderProgramJitter_->setUniform("entryPoints_", entryUnit.getUnitNumber());

    entryPort.bindDepthTexture(entryDepthUnit.getEnum());
    shaderProgramJitter_->setUniform("entryPointsDepth_", entryDepthUnit.getUnitNumber());
    entryPort.setTextureParameters(shaderProgramJitter_, "entryParameters_");

    // bind exit points texture
    exitPort.bindColorTexture(exitUnit.getEnum());
    shaderProgramJitter_->setUniform("exitPoints_", exitUnit.getUnitNumber());
    exitPort.setTextureParameters(shaderProgramJitter_, "exitParameters_");

    shaderProgramJitter_->setUniform("stepLength_", jitterStepLength_.get());

    outport.activateTarget("jitteredEntryParams");
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // render screen aligned quad
    renderQuad();

    shaderProgramJitter_->deactivate();
    outport.deactivateTarget();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;
}

void MeshEntryExitPoints::fillEntryPoints(RenderPort& firstBackPort, RenderPort& firstFrontPort, RenderPort& outport, const Geometry* geometry) {
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    shaderProgramInsideVolume_->activate();
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(shaderProgramInsideVolume_, &cam);

    // bind texture
    TextureUnit firstFrontUnit, firstFrontDepthUnit, firstBackUnit, firstBackDepthUnit;

    // bind firstFront points texture and depth texture
    firstFrontPort.bindColorTexture(firstFrontUnit.getEnum());
    shaderProgramInsideVolume_->setUniform("firstFront_", firstFrontUnit.getUnitNumber());
    firstFrontPort.bindDepthTexture(firstFrontDepthUnit.getEnum());
    shaderProgramInsideVolume_->setUniform("firstFrontDepth_", firstFrontDepthUnit.getUnitNumber());
    firstFrontPort.setTextureParameters(shaderProgramInsideVolume_, "firstFrontParameters_");

    // bind firstBack points texture
    firstBackPort.bindColorTexture(firstBackUnit.getEnum());
    shaderProgramInsideVolume_->setUniform("firstBack_", firstBackUnit.getUnitNumber());
    //firstBackPort.bindDepthTexture(firstBackDepthUnit.getEnum());
    //shaderProgramInsideVolume_->setUniform("firstBackDepth_", firstBackDepthUnit.getUnitNumber());
    firstBackPort.setTextureParameters(shaderProgramInsideVolume_, "firstBackParameters_");

    tgt::Bounds bounds = geometry->getBoundingBox(outputCoordinateSystem_.isSelected("world"));
    shaderProgramInsideVolume_->setUniform("llf_", bounds.getLLF());
    shaderProgramInsideVolume_->setUniform("urb_", bounds.getURB());

    if(outputCoordinateSystem_.isSelected("texture")) {
        shaderProgramInsideVolume_->setUniform("useTextureCoordinates_", true);
        mat4 worldToTextureMatrix;
        geometry->getTransformationMatrix().invert(worldToTextureMatrix);
        shaderProgramInsideVolume_->setUniform("worldToTexture_", worldToTextureMatrix);
        // TODO: does not work for GeometrySequences, which are not to be expected when using texture coordinates (=> single volume raycasting)
    }
    else {
        shaderProgramInsideVolume_->setUniform("useTextureCoordinates_", false);
        shaderProgramInsideVolume_->setUniform("worldToTexture_", tgt::mat4::identity);
    }

    outport.activateTarget("fillEntryParams");
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // render screen aligned quad
    renderQuad();

    shaderProgramInsideVolume_->deactivate();
    outport.deactivateTarget();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;
}

void MeshEntryExitPoints::onJitterEntryPointsChanged() {
    jitterStepLength_.setVisibleFlag(jitterEntryPoints_.get());
}

void MeshEntryExitPoints::adjustRenderOutportSizes() {
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
        tmpPort_.resize(size);
    }
}

} // namespace voreen
