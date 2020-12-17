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

#include "geometryprocessor.h"
#include "voreen/core/interaction/idmanager.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/processors/geometryrendererbase.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "tgt/textureunit.h"
#include "tgt/glmath.h"
#include "tgt/vector.h"

using tgt::vec4;
using tgt::vec3;
using tgt::TextureUnit;

namespace voreen {

GeometryProcessor::GeometryProcessor()
    : RenderProcessor()
    , renderGeometries_("renderGeometries", "Render Geometries", true)
    , adaptToScene_("adaptToScene", "Adapt Camera to Scene", false)
    , forceAdaptToScene_("forceAdaptToScene", "Force Camera adaption")
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)))
    , applyOrderIndependentTransparency_("applyOrderIndependentTransparency_", "Enabled", true)
    , maxFragmentsPerPixel_("maxFragmentsPerPixelAvgAbs", "Global per Pixel Buffer Size\nLocal per Pixel Buffer Size", tgt::ivec2(15,20), 1, 50, 0, INT_MAX, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , inport_(Port::INPORT, "image.input", "Image Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , outport_(Port::OUTPORT, "image.output", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , tempPort_(Port::OUTPORT, "image.temp")
    , oitImageAddTargetPort_(Port::OUTPORT, "blank.temp")
    , pickingPort_(Port::OUTPORT, "pickingTarget")
    , composeShader_("composeShader", "Compose Shader", "image/compositor.frag", "passthrough.vert", "", Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , cpPort_(Port::INPORT, "coprocessor.geometryrenderers", "GeometryRenderers", true)
    , oirAddImageShader_("oirAddImageShader", "OIT Add Image Shader", "oit_addimage.frag", "oit_passthrough.vert", "", Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , oirBlendShader_("oirBlendShader", "OIT Blend Shader", "oit_blend.frag", "oit_passthrough.vert", "", Processor::INVALID_PROGRAM,Property::LOD_DEBUG)
    , headPointerTexture_(0)
    , headPointerInitializer_(0)
    , atomicFragmentCounterBuffer_(0)
    , fragmentStorageBuffer_(0)
    , fragmentStorageTexture_(0)
    , updateSceneAdaptation_(false)
{
    addProperty(renderGeometries_);
    addProperty(adaptToScene_);
        ON_CHANGE_LAMBDA(adaptToScene_, [this] {
            forceAdaptToScene_.setVisibleFlag(adaptToScene_.get());
        });
        adaptToScene_.setGroupID("camera");
    addProperty(forceAdaptToScene_);
        ON_CHANGE_LAMBDA(forceAdaptToScene_, [this] {
            adaptToScene(true);
        });
        forceAdaptToScene_.setVisibleFlag(false);
        forceAdaptToScene_.setGroupID("camera");
    addProperty(camera_);
        camera_.setGroupID("camera");
    setPropertyGroupGuiName("camera", "Camera Settings");

    addProperty(applyOrderIndependentTransparency_);
        applyOrderIndependentTransparency_.setGroupID("OIT");
        ON_CHANGE(applyOrderIndependentTransparency_, GeometryProcessor, adjustOITProperties);
    addProperty(maxFragmentsPerPixel_);
        maxFragmentsPerPixel_.setGroupID("OIT");
        ON_CHANGE(maxFragmentsPerPixel_, GeometryProcessor, maxFragmentsPerPixelChanged);
    setPropertyGroupGuiName("OIT", "Order independent transparency");

    addProperty(composeShader_);
    addProperty(oirAddImageShader_);
    addProperty(oirBlendShader_);

    cameraHandler_ = new CameraInteractionHandler("cameraHandler", "Camera", &camera_);
    addInteractionHandler(cameraHandler_);

    addPort(inport_);
    addPort(outport_);
    outport_.onSizeReceiveChange<GeometryProcessor>(this, &GeometryProcessor::passThroughSizeRequest);
    outport_.onSizeReceiveChange<GeometryProcessor>(this, &GeometryProcessor::setupOITBuffers);

    addPrivateRenderPort(tempPort_);
    addPrivateRenderPort(oitImageAddTargetPort_);
    addPrivateRenderPort(pickingPort_);
    addPort(cpPort_);
        cpPort_.onChange(LambdaFunctionCallback([this] { updateSceneAdaptation_ = true; }));
}

GeometryProcessor::~GeometryProcessor() {
    delete cameraHandler_;
}

void GeometryProcessor::initialize() {
    RenderProcessor::initialize();

    idManager_.setRenderTarget(pickingPort_.getRenderTarget());
    idManager_.initializeTarget();

    composeShader_.setHeader(generateHeader());
    composeShader_.rebuild();

    setupOITShaders();

    setupOITBuffers();
}

void GeometryProcessor::deinitialize() {
    // Delete OIT buffers
    if(headPointerTexture_) {
        glDeleteTextures(1, &headPointerTexture_);
    }
    if(headPointerInitializer_) {
        glDeleteBuffers(1, &headPointerInitializer_);
    }
    if(atomicFragmentCounterBuffer_) {
        glDeleteBuffers(1, &atomicFragmentCounterBuffer_);
    }
    if(fragmentStorageBuffer_) {
        glDeleteBuffers(1, &fragmentStorageBuffer_);
    }
    if(fragmentStorageTexture_) {
        glDeleteTextures(1, &fragmentStorageTexture_);
    }

    RenderProcessor::deinitialize();
}

Processor* GeometryProcessor::create() const {
    return new GeometryProcessor();
}

bool GeometryProcessor::isReady() const {
    return outport_.isReady();
}
void GeometryProcessor::process() {

    // Update enclosing Bounding Box and adapt to scene.
    if (updateSceneAdaptation_ && adaptToScene_.get())
        adaptToScene();

    // set modelview and projection matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(camera_.get().getProjectionMatrix(outport_.getSize()));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadMatrix(camera_.get().getViewMatrix());
    LGL_ERROR;

    //
    // Render geometry
    //

    // Decide which port to use for rendering
    RenderPort& geometryTargetPort = inport_.isReady() ? tempPort_ : outport_;

    std::vector<GeometryRendererBase*> portData = cpPort_.getConnectedProcessors();

    //
    // First render all opaque geometry
    //
    geometryTargetPort.activateTarget();
    geometryTargetPort.clearTarget();
    if(renderGeometries_.get()) {
        for (GeometryRendererBase* geomRenderer : portData) {
            if(geomRenderer->isReady() && (!geomRenderer->usesTransparency()
                        // If we don't use OIT: render all Geometry now.
                        || !applyOrderIndependentTransparency_.get())) {
                geomRenderer->setCamera(camera_.get());
                geomRenderer->setViewport(outport_.getSize());
                geomRenderer->render();
                LGL_ERROR;

                tgtAssert(tgt::Shader::getCurrentProgram() == 0, "Shader active after rendering geometry.");
            }
        }
    }
    geometryTargetPort.deactivateTarget();

    //
    // If enabled: Apply OIT-rendering techniques.
    //
    if(applyOrderIndependentTransparency_.get() && renderGeometries_.get()) {
        // First clear all fragments from OIT buffers:
        clearOITBuffers();

        // If inport is ready, write it to OITBuffer
        renderIntoOITBuffer(geometryTargetPort);

        if(inport_.isReady()) {
            renderIntoOITBuffer(inport_);
        }

        // Now render all transparent geometry
        // The geometryTargetPort still has the depth values of the opaque geometry stored, so that
        // we even omit hidden pixels into the OIT buffers.
        geometryTargetPort.activateTarget();

        // Bind OIT buffers
        glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, atomicFragmentCounterBuffer_);
        glBindImageTexture(0, headPointerTexture_, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32UI);
        glBindImageTexture(1, fragmentStorageTexture_, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32UI);

        // We want to use the depth values generated by rendering opaque geometry, but must not
        // overwrite them.
        glDepthMask(GL_FALSE);

        for (GeometryRendererBase* geomRenderer : portData) {
            if(geomRenderer->isReady() && geomRenderer->usesTransparency()) {
                geomRenderer->setCamera(camera_.get());
                geomRenderer->setViewport(outport_.getSize());
                geomRenderer->renderTransparent();
                LGL_ERROR;

                tgtAssert(tgt::Shader::getCurrentProgram() == 0, "Shader active after rendering geometry.");
            }
        }
        glDepthMask(GL_TRUE);

        geometryTargetPort.deactivateTarget();
    }

    //
    // render picking objects
    //
    idManager_.activateTarget(getID());
    idManager_.clearTarget();
    if (renderGeometries_.get()) {
        for (size_t i=0; i<portData.size(); i++) {
            GeometryRendererBase* geomRenderer = portData.at(i);
            if (geomRenderer->isReady()) {
                if (geomRenderer->getIDManager() != &idManager_)
                    geomRenderer->setIDManager(&idManager_);
                geomRenderer->renderPicking();
                LGL_ERROR;

                tgtAssert(tgt::Shader::getCurrentProgram() == 0, "Shader active after rendering picking.");
            }
        }
    }
    idManager_.deactivateTarget();

    // restore matrices
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    LGL_ERROR;

    if(applyOrderIndependentTransparency_.get() && renderGeometries_.get()) {
        // If we use OIT and have geometries to render (and thus cleared the and rendered into the OIT buffers),
        // we always have to generate the output image from OIT buffers:
        outport_.activateTarget();
        outport_.clearTarget();

        blendOITBuffer();

        outport_.deactivateTarget();

    }
    else if(inport_.isReady()) {
        // If we do not use OIT or do not display geometries, and have an input image,
        // we have to blend geometry and input image
        outport_.activateTarget();
        outport_.clearTarget();

        // In this case we used tempPort_ before as geometryTargetPort.
        compose(inport_, tempPort_);

        outport_.deactivateTarget();

    } // Otherwise we rendered all geometry into outport_ and are done now.


    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

std::string GeometryProcessor::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = RenderProcessor::generateHeader();
    header += "#define MODE_ALPHA_COMPOSITING\n";
    return header;
}

void GeometryProcessor::adjustOITProperties() {
    maxFragmentsPerPixel_.setVisibleFlag(applyOrderIndependentTransparency_.get());
}

void GeometryProcessor::blendOITBuffer() {
    tgt::Shader* shader = oirBlendShader_.getShader();
    if(!shader) {
        LWARNING("No blend shader");
        return;
    }


    shader->activate();

    glBindImageTexture(0, headPointerTexture_, 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32UI);
    glBindImageTexture(1, fragmentStorageTexture_, 0, GL_FALSE, 0, GL_READ_ONLY, GL_RGBA32UI);
    LGL_ERROR;

    renderQuad();
    LGL_ERROR;

    shader->deactivate();
    LGL_ERROR;

    // Get the number of fragments written to fragmentStorageBuffer
    GLuint fragmentsWritten = 0;
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicFragmentCounterBuffer_);
    glGetBufferSubData(GL_ATOMIC_COUNTER_BUFFER, 0, sizeof(fragmentsWritten), &fragmentsWritten);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);

    // Calculate the size of the fragmentStorageBuffer (number of fragments)
    size_t maxFragments = maxFragmentsPerPixel_.get().x * tgt::hmul(outport_.getSize());
    // Warn during this rendering pass the buffer has almost been used completely
    if(fragmentsWritten >= maxFragments) {
        LWARNING("Fragment list buffer utilization is at or above 100%!");
        LDEBUG("Written: " << fragmentsWritten << ", available: " << maxFragments);
    }
}


void GeometryProcessor::renderIntoOITBuffer(RenderPort& p) {
    tgt::Shader* shader = oirAddImageShader_.getShader();
    if(!shader) {
        LWARNING("No Image add shader");
        return;
    }
    oitImageAddTargetPort_.activateTarget();
    oitImageAddTargetPort_.clearTarget();
    glDisable(GL_DEPTH_TEST);
    glDepthMask(false);

    // Set texture input
    shader->activate();
    p.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
    shader->setUniform("imageToRender_", 0);
    shader->setUniform("imageToRenderDepth_", 1);

    // Set OIT buffers
    glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, atomicFragmentCounterBuffer_);
    glBindImageTexture(0, headPointerTexture_, 0, GL_FALSE, 0, GL_READ_WRITE, GL_R32UI);
    glBindImageTexture(1, fragmentStorageTexture_, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32UI);
    shader->setIgnoreUnsetUniform("headPointerImage_"); // disable warning for unset uniform.

    renderQuad();

    shader->deactivate();

    glDepthMask(true);
    glEnable(GL_DEPTH_TEST);
    oitImageAddTargetPort_.deactivateTarget();
    LGL_ERROR;

    //glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    //glMemoryBarrier(GL_ALL_BARRIER_BITS);
}


void GeometryProcessor::clearOITBuffers() {
    if(!outport_.isReady()) {
        //LWARNING("Clear: Outport not ready");
        return;
    }
    // Reset headPointerTexture with values from headPointerInitializer
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, headPointerInitializer_);
    glBindTexture(GL_TEXTURE_2D, headPointerTexture_);
    glTexImage2D(GL_TEXTURE_2D, 0, // 2D texture, first level
            GL_R32UI, // 32-bit GLuint per texel
            outport_.getSize().x, // Width
            outport_.getSize().y, // Height
            0, // No border
            GL_RED_INTEGER, // Single channel
            GL_UNSIGNED_INT, // Unsigned int
            NULL); // Consume data from PBO
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

    // Reset atomicFragmentCounterBuffer_
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicFragmentCounterBuffer_);
    const GLuint zero = 0;
    glBufferSubData(GL_ATOMIC_COUNTER_BUFFER, 0, sizeof(zero), &zero);
    LGL_ERROR;
}

void GeometryProcessor::setupOITBuffers() {
    if(!outport_.isReady()) {
        //LWARNING("Setup: Outport not ready");
        return;
    }
    GLuint* data;
    size_t total_pixels = tgt::hmul(outport_.getReceivedSize());
    if(total_pixels == 0) {
        LWARNING("Outport has zero size");
        return;
    }

    // Create the 2D image that will be used to store the head pointers
    // for the per-pixel linked lists.
    if(headPointerTexture_) {
        glDeleteTextures(1, &headPointerTexture_);
    }
    glGenTextures(1, &headPointerTexture_);
    glBindTexture(GL_TEXTURE_2D, headPointerTexture_);
    glTexImage2D(GL_TEXTURE_2D, 0, // 2D texture, level 0
        GL_R32UI, // 32-bit GLuint per texel
        outport_.getReceivedSize().x, // Width
        outport_.getReceivedSize().y, // Height
        0, //No border
        GL_RED_INTEGER, //Single channel
        GL_UNSIGNED_INT, //Unsigned int
        NULL); //No data... yet
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
    glBindTexture(GL_TEXTURE_2D, 0);

    // Create initialization buffer for headPointerTexture.
    // This way clearing headPointerTexture can be done without copying
    // values between CPU and GPU later.
    if(headPointerInitializer_) {
        glDeleteBuffers(1, &headPointerInitializer_);
    }
    glGenBuffers(1, &headPointerInitializer_);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, headPointerInitializer_);
    glBufferData(GL_PIXEL_UNPACK_BUFFER,
            total_pixels * sizeof(GLuint), // 1 uint per pixel
            NULL, // No data - we'll map it
            GL_STATIC_DRAW); // Never going to change
    data = (GLuint*)glMapBuffer(GL_PIXEL_UNPACK_BUFFER, GL_WRITE_ONLY);
    // 0xFFFFFFFF will be our "end of list" marker.
    memset(data, 0xFF, total_pixels * sizeof(GLuint));
    glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    LGL_ERROR;

    // Create the atomicFragmentCounterBuffer_
    if(atomicFragmentCounterBuffer_) {
        glDeleteBuffers(1, &atomicFragmentCounterBuffer_);
    }
    glGenBuffers(1, &atomicFragmentCounterBuffer_);
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, atomicFragmentCounterBuffer_);
    glBufferData(GL_ATOMIC_COUNTER_BUFFER, // Allocate buffer...
            sizeof(GLuint), NULL, // with space for 1 GLuint
            GL_DYNAMIC_COPY); // written to by GPU
    glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);
    LGL_ERROR;

    // Create the buffer for the fragment list.
    const size_t maxFragmentsPerPixelAvg = maxFragmentsPerPixel_.get().x;
    if(fragmentStorageBuffer_) {
        glDeleteBuffers(1, &fragmentStorageBuffer_);
    }
    glGenBuffers(1, &fragmentStorageBuffer_);
    glBindBuffer(GL_TEXTURE_BUFFER, fragmentStorageBuffer_);
    glBufferData(GL_TEXTURE_BUFFER,
            maxFragmentsPerPixelAvg * total_pixels * // N times the maximum number of pixels
            sizeof(vec4), // Times vec4
            NULL, // No data
            GL_DYNAMIC_COPY); // Updated often by GPU
    glBindBuffer(GL_TEXTURE_BUFFER, 0);
    LGL_ERROR;

    // Create the texture handle for fragmentStorageBuffer_
    if(fragmentStorageTexture_) {
        glDeleteTextures(1, &fragmentStorageTexture_);
    }
    glGenTextures(1, &fragmentStorageTexture_);
    glBindTexture(GL_TEXTURE_BUFFER, fragmentStorageTexture_);
    glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32UI, fragmentStorageBuffer_);
    glBindTexture(GL_TEXTURE_BUFFER, 0);
    LGL_ERROR;
}
void GeometryProcessor::setupOITShaders() {
    std::string maxFragmentsDefine =
        "#define MAX_FRAGMENTS " + std::to_string(maxFragmentsPerPixel_.get().y);

    oirBlendShader_.setHeader(maxFragmentsDefine);
    oirBlendShader_.rebuild();
    oirAddImageShader_.rebuild();
}

void GeometryProcessor::compose(RenderPort& r1, RenderPort& r2) {
    tgt::Shader* shader = composeShader_.getShader();
    if(!shader) {
        LWARNING("No compose shader");
        return;
    }

    r1.bindTextures(GL_TEXTURE0, GL_TEXTURE1);
    r2.bindTextures(GL_TEXTURE2, GL_TEXTURE3);

    shader->activate();

    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(shader, &cam);
    shader->setUniform("colorTex0_", 0);
    shader->setUniform("depthTex0_", 1);
    r1.setTextureParameters(shader, "textureParameters0_");

    shader->setUniform("colorTex1_", 2);
    shader->setUniform("depthTex1_", 3);
    r2.setTextureParameters(shader, "textureParameters1_");

    renderQuad();
    shader->deactivate();

    LGL_ERROR;
}
void GeometryProcessor::maxFragmentsPerPixelChanged() {
    setupOITBuffers();
    setupOITShaders();
}

void GeometryProcessor::passThroughSizeRequest() {
    inport_.requestSize(outport_.getReceivedSize());
}

void GeometryProcessor::adaptToScene(bool forced) {

    tgt::Bounds boundingBox; // Undefined as per definition.
    for (GeometryRendererBase* geomRenderer : cpPort_.getConnectedProcessors()) {
        const tgt::Bounds& box = geomRenderer->getBoundingBox();
        if (geomRenderer->isReady() && box.isDefined())
            boundingBox.addVolume(box);
    }

    // In case the bounding box is defined, adapt camera to scene.
    if (boundingBox.isDefined())
        camera_.adaptInteractionToScene(boundingBox, 0.0f, forced);

    // Update done.
    updateSceneAdaptation_ = false;
}

} // namespace voreen
