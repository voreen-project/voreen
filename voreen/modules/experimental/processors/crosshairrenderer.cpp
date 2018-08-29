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

#include "crosshairrenderer.h"

#include "tgt/textureunit.h"

using tgt::TextureUnit;

namespace voreen {

CrosshairRenderer::CrosshairRenderer()
    : ImageProcessor("image/compositor")
    , inport_(Port::INPORT, "image.in")
    , outport_(Port::OUTPORT, "image.out")
    , render_("render", "Render", true)
    , width_("width", "Width", 1.f, 0.1f, 10.f)
    , opacity_("opacity", "Opacity", 1.0f, 0.0f, 1.f)
    , color_("color", "Color", tgt::Color(1.f, 1.f, 1.f, 1.f))
    , position_("position", "Position", tgt::vec3(0.0f), tgt::vec3(-10.0f), tgt::vec3(10.0f))
    , camera_("camera", "Camera", tgt::Camera(tgt::vec3(0.f, 0.f, 3.5f), tgt::vec3(0.f, 0.f, 0.f), tgt::vec3(0.f, 1.f, 0.f)))
    , copyShader_(0)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(opacity_);
    addProperty(render_);
    addProperty(width_);
    addProperty(color_);
    addProperty(position_);
    addProperty(camera_);
}

Processor* CrosshairRenderer::create() const {
    return new CrosshairRenderer();
}

void CrosshairRenderer::initialize() {
    ImageProcessor::initialize();

    try {
        copyShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
            ImageProcessor::generateHeader(), false);
        copyShader_->deactivate();
    } catch(tgt::Exception) {
        processorState_ = PROCESSOR_STATE_NOT_INITIALIZED;
        throw VoreenException("Failed to load shader: passthrough.vert/copyimage.frag");
    }
}

void CrosshairRenderer::deinitialize() {
    ShdrMgr.dispose(copyShader_);
    copyShader_ = 0;

    ImageProcessor::deinitialize();
}

void CrosshairRenderer::process() {
    tgtAssert(outport_.isReady(), "Outport not ready");
    tgtAssert(inport_.isReady(), "Inport not ready");

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // bind input image to tex unit
    TextureUnit imageUnit, imageUnitDepth;
    inport_.bindTextures(imageUnit.getEnum(), imageUnitDepth.getEnum());

    // 1. copy input image to outport
    copyShader_->activate();
    setGlobalShaderParameters(copyShader_);
    inport_.setTextureParameters(copyShader_, "texParams_");
    copyShader_->setUniform("colorTex_", imageUnit.getUnitNumber());
    copyShader_->setUniform("depthTex_", imageUnitDepth.getUnitNumber());
    renderQuad();
    copyShader_->deactivate();
    LGL_ERROR;

    // render border around overlay
    if (render_.get()) {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glColor4f(color_.get().x, color_.get().y, color_.get().z, opacity_.get());
        glLineWidth(width_.get());
        glDepthFunc(GL_ALWAYS);
        glEnable(GL_BLEND);
        tgt::vec2 c = camera_.get().project(outport_.getSize(), position_.get()).xy();
        c /= tgt::vec2(static_cast<float>(outport_.getSize().x), static_cast<float>(outport_.getSize().y));
        c *= 2.0f;
        c -= 1.0f;

        glBegin(GL_LINES);
        glVertex2f(-1.0f, c.y);
        glVertex2f(1.0f, c.y);
        glVertex2f(c.x, 1.0f);
        glVertex2f(c.x, -1.0f);
        glEnd();
        glPopAttrib();
        glDisable(GL_BLEND);
        LGL_ERROR;
    }

    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

} // namespace voreen
