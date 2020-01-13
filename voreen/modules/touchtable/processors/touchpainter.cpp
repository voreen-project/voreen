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

#include "touchpainter.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/textureunit.h"
#include "tgt/event/timeevent.h"
#include "tgt/timer.h"

#include <deque>
#include <iterator>
#include <algorithm>
#include <time.h>

using tgt::TextureUnit;
using tgt::vec2;
using tgt::TouchPoint;

namespace voreen {

TouchPainter::TouchPainter()
    : ImageProcessor("image/compositor")
    , inport_(Port::INPORT, "image.in")
    , outport_(Port::OUTPORT, "image.out", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    , randomizeColors_("randomColors", "Randomize Colors", false)
    , touchRadius_("touchRadius", "Radius", 1.f, 0.1f, 100.f)
    , color_("color", "Color", tgt::Color(1.f, 1.f, 1.f, 1.f))
    , opacity_("opacity", "Opacity", 1.f, 0.f, 1.f)
    , clearImage_("clearImage", "Clear Image")
    , copyShader_(0)
{
    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(privatePort_);

    addProperty(randomizeColors_);
    addProperty(touchRadius_);
    addProperty(color_);
    addProperty(opacity_);
    addProperty(clearImage_);

    clearImage_.onChange(MemberFunctionCallback<TouchPainter>(this, &TouchPainter::onClear));
    std::srand((unsigned)time(NULL));
}

bool TouchPainter::isReady() const {
    return outport_.isReady();
}

void TouchPainter::onClear() {
    privatePort_.activateTarget();
    privatePort_.clearTarget();
    privatePort_.deactivateTarget();
    invalidate();
}

void TouchPainter::initialize() {
    ImageProcessor::initialize();

    try {
        copyShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
            ImageProcessor::generateHeader(), false);
        copyShader_->deactivate();
    } catch(tgt::Exception) {
        processorState_ = PROCESSOR_STATE_NOT_INITIALIZED;
        throw VoreenException("Failed to load shader: passthrough.vert/copyimage.frag");
    }

    onClear();
}

void TouchPainter::deinitialize() {
    ShdrMgr.dispose(copyShader_);
    ImageProcessor::deinitialize();
}

void TouchPainter::process() {

    if(privatePort_.getSize() != outport_.getSize())
        privatePort_.resize(outport_.getSize());

    privatePort_.activateTarget();

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    if(!randomizeColors_.get())
        glColor4f(color_.get().x, color_.get().y, color_.get().z, opacity_.get());
    glPointSize(touchRadius_.get());
    glDepthFunc(GL_ALWAYS);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH );
    glBegin(GL_POINTS);
    while(!toPaint_.empty()) {
        tgt::TouchPoint tp = toPaint_.front();
        toPaint_.pop_front();

        if(randomizeColors_.get()) {
            tgt::vec3 col = colorMap_[tp.id()];
            glColor4f(col.x, col.y, col.z, opacity_.get());
        }

        tgt::vec2 pos = tp.pos();
        pos.y = outport_.getSize().y - pos.y;
        pos /= tgt::vec2(static_cast<float>(outport_.getSize().x), static_cast<float>(outport_.getSize().y));
        pos *= 2.0f;
        pos -= 1.0f;
        glVertex2f(pos.x, pos.y);

        if(tp.state() & tgt::TouchPoint::TouchPointReleased)
            colorMap_.erase(tp.id());
    }
    glEnd();
    glPopAttrib();

    privatePort_.deactivateTarget();

    outport_.activateTarget();
    outport_.clearTarget();

    if(inport_.isReady()) {
        // blend painted image over incoming image
        TextureUnit imageUnit, imageUnitDepth;
        inport_.bindTextures(imageUnit.getEnum(), imageUnitDepth.getEnum());

        TextureUnit overlayUnit;
        tgt::Texture* overlayTex = privatePort_.getColorTexture();
        overlayUnit.activate();
        overlayTex->bind();

        program_->activate();
        setGlobalShaderParameters(program_);

        // image texture parameters
        outport_.setTextureParameters(program_, "textureParameters0_");
        outport_.setTextureParameters(program_, "textureParameters1_");
        program_->setUniform("colorTex0_", overlayUnit.getUnitNumber());
        program_->setUniform("colorTex1_", imageUnit.getUnitNumber());
        program_->setUniform("depthTex1_", imageUnitDepth.getUnitNumber());

        glDepthFunc(GL_ALWAYS);
        renderQuad();
        glDepthFunc(GL_LESS);

        program_->deactivate();
        LGL_ERROR;
    } else {
        // copy painted image to outport
        TextureUnit imageUnit;
        privatePort_.bindColorTexture(imageUnit.getEnum());

        copyShader_->activate();
        setGlobalShaderParameters(copyShader_);
        privatePort_.setTextureParameters(copyShader_, "texParams_");
        copyShader_->setUniform("colorTex_", imageUnit.getUnitNumber());
        renderQuad();
        copyShader_->deactivate();
        LGL_ERROR;
    }

    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

void TouchPainter::touchEvent(tgt::TouchEvent* e) {

    e->accept();

    for(std::deque<tgt::TouchPoint>::const_iterator it = e->touchPoints().begin(); it != e->touchPoints().end(); it++) {
        toPaint_.push_back(*it);

        if(it->state() == tgt::TouchPoint::TouchPointPressed)
            colorMap_[it->id()] = tgt::vec3(getRandomFloat(), getRandomFloat(), getRandomFloat());
    }

    invalidate();
}

void TouchPainter::mouseDoubleClickEvent(tgt::MouseEvent* e) {
    onClear();
    e->accept();
}

std::string TouchPainter::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = ImageProcessor::generateHeader(version);
    header += "#define MODE_ALPHA_BLENDING\n";
    return header;
}

} // namespace voreen
