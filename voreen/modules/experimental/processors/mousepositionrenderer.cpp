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

#include "mousepositionrenderer.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/textureunit.h"
#include "tgt/event/timeevent.h"
#include "tgt/event/touchevent.h"
#include "tgt/event/touchpoint.h"
#include "tgt/timer.h"

using tgt::TextureUnit;
using tgt::vec2;
using tgt::TouchPoint;

namespace voreen {

MousePositionRenderer::MousePositionRenderer()
    : ImageProcessor("image/compositor")
    , inport_(Port::INPORT, "image.in")
    , outport_(Port::OUTPORT, "image.out")
    , renderCursor_("renderCursor", "Render Cursor", true)
    , renderTrail_("renderTrail", "Render Trail", true)
    , width_("width", "Width", 1.f, 0.1f, 10.f)
    , opacity_("opacity", "Opacity", 1.0f, 0.0f, 1.f)
    , color_("color", "Color", tgt::Color(1.f, 1.f, 1.f, 1.f))
    , cursorSize_("cursorsize", "Cursor Size", 0.05f, 0.001f, 1.0f)
    , isPressed_(false)
    , isInside_(false)
    , timer_(0)
    , eventHandler_()
    , copyShader_(0)
{
    eventHandler_.addListenerToBack(this);

    addPort(inport_);
    addPort(outport_);

    addProperty(opacity_);
    addProperty(renderCursor_);
    addProperty(renderTrail_);
    addProperty(width_);
    addProperty(color_);
    addProperty(cursorSize_);
}

Processor* MousePositionRenderer::create() const {
    return new MousePositionRenderer();
}

void MousePositionRenderer::initialize() {
    ImageProcessor::initialize();

    try {
        copyShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
            ImageProcessor::generateHeader(), false);
        copyShader_->deactivate();
    } catch(tgt::Exception) {
        processorState_ = PROCESSOR_STATE_NOT_INITIALIZED;
        throw VoreenException("Failed to load shader: passthrough.vert/copyimage.frag");
    }

    timer_ = VoreenApplication::app()->createTimer(&eventHandler_);
    if (timer_)
        timer_->start(30);
}

void MousePositionRenderer::deinitialize() {
    ShdrMgr.dispose(copyShader_);
    copyShader_ = 0;

    delete timer_;
    timer_ = 0;

    ImageProcessor::deinitialize();
}

void MousePositionRenderer::process() {
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

    if(renderCursor_.get() || renderTrail_.get()) {
        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glColor4f(color_.get().x, color_.get().y, color_.get().z, opacity_.get());
        glLineWidth(width_.get());
        glDepthFunc(GL_ALWAYS);
        glEnable(GL_BLEND);

        if(renderTrail_.get()) {
            //render mousetrail or trail for first touch point
            glColor4f(0.0f, 0.0f, 1.0f, 0.7f);
            glLineWidth(2.0f);
            glBegin(GL_LINE_STRIP);
            for(std::deque< std::pair<tgt::vec2, int> >::iterator i = history_.begin(); i != history_.end(); i++) {
                vec2 c = (*i).first;
                glVertex2f(c.x, c.y);
            }
            glEnd();
            //render trail for second touch point
            glBegin(GL_LINE_STRIP);
            for(std::deque< std::pair<tgt::vec2, int> >::iterator i = history2_.begin(); i != history2_.end(); i++) {
                vec2 c = (*i).first;
                glVertex2f(c.x, c.y);
            }
            glEnd();
            glLineWidth(1.0f);
        }

        if(isInside_) {
            if (renderCursor_.get()) {
                //render mouse cursor:
                float ratio = (float) outport_.getSize().x / (float) outport_.getSize().y;
                float wmul;
                if(ratio < 0.0f)
                    wmul = 1.0f / ratio;
                else
                    wmul = ratio;

                MatStack.pushMatrix();
                MatStack.translate(mpos_.x, mpos_.y, 0.0f);
                MatStack.scale(cursorSize_.get(), cursorSize_.get() * wmul, 1.0f);

                vec2 a(0.7f, -0.7f);
                vec2 b(0.0f, -1.0f);

                float os = 0.1f;
                vec2 c = (0.5f + os) * a + (1.0f - (0.5f + os)) * b;
                vec2 d = (0.5f - os) * a + (1.0f - (0.5f - os)) * b;

                vec2 e = c + vec2(0.1f, -0.2f);
                vec2 f = d + vec2(0.1f, -0.2f);

                if(isPressed_)
                    glColor3ub(255,0,0);
                else
                    glColor3ub(255,255,255);
                glBegin(GL_POLYGON);
                glVertex2f(0.0f, 0.0f);

                glVertex2f(b.x, b.y);

                glVertex2f(d.x, d.y);
                glVertex2f(f.x, f.y);

                glVertex2f(e.x, e.y);
                glVertex2f(c.x, c.y);

                glVertex2f(a.x, a.y);
                glEnd();

                glLineWidth(2.0f);
                glColor3ub(0,0,0);
                glBegin(GL_LINE_LOOP);
                glVertex2f(a.x, a.y);

                glVertex2f(c.x, c.y);
                glVertex2f(e.x, e.y);

                glVertex2f(f.x, f.y);
                glVertex2f(d.x, d.y);

                glVertex2f(b.x, b.y);
                glVertex2f(0.0f, 0.0f);
                glEnd();
                glLineWidth(1.0f);

                MatStack.popMatrix();
            }
        }

        glPopAttrib();
        glDisable(GL_BLEND);
        LGL_ERROR;
    }

    outport_.deactivateTarget();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

void MousePositionRenderer::onEvent(tgt::Event* e) {

    //react on event if it is a touch event
    tgt::TouchEvent* touche = dynamic_cast<tgt::TouchEvent*>(e);

    if(touche) {
        if(touche->touchPointStates() & TouchPoint::TouchPointPressed) {
            isPressed_ = true;
            //history_.clear();
            //history2_.clear();
        }
        if (touche->touchPointStates() & TouchPoint::TouchPointReleased) {
            isPressed_ = false;
        }
        if (touche->touchPointStates() & TouchPoint::TouchPointMoved) {

            mpos_ = touche->touchPoints()[0].pos();
            mpos2_ = touche->touchPoints()[1].pos();

            mpos_.y = outport_.getSize().y - mpos_.y;
            mpos2_.y = outport_.getSize().y - mpos2_.y;

            mpos_ /= tgt::vec2(static_cast<float>(outport_.getSize().x), static_cast<float>(outport_.getSize().y));
            mpos2_ /= tgt::vec2(static_cast<float>(outport_.getSize().x), static_cast<float>(outport_.getSize().y));
            mpos_ *= 2.0f;
            mpos2_ *= 2.f;
            mpos2_ -= 1.f;
            mpos_ -= 1.0f;

            if(timer_) {
                if(isPressed_) {
                    history_.push_front(std::pair< tgt::vec2, int>(mpos_, timer_->getCount()));
                    history2_.push_front(std::pair< tgt::vec2, int>(mpos2_, timer_->getCount()));
                }
            }
            invalidate();
        }
        e->ignore();
        RenderProcessor::onEvent(e);
        return;
    }

    //react on event if it is a mouse event
    tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);

    if(me) {
        if((me->action() & tgt::MouseEvent::EXIT) == tgt::MouseEvent::EXIT)
            isInside_ = false;
        if((me->action() & tgt::MouseEvent::ENTER) == tgt::MouseEvent::ENTER)
            isInside_ = true;

        if((me->action() & tgt::MouseEvent::PRESSED) == tgt::MouseEvent::PRESSED) {
            isPressed_ = true;
            history_.clear();
            history2_.clear();
        }
        if((me->action() & tgt::MouseEvent::RELEASED) == tgt::MouseEvent::RELEASED)
            isPressed_ = false;

        if((me->action() & tgt::MouseEvent::MOTION) == tgt::MouseEvent::MOTION) {
            mpos_ = me->coord();
            mpos_.y = outport_.getSize().y - mpos_.y;

            mpos_ /= tgt::vec2(static_cast<float>(outport_.getSize().x), static_cast<float>(outport_.getSize().y));
            mpos_ *= 2.0f;
            mpos_ -= 1.0f;

            if(timer_) {
                if(isPressed_)
                    history_.push_front(std::pair< tgt::vec2, int>(mpos_, timer_->getCount()));
            }
            invalidate();
        }
        e->ignore();
        RenderProcessor::onEvent(e);
        return;
    }

    //react on event if it is a time event
    tgt::TimeEvent* te = dynamic_cast<tgt::TimeEvent*>(e);
    if(te) {
        if(te->getTimer() == timer_) {
            te->accept();
            while(!history_.empty() && ((timer_->getCount() - history_.back().second) > 20))
                history_.pop_back();
            while(!history2_.empty() && ((timer_->getCount() - history2_.back().second) > 20))
                history2_.pop_back();


            invalidate();
            return;
        }
    }
    RenderProcessor::onEvent(e);
}
} // namespace voreen
