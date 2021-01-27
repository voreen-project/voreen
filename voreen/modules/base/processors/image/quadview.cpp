/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "quadview.h"

namespace voreen {

QuadView::QuadView()
    : RenderProcessor()
    , showGrid_("showGrid", "Show grid", true)
    , gridColor_("gridColor", "Grid color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    , maximized_("maximized", "Maximized sub-view", 0, 0, 4)
    , maximizeOnDoubleClick_("maximizeOnDoubleClick", "Maximize on double click", true)
    , maximizeEventProp_("mouseEvent.maximize", "Maximize Event", this, &QuadView::toggleMaximization,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT, tgt::MouseEvent::DOUBLECLICK, tgt::MouseEvent::MODIFIER_NONE)
    , mouseMoveEventProp_("mouseEvent.move", "Move Event", this, &QuadView::mouseMove,
        tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::MOTION | tgt::MouseEvent::CLICK | tgt::MouseEvent::ENTER_EXIT, tgt::MouseEvent::MODIFIER_NONE)
    , forwardDoubleClickEvent_("forwardDoubleClickEvent", "Forward double click events", false, Processor::INVALID_RESULT, Property::LOD_DEBUG)
    , outport_(Port::OUTPORT, "outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport1_(Port::INPORT, "inport1", "Image1 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , inport2_(Port::INPORT, "inport2", "Image2 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , inport3_(Port::INPORT, "inport3", "Image3 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , inport4_(Port::INPORT, "inport4", "Image4 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , currentPort_(-1)
    , isDragging_(false)
{
    addProperty(showGrid_);
    addProperty(gridColor_);
    addProperty(maximized_);
    maximized_.setVisibleFlag(false);
    addProperty(maximizeOnDoubleClick_);
    addEventProperty(maximizeEventProp_);
    addEventProperty(mouseMoveEventProp_);
    addProperty(forwardDoubleClickEvent_);

    addPort(outport_);
    addPort(inport1_);
    addPort(inport2_);
    addPort(inport3_);
    addPort(inport4_);

    outport_.onSizeReceiveChange<QuadView>(this, &QuadView::portSizeReceiveChanged);
}

QuadView::~QuadView() {
}

Processor* QuadView::create() const {
    return new QuadView();
}

bool QuadView::isReady() const {
    if (!outport_.isReady())
        return false;

    if (!inport1_.isReady() && !inport2_.isReady() && !inport3_.isReady() && !inport4_.isReady())
        return false;

    if(maximized_.get() != 0) {
        switch(maximized_.get()) {
            case 1: if(!inport1_.isReady())
                        return false;
                    break;
            case 2: if(!inport2_.isReady())
                        return false;
                    break;
            case 3: if(!inport3_.isReady())
                        return false;
                    break;
            case 4: if(!inport4_.isReady())
                        return false;
                    break;
        }
    }
    return true;
}

void QuadView::process() {
    shader_->activate();
    shader_->setUniform("useTexcoords", true);

    // Set default values to unused textures to prevent undefined behavior.
    shader_->setUniform("colorTex_", 0);
    shader_->setUniform("depthTex_", 0);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    if (maximized_.get() == 0) {
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        outport_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (inport1_.isReady()) {
            inport1_.bindColorTexture(GL_TEXTURE0);
            inport1_.getColorTexture()->enable();
            MatStack.translate(-0.5f, 0.5f, 0.0f);
            MatStack.scale(0.5f, 0.5f, 1.0f);

            glDepthFunc(GL_ALWAYS);
            setGlobalShaderParameters(shader_);
            inport1_.setTextureParameters(shader_, "texParams_");
            renderQuad();
            glDepthFunc(GL_LESS);

            MatStack.loadIdentity();
            inport1_.getColorTexture()->disable();
        }

        if (inport2_.isReady()) {
            inport2_.bindColorTexture(GL_TEXTURE0);
            inport2_.getColorTexture()->enable();
            MatStack.translate(0.5f, 0.5f, 0.0f);
            MatStack.scale(0.5f, 0.5f, 1.0f);

            glDepthFunc(GL_ALWAYS);
            setGlobalShaderParameters(shader_);
            inport2_.setTextureParameters(shader_, "texParams_");
            renderQuad();
            glDepthFunc(GL_LESS);

            MatStack.loadIdentity();
            inport2_.getColorTexture()->disable();
        }

        if (inport3_.isReady()) {
            inport3_.bindColorTexture(GL_TEXTURE0);
            inport3_.getColorTexture()->enable();
            MatStack.translate(-0.5f, -0.5f, 0.0f);
            MatStack.scale(0.5f, 0.5f, 1.0f);

            glDepthFunc(GL_ALWAYS);
            setGlobalShaderParameters(shader_);
            inport3_.setTextureParameters(shader_, "texParams_");
            renderQuad();
            glDepthFunc(GL_LESS);

            MatStack.loadIdentity();
            inport3_.getColorTexture()->disable();
        }

        if (inport4_.isReady()) {
            inport4_.bindColorTexture(GL_TEXTURE0);
            inport4_.getColorTexture()->enable();
            MatStack.translate(0.5f, -0.5f, 0.0f);
            MatStack.scale(0.5f, 0.5f, 1.0f);

            glDepthFunc(GL_ALWAYS);
            setGlobalShaderParameters(shader_);
            inport4_.setTextureParameters(shader_, "texParams_");
            renderQuad();
            glDepthFunc(GL_LESS);

            MatStack.loadIdentity();
            inport4_.getColorTexture()->disable();
        }

        glActiveTexture(GL_TEXTURE0);

        if(showGrid_.get()) {
            glDepthFunc(GL_ALWAYS);
            glColor4f(gridColor_.get().r, gridColor_.get().g, gridColor_.get().b, gridColor_.get().a);
            setGlobalShaderParameters(shader_);
            shader_->deactivate();
            glBegin(GL_LINES);
            glVertex2f(-1.0f, 0.0f);
            glVertex2f(1.0f, 0.0f);

            glVertex2f(0.0f, 1.0f);
            glVertex2f(0.0f, -1.0f);
            glEnd();
            glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
            shader_->activate();
            glDepthFunc(GL_LESS);
        }

        outport_.deactivateTarget();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadIdentity();
        LGL_ERROR;
    }
    else {
        //maximized:
        switch(maximized_.get()) {
            case 1: if(!inport1_.isReady())
                        return;

                    outport_.activateTarget();
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    inport1_.bindColorTexture(GL_TEXTURE0);
                    inport1_.getColorTexture()->enable();

                    glDepthFunc(GL_ALWAYS);
                    setGlobalShaderParameters(shader_);
                    inport1_.setTextureParameters(shader_, "texParams_");
                    renderQuad();
                    glDepthFunc(GL_LESS);

                    MatStack.loadIdentity();
                    inport1_.getColorTexture()->disable();
                    outport_.deactivateTarget();
                    break;
            case 2: if(!inport2_.isReady())
                        return;

                    outport_.activateTarget();
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    inport2_.bindColorTexture(GL_TEXTURE0);
                    inport2_.getColorTexture()->enable();

                    glDepthFunc(GL_ALWAYS);
                    setGlobalShaderParameters(shader_);
                    inport2_.setTextureParameters(shader_, "texParams_");
                    renderQuad();
                    glDepthFunc(GL_LESS);

                    MatStack.loadIdentity();
                    inport2_.getColorTexture()->disable();
                    outport_.deactivateTarget();
                    break;
            case 3: if(!inport3_.isReady())
                        return;

                    outport_.activateTarget();
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    inport3_.bindColorTexture(GL_TEXTURE0);
                    inport3_.getColorTexture()->enable();

                    glDepthFunc(GL_ALWAYS);
                    setGlobalShaderParameters(shader_);
                    inport3_.setTextureParameters(shader_, "texParams_");
                    renderQuad();
                    glDepthFunc(GL_LESS);

                    MatStack.loadIdentity();
                    inport3_.getColorTexture()->disable();
                    outport_.deactivateTarget();
                    break;
            case 4: if(!inport4_.isReady())
                        return;

                    outport_.activateTarget();
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    inport4_.bindColorTexture(GL_TEXTURE0);
                    inport4_.getColorTexture()->enable();

                    glDepthFunc(GL_ALWAYS);
                    setGlobalShaderParameters(shader_);
                    inport4_.setTextureParameters(shader_, "texParams_");
                    renderQuad();
                    glDepthFunc(GL_LESS);

                    MatStack.loadIdentity();
                    inport4_.getColorTexture()->disable();
                    outport_.deactivateTarget();
                    break;
        }
    }

    shader_->deactivate();
}

void QuadView::initialize() {
    RenderProcessor::initialize();
    updateSizes();
    shader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag", generateHeader(), false);
}

void QuadView::deinitialize() {

    ShdrMgr.dispose(shader_);

    RenderProcessor::deinitialize();
}

void QuadView::portSizeReceiveChanged() {
    updateSizes();
}

void QuadView::updateSizes() {
    if (outport_.getSize() == tgt::ivec2(0))
        return;

    if(maximized_.get() == 0) {
        tgt::ivec2 subsize = outport_.getReceivedSize() / 2;
        inport1_.requestSize(subsize);
        inport2_.requestSize(subsize);
        inport3_.requestSize(subsize);
        inport4_.requestSize(subsize);
    }
    else {
        switch(maximized_.get()) {
            case 1: inport1_.requestSize(outport_.getReceivedSize());
                    break;
            case 2: inport2_.requestSize(outport_.getReceivedSize());
                    break;
            case 3: inport3_.requestSize(outport_.getReceivedSize());
                    break;
            case 4: inport4_.requestSize(outport_.getReceivedSize());
                    break;
            default:;
        }
    }
}

void QuadView::mouseMove(tgt::MouseEvent* e) {
    e->accept();
    int prevCurrenPort = currentPort_;

    if (maximized_.get() == 0) {
        if ((e->action() & tgt::MouseEvent::EXIT) == tgt::MouseEvent::EXIT)
            currentPort_ = -1;

        if ((e->action() & tgt::MouseEvent::PRESSED) == tgt::MouseEvent::PRESSED)
            isDragging_ = true;
        if ((e->action() & tgt::MouseEvent::RELEASED) == tgt::MouseEvent::RELEASED)
            isDragging_ = false;

        if (!isDragging_) {
            if (e->y() < (e->viewport().y/2)) {
                if(e->x() < (e->viewport().x/2)) {
                    currentPort_ = 1;
                }
                else {
                    currentPort_ = 2;
                }
            }
            else {
                if(e->x() < (e->viewport().x/2)) {
                    currentPort_ = 3;
                }
                else {
                    currentPort_ = 4;
                }
            }
        }

        if (currentPort_ != prevCurrenPort) {
            tgt::MouseEvent leaveEvent(1, 1, tgt::MouseEvent::EXIT, e->modifiers(), e->button(), e->viewport() / 2);
            tgt::MouseEvent enterEvent(1, 1, tgt::MouseEvent::ENTER, e->modifiers(), e->button(), e->viewport() / 2);
            leaveEvent.ignore();
            enterEvent.ignore();
            switch(prevCurrenPort) {
                case 1:
                    inport1_.distributeEvent(&leaveEvent);
                    break;
                case 2:
                    inport2_.distributeEvent(&leaveEvent);
                    break;
                case 3:
                    inport3_.distributeEvent(&leaveEvent);
                    break;
                case 4:
                    inport4_.distributeEvent(&leaveEvent);
                    break;
            }
            switch(currentPort_) {
                case 1:
                    inport1_.distributeEvent(&enterEvent);
                    break;
                case 2:
                    inport2_.distributeEvent(&enterEvent);
                    break;
                case 3:
                    inport3_.distributeEvent(&enterEvent);
                    break;
                case 4:
                    inport4_.distributeEvent(&enterEvent);
                    break;
            }
        }
        tgt::MouseEvent moveEvent(e->x() % (e->viewport().x/2), e->y() % (e->viewport().y/2), tgt::MouseEvent::MOTION, e->modifiers(), e->button(), e->viewport() / 2);
        moveEvent.ignore();
        switch(currentPort_) {
            case 1:
                inport1_.distributeEvent(&moveEvent);
                break;
            case 2:
                inport2_.distributeEvent(&moveEvent);
                break;
            case 3:
                inport3_.distributeEvent(&moveEvent);
                break;
            case 4:
                inport4_.distributeEvent(&moveEvent);
                break;
        }
    }
    else {
        switch(maximized_.get()) {
            case 1: inport1_.distributeEvent(e);
                break;
            case 2: inport2_.distributeEvent(e);
                break;
            case 3: inport3_.distributeEvent(e);
                break;
            case 4: inport4_.distributeEvent(e);
                break;
            default:;
        }
    }
}

void QuadView::invalidate(int inv) {
    RenderProcessor::invalidate(inv);
}

void QuadView::onEvent(tgt::Event* e) {

    tgt::TouchEvent* te = dynamic_cast<tgt::TouchEvent*>(e);
    if(te) {
        distributeTouchEvent(te);
        return;
    }

    tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);
    if (me && !(maximized_.get() == 0 && mouseMoveEventProp_.accepts(me)) && !(maximizeEventProp_.accepts(me) && maximizeOnDoubleClick_.get())) {
        distributeMouseEvent(me);
        return;
    }

    // distributes key events according to the mouse (!) position
    tgt::KeyEvent* ke = dynamic_cast<tgt::KeyEvent*>(e);
    if (ke) {
        distributeKeyEvent(ke);
        return;
    }

    RenderProcessor::onEvent(e);
}

// TODO For now, Touch Events are sent to inports only if all touchpoints are located within a single subsection of the viewport
void QuadView::distributeTouchEvent(tgt::TouchEvent* te) {

    if (maximized_.get() == 0) {
        const std::deque<tgt::TouchPoint>& tps = te->touchPoints();
        const tgt::TouchPoint& first = tps.front();
        tgt::vec2 outSize = outport_.getSize();
        int section = first.pos().y < outSize.y / 2 ? (first.pos().x < outSize.x / 2 ? 0 : 1) : (first.pos().x < outSize.x / 2 ? 2 : 3);
        for(std::deque<tgt::TouchPoint>::const_iterator it = tps.begin() + 1; it != tps.end(); it++) {
            const tgt::TouchPoint& tp = *it;
            // TODO different sections -> handle instead of doing nothing
            if(section != (tp.pos().y < outSize.y / 2 ? (tp.pos().x < outSize.x / 2 ? 0 : 1) : (tp.pos().x < outSize.x / 2 ? 2 : 3)))
                return;
        }

        RenderPort* inport = (section == 0 ? &inport1_ : (section == 1 ? &inport2_ : (section == 2 ? &inport3_ : &inport4_)));
        tgt::vec2 offset = section == 0 ? tgt::vec2(0.f) : (section == 1 ? tgt::vec2(-outSize.x / 2.f, 0.f) : (section == 2 ? tgt::vec2(0.f, -outSize.y / 2.f) : -outSize / 2.f));
        std::deque<tgt::TouchPoint> tpsTrafo;
        for(std::deque<tgt::TouchPoint>::const_iterator it = tps.begin(); it != tps.end(); it++) {
            const tgt::TouchPoint& tp = *it;
            tgt::TouchPoint newTP = tp;
            newTP.setPos(tp.pos() + offset);
            tpsTrafo.push_back(newTP);
        }

        tgt::TouchEvent nte = tgt::TouchEvent(tgt::Event::MODIFIER_NONE, te->touchPointStates(), te->deviceType(), tpsTrafo);
        nte.ignore();  // accepted is set to true by default
        inport->distributeEvent(&nte);
        if(nte.isAccepted())
            te->accept();
    }
    else {
        switch(maximized_.get()) {
            case 1: inport1_.distributeEvent(te);
                    break;
            case 2: inport2_.distributeEvent(te);
                    break;
            case 3: inport3_.distributeEvent(te);
                    break;
            case 4: inport4_.distributeEvent(te);
                    break;
            default:;
        }
    }
}

void QuadView::distributeKeyEvent(tgt::KeyEvent* ke) {
    if (maximized_.get() == 0) {
        // not maximized -> find out the right quad
        switch (currentPort_) {
            case 1: inport1_.distributeEvent(ke);
                    break;
            case 2: inport2_.distributeEvent(ke);
                    break;
            case 3: inport3_.distributeEvent(ke);
                    break;
            case 4: inport4_.distributeEvent(ke);
                    break;
            default:;
        }
    }
    else {
        // switch according to maximized inport
        switch(maximized_.get()) {
            case 1: inport1_.distributeEvent(ke);
                    break;
            case 2: inport2_.distributeEvent(ke);
                    break;
            case 3: inport3_.distributeEvent(ke);
                    break;
            case 4: inport4_.distributeEvent(ke);
                    break;
            default:;
        }
    }
}


void QuadView::distributeMouseEvent(tgt::MouseEvent* me) {
    if (maximized_.get() == 0) {
        if (me->y() < (me->viewport().y / 2)) {
            if (me->x() < (me->viewport().x / 2)) {
                tgt::MouseEvent newme(me->x(), me->y(), me->action(), me->modifiers(), me->button(), me->viewport() / 2);
                newme.ignore();  // accepted is set to true by default
                inport1_.distributeEvent(&newme);
                if (newme.isAccepted())
                    me->accept();
            }
            else {
                tgt::MouseEvent newme(me->x() - (me->viewport().x / 2), me->y(), me->action(), me->modifiers(), me->button(), me->viewport() / 2);
                newme.ignore();
                inport2_.distributeEvent(&newme);
                if (newme.isAccepted())
                    me->accept();
            }
        }
        else {
            if (me->x() < (me->viewport().x / 2)) {
                tgt::MouseEvent newme(me->x(), me->y() - (me->viewport().y / 2), me->action(), me->modifiers(), me->button(), me->viewport() / 2);
                newme.ignore();
                inport3_.distributeEvent(&newme);
                if (newme.isAccepted())
                    me->accept();
            }
            else {
                tgt::MouseEvent newme(me->x() - (me->viewport().x / 2), me->y() - (me->viewport().y / 2), me->action(), me->modifiers(), me->button(), me->viewport() / 2);
                newme.ignore();
                inport4_.distributeEvent(&newme);
                if (newme.isAccepted())
                    me->accept();
            }
        }
    }
    else {
        switch(maximized_.get()) {
            case 1: inport1_.distributeEvent(me);
                    break;
            case 2: inport2_.distributeEvent(me);
                    break;
            case 3: inport3_.distributeEvent(me);
                    break;
            case 4: inport4_.distributeEvent(me);
                    break;
            default:;
        }
    }
}

void QuadView::toggleMaximization(tgt::MouseEvent* me) {

    if (maximizeOnDoubleClick_.get()) {
        if (maximized_.get() == 0) {
            if (me->y() < (me->viewport().y / 2)) {
                if (me->x() < (me->viewport().x / 2)) {
                    maximized_.set(1);
                }
                else {
                    maximized_.set(2);
                }
            }
            else {
                if (me->x() < (me->viewport().x / 2)) {
                    maximized_.set(3);
                }
                else {
                    maximized_.set(4);
                }
            }
        }
        else {
            maximized_.set(0);
        }
        updateSizes();


        if (forwardDoubleClickEvent_.get())
            distributeMouseEvent(me);
        else
            me->accept();
    }
}

} // namespace voreen
