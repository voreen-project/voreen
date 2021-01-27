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

#include "tripleview.h"
#include "tgt/immediatemode/immediatemode.h"
namespace voreen {

using tgt::ivec2;
using tgt::vec3;

TripleView::TripleView()
    : RenderProcessor()
    , showGrid_("showGrid", "Show grid", true)
    , gridColor_("gridColor", "Grid color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    , configuration_("configuration", "Configuration")
    , maximized_("maximized", "Maximized sub-view", 0, 0, 4)
    , maximizeOnDoubleClick_("maximizeOnDoubleClick", "Maximize on double click", true)
    , maximizeEventProp_("mouseEvent.maximize", "Maximize Event", this, &TripleView::toggleMaximization,
    tgt::MouseEvent::MOUSE_BUTTON_LEFT, tgt::MouseEvent::DOUBLECLICK, tgt::MouseEvent::MODIFIER_NONE)
    , outport_(Port::OUTPORT, "outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , inport1_(Port::INPORT, "inport1", "Image1 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , inport2_(Port::INPORT, "inport2", "Image2 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , inport3_(Port::INPORT, "inport3", "Image3 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , currentPort_(-1)
    , isDragging_(false)
    , lastWindow_(1)
{
    addProperty(showGrid_);
    addProperty(gridColor_);
    configuration_.addOption("abc", "abc", abc);
    configuration_.addOption("Abc", "Abc", Abc);
    configuration_.addOption("Bac", "Bac", Bac);
    configuration_.addOption("Cab", "Cab", Cab);
    //configuration_.addOption("A", "A", A);
    //configuration_.addOption("B", "B", B);
    //configuration_.addOption("C", "C", C);
    addProperty(configuration_);
    configuration_.onChange(MemberFunctionCallback<TripleView>(this, &TripleView::updateSizes));
    maximized_.onChange(MemberFunctionCallback<TripleView>(this, &TripleView::updateSizes));

    addPort(outport_);
    addPort(inport1_);
    addPort(inport2_);
    addPort(inport3_);


    addProperty(maximized_);
    maximized_.setVisibleFlag(false);
    addProperty(maximizeOnDoubleClick_);
    addEventProperty(maximizeEventProp_);

    outport_.onSizeReceiveChange<TripleView>(this, &TripleView::updateSizes);
}

TripleView::~TripleView() {
}

Processor* TripleView::create() const {
    return new TripleView();
}

bool TripleView::isReady() const {
    if (!outport_.isReady())
        return false;

    if (!inport1_.isReady() && !inport2_.isReady() && !inport3_.isReady())
        return false;
    switch(getWindowConfiguration()) {
        case A: if(!inport1_.isReady())
                    return false;
                break;
        case B: if(!inport2_.isReady())
                    return false;
                break;
        case C: if(!inport3_.isReady())
                    return false;
                break;
        default: // suppress compiler warning
                ;
    }
    return true;
}

void TripleView::renderPortQuad(RenderPort& rp, tgt::vec3 translate, tgt::vec3 scale) {
    rp.bindColorTexture(GL_TEXTURE0);
    rp.getColorTexture()->enable();
    MatStack.translate(translate.x, translate.y, translate.z);
    MatStack.scale(scale.x, scale.y, scale.z);

    glDepthFunc(GL_ALWAYS);
    setGlobalShaderParameters(shader_);
    rp.setTextureParameters(shader_, "texParams_");
    renderQuad();
    glDepthFunc(GL_LESS);

    MatStack.loadIdentity();
    rp.getColorTexture()->disable();
}

void TripleView::renderLargeSmallSmall(RenderPort& large, RenderPort& small1, RenderPort& small2) {
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);

    if (large.isReady())
        renderPortQuad(large, vec3(-0.333333f, 0.0f, 0.0f), vec3(0.666666, 1.0f, 1.0f));

    if (small1.isReady())
        renderPortQuad(small1, vec3(0.666666f, 0.5f, 0.0f), vec3(0.333333f, 0.5f, 1.0f));

    if (small2.isReady())
        renderPortQuad(small2, vec3(0.666666f, -0.5f, 0.0f), vec3(0.333333f, 0.5f, 1.0f));

    glActiveTexture(GL_TEXTURE0);

    if(showGrid_.get()) {
        shader_->deactivate(); //deactivate shader to render lines
        glDepthFunc(GL_ALWAYS);
        IMode.color(gridColor_.get().r, gridColor_.get().g, gridColor_.get().b, gridColor_.get().a);
        IMode.begin(tgt::ImmediateMode::LINES);
        IMode.vertex(0.333333f, -1.0f);
        IMode.vertex(0.333333f, 1.0f);

        IMode.vertex(0.333333f, 0.0f);
        IMode.vertex(1.f, 0.0f);

        IMode.end();
        IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
        glDepthFunc(GL_LESS);
    }

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    LGL_ERROR;
}

void TripleView::process() {
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //shader should be always active
    shader_->activate();
    shader_->setUniform("useTexcoords", true);
    shader_->setUniform("colorTex_", 0);
    shader_->setUniform("depthTex_", 0);
    switch(getWindowConfiguration()) {
        case abc: {
                      MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
                      if (inport1_.isReady())
                          renderPortQuad(inport1_, vec3(-0.66666f, 0.0f, 0.0f), vec3(0.333333f, 1.0f, 1.0f));

                      if (inport2_.isReady())
                          renderPortQuad(inport2_, vec3(0.0f, 0.0f, 0.0f), vec3(0.333333f, 1.0f, 1.0f));

                      if (inport3_.isReady())
                          renderPortQuad(inport3_, vec3(0.666666f, 0.0f, 0.0f), vec3(0.333333f, 1.0f, 1.0f));

                      glActiveTexture(GL_TEXTURE0);

                      if(showGrid_.get()) {
                          shader_->deactivate(); //deactivate shader to render lines
                          glDepthFunc(GL_ALWAYS);
                          IMode.color(gridColor_.get().r, gridColor_.get().g, gridColor_.get().b, gridColor_.get().a);
                          IMode.begin(tgt::ImmediateMode::LINES);
                          IMode.vertex(-0.333333f, -1.0f);
                          IMode.vertex(-0.333333f, 1.0f);

                          IMode.vertex(0.333333f, -1.0f);
                          IMode.vertex(0.333333f, 1.0f);
                          IMode.end();
                          IMode.color(1.0f, 1.0f, 1.0f, 1.0f);
                          glDepthFunc(GL_LESS);
                      }

                      MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
                      MatStack.loadIdentity();
                      LGL_ERROR;
                  }
                  break;
        case Abc: renderLargeSmallSmall(inport1_, inport2_, inport3_);
                  break;
        case Bac: renderLargeSmallSmall(inport2_, inport1_, inport3_);
                  break;
        case Cab: renderLargeSmallSmall(inport3_, inport1_, inport2_);
                  break;
        //maximized:
        case A: if(!inport1_.isReady())
                    return;
                renderPortQuad(inport1_, vec3(0.0f, 0.0f, 0.0f), vec3(1.0f, 1.0f, 1.0f));
                break;
        case B: if(!inport2_.isReady())
                    return;
                renderPortQuad(inport2_, vec3(0.0f, 0.0f, 0.0f), vec3(1.0f, 1.0f, 1.0f));
                break;
        case C: if(!inport3_.isReady())
                    return;
                renderPortQuad(inport3_, vec3(0.0f, 0.0f, 0.0f), vec3(1.0f, 1.0f, 1.0f));
                break;
    }
    //deactivate shader if it is active
    if(shader_->isActivated())
        shader_->deactivate();
    outport_.deactivateTarget();
    LGL_ERROR;
}

void TripleView::initialize() {
    RenderProcessor::initialize();
    updateSizes();
    shader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag", generateHeader(), false);
}

void TripleView::deinitialize() {

    ShdrMgr.dispose(shader_);

    RenderProcessor::deinitialize();
}

void TripleView::updateSizes() {
    if (outport_.getReceivedSize() == tgt::ivec2(0))
        return;
    else {
        tgt::ivec2 outsize = outport_.getReceivedSize();
        switch(getWindowConfiguration()) {
            case abc: inport1_.requestSize(ivec2(static_cast<int>(outsize.x * 0.333333f), outsize.y));
                      inport2_.requestSize(ivec2(static_cast<int>(outsize.x * 0.333333f), outsize.y));
                      inport3_.requestSize(ivec2(static_cast<int>(outsize.x * 0.333333f), outsize.y));
                      break;
            case Abc: inport1_.requestSize(ivec2(static_cast<int>(outsize.x * 0.666666f), outsize.y));
                      inport2_.requestSize(ivec2(static_cast<int>(outsize.x * 0.333333f), outsize.y / 2));
                      inport3_.requestSize(ivec2(static_cast<int>(outsize.x * 0.333333f), outsize.y / 2));
                      break;
            case Bac: inport2_.requestSize(ivec2(static_cast<int>(outsize.x * 0.666666f), outsize.y));
                      inport1_.requestSize(ivec2(static_cast<int>(outsize.x * 0.333333f), outsize.y / 2));
                      inport3_.requestSize(ivec2(static_cast<int>(outsize.x * 0.333333f), outsize.y / 2));
                      break;
            case Cab: inport3_.requestSize(ivec2(static_cast<int>(outsize.x * 0.666666f), outsize.y));
                      inport1_.requestSize(ivec2(static_cast<int>(outsize.x * 0.333333f), outsize.y / 2));
                      inport2_.requestSize(ivec2(static_cast<int>(outsize.x * 0.333333f), outsize.y / 2));
                      break;
            case A: inport1_.requestSize(outport_.getReceivedSize());
                    break;
            case B: inport2_.requestSize(outport_.getReceivedSize());
                    break;
            case C: inport3_.requestSize(outport_.getReceivedSize());
                    break;
            default:;
        }
    }
}


void TripleView::onEvent(tgt::Event* e) {
    tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);

    if (!me || maximizeEventProp_.accepts(me)) {
        RenderProcessor::onEvent(e);
        return;
    }

    tgt::MouseEvent newme(0, 0, tgt::MouseEvent::ACTION_NONE, tgt::Event::Modifier::MODIFIER_NONE); // no empty constructor
    int window = getWindowForEvent(*me, &newme);
    if (window != lastWindow_ && window < 3){
        tgt::MouseEvent leaveEvent(1, 1, tgt::MouseEvent::EXIT, me->modifiers(), me->button(), getWindowViewport(getWindowConfiguration(), lastWindow_));
        tgt::MouseEvent enterEvent(1, 1, tgt::MouseEvent::ENTER, me->modifiers(), me->button(), getWindowViewport(getWindowConfiguration(), window));
        leaveEvent.ignore();
        enterEvent.ignore();
        distributeMouseEvent(lastWindow_, &leaveEvent);
        distributeMouseEvent(window, &enterEvent);
    }
    lastWindow_ = window;
    distributeMouseEvent(window, &newme);

    if (newme.isAccepted())
        me->accept();
}

void TripleView::toggleMaximization(tgt::MouseEvent* me){
    if (!maximizeOnDoubleClick_.get())
        return;
    if (maximized_.get() != 0){
        maximized_.set(0);
        return;
    }
    int window = getWindowForEvent(*me, nullptr);
    maximized_.set(window);
}

TripleView::WindowConfiguration TripleView::getWindowConfiguration() const{
    WindowConfiguration c;
    switch (maximized_.get())
    {
    case 1:
        c = WindowConfiguration::A;
        break;
    case 2:
        c = WindowConfiguration::B;
        break;
    case 3:
        c = WindowConfiguration::C;
        break;
    default:
        c = static_cast<WindowConfiguration>(configuration_.getValue());
    }
    return c;
}

int TripleView::getWindowForEvent(tgt::MouseEvent me, tgt::MouseEvent* translatedMouseEvent){
    switch(getWindowConfiguration()) {
    case abc:
        if (me.x() < (me.viewport().x * 0.33333f)) {
            if (translatedMouseEvent)
                *translatedMouseEvent = tgt::MouseEvent(me.x(), me.y(), me.action(), me.modifiers(), me.button(), ivec2(me.viewport().x / 3, me.viewport().y));
            return 1;
        }
        else if (me.x() < (me.viewport().x * 0.66666f)) {
            if (translatedMouseEvent)
                *translatedMouseEvent = tgt::MouseEvent(me.x() - (me.viewport().x / 3), me.y(), me.action(), me.modifiers(), me.button(), ivec2(me.viewport().x / 3, me.viewport().y));
            return 2;
        }
        else {
            if (translatedMouseEvent)
                *translatedMouseEvent = tgt::MouseEvent(me.x() - (me.viewport().x * 2 / 3), me.y(), me.action(), me.modifiers(), me.button(), ivec2(me.viewport().x / 3, me.viewport().y));
            return 3;
        }
        tgtAssert(false, "Base.TrippleView: could not find window!");
        break;
    case Abc:
        return getWindowForEventLargeSmallSmall(me, 1, 2, 3, translatedMouseEvent);
        break;
    case Bac:
        return getWindowForEventLargeSmallSmall(me, 2, 1, 3, translatedMouseEvent);
        break;
    case Cab:
        return getWindowForEventLargeSmallSmall(me, 3, 1, 2, translatedMouseEvent);
        break;
    case A:
        if (translatedMouseEvent)
            *translatedMouseEvent = me;
        return 1;
    case B:
        if (translatedMouseEvent)
            *translatedMouseEvent = me;
        return 2;
    case C:
        if (translatedMouseEvent)
            *translatedMouseEvent = me;
        return 3;
    default:
        tgtAssert(false, "Base.TrippleView: could not find window!");
        return 0;
    }
}

int TripleView::getWindowForEventLargeSmallSmall(tgt::MouseEvent me, int large, int small1, int small2, tgt::MouseEvent* translatedMouseEvent){
    if (me.x() < (me.viewport().x * 0.666666f)) {
        if (translatedMouseEvent)
            *translatedMouseEvent = tgt::MouseEvent(me.x(), me.y(), me.action(), me.modifiers(), me.button(), ivec2(static_cast<int>(me.viewport().x * 0.666666f), me.viewport().y));
        return large;
    }
    else if (me.y() < (me.viewport().y * 0.5)) {
        if (translatedMouseEvent)
            *translatedMouseEvent = tgt::MouseEvent(me.x() - static_cast<int>(me.viewport().x * 0.666666f), me.y(), me.action(), me.modifiers(), me.button(), ivec2(me.viewport().x / 3, me.viewport().y / 2));
        return small1;
    }
    else {
        if (translatedMouseEvent)
            *translatedMouseEvent = tgt::MouseEvent(me.x() - static_cast<int>(me.viewport().x * 0.666666f), me.y() - static_cast<int>(me.viewport().y * 0.5f), me.action(), me.modifiers(), me.button(), ivec2(me.viewport().x / 3, me.viewport().y / 2));
        return small2;
    }
}

void TripleView::distributeMouseEvent( int Window, tgt::MouseEvent *newme )
{
    RenderPort* distributionPort = 0;
    switch (Window)
    {
    case 1:
        distributionPort = &inport1_;
        break;
    case 2:
        distributionPort = &inport2_;
        break;
    case 3:
        distributionPort = &inport3_;
        break;
    default:
        tgtAssert(false, "TripleView: Should not happen!")
        break;
    }
    newme->ignore();
    distributionPort->distributeEvent(newme);
}

tgt::ivec2 TripleView::getWindowViewportLargeSmallSmall( bool isFirst )
{
    tgt::ivec2 size = outport_.getReceivedSize();
    if (isFirst){
        return tgt::ivec2(size.x*2/3, size.y);
    }else{
        return tgt::ivec2(size.x/3, size.y/2);
    }
}

tgt::ivec2 TripleView::getWindowViewport( WindowConfiguration configuration, int window )
{
    tgt::ivec2 size = outport_.getReceivedSize();
    switch (configuration)
    {
    case voreen::TripleView::abc:
        return tgt::ivec2(size.x/3, size.y);
    case voreen::TripleView::Abc:
        return getWindowViewportLargeSmallSmall(window == 1);
        break;
    case voreen::TripleView::Bac:
        return getWindowViewportLargeSmallSmall(window == 2);
    case voreen::TripleView::Cab:
        return getWindowViewportLargeSmallSmall(window == 3);
        break;
    case voreen::TripleView::A:
        return (window==1)?size:tgt::ivec2::zero;
    case voreen::TripleView::B:
        return (window==2)?size:tgt::ivec2::zero;
    case voreen::TripleView::C:
        return (window==3)?size:tgt::ivec2::zero;
    default:
        return size;
    }
}

} // namespace voreen
