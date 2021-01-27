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

#include "tabbedview.h"

#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

using tgt::ivec2;

TabbedView::TabbedView()
    : RenderProcessor()
    //ports
    , inport1_(Port::INPORT, "inport1", "Image1 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , inport2_(Port::INPORT, "inport2", "Image2 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , inport3_(Port::INPORT, "inport3", "Image3 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , inport4_(Port::INPORT, "inport4", "Image4 Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , outport_(Port::OUTPORT, "outport", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    //properties
    , currentView_("currentView", "Current View")
        //tab-bar
    , hideTabbar_("hideTabbar_", "Hide Tab-bar", false)
    , renderAtBottom_("renderAtBottom", "Render Tab-bar at Bottom", false)
    , borderColor_("borderColor", "Border Color", tgt::vec4(0.0f, 0.0f, 0.0f, 1.0f))
    , buttonColor_("buttonColor", "Button Color", tgt::vec4(1.0f, 1.0f, 1.0f, 1.0f))
    , buttonHoverColor_("buttonHoverColor", "Button Hover Color", tgt::vec4(0.7f, 0.7f, 0.7f, 1.0f))
    , buttonActiveColor_("buttonActiveColor", "Button Active Color", tgt::vec4(0.3f, 0.3f, 0.3f, 1.0f))
    , textColor_("textColor", "Text Color", tgt::vec4(0.0f, 0.0f, 0.0f, 1.0f))
    , textHoverColor_("textHoverColor", "Text Hover Color", tgt::vec4(0.0f, 0.0f, 0.0f, 1.0f))
    , textActiveColor_("textActiveColor", "Text Active Color", tgt::vec4(0.9f, 0.9f, 0.9f, 1.0f))
    , fontProp_("voreen.fontprop", "Font:")
        //tab configuration
    , tabText1_("tabText1", "Tab 1 Label:", "Tab 1")
    , tabText2_("tabText2", "Tab 2 Label:", "Tab 2")
    , enableTab3_("enableTab3", "Enable Tab 3", true)
    , tabText3_("tabText3", "Tab 3 Label:", "Tab 3")
    , enableTab4_("enableTab4", "Enable Tab 4", true)
    , tabText4_("tabText4", "Tab 4 Label:", "Tab 4")
        //event helper
    , insideViewPort_(false)
    , mouseOverButton_(-1)
    , isDragging_(false)
{
    //ports
    addPort(outport_);
    addPort(inport1_);
    addPort(inport2_);
    addPort(inport3_);
    addPort(inport4_);
    outport_.onSizeReceiveChange<TabbedView>(this, &TabbedView::updateSizes);

    //properties
    addProperty(currentView_);
    currentView_.addOption("first","First Inport",1);
    currentView_.addOption("second","Second Inport",2);
    currentView_.addOption("third","Third Inport",3);
    currentView_.addOption("forth","Forth Inport",4);
        //tab-bar settings
    addProperty(hideTabbar_);
    hideTabbar_.onChange(MemberFunctionCallback<TabbedView>(this, &TabbedView::updateSizes));
    hideTabbar_.setGroupID("bar");
    addProperty(renderAtBottom_);
    renderAtBottom_.setGroupID("bar");
    addProperty(borderColor_);
    borderColor_.setGroupID("bar");
    addProperty(buttonColor_);
    buttonColor_.setGroupID("bar");
    addProperty(buttonHoverColor_);
    buttonHoverColor_.setGroupID("bar");
    addProperty(buttonActiveColor_);
    buttonActiveColor_.setGroupID("bar");
    addProperty(textColor_);
    textColor_.setGroupID("bar");
    addProperty(textHoverColor_);
    textHoverColor_.setGroupID("bar");
    addProperty(textActiveColor_);
    textActiveColor_.setGroupID("bar");
    addProperty(fontProp_);
    fontProp_.onChange(MemberFunctionCallback<TabbedView>(this, &TabbedView::updateSizes));
    fontProp_.setGroupID("bar");
    setPropertyGroupGuiName("bar","Tab-Bar Settings");
        //tab configuration
    addProperty(tabText1_);
    tabText1_.setGroupID("tab");
    addProperty(tabText2_);
    tabText2_.setGroupID("tab");
    addProperty(enableTab3_);
    enableTab3_.setGroupID("tab");
    addProperty(tabText3_);
    enableTab3_.onChange(MemberFunctionCallback<TabbedView>(this, &TabbedView::updateNumOptions));
    tabText3_.setGroupID("tab");
    addProperty(enableTab4_);
    enableTab4_.onChange(MemberFunctionCallback<TabbedView>(this, &TabbedView::updateNumOptions));
    enableTab4_.setGroupID("tab");
    addProperty(tabText4_);
    tabText4_.setGroupID("tab");
    setPropertyGroupGuiName("tab","Tab Configuration");
}

TabbedView::~TabbedView() {
}

Processor* TabbedView::create() const {
    return new TabbedView();
}

bool TabbedView::isReady() const {
    if (!outport_.isReady()) {
        setNotReadyErrorMessage("Outport not connected!");
        return false;
    }

    return true;
}

void TabbedView::process() {
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    float scale = 1.0f;
    float offset = 0.0f;
    if(!hideTabbar_.get()) {
        scale = 1.0f - ((float) getBorderWidth() / (float)outport_.getSize().y);
        offset = (float) getBorderWidth() / (float) outport_.getSize().y;
    }

    if(renderAtBottom_.get())
        MatStack.translate(0.0f, +offset, 0.0f);
    else
        MatStack.translate(0.0f, -offset, 0.0f);

    MatStack.scale(1.0f, scale, 1.0f);

    RenderPort* currentPort = 0;
    switch(currentView_.getValue()) {
        case 1: currentPort = &inport1_;
                break;
        case 2: currentPort = &inport2_;
                break;
        case 3: currentPort = &inport3_;
                break;
        case 4: currentPort = &inport4_;
                break;
    }

    if(currentPort && currentPort->isReady()) {
        currentPort->bindColorTexture(GL_TEXTURE0);
        currentPort->getColorTexture()->enable();

        glDepthFunc(GL_ALWAYS);
        renderQuad();
        glDepthFunc(GL_LESS);

        currentPort->getColorTexture()->disable();
    }

    if(!hideTabbar_.get()) {
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadIdentity();
        // render buttons with fonts in screenspace
        MatStack.translate(-1.0f, -1.0f, 0.0f);
        float scaleFactorX = 2.0f / (float)outport_.getSize().x;
        float scaleFactorY = 2.0f / (float)outport_.getSize().y;
        MatStack.scale(scaleFactorX, scaleFactorY, 1.0f);

        //float yStart = 0.0f;
        float yStart = (float)outport_.getSize().y;
        float yEnd = yStart - getBorderWidth();
        if(renderAtBottom_.get()) {
            yStart = (float)getBorderWidth();
            yEnd = 0;
        }
        glLineWidth(2.0f);
        int numOptions = currentView_.getOptions().size();
        int buttonWidth = tgt::iround(outport_.getSize().x / (float)numOptions);

        fontProp_.get()->setLineWidth(buttonWidth - 10.f); //TODO: use linewidth
        fontProp_.get()->setTextAlignment(tgt::Font::MiddleLeft);

        for(int i=0; i<numOptions; i++) {
            if (mouseOverButton_ == i)
                IMode.color(buttonHoverColor_.get().r, buttonHoverColor_.get().g, buttonHoverColor_.get().b, buttonHoverColor_.get().a);
            else if (currentView_.getValue() - 1 == i)
                IMode.color(buttonActiveColor_.get().r, buttonActiveColor_.get().g, buttonActiveColor_.get().b, buttonActiveColor_.get().a);
            else
                IMode.color(buttonColor_.get().r, buttonColor_.get().g, buttonColor_.get().b, buttonColor_.get().a);

            float xStart = 0.f + (i*buttonWidth);
            float xEnd = xStart + buttonWidth;;

            IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(xStart, yStart, 0.1f);
            IMode.vertex(xEnd, yStart, 0.1f);

            IMode.vertex(xEnd, yEnd, 0.1f);
            IMode.vertex(xStart, yEnd, 0.1f);
            IMode.end();

            glDepthFunc(GL_ALWAYS);
            IMode.color(borderColor_.get().r, borderColor_.get().g, borderColor_.get().b, borderColor_.get().a);

            IMode.begin(tgt::ImmediateMode::LINE_LOOP);
            IMode.vertex(xStart, yStart);
            IMode.vertex(xEnd, yStart);

            IMode.vertex(xEnd, yEnd);
            IMode.vertex(xStart, yEnd);
            IMode.end();
            glDepthFunc(GL_LESS);

            std::string label;
            switch(i) {
                case 0: label = tabText1_.get();
                        break;
                case 1: label = tabText2_.get();
                        break;
                case 2: if(enableTab3_.get())
                            label = tabText3_.get();
                        else
                            label = tabText4_.get();
                        break;
                case 3: label = tabText4_.get();
                        break;
            }

            if (mouseOverButton_ == i)
                fontProp_.get()->setFontColor(textHoverColor_.get());
            else if (currentView_.getValue() - 1 == i)
                fontProp_.get()->setFontColor(textActiveColor_.get());
            else
                fontProp_.get()->setFontColor(textColor_.get());
            

            float yMid = (yStart + yEnd) / 2.0f;
            fontProp_.get()->setTextAlignment(tgt::Font::MiddleLeft);
            fontProp_.get()->render(tgt::vec3(xStart + 5, yMid, 0), label, outport_.getSize());
        }

        glLineWidth(1.0f);
    }
    MatStack.loadIdentity();

    outport_.deactivateTarget();

    IMode.color(1.f, 1.f, 1.f, 1.f);

    LGL_ERROR;
}

void TabbedView::initialize() {
    RenderProcessor::initialize();
    updateSizes();
}

int TabbedView::getBorderWidth() {
    if(hideTabbar_.get())
        return 0;

    if (fontProp_.get()) {
        return tgt::iround(fontProp_.get()->getLineHeight() + 5);
    }
    else
        return 20;
}

ivec2 TabbedView::getInternalSize() {
    if (outport_.getReceivedSize() == tgt::ivec2(0))
        return ivec2(0, 0);
    else {
        if(!hideTabbar_.get())
            return outport_.getReceivedSize() - ivec2(0, getBorderWidth());
        else
            return outport_.getReceivedSize();
    }
}



void TabbedView::handleMouseEvent(tgt::MouseEvent* e) {
    e->accept();
    bool wasInsideViewPort = insideViewPort_;

    if ((e->action() & tgt::MouseEvent::EXIT) == tgt::MouseEvent::EXIT) {
        insideViewPort_ = false;

        if(mouseOverButton_ != -1)
            invalidate();
        mouseOverButton_ = -1;
    }

    if ((e->action() & tgt::MouseEvent::PRESSED) == tgt::MouseEvent::PRESSED)
        isDragging_ = true;
    if ((e->action() & tgt::MouseEvent::RELEASED) == tgt::MouseEvent::RELEASED)
        isDragging_ = false;

    if (!isDragging_) {
        bool insideY;
        if(!hideTabbar_.get()) {
            if(renderAtBottom_.get()) {
                insideY = e->y() < (outport_.getSize().y - getBorderWidth());
                insideY &= e->y() > 0;
            }
            else {
                insideY = e->y() < outport_.getSize().y;
                insideY &= e->y() > getBorderWidth();
            }

            if(insideY && (e->x() >= 0) && (e->x() < outport_.getSize().x)) {
                insideViewPort_ = true;
            }
            else {
                insideViewPort_ = false;
            }
        }
        else
            insideViewPort_ = true;
    }

    if (wasInsideViewPort != insideViewPort_) {
        if(wasInsideViewPort) {
            tgt::MouseEvent leaveEvent(1, 1, tgt::MouseEvent::EXIT, e->modifiers(), e->button(), e->viewport() - ivec2(0, getBorderWidth()));
            leaveEvent.ignore();
            switch(currentView_.getValue()) {
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
        }
        else {
            tgt::MouseEvent enterEvent(1, 1, tgt::MouseEvent::ENTER, e->modifiers(), e->button(), e->viewport() - ivec2(0, getBorderWidth()));
            enterEvent.ignore();
            switch(currentView_.getValue()) {
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
    }

    if(insideViewPort_) {
        //hidden tabbar is implicitely handled through getBorderWidth
        int newY =  e->y() - getBorderWidth();
        if(renderAtBottom_.get())
            newY = e->y();

        tgt::MouseEvent moveEvent(e->x(), newY, e->action(), e->modifiers(), e->button(), e->viewport() - ivec2(0, getBorderWidth()));
        moveEvent.ignore();

        switch(currentView_.getValue()) {
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

        if(mouseOverButton_ != -1)
            invalidate();
        mouseOverButton_ = -1;
    }
    else {
        if(hideTabbar_.get()) {
            mouseOverButton_ = -1;
        }
        else {
            //in tab-bar:
            if((e->y() > 0) && (e->y() < outport_.getSize().y) && (e->x() > 0) && (e->x() < outport_.getSize().x)) {
                int div = outport_.getSize().x / currentView_.getOptions().size();
                int n = e->x() / div;

                if(e->action() == tgt::MouseEvent::PRESSED) {
                    if(e->button() == tgt::MouseEvent::MOUSE_BUTTON_LEFT) {
                        currentView_.selectByIndex(n);
                    }
                }

                if(mouseOverButton_ != n)
                    invalidate();
                mouseOverButton_ = n;
            }
            else {
                if(mouseOverButton_ != -1)
                    invalidate();
                mouseOverButton_ = -1;
            }
        }
    }
}

void TabbedView::invalidate(int inv) {
    RenderProcessor::invalidate(inv);
}

void TabbedView::onEvent(tgt::Event* e) {
    tgt::MouseEvent* me = dynamic_cast<tgt::MouseEvent*>(e);

    if(me) {
        handleMouseEvent(me);
    }
    else {
        switch(currentView_.getValue()) {
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

//-----------------------------------------------------------------------------------------
//  Callbacks
//-----------------------------------------------------------------------------------------
void TabbedView::updateSizes() {
    if (outport_.getReceivedSize() == tgt::ivec2(0))
        return;

    tgt::ivec2 subsize = getInternalSize();
    inport1_.requestSize(subsize);
    inport2_.requestSize(subsize);
    inport3_.requestSize(subsize);
    inport4_.requestSize(subsize);
}

void TabbedView::updateNumOptions() {
    if(enableTab3_.get()) {
        if(!currentView_.hasKey("third"))
            currentView_.addOption("third","Third Inport",3,tgt::vec4::zero,2);
    } else {
        //test, if "third" exists is done in removeOption
        currentView_.removeOption("third");
    }

    if(enableTab4_.get()) {
        if(!currentView_.hasKey("forth"))
            currentView_.addOption("forth","Forth Inport",4,tgt::vec4::zero,3);
    } else {
        //test, if "forth" exists is done in removeOption
        currentView_.removeOption("forth");
    }

}

} // namespace voreen
