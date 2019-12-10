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

#include "buttonoverlayprocessor.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "tgt/textureunit.h"

#include <sstream>

using tgt::vec3;
using tgt::TextureUnit;

namespace voreen {

ButtonOverlayProcessor::ButtonOverlayProcessor()
    : ImageProcessor("image/compositor")
    , inport_(Port::INPORT, "image.input")
    , outport_(Port::OUTPORT, "image.output", "image.output", true, Processor::INVALID_PROGRAM, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA)
    , overlayPort_(Port::OUTPORT, "image.overlay", "image.overlay", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA)
    , pickingPort_(Port::OUTPORT, "image.picking", "image.picking", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_DEFAULT, GL_RGBA)
    , buttonRadiusX_("buttonRadiusX", "Button Size X", 20, 1, 100)
    , buttonRadiusY_("buttonRadiusY", "Button Size Y", 20, 1, 100)
    , groupBorder_("groupBorder", "Border between Buttons", 5, 0, 30)
    , horzBorder_("horzBorder", "Distance from left/right Border", 30, 0, 1000)
    , vertBorder_("vertBorder", "Distance from top/bottom Border", 100, 0, 1000)
    , fontSize_("fontSize", "Font size", 12, 1, 48)
    , opacity_("opacity", "Opacity", 1.0f)
    , pushButtonMode_("pushButtonMode", "Use pushbutton mode", false)
#ifdef _MSC_VER
#pragma warning(disable:4355)  // passing 'this' is safe here
#endif
    , buttonChecker_("mouseEvent.clickButton", "Clicking Button", this,
                       &ButtonOverlayProcessor::onMouseClick,
                       tgt::MouseEvent::MOUSE_BUTTON_LEFT, tgt::MouseEvent::PRESSED | tgt::MouseEvent::MOTION | tgt::MouseEvent::RELEASED, tgt::Event::MODIFIER_NONE)
    , buttonTexturePressed_(0)
    , buttonTextureReleased_(0)
    , buttonProps_()
    , pushButtonProps_()
    , isVisibleProps_()
    , nameProps_()
    , regenerateOverlay_(false)
    , oldSize_(tgt::ivec2(0))
{
    addEventProperty(buttonChecker_);
    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(&overlayPort_);
    addPrivateRenderPort(&pickingPort_);

    addProperty(horzBorder_);
    addProperty(vertBorder_);
    addProperty(buttonRadiusX_);
    addProperty(buttonRadiusY_);
    addProperty(groupBorder_);
    addProperty(opacity_);
    addProperty(fontSize_);
    addProperty(pushButtonMode_);

    buttonRadiusX_.onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));
    buttonRadiusY_.onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));
    opacity_.onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));
    groupBorder_.onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));
    horzBorder_.onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));
    vertBorder_.onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));
    fontSize_.onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));

    for(size_t i = 0; i < 5; i++) {
        std::ostringstream os1, os2, os3, os4, os5, os6, os7, os8, os9;
        os1 << "showProp" << i;
        os2 << "Show Button No. " << i + 1;
        isVisibleProps_.push_back(new BoolProperty(os1.str(), os2.str(), false));
        isVisibleProps_.back()->onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));
        addProperty(isVisibleProps_.back());

        os3 << "buttonProp" << i;
        os4 << "Button Property " << i + 1;
        buttonProps_.push_back(new BoolProperty(os3.str(), os4.str(), false));
        buttonProps_.back()->onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));
        addProperty(buttonProps_.back());

        os5 << "nameProp" << i;
        os6 << "Button Name Property " << i + 1;
        os7 << "Button " << i + 1;
        nameProps_.push_back(new StringProperty(os5.str(), os6.str(), os7.str()));
        nameProps_.back()->onChange(MemberFunctionCallback<ButtonOverlayProcessor>(this, &ButtonOverlayProcessor::onOverlayChange));
        addProperty(nameProps_.back());

        os8 << "pushButtonProp" << i;
        os9 << "Pushbutton Property " << i + 1;
        pushButtonProps_.push_back(new ButtonProperty(os8.str(), os9.str()));
        addProperty(pushButtonProps_.back());
    }

    font_ = new tgt::Font(VoreenApplication::app()->getFontPath("Vera.ttf"), 12);
    //font_ = new tgt::Font(VoreenApplication::app()->getFontPath("VeraMono.ttf"), 12, tgt::BitmapFont);
}

void ButtonOverlayProcessor::deinitialize() {
    TexMgr.dispose(buttonTexturePressed_);
    TexMgr.dispose(buttonTextureReleased_);
    ImageProcessor::deinitialize();
}

ButtonOverlayProcessor::~ButtonOverlayProcessor() {
    for(size_t i = 0; i < isVisibleProps_.size(); i++) {
        delete isVisibleProps_.at(i);
        delete buttonProps_.at(i);
        delete nameProps_.at(i);
        delete pushButtonProps_.at(i);
    }
    delete font_;
}

std::string ButtonOverlayProcessor::getProcessorInfo() const {
    return "";
}

void ButtonOverlayProcessor::initialize() {
    //TexMgr.addPath(VoreenApplication::app()->getModulePath("viscontest/icons"));
    buttonTexturePressed_ = TexMgr.load("pressedButton.png");
    buttonTextureReleased_ = TexMgr.load("releasedButton.png");

    ImageProcessor::initialize();
}

bool ButtonOverlayProcessor::isReady() const {
    //check if at least one outport is connected:
    if (!outport_.isReady() && !inport_.isReady())
        return false;

    return true;
}

void ButtonOverlayProcessor::process() {
    if (!outport_.isReady())
        return;

    if (!inport_.isReady())
        return;

    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();

    if(outport_.getSize() != oldSize_) {
        regenerateOverlay_ = true;
        oldSize_ = outport_.getSize();
    }

    if(regenerateOverlay_)
        renderOverlayImage();

    outport_.activateTarget();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    TextureUnit colorUnit0, depthUnit0, colorUnit1, depthUnit1;
    overlayPort_.bindTextures(colorUnit0.getEnum(), depthUnit0.getEnum());
    inport_.bindTextures(colorUnit1.getEnum(), depthUnit1.getEnum());

    // initialize shader
    program_->activate();
    setGlobalShaderParameters(program_);
    program_->setUniform("colorTex0_", colorUnit0.getUnitNumber());
    program_->setUniform("depthTex0_", depthUnit0.getUnitNumber());
    program_->setUniform("colorTex1_", colorUnit1.getUnitNumber());
    program_->setUniform("depthTex1_", depthUnit1.getUnitNumber());
    overlayPort_.setTextureParameters(program_, "textureParameters0_");
    inport_.setTextureParameters(program_, "textureParameters1_");

    glDepthFunc(GL_ALWAYS);
    renderQuad();
    glDepthFunc(GL_LESS);
    outport_.deactivateTarget();

    program_->deactivate();
    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

void ButtonOverlayProcessor::renderOverlayImage() {

    regenerateOverlay_ = false;

    tgt::Texture* pickingTex = pickingPort_.getColorTexture();
    if(pickingTex->getCpuTextureData())
        pickingTex->setCpuTextureData(0,false);

    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);

    overlayPort_.activateTarget();
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::TEXTURE);
    MatStack.loadIdentity();

    glClearDepth(1.0);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    pickingPort_.activateTarget();
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    glClearDepth(1.0);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    TextureUnit overlayUnit;
    tgt::ivec2 groupCoords = inport_.getSize() - tgt::ivec2(horzBorder_.get(), vertBorder_.get());

    int count = 0;

    for(size_t i = 0; i < buttonProps_.size(); i++) {

        if(!isVisibleProps_.at(i)->get())
            continue;

        overlayPort_.activateTarget();
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_COLOR, GL_ONE_MINUS_SRC_COLOR);

        // add support for transparent buttons in combination with user controlled transparency, otherwise
        // the texture transparency always has priority
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

        overlayUnit.activate();
        if(buttonProps_.at(i)->get()) {
            //buttonTexturesPressed_.at(i)->bind();
            //buttonTexturesPressed_.at(i)->enable();
            buttonTexturePressed_->bind();
            buttonTexturePressed_->enable();
        } else {
            //buttonTexturesNonPressed_.at(i)->bind();
            //buttonTexturesNonPressed_.at(i)->enable();
            buttonTextureReleased_->bind();
            buttonTextureReleased_->enable();
        }
        LGL_ERROR;

        tgt::ivec2 radius = tgt::ivec2(buttonRadiusX_.get(), buttonRadiusY_.get());
        tgt::ivec2 centerCoords = groupCoords - tgt::ivec2(0, count*(groupBorder_.get() + 2*radius.y)) - radius;

        //radius.x = (float(buttonTexturePressed_->getWidth()) / float(buttonTexturePressed_->getHeight())) * radius.y;

        glColor4f(1.0f, 1.0f, 1.0f, opacity_.get());
        glBegin(GL_QUADS);

        tgt::vec2 ll = tgt::vec2(float(centerCoords.x - radius.x) / inport_.getSize().x, float(centerCoords.y - radius.y) / inport_.getSize().y);
        tgt::vec2 ur = tgt::vec2(float(centerCoords.x + radius.x) / inport_.getSize().x, float(centerCoords.y + radius.y) / inport_.getSize().y);
        ll = 2.0f * ll - 1.0f;
        ur = 2.0f * ur - 1.0f;

        glTexCoord2f(0.0f, 0.0f);
        glVertex2f(ll.x, ll.y);
        glTexCoord2f(1.0f, 0.0f);
        glVertex2f(ur.x, ll.y);
        glTexCoord2f(1.0f, 1.0f);
        glVertex2f(ur.x, ur.y);
        glTexCoord2f(0.0f, 1.0f);
        glVertex2f(ll.x, ur.y);

        glEnd();

        if(buttonProps_.at(i)->get())
            buttonTexturePressed_->disable();
        else
            buttonTextureReleased_->disable();

        renderFont(static_cast<int>(i), centerCoords);

        overlayPort_.deactivateTarget();
        glBlendFunc(GL_ONE, GL_ZERO);
        glDisable(GL_BLEND);

        pickingPort_.activateTarget();
        float col = float(i+1) / 255.f;
        glColor4f(col, col, col, 1.f);

        glRectf(ll.x, ll.y, ur.x, ur.y);
        pickingPort_.deactivateTarget();

        count++;
    }

    // download the picking texture to a buffer so we can perform picking on the cpu
    pickingPort_.getColorTexture()->downloadTexture();
    TextureUnit::setZeroUnit();
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
}

void ButtonOverlayProcessor::renderFont(int i, tgt::ivec2& center) {

    MatStack.pushMatrix();
    MatStack.loadIdentity();
    MatStack.translate(-1,-1,0);
    float scaleFactorX = 2.0f / (float)inport_.getSize().x;
    float scaleFactorY = 2.0f / (float)inport_.getSize().y;
    MatStack.scale(scaleFactorX, scaleFactorY, 1);

    glColor4f(1.0f, 1.0f, 1.0f, opacity_.get());

    font_->setFontSize(fontSize_.get());
    font_->setLineWidth(2.f * float(buttonRadiusX_.get()));

    font_->setTextAlignment(tgt::Font::Centered);
    font_->render(tgt::vec3((float)center.x, (float)center.y, 0),
                  nameProps_.at(i)->get(), outport_.getSize());
    MatStack.popMatrix();
}


std::string ButtonOverlayProcessor::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = ImageProcessor::generateHeader(version);
    header += "#define MODE_ALPHA_BLENDING\n";
    return header;
}


void ButtonOverlayProcessor::onMouseClick(tgt::MouseEvent* e) {

    if(e->getEventType() == tgt::MouseEvent::MOUSEPRESSEVENT) {
        tgt::col4 color = tgt::col4(0, 0, 0, 0);

        tgt::Texture* pickingTex = pickingPort_.getColorTexture();
        if (!pickingTex->getCpuTextureData()) {
            LWARNING("Picking texture has no cpu data");
        }
        else {
            color = pickingTex->texel<tgt::col4>(e->coord().x, pickingTex->getHeight() - e->coord().y);
        }

        if (color.x > 0) {
            e->accept();
            bool invert = !(buttonProps_.at(color.x - 1)->get());
            buttonProps_.at(color.x - 1)->set(invert);
            regenerateOverlay_ = true;
        }
    }

    if(e->getEventType() == tgt::MouseEvent::MOUSERELEASEEVENT && pushButtonMode_.get()) {
        tgt::col4 color = tgt::col4(0, 0, 0, 0);

        tgt::Texture* pickingTex = pickingPort_.getColorTexture();
        if (!pickingTex->getCpuTextureData()) {
            LWARNING("Picking texture has no cpu data");
        }
        else {
            color = pickingTex->texel<tgt::col4>(e->coord().x, pickingTex->getHeight() - e->coord().y);
        }

        if (color.x > 0) {
            e->accept();
            buttonProps_.at(color.x - 1)->set(false);
            regenerateOverlay_ = true;
            pushButtonProps_.at(color.x - 1)->clicked();
        }
    }
}

void ButtonOverlayProcessor::onOverlayChange() {
    regenerateOverlay_ = true;
}

} // namespace
