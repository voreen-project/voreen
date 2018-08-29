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

#include "textoverlay.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/textureunit.h"

#include <sstream>

using tgt::TextureUnit;

namespace {
#ifdef __APPLE__
    const int fontSize = 11;
#else
    const int fontSize = 12;
#endif
}

namespace voreen {

const std::string TextOverlay::loggerCat_("voreen.TextOverlay");

TextOverlay::TextOverlay()
    : ImageProcessor("textoverlay")
    , inport_(Port::INPORT, "image.input", "Image Input")
    , outport_(Port::OUTPORT, "image.output", "Image Output")
    , privatePort_(Port::OUTPORT, "image.tmp")
    , text0_(Port::INPORT, "text.text0", "Text0 Input", true)
    , text1_(Port::INPORT, "text.text1", "Text1 Input", true)
    , text2_(Port::INPORT, "text.text2", "Text2 Input", true)
    , text3_(Port::INPORT, "text.text3", "Text3 Input", true)
    , enabled_("enabled", "Enabled", true)
    , layout0_("layout0", "Text position 1:")
    , layout1_("layout1", "Text position 2:")
    , layout2_("layout2", "Text position 3:")
    , layout3_("layout3", "Text position 4:")
#ifdef _MSC_VER
#pragma warning(disable:4355)  // passing 'this' is safe here
#endif
    , mouseMoveEventProp_("mouseEvent.move", "Move Event", this, &TextOverlay::mouseMove,
        tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::MOTION, tgt::MouseEvent::MODIFIER_NONE)
#ifdef _MSC_VER
#pragma warning(disable:4355)  // passing 'this' is safe here
#endif
    , mouseEnterExitEventProp_("mouseEvent.EnterExit", "Enter/Exit Event", this, &TextOverlay::mouseEnterExit,
        tgt::MouseEvent::MOUSE_BUTTON_NONE, tgt::MouseEvent::ENTER_EXIT, tgt::MouseEvent::MODIFIER_NONE)
    , viewportSize_(0,0)
    , mousePos_(0,0)
    , fontProp_("voreen.fontprop", "Font:")
    , blendMode_("blendMode","Color Mode:")
    , fontColor_("fontColor","Color of Font:",tgt::vec4(1.0f))
    , renderFollowMouseText_(false)
{
    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(&privatePort_);
    addPort(text0_);
    addPort(text1_);
    addPort(text2_);
    addPort(text3_);

    std::deque<Option<std::string> > options;
    options.push_back(Option<std::string>("FOLLOW", "Follow mouse", ""));
    options.push_back(Option<std::string>("N", "North", ""));
    options.push_back(Option<std::string>("NE", "North-East", ""));
    options.push_back(Option<std::string>("E", "East", ""));
    options.push_back(Option<std::string>("SE", "South-East", ""));
    options.push_back(Option<std::string>("S", "South", ""));
    options.push_back(Option<std::string>("SW", "South-West", ""));
    options.push_back(Option<std::string>("W", "West", ""));
    options.push_back(Option<std::string>("NW", "North-West", ""));
    options.push_back(Option<std::string>("CENTER", "Center", ""));

    layout0_.setOptions(options);
    layout0_.select("NW");
    layout1_.setOptions(options);
    layout1_.select("NW");
    layout2_.setOptions(options);
    layout2_.select("NW");
    layout3_.setOptions(options);
    layout3_.select("NW");

    blendMode_.addOption("auto","Inverse Color of Background",0);
    blendMode_.addOption("fix","Chosen Color",1);
    blendMode_.select("auto");
    blendMode_.onChange(MemberFunctionCallback<TextOverlay>(this, &TextOverlay::blendModeOnChange));

    addProperty(enabled_);
    addProperty(layout0_);
    addProperty(layout1_);
    addProperty(layout2_);
    addProperty(layout3_);

    addEventProperty(mouseMoveEventProp_);
    addEventProperty(mouseEnterExitEventProp_);

    addProperty(fontProp_);
    addProperty(blendMode_);
    addProperty(fontColor_);

    blendModeOnChange();
}

TextOverlay::~TextOverlay() {}

Processor* TextOverlay::create() const {
    return new TextOverlay();
}

void TextOverlay::initialize() {
    ImageProcessor::initialize();
}

void TextOverlay::deinitialize() {
    ImageProcessor::deinitialize();
}

bool TextOverlay::isReady() const {
    if (!isInitialized())
        return false;

    if (!inport_.isReady())
        return false;

    if (!outport_.isReady())
        return false;

    return true;
}

void TextOverlay::blendModeOnChange(){
    switch(blendMode_.getValue()){
    case 1:
        fontColor_.setVisibleFlag(true);
        break;
    case 0:
        fontColor_.setVisibleFlag(false);
        break;
    default:
        break;
    }
}

void TextOverlay::process() {
    if (!inport_.isReady())
        outport_.activateTarget();
    else
        privatePort_.activateTarget();

    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    glClearDepth(1.0);

    glDisable(GL_DEPTH_TEST);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    if(enabled_.get())
        renderOverlay();

    if (!inport_.isReady())
        outport_.deactivateTarget();
    else
        privatePort_.deactivateTarget();

    glEnable(GL_DEPTH_TEST);

    // skip if there's nothing to render above
    if (inport_.isReady()) {
        outport_.activateTarget();
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        TextureUnit colorUnit0, colorUnit1, depthUnit0, depthUnit1;
        privatePort_.bindTextures(colorUnit0.getEnum(), depthUnit0.getEnum());
        inport_.bindTextures(colorUnit1.getEnum(), depthUnit1.getEnum(), GL_NEAREST);

        // initialize shader
        program_->activate();
        setGlobalShaderParameters(program_);
        program_->setUniform("colorTex0_", colorUnit0.getUnitNumber());
        program_->setUniform("colorTex1_", colorUnit1.getUnitNumber());
        program_->setUniform("depthTex1_", depthUnit1.getUnitNumber());
        privatePort_.setTextureParameters(program_, "textureParameters0_");
        inport_.setTextureParameters(program_, "textureParameters1_");
        program_->setUniform("option_", blendMode_.getValue());
        program_->setUniform("threshold_", 0.2f); // default set within shader file

        ImageProcessor::renderQuad();

        outport_.deactivateTarget();

        TextureUnit::setZeroUnit();
        program_->deactivate();
    }
    LGL_ERROR;
}

void TextOverlay::renderOverlay() {
    switch(blendMode_.getValue()){
    case 1:
        glClear(GL_COLOR_BUFFER_BIT);
        break;
    case 0:
    default:
        glClear(GL_COLOR_BUFFER_BIT);
        break;
    }
    // render font(s)

    tgt::ivec2 screensize = inport_.getSize();

    fontProp_.get()->setFontColor(fontColor_.get());
    fontProp_.get()->setLineWidth((float)inport_.getSize().x);


    // left
    fontProp_.get()->setTextAlignment(tgt::Font::BottomLeft);
    fontProp_.get()->render(tgt::vec3(0, (float)inport_.getSize().y, 0), collectText("NW"), screensize);

    fontProp_.get()->setTextAlignment(tgt::Font::MiddleLeft);
    fontProp_.get()->render(tgt::vec3(0, (float)inport_.getSize().y * 0.5f, 0), collectText("W"), screensize);

    fontProp_.get()->setTextAlignment(tgt::Font::TopLeft);
    fontProp_.get()->render(tgt::vec3(0, 0, 0), collectText("SW"), screensize);

    // center
    fontProp_.get()->setTextAlignment(tgt::Font::BottomCentered);
    fontProp_.get()->render(tgt::vec3(0, (float)inport_.getSize().y , 0), collectText("N"), screensize);

    fontProp_.get()->setTextAlignment(tgt::Font::Centered);
    fontProp_.get()->render(tgt::vec3(0, (float)inport_.getSize().y * 0.5f, 0), collectText("CENTER"), screensize);

    fontProp_.get()->setTextAlignment(tgt::Font::TopCentered);
    fontProp_.get()->render(tgt::vec3(0, 0, 0), collectText("S"), screensize);

    // right
    fontProp_.get()->setTextAlignment(tgt::Font::BottomRight);
    fontProp_.get()->render(tgt::vec3(screensize.x, (float)inport_.getSize().y, 0), collectText("NE"), screensize);

    fontProp_.get()->setTextAlignment(tgt::Font::MiddleRight);
    fontProp_.get()->render(tgt::vec3(screensize.x, (float)inport_.getSize().y * 0.5f, 0), collectText("E"), screensize);

    fontProp_.get()->setTextAlignment(tgt::Font::TopRight);
    fontProp_.get()->render(tgt::vec3(screensize.x, 0, 0), collectText("SE"), screensize);


    // Follow mouse
    if(renderFollowMouseText_) {
        fontProp_.get()->setPadding(tgt::vec2(10, 0));
        if(mousePos_.x < (float)viewportSize_.x / 2.0f) {
            if(mousePos_.y >= (float)viewportSize_.y / 2.0f)
                fontProp_.get()->setTextAlignment(tgt::Font::BottomLeft);
            else
                fontProp_.get()->setTextAlignment(tgt::Font::TopLeft);
        } else {
            if(mousePos_.y >= (float)viewportSize_.y / 2.0f)
                fontProp_.get()->setTextAlignment(tgt::Font::BottomRight);
            else
                fontProp_.get()->setTextAlignment(tgt::Font::TopRight);
        }
        fontProp_.get()->render(tgt::vec3((float)mousePos_.x, (float)mousePos_.y, 0), collectText("FOLLOW"), screensize);
    }
    fontProp_.get()->setPadding(0);
}

std::string TextOverlay::collectText(std::string key) {
    std::stringstream ss("");

    if (text0_.isReady() && layout0_.getKey() == key) {
        for(std::string str: text0_.getAllData()) {
            if (ss.str().length() > 0)
                ss << "\n";
            ss << str;
        }
    }
    if (text1_.isReady() && layout1_.getKey() == key) {
        for(std::string str: text1_.getAllData())  {
            if (ss.str().length() > 0)
                ss << "\n";
            ss << str;
        }
    }
    if (text2_.isReady() && layout2_.getKey() == key) {
        for(std::string str: text2_.getAllData()) {
            if (ss.str().length() > 0)
                ss << "\n";
            ss << str;
        }
    }
    if (text3_.isReady() && layout3_.getKey() == key) {
        for(std::string str: text3_.getAllData()){
            if (ss.str().length() > 0)
                ss << "\n";
            ss << str;
        }
    }
    return ss.str();
}

int TextOverlay::getNumberOfLines(std::string s) {
    int numLines = 0;
    std::string line;
    std::stringstream ss(s);
    while (std::getline(ss, line))
        numLines++;

    return numLines;
}

void TextOverlay::mouseMove(tgt::MouseEvent* e) {
    viewportSize_ = tgt::ivec2(e->viewport().x, e->viewport().y);
    mousePos_ = tgt::ivec2(e->coord().x, e->viewport().y - e->coord().y);
    invalidate();
    e->ignore();
}

void TextOverlay::mouseEnterExit(tgt::MouseEvent* e) {
    renderFollowMouseText_ = e->action() == tgt::MouseEvent::ENTER;
    invalidate();
    e->ignore();
}

} // namespace
