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

#include "touchtableoverlay.h"

#include "tgt/textureunit.h"
#include "tgt/material.h" //bp3d
#include "tgt/immediatemode/immediatemode.h"

#include "widgets/touchtablemenuwidget.h"

#define _USE_MATH_DEFINES
#include <math.h>


namespace voreen {

const std::string TouchTableOverlay::loggerCat_("voreen.touchtable.TouchTableOverlay");

#define M_PI_2 1.57079632679489661923

using tgt::vec2;
using tgt::ivec2;
using tgt::vec4;
using tgt::TouchPoint;
using tgt::TextureUnit;


TouchTableOverlay::TouchTableOverlay()
    : ImageProcessor("touchinteraction/touchtableoverlay")
    , inport_(Port::INPORT, "image.in", "Image Input")
    , outport_(Port::OUTPORT, "image.out", "Image Output")
    , privateRenderPort_(Port::OUTPORT,"PrivateRenderPort","Private Renderport") //bp3d
    , widgetPort_(Port::INPORT, "widgetport", "Widget Port", true)
    , widgetRadius_("widgetradius", "Widget Radius", 40, 5, 400)
    , menuRadius_("menuradius", "Menu Radius", 60, 5, 400)
    , movementThreshold_("threshold", "Movement Threshold", 40, 5, 400)
    , menuOpen_("menuopen", "Menu Open", false)
    , moveView_("moveview", "Move View", false)
    , isInExclusiveMode_("exclusivemode", "Exclusive Mode", false)
    , quadShader_(0)
    , controlElementShader_(0)
    , currentMenuWidget_(0)
    , overlayMenu_(tgt::ivec2(0,0), tgt::vec4(0.4f, 0.4f, 0.4f, 0.5f))
    , shiftDiameter_("shiftdiameter", "Shift Diameter", 400.f, FLT_MIN, FLT_MAX)
{
    addPort(inport_);
    addPort(outport_);
    addPort(widgetPort_);
    addPrivateRenderPort(privateRenderPort_); //bp3d

    addProperty(widgetRadius_);
    addProperty(menuRadius_);
    addProperty(movementThreshold_);
    addProperty(menuOpen_);
    addProperty(moveView_);
    addProperty(isInExclusiveMode_);

    addProperty(shiftDiameter_);

    shaderProp_.setVisibleFlag(false);

    //initialize Menu
    overlayMenu_.setID(-1);
    overlayMenu_.setLL(tgt::ivec2(0,0));
    overlayMenu_.setUR(tgt::ivec2(100, 50));
}

const tgt::vec4 TouchTableOverlay::getGlobalMenuColor() {
    return tgt::vec4(0.4f, 0.4f, 0.4f, 0.5f);
}

Processor* TouchTableOverlay::create() const {
    return new TouchTableOverlay();
}

void TouchTableOverlay::initialize() {
    ImageProcessor::initialize();

    // load all shader programs
    controlElementShader_ = ShdrMgr.loadSeparate("passthrough.vert", "touchinteraction/touchtablecontrolelement.frag", ImageProcessor::generateHeader(), false);
    quadShader_ = ShdrMgr.loadSeparate("passthrough.vert", "touchinteraction/touchtablequad.frag", ImageProcessor::generateHeader(), false);
    quadShader1D_ = ShdrMgr.loadSeparate("passthrough.vert", "touchinteraction/touchtablequad1D.frag", ImageProcessor::generateHeader(), false);

    // load image files as textures
    TexMgr.addPath(VoreenApplication::app()->getModulePath("touchtable") + "/textures");
    puckOnTex_ = TexMgr.load("newPuckOn2.png");
    puckOffTex_ = TexMgr.load("newPuckOff2.png");
    alternateControlTex_ = TexMgr.load("alternate.png");
    menuButtonTex_ = TexMgr.load("menuButton.png");
    innerTex_ = TexMgr.load("innerTex.png");

    sliderMidTex_ = TexMgr.load("slidermid.png");
    sliderIndicatorTex_ = TexMgr.load("sliderbutton.png");
    sliderLeftEndTex_ = TexMgr.load("sliderleftend.png");
    sliderRightEndTex_ = TexMgr.load("sliderrightend.png");

    checkerBoardTex_ = TexMgr.load("checkerboard.png");

    privateRenderPort_.changeFormat(GL_RGBA); //bp3d
    privateRenderPort_.resize(250,250);
}

void TouchTableOverlay::deinitialize() {

    // free all ressources (shaders, textures)
    ShdrMgr.dispose(program_);
    ShdrMgr.dispose(controlElementShader_);
    ShdrMgr.dispose(quadShader_);
    ShdrMgr.dispose(quadShader1D_);
    program_ = 0;
    controlElementShader_ = 0;
    quadShader_ = 0;
    quadShader1D_ = 0;

    TexMgr.dispose(puckOnTex_);
    TexMgr.dispose(puckOffTex_);
    TexMgr.dispose(menuButtonTex_);
    TexMgr.dispose(innerTex_);
    TexMgr.dispose(alternateControlTex_);
    TexMgr.dispose(sliderMidTex_);
    TexMgr.dispose(sliderIndicatorTex_);
    TexMgr.dispose(sliderLeftEndTex_);
    TexMgr.dispose(sliderRightEndTex_);

    TexMgr.dispose(checkerBoardTex_);

    // clear widget lists
    widgetsOnTable_.clear();
    widgetsInMenu_.clear();

    // call parent class deinitialize
    ImageProcessor::deinitialize();
}

bool TouchTableOverlay::isReady() const {
    return (inport_.isReady() && outport_.isReady());
}

void TouchTableOverlay::beforeProcess() {

    ImageProcessor::beforeProcess();

    //set positions for menu buttons
    //TODO: find a better solution for this
    overlayMenu_.clearButtons();
    overlayMenu_.addButton(vec2(0,0));
    overlayMenu_.addButton(outport_.getSize());

    ImageProcessor::beforeProcess();
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();

    if (widgetPort_.hasChanged()) {

        // if a widget has been removed or added: update lists of connected widgets
        if (widgetPort_.getConnectedProcessors().size() != widgetsOnTable_.size() + widgetsInMenu_.size()) {
            widgetsOnTable_.clear();
            widgetsInMenu_.clear();
            //TODO: find better solution to avoid crashes when widget is deleted
            currentMenuWidget_ = 0;

            //get list of connected widgets
            std::vector<TouchTableWidget*> widgetVector = widgetPort_.getConnectedProcessors();

            std::vector<TouchTableWidget*>::iterator iter;

            for (iter = widgetVector.begin(); iter != widgetVector.end(); ++iter) {
                // put widgets in
                if ((*iter)->isOnTable())
                    widgetsOnTable_.push_back(*iter);
                else
                    widgetsInMenu_.push_back(*iter);

                (*iter)->setOverlay(this);

                //this fixes two problems:
                // 1) an activated menu widget cannot set itself as the current menu widget at the beginning (ie. when deserializing a workspace and the widget is activated)
                // 2) if a widget is removed this leads to setting the current mode widget to 0, even if the activated menu widget is not the one that has been removed
                //since we are here, the number of widgets has changed and the current mode widget is set to 0
                //-> look for a menu widget that is activated and set it as current menu widget
                //TODO: implement better handling
                TouchTableMenuWidget* menuWidget = dynamic_cast<TouchTableMenuWidget*>(*iter);
                if (menuWidget && menuWidget->isActive())
                    currentMenuWidget_ = menuWidget;
            }
        }
    }
}

void TouchTableOverlay::process() {

    tgtAssert(outport_.isReady(), "Outport not ready");
    tgtAssert(inport_.isReady(), "Inport not ready");
    tgtAssert(program_, "Shader missing");
    tgtAssert(quadShader_, "touchtablequad.frag shader is missing");

    outport_.activateTarget();
    outport_.clearTarget();

    // pass through input rendering
    renderInputImage();

    LGL_ERROR;

    //render menu buttons
    renderMenuButtons();

    //render widgets on table
    renderWidgetsOnTable();

    if (currentMenuWidget_ != 0)
       currentMenuWidget_->render();

    if (menuOpen_.get()) {
        renderMenu();
        renderWidgetsInMenu();
    }

    outport_.deactivateTarget();

    LGL_ERROR;

}

void TouchTableOverlay::renderInputImage() {

    // bind input rendering to texture units
    TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());

    // copy input image to outport
    program_->activate();
    setGlobalShaderParameters(program_);
    inport_.setTextureParameters(program_, "texParams_");

    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());

    renderQuad();

    program_->deactivate();

    tgt::TextureUnit::setZeroUnit();

    inport_.getColorTexture()->disable();
    inport_.getDepthTexture()->disable();
    LGL_ERROR;

}

void TouchTableOverlay::renderMenuButtons() {

    std::vector<ivec2>::iterator it;
    for (it = overlayMenu_.getButtons()->begin(); it != overlayMenu_.getButtons()->end(); ++it) {

        //make menu button slightly transparent
        float opacity = 0.5f;
        tgt::vec4 colorMod = tgt::vec4(0.f);

        //set coordinates
        tgt::ivec2 ll = tgt::ivec2(*it) - tgt::ivec2(menuRadius_.get());
        tgt::ivec2 ur = tgt::ivec2(*it) + tgt::ivec2(menuRadius_.get());

        //render button
        renderTexturedQuad(ll, ur, menuButtonTex_, colorMod, opacity, tgt::mat4::identity, 0.f, false);
    }

    LGL_ERROR;
}

void TouchTableOverlay::renderWidgetsInMenu() {

    if (!menuOpen_.get())
        return;

    int i = 0; //counter for current widget position

    std::vector<TouchTableWidget*>::iterator iter;
    for (iter = widgetsInMenu_.begin(); iter != widgetsInMenu_.end(); ++iter) {

        TouchTableWidget* widget = *iter;

        //TODO: change this...
        int xMod = (overlayMenu_.getLL().x == 0) ? menuRadius_.get() : 0;

        widget->setPosition(vec2(static_cast<float>(xMod + overlayMenu_.getLL().x + 5 + widgetRadius_.get() + i*(2*widgetRadius_.get() + 10)), static_cast<float>((overlayMenu_.getUR().y - overlayMenu_.getLL().y) / 2 + overlayMenu_.getLL().y)));
        renderWidget(widget);

        i++;
    }
}

void TouchTableOverlay::renderWidgetsOnTable() {

    //vector for widgets in motion -> render at end
    std::vector<TouchTableWidget*> motionWidgets;

    std::vector<TouchTableWidget*>::iterator iter;
    for (iter = widgetsOnTable_.begin(); iter != widgetsOnTable_.end(); ++iter) {

        TouchTableWidget* widget = *iter;

        // if the widget is in motion: render later
        if (widget->isInMotion()) {
            motionWidgets.push_back(widget);
            continue;
        }

        renderWidget(widget);

    }

    for (iter = motionWidgets.begin(); iter != motionWidgets.end(); ++iter) {
        TouchTableWidget* widget = *iter;
        renderWidget(widget);
    }
}

void TouchTableOverlay::renderTexturedQuad(tgt::ivec2 ll, tgt::ivec2 ur, const tgt::Texture* tex, vec4 colorMod, float opacity, const tgt::mat4& textureMatrix, float greyOutFactor, bool useMenuWidgetOrientationInvert) {

    // set OpenGL status: depth func and blending
    glDepthFunc(GL_ALWAYS);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //set transformation to use pixel coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, 1.0f*outport_.getSize().x, 0.0f, 1.0f*outport_.getSize().y, -1.0f, 1.0f));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    tgt::ivec2 shaderLL = ll;
    tgt::ivec2 shaderUR = ur;
    //if menu of the mode widget is inverted and this should be taken into account: set transformation
    if (currentMenuWidget_ && currentMenuWidget_->isMenuInverted() && useMenuWidgetOrientationInvert) {

        //set OpenGL traansformation
        float xt = - (static_cast<float>(currentMenuWidget_->getMenuLL().x) + static_cast<float>(currentMenuWidget_->getMenuDimensions().x) / 2.f);
        float yt = - (static_cast<float>(currentMenuWidget_->getMenuLL().y) + static_cast<float>(currentMenuWidget_->getMenuDimensions().y) / 2.f);

        MatStack.translate(-xt,-yt, 0.f);
        MatStack.rotate(180.f, 0.f, 0.f, 1.f);
        MatStack.translate(xt, yt, 0.f);

        //compute shader parameters
        shaderLL -= currentMenuWidget_->getMenuLL();
        shaderUR -= currentMenuWidget_->getMenuLL();
        shaderLL = currentMenuWidget_->getMenuDimensions() - shaderLL;
        shaderUR = currentMenuWidget_->getMenuDimensions() - shaderUR;
        shaderLL += currentMenuWidget_->getMenuLL();
        shaderUR += currentMenuWidget_->getMenuLL();
    }

    tgt::TextureUnit texUnit;
    texUnit.activate();
    tex->bind();

    quadShader_->activate();
    setGlobalShaderParameters(quadShader_);

    //set the texture
    quadShader_->setUniform("tex_", texUnit.getUnitNumber());

    //set texture parameters
    bool oldIgnoreError = quadShader_->getIgnoreUniformLocationError();
    quadShader_->setIgnoreUniformLocationError(true);
    quadShader_->setUniform("texParams_.dimensions_", vec2(tex->getDimensions().xy()));
    quadShader_->setUniform("texParams_.dimensionsRCP_", vec2(1.f) / vec2(tex->getDimensions().xy()));
    quadShader_->setUniform("texParams_.matrix_", textureMatrix);
    quadShader_->setIgnoreUniformLocationError(oldIgnoreError);

    //set uniforms
    quadShader_->setUniform("ll_", tgt::vec2(shaderLL));
    quadShader_->setUniform("ur_", tgt::vec2(shaderUR));
    quadShader_->setUniform("opacity_", opacity);
    quadShader_->setUniform("colorMod_", colorMod);
    quadShader_->setUniform("greyOutFactor_", tgt::clamp(greyOutFactor, 0.f, 1.f));

    //render quad
    glBegin(GL_QUADS);
        glVertex2f(static_cast<float>(ll.x), static_cast<float>(ll.y));
        glVertex2f(static_cast<float>(ur.x), static_cast<float>(ll.y));
        glVertex2f(static_cast<float>(ur.x), static_cast<float>(ur.y));
        glVertex2f(static_cast<float>(ll.x), static_cast<float>(ur.y));
    glEnd();

    tgt::TextureUnit::setZeroUnit();
    tex->disable();

    quadShader_->deactivate();

    // set back OpenGL status
    glDepthFunc(GL_LESS);
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);

    //set transformation to default
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
}

void TouchTableOverlay::renderTexturedQuad1DTexture(tgt::ivec2 ll, tgt::ivec2 ur, const tgt::Texture* tex, tgt::vec4 colorMod, float opacity, float greyOutFactor, bool useMenuWidgetOrientationInvert) {

    // set OpenGL status: depth func and blending
    glDepthFunc(GL_ALWAYS);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //set transformation to use pixel coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, 1.0f*outport_.getSize().x, 0.0f, 1.0f*outport_.getSize().y, -1.0f, 1.0f));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    //if menu of the mode widget is inverted and this should be taken into account: set transformation
    tgt::ivec2 shaderLL = ll;
    tgt::ivec2 shaderUR = ur;
    if (currentMenuWidget_ && currentMenuWidget_->isMenuInverted() && useMenuWidgetOrientationInvert) {
        float xt = - (static_cast<float>(currentMenuWidget_->getMenuLL().x) + static_cast<float>(currentMenuWidget_->getMenuDimensions().x) / 2.f);
        float yt = - (static_cast<float>(currentMenuWidget_->getMenuLL().y) + static_cast<float>(currentMenuWidget_->getMenuDimensions().y) / 2.f);

        MatStack.translate(-xt,-yt, 0.f);
        MatStack.rotate(180.f, 0.f, 0.f, 1.f);
        MatStack.translate(xt,yt, 0.f);

        //compute shader parameters
        shaderLL -= currentMenuWidget_->getMenuLL();
        shaderUR -= currentMenuWidget_->getMenuLL();
        shaderLL = currentMenuWidget_->getMenuDimensions() - shaderLL;
        shaderUR = currentMenuWidget_->getMenuDimensions() - shaderUR;
        shaderLL += currentMenuWidget_->getMenuLL();
        shaderUR += currentMenuWidget_->getMenuLL();
    }

    tgt::TextureUnit texUnit;
    texUnit.activate();
    tex->bind();

    quadShader1D_->activate();
    setGlobalShaderParameters(quadShader1D_);

    //set the texture
    quadShader1D_->setUniform("tex_", texUnit.getUnitNumber());

    //set uniforms
    quadShader1D_->setUniform("ll_", tgt::vec2(shaderLL));
    quadShader1D_->setUniform("ur_", tgt::vec2(shaderUR));
    quadShader1D_->setUniform("opacity_", opacity);
    quadShader1D_->setUniform("colorMod_", colorMod);
    quadShader1D_->setUniform("greyOutFactor_", tgt::clamp(greyOutFactor, 0.f, 1.f));

    //render quad
    glBegin(GL_QUADS);
        glVertex2f(static_cast<float>(ll.x), static_cast<float>(ll.y));
        glVertex2f(static_cast<float>(ur.x), static_cast<float>(ll.y));
        glVertex2f(static_cast<float>(ur.x), static_cast<float>(ur.y));
        glVertex2f(static_cast<float>(ll.x), static_cast<float>(ur.y));
    glEnd();

    tgt::TextureUnit::setZeroUnit();
    tex->disable();

    quadShader1D_->deactivate();

    // set back OpenGL status
    glDepthFunc(GL_LESS);
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);

    //set transformation to default
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
}

void TouchTableOverlay::renderColoredQuad(tgt::ivec2 ll, tgt::ivec2 ur, vec4 color, bool useMenuWidgetOrientationInvert) const {

    // set OpenGL status: depth func and blending
    glDepthFunc(GL_ALWAYS);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //set transformation to use pixel coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, 1.0f*outport_.getSize().x, 0.0f, 1.0f*outport_.getSize().y, -1.0f, 1.0f));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    //if menu of the mode widget is inverted and this should be taken into account: set transformation
    if (currentMenuWidget_ && currentMenuWidget_->isMenuInverted() && useMenuWidgetOrientationInvert) {
        float xt = - (static_cast<float>(currentMenuWidget_->getMenuLL().x) + static_cast<float>(currentMenuWidget_->getMenuDimensions().x) / 2.f);
        float yt = - (static_cast<float>(currentMenuWidget_->getMenuLL().y) + static_cast<float>(currentMenuWidget_->getMenuDimensions().y) / 2.f);

        MatStack.translate(-xt,-yt, 0.f);
        MatStack.rotate(180.f, 0.f, 0.f, 1.f);
        MatStack.translate(xt, yt, 0.f);
    }

    glColor4f(color.r, color.g, color.b, color.a);

    //render quad
    glBegin(GL_QUADS);
        glVertex2f(static_cast<float>(ll.x), static_cast<float>(ll.y));
        glVertex2f(static_cast<float>(ur.x), static_cast<float>(ll.y));
        glVertex2f(static_cast<float>(ur.x), static_cast<float>(ur.y));
        glVertex2f(static_cast<float>(ll.x), static_cast<float>(ur.y));
    glEnd();

    // set back OpenGL status
    glColor4f(1.f,1.f,1.f, 1.f);
    glDepthFunc(GL_LESS);
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);

    //set transformation to default
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
}

void TouchTableOverlay::renderWidget(TouchTableWidget* widget) {

    //select basic texture depending on the status of the widget
    tgt::Texture* puckTex;
    //TODO: load separate texture if widget is not activatable
    if (widget->isActive())
       puckTex = puckOnTex_;
    else
       puckTex = puckOffTex_;

    //if the widget is currently moved: render transparent
    float opacity = (widget->isInMotion()) ? 0.5f : 1.f;

    vec4 colorMod = vec4(0.f);
    //if position is invalid: set colorMod to red
    if(widget->isPositionInvalid())
        colorMod = vec4(0.5f, -0.5f, -0.5f, 0.f);

    // if widget is in motion and over menu: set color mod to blue
    if(menuOpen_.get() && (isWidgetInMenu(widget->getPosition(), 1.0f*widgetRadius_.get())) && (widget->isInMotion())){
        colorMod = vec4(-0.3f, -0.3f, 0.5f, 0.0f);
    }

    // if widget is in menu: set color mod to slighty blue / gray
    if (!widget->isOnTable())
        colorMod = vec4(-0.25f, -0.25f, 0.1f, -0.1f);

    int radius = widgetRadius_.get();
    tgt::ivec2 pos = widget->getPosition();

    //render puck
    renderTexturedQuad(pos - tgt::ivec2(radius), pos + tgt::ivec2(radius), puckTex, colorMod, opacity, tgt::mat4::identity, 0.f, false);

    //render symbol on top
    if (!widget->getSymbolTexture())
        LERROR("Could not render symbol for widget: symbol texture not set");
    else
        renderTexturedQuad(pos - tgt::ivec2(radius/2), pos + tgt::ivec2(radius/2), widget->getSymbolTexture(), colorMod, opacity, tgt::mat4::identity, 0.f, false);

}

void TouchTableOverlay::renderMenu() {
    renderColoredQuad(overlayMenu_.getLL(), overlayMenu_.getUR(), overlayMenu_.getColor(), false);
}

void TouchTableOverlay::renderMenuFrame(const  TouchTableMenuFrame& menu, bool useMenuWidgetOrientationInvert, ivec2 offset){
    renderColoredQuad(menu.getLL() + offset, menu.getUR() + offset, menu.getColor(), useMenuWidgetOrientationInvert);
}

void TouchTableOverlay::renderControlElement(TouchTableControlElement* controlElement, tgt::ivec2 menuLL, float greyOutFactor, bool useMenuWidgetOrientationInvert)  {

    // set OpenGL status: depth func and blending
    glDepthFunc(GL_ALWAYS);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    //set transformation to use pixel coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, 1.0f*outport_.getSize().x, 0.0f, 1.0f*outport_.getSize().y, -1.0f, 1.0f));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    //if menu of the mode widget is inverted and this should be taken into account: set transformation
    if (currentMenuWidget_ && currentMenuWidget_->isMenuInverted() && useMenuWidgetOrientationInvert) {
        float xt = - (static_cast<float>(currentMenuWidget_->getMenuLL().x) + static_cast<float>(currentMenuWidget_->getMenuDimensions().x) / 2.f);
        float yt = - (static_cast<float>(currentMenuWidget_->getMenuLL().y) + static_cast<float>(currentMenuWidget_->getMenuDimensions().y) / 2.f);

        MatStack.translate(- xt,- yt, 0.f);
        MatStack.rotate(180.f, 0.f, 0.f, 1.f);
        MatStack.translate(xt, yt, 0.f);
    }

    //activate shader
    controlElementShader_->activate();
    setGlobalShaderParameters(controlElementShader_);

    //bind textures
    tgt::TextureUnit texUnit, symbUnit, innerUnit;

    //get background texture
    texUnit.activate();

    tgt::Texture* puckTex;
    if (!controlElement->isActivatable())
        puckTex = alternateControlTex_;
    else if (controlElement->isActive())
        puckTex = puckOnTex_;
    else
        puckTex = puckOffTex_;

    puckTex->bind();
    //puckTex->enable();

    //get the inner Texture
    innerUnit.activate();
    innerTex_->bind();
    //innerTex_->enable();

    float opacity = controlElement->getOpacity();

    vec4 colorMod = controlElement->getColorMod();

    //set the uniforms
    controlElementShader_->setUniform("opacity_", opacity);
    controlElementShader_->setUniform("puckTex_", texUnit.getUnitNumber());
    controlElementShader_->setUniform("innerTex_", innerUnit.getUnitNumber());
    controlElementShader_->setUniform("colorMod_", colorMod);
    controlElementShader_->setUniform("greyOutFactor_", tgt::clamp(greyOutFactor, 0.f, 1.f));

    //set texture parameters
    bool oldIgnoreError = controlElementShader_->getIgnoreUniformLocationError();
    controlElementShader_->setIgnoreUniformLocationError(true);
    controlElementShader_->setUniform("texParams_.dimensions_", vec2(puckTex->getDimensions().xy()));
    controlElementShader_->setUniform("texParams_.dimensionsRCP_", vec2(1.f) / vec2(puckTex->getDimensions().xy()));
    controlElementShader_->setUniform("texParams_.matrix_", tgt::mat4::identity);

    //set position and radius of the control element as uniform
    vec2 pos = tgt::vec2(controlElement->getPos());
    float radius = static_cast<float>(controlElement->getRadius());

    //calculate position relative to the ModeMenu
    vec2 posInMenu = pos + tgt::vec2(menuLL);

    tgt::vec2 ll(posInMenu - vec2(radius));
    tgt::vec2 ur (posInMenu + vec2(radius));

    if (currentMenuWidget_ && currentMenuWidget_->isMenuInverted() && useMenuWidgetOrientationInvert) {
        ll -= vec2(currentMenuWidget_->getMenuLL());
        ur -= vec2(currentMenuWidget_->getMenuLL());
        ll = vec2(currentMenuWidget_->getMenuDimensions()) - ll;
        ur = vec2(currentMenuWidget_->getMenuDimensions()) - ur;
        ll += vec2(currentMenuWidget_->getMenuLL());
        ur += vec2(currentMenuWidget_->getMenuLL());
    }

    controlElementShader_->setUniform("ll_", ll);
    controlElementShader_->setUniform("ur_", ur);
    controlElementShader_->setIgnoreUniformLocationError(oldIgnoreError);

    //render quad for overlay
    glBegin(GL_QUADS);
        glVertex2f(posInMenu.x - radius, posInMenu.y - radius);
        glVertex2f(posInMenu.x + radius, posInMenu.y - radius);
        glVertex2f(posInMenu.x + radius, posInMenu.y + radius);
        glVertex2f(posInMenu.x - radius, posInMenu.y + radius);
    glEnd();

    controlElementShader_->deactivate();

    tgt::TextureUnit::setZeroUnit();

    puckTex->disable();
    innerTex_->disable();
    alternateControlTex_->disable();


    // set back OpenGL status
    glDepthFunc(GL_LESS);
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);

    //set transformation back
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    //render symbol on top
    if (!controlElement->getSymbolTexture())
        LERROR("Could not render control element symbol: no symbol texture set");
    else {
        tgt::ivec2 ll = tgt::ivec2(static_cast<int>(posInMenu.x - radius/2.0f), static_cast<int>(posInMenu.y - radius/2.0f));
        tgt::ivec2 ur = tgt::ivec2(static_cast<int>(posInMenu.x + radius/2.0f), static_cast<int>(posInMenu.y + radius/2.0f));
        renderTexturedQuad(ll, ur, controlElement->getSymbolTexture(), tgt::vec4(0.f), opacity, tgt::mat4::identity, greyOutFactor, useMenuWidgetOrientationInvert);
    }
}

void TouchTableOverlay::renderSlider(TouchTableTwoIndicatorSlider* slider, tgt::ivec2 menuLL, float greyOutFactor) {

    //get slider orientation
    Orientation orientation = slider->getOrientation();

    //render bar
    if (orientation == HORIZONTAL)
        renderTexturedQuad(slider->getBarLL() + menuLL, slider->getBarUR() + menuLL, sliderMidTex_, tgt::vec4(0.f), 1.f, tgt::mat4::identity, greyOutFactor);
    else
        renderTexturedQuad(slider->getBarLL() + menuLL, slider->getBarUR() + menuLL, sliderMidTex_, tgt::vec4(0.f), 1.f, tgt::mat4::createRotationZ(static_cast<float>(M_PI_2)), greyOutFactor);

    //get indicator positions
    std::deque<tgt::ivec2> positions = slider->getIndicatorPositions();

    //render colored line between the indicators
    if (orientation == HORIZONTAL)
        renderColoredQuad(tgt::ivec2(positions.at(1).x + menuLL.x, slider->getBarLL().y + menuLL.y) , tgt::ivec2(positions.at(2).x + menuLL.x, slider->getBarUR().y + menuLL.y), greyOutFactor * tgt::vec4(tgt::vec3(0.3f), 0.7f) + (1.f - greyOutFactor) * tgt::vec4(0.f, 0.9f, 0.2f, 0.7f));
    else
        renderColoredQuad(tgt::ivec2(slider->getBarLL().x + menuLL.x, positions.at(1).y + menuLL.y), tgt::ivec2(slider->getBarUR().x + menuLL.x, positions.at(2).y + menuLL.y), greyOutFactor * tgt::vec4(tgt::vec3(0.3f), 0.7f) + (1.f - greyOutFactor) * tgt::vec4(0.f, 0.9f, 0.2f, 0.7f));

    //render end pieces
    if (orientation == HORIZONTAL) {
        renderTexturedQuad(slider->getBarLL() + menuLL - tgt::ivec2(slider->getIndicatorWidth()-1, 0), tgt::ivec2(slider->getBarLL().x, slider->getBarUR().y) + menuLL, sliderLeftEndTex_, tgt::vec4(0.f), 1.f, tgt::mat4::identity, greyOutFactor);
        renderTexturedQuad(tgt::ivec2(slider->getBarUR().x, slider->getBarLL().y) + menuLL, slider->getBarUR() + menuLL + tgt::ivec2(slider->getIndicatorWidth()-1, 0), sliderRightEndTex_, tgt::vec4(0.f), 1.f, tgt::mat4::identity, greyOutFactor);
    } else {
        renderTexturedQuad(slider->getBarLL() + menuLL - tgt::ivec2(0, slider->getIndicatorWidth()-1), tgt::ivec2(slider->getBarUR().x, slider->getBarLL().y) + menuLL, sliderLeftEndTex_, tgt::vec4(0.f), 1.f, tgt::mat4::createRotationZ(static_cast<float>(M_PI_2)), greyOutFactor);
        renderTexturedQuad(tgt::ivec2(slider->getBarLL().x, slider->getBarUR().y) + menuLL, slider->getBarUR() + menuLL + tgt::ivec2(0, slider->getIndicatorWidth()-1), sliderRightEndTex_, tgt::vec4(0.f), 1.f, tgt::mat4::createRotationZ(static_cast<float>(M_PI_2)), greyOutFactor);
    }

    //render indicators
    renderTexturedQuad(positions.at(0) + menuLL, positions.at(1) + menuLL, sliderIndicatorTex_, tgt::vec4(0.f), 1.f, tgt::mat4::identity, greyOutFactor);
    renderTexturedQuad(positions.at(2) + menuLL, positions.at(3) + menuLL, sliderIndicatorTex_, tgt::vec4(0.f), 1.f, tgt::mat4::identity, greyOutFactor);
}

void TouchTableOverlay::renderSlider(TouchTableSlider* slider, tgt::ivec2 menuLL, float greyOutFactor){

    //get slider orientation
    Orientation orientation = slider->getOrientation();

    //render bar
    if (orientation == HORIZONTAL)
        renderTexturedQuad(slider->getBarLL() + menuLL + tgt::ivec2(20, 0), slider->getBarUR() + menuLL + tgt::ivec2(-20, 0), sliderMidTex_, tgt::vec4(0.f), 1.f, tgt::mat4::identity, greyOutFactor);
    else
        renderTexturedQuad(slider->getBarLL() + menuLL + tgt::ivec2(0, 20), slider->getBarUR() + menuLL  + tgt::ivec2(0, -20), sliderMidTex_, tgt::vec4(0.f), 1.f, tgt::mat4::createRotationZ(static_cast<float>(M_PI_2)), greyOutFactor);


    //render end pieces
    if (orientation == HORIZONTAL) {
           renderTexturedQuad(slider->getBarLL() + menuLL, tgt::ivec2(slider->getBarLL().x, slider->getBarUR().y) + menuLL + tgt::ivec2(20, 0), sliderLeftEndTex_, tgt::vec4(0.f), 1.f, tgt::mat4::identity, greyOutFactor);
        renderTexturedQuad(tgt::ivec2(slider->getBarUR().x, slider->getBarLL().y) + menuLL - tgt::ivec2(20, 0), tgt::ivec2(slider->getBarUR()) + menuLL,sliderRightEndTex_, tgt::vec4(0.f), 1.f, tgt::mat4::identity, greyOutFactor);
    } else {
        renderTexturedQuad(slider->getBarLL() + menuLL - tgt::ivec2(0, -20), tgt::ivec2(slider->getBarUR().x, slider->getBarLL().y) + menuLL, sliderLeftEndTex_, tgt::vec4(0.f), 1.f, tgt::mat4::createRotationZ(static_cast<float>(M_PI_2)), greyOutFactor);
        renderTexturedQuad(tgt::ivec2(slider->getBarLL().x, slider->getBarUR().y -20) + menuLL, slider->getBarUR() + menuLL, sliderLeftEndTex_, tgt::vec4(0.f), 1.f, tgt::mat4::createRotationZ(static_cast<float>(M_PI_2)), greyOutFactor);
    }

    //render indicator
    renderTexturedQuad(slider->getIndicatorLL() + menuLL, slider->getIndicatorUR() + menuLL, sliderIndicatorTex_, tgt::vec4(0.f), 1.f, tgt::mat4::identity, greyOutFactor);
}

void TouchTableOverlay::renderScrollableMenu(TouchTableScrollableMenu* menu, ivec2 offset){

    renderMenuFrame(*menu, true, offset);
    renderMenuFrame(menu->getContentMenu(), true, offset);
    float grayOutFactor = 0.f;
    if(menu->getContent().size() <= menu->getMaxElementsInMenu())
        grayOutFactor = 0.8;
    renderSlider(&menu->getSlider(), offset, grayOutFactor);

    //render volume file names

    int j=0;

    tgt::Font* font = menu->getFont();
    font->setTextAlignment(tgt::Font::Centered);

    for (size_t i = menu->getContentOffset(); i < tgt::clamp(menu->getContentOffset() + menu->getMaxElementsInMenu(),0, static_cast<int>(menu->getContent().size())); ++i, ++j) {
        tgt::vec3 position(static_cast<float>(menu->getContentMenu().getLL().x + offset.x), static_cast<float>(menu->getContentMenu().getUR().y - j * menu->getRowHeight() + offset.y - menu->getRowHeight()/2), 0.f);
        //set transformation to use pixel coordinates
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, 1.0f*outport_.getSize().x, 0.0f, 1.0f*outport_.getSize().y, -1.0f, 1.0f));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadIdentity();
        //if menu of the mode widget is inverted and this should be taken into account: set transformation
        if (currentMenuWidget_ && currentMenuWidget_->isMenuInverted()) {
            float xt = - (static_cast<float>(currentMenuWidget_->getMenuLL().x) + static_cast<float>(currentMenuWidget_->getMenuDimensions().x) / 2.f);
            float yt = - (static_cast<float>(currentMenuWidget_->getMenuLL().y) + static_cast<float>(currentMenuWidget_->getMenuDimensions().y) / 2.f);

            MatStack.translate(-xt,-yt, 0.f);
            MatStack.rotate(180.f, 0.f, 0.f, 1.f);
            MatStack.translate(xt,yt, 0.f);
        }
        glDisable(GL_DEPTH_TEST);

        std::string stringToRender =  menu->getContent().at(i).first;

        // truncate text if it does not fix
        float maxContentMenuLength = static_cast<float>(menu->getContentMenu().getUR().x - menu->getContentMenu().getLL().x);
        float textLength = font->getSize(position,stringToRender,tgt::ivec2(outport_.getSize())).x;
        bool textTruncated = false;
        while(textLength > maxContentMenuLength && stringToRender.size() > 1){
            stringToRender = stringToRender.substr(0, static_cast<int>(stringToRender.size()-1));
            textLength = font->getSize(position, stringToRender+"...",tgt::ivec2(outport_.getSize())).x;
            textTruncated = true;
        }
        if (textTruncated)
            stringToRender.append("...");




        // get text into texture
        tgt::Texture* tex = font->renderToTexture(stringToRender);
        tex->bind();

        // blending
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

        tgt::vec4 color = font->getFontColor();

        // layout and render text
        tgt::ivec2 dim = tex->getDimensions().xy();
        IMode.setTextureMode(tgt::ImmediateMode::TEX2D);
        float lineWidth = menu->getFont()->getLineWidth();
        float offsetx = (lineWidth - dim.x) / 2.0f;
        float offsety = -dim.y / 2.0f;

        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.texcoord(0, 0);
            IMode.vertex(position+tgt::vec3(offsetx, offsety, 0));
            IMode.texcoord(1, 0);
            IMode.vertex(position+tgt::vec3(1.0f*dim.x, 0, 0)+tgt::vec3(offsetx, offsety, 0));
            IMode.texcoord(1, 1);
            IMode.vertex(position+tgt::vec3(1.0f*dim.x, 1.0f*dim.y, 0)+tgt::vec3(offsetx, offsety, 0));
            IMode.texcoord(0, 1);
            IMode.vertex(position+tgt::vec3(0, 1.0f*dim.y, 0)+tgt::vec3(offsetx, offsety, 0));
        IMode.end();
        IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);

        // cleanup
        glBindTexture(GL_TEXTURE_2D, 0);
        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadIdentity();
        delete tex;


        glEnable(GL_DEPTH_TEST);
        if(menu->isVisualizingSelection() && (menu->getSelectedRow() == i))
             renderColoredQuad(tgt::ivec2(menu->getContentMenu().getLL().x + offset.x, static_cast<int>(position.y) - menu->getRowHeight() / 2),
                               tgt::ivec2(menu->getContentMenu().getUR().x + offset.x, static_cast<int>(position.y) + menu->getRowHeight() / 2),
                               tgt::vec4(0.4f, 0.4f, 0.4f, 0.5f));

        //render texture
        if (menu->getContent().at(i).second) {
            //render "frame"
            renderColoredQuad(tgt::ivec2(menu->getContentMenu().getLL().x + offset.x + menu->getTextColumnWidth() - 1 ,
                        static_cast<int>(position.y) - menu->getRowHeight() / 2 + 1)
                    , tgt::ivec2(menu->getContentMenu().getUR().x + offset.x - 1,
                        static_cast<int>(position.y) + menu->getRowHeight() / 2 - 1)
                    , tgt::vec4(tgt::vec3(0.f), 1.f));
            //render checkerboard background
            renderTexturedQuad(
                    tgt::ivec2(menu->getContentMenu().getLL().x + offset.x + menu->getTextColumnWidth(),
                        static_cast<int>(position.y) - menu->getRowHeight() / 2 + 2)
                    , tgt::ivec2(menu->getContentMenu().getUR().x + offset.x - 2,
                        static_cast<int>(position.y) + menu->getRowHeight() / 2 - 2)
                    , checkerBoardTex_);

            //render texture
            if (menu->getContent().at(i).second->getType() == GL_TEXTURE_2D)
                renderTexturedQuad(
                    tgt::ivec2(menu->getContentMenu().getLL().x + offset.x + menu->getTextColumnWidth(),
                        static_cast<int>(position.y) - menu->getRowHeight() / 2 + 2)
                    , tgt::ivec2(menu->getContentMenu().getUR().x + offset.x - 2,
                        static_cast<int>(position.y) + menu->getRowHeight() / 2 - 2)
                    , menu->getContent().at(i).second);
            else if (menu->getContent().at(i).second->getType() == GL_TEXTURE_1D)
                renderTexturedQuad1DTexture(
                    tgt::ivec2(menu->getContentMenu().getLL().x + offset.x + menu->getTextColumnWidth(),
                        static_cast<int>(position.y) - menu->getRowHeight() / 2 + 2)
                    , tgt::ivec2(menu->getContentMenu().getUR().x + offset.x - 2 ,
                        static_cast<int>(position.y) + menu->getRowHeight() / 2 - 2)
                    , menu->getContent().at(i).second);
        }

        if (j < menu->getMaxElementsInMenu()-1) {
            renderColoredQuad(tgt::ivec2(menu->getContentMenu().getLL().x + offset.x, static_cast<int>(position.y) - menu->getRowHeight() / 2 - 1),
                                tgt::ivec2(menu->getContentMenu().getUR().x + offset.x, static_cast<int>(position.y) - menu->getRowHeight() / 2),
                                tgt::vec4(tgt::vec3(0.f), 0.3f));
        }
    }
}

void TouchTableOverlay::renderMenuWidgetMessage(tgt::ivec2 position, tgt::ivec2 menuLL, std::string message, tgt::vec4 color, tgt::Font* font) {

    //set transformation to use pixel coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, 1.0f*outport_.getSize().x, 0.0f, 1.0f*outport_.getSize().y, -1.0f, 1.0f));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    //if menu of the mode widget is inverted and this should be taken into account: set transformation
    if (currentMenuWidget_ && currentMenuWidget_->isMenuInverted()) {
        float xt = - (static_cast<float>(currentMenuWidget_->getMenuLL().x) + static_cast<float>(currentMenuWidget_->getMenuDimensions().x) / 2.f);
        float yt = - (static_cast<float>(currentMenuWidget_->getMenuLL().y) + static_cast<float>(currentMenuWidget_->getMenuDimensions().y) / 2.f);

        MatStack.translate(-xt,-yt, 0.f);
        MatStack.rotate(180.f, 0.f, 0.f, 1.f);
        MatStack.translate(xt,yt, 0.f);
    }
    glDisable(GL_DEPTH_TEST);

    glColor4f(color.x, color.y, color.z, color.w);

    tgt::vec3 fontPos = tgt::vec3(tgt::vec2(position + menuLL), 1.f);



    if (font){
        //font->render(fontPos, message, outport_.getSize());
        // get text into texture
        tgt::Texture* tex = font->renderToTexture(message);
        tex->bind();

        // blending
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

        tgt::vec4 color = font->getFontColor();

        // layout and render text
        tgt::ivec2 dim = tex->getDimensions().xy();
        IMode.setTextureMode(tgt::ImmediateMode::TEX2D);
        float lineWidth = font->getLineWidth();
        float offsetx = (lineWidth - dim.x) / 2.0f;
        float offsety = -dim.y / 2.0f;

        IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.texcoord(0, 0);
        IMode.vertex(fontPos+tgt::vec3(offsetx, offsety, 0));
        IMode.texcoord(1, 0);
        IMode.vertex(fontPos+tgt::vec3(1.0f*dim.x, 0, 0)+tgt::vec3(offsetx, offsety, 0));
        IMode.texcoord(1, 1);
        IMode.vertex(fontPos+tgt::vec3(1.0f*dim.x, 1.0f*dim.y, 0)+tgt::vec3(offsetx, offsety, 0));
        IMode.texcoord(0, 1);
        IMode.vertex(fontPos+tgt::vec3(0, 1.0f*dim.y, 0)+tgt::vec3(offsetx, offsety, 0));
        IMode.end();
        IMode.setTextureMode(tgt::ImmediateMode::TEXNONE);

        // cleanup
        glBindTexture(GL_TEXTURE_2D, 0);
        glDisable(GL_BLEND);
        glEnable(GL_DEPTH_TEST);
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadIdentity();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadIdentity();
        delete tex;
    }


    //clean-up
    glColor4f(1.f,1.f,1.f,1.f);
    glEnable(GL_DEPTH_TEST);
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    glEnable(GL_DEPTH_TEST);
}

void TouchTableOverlay::renderSphere(float radius, float shininess, const tgt::vec4& lightPosition, const tgt::vec4& diffuseLight, const tgt::vec4& specularLight, const tgt::vec4& ambientLight, const tgt::ivec2& position){
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glShadeModel(GL_SMOOTH);

    GLUquadricObj* quadric = gluNewQuadric();
    LGL_ERROR;

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight.elem);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight.elem);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight.elem);
    glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 128.f);

    tgt::vec4 color = tgt::vec4(0.5f,0.5f,0.5f,1.0f);

    tgt::Material material(color, color, color, shininess);
    material.activate();


    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, 1.0f*outport_.getSize().x, 0.0f, 1.0f*outport_.getSize().y, -radius, radius));
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition.elem);
    LGL_ERROR;

    //replace with translation?
    if (currentMenuWidget_ && currentMenuWidget_->isMenuInverted()) {
        float xt = - (static_cast<float>(currentMenuWidget_->getMenuLL().x) + static_cast<float>(currentMenuWidget_->getMenuDimensions().x) / 2.f);
        float yt = - (static_cast<float>(currentMenuWidget_->getMenuLL().y) + static_cast<float>(currentMenuWidget_->getMenuDimensions().y) / 2.f);

        MatStack.translate(-xt,-yt, 0.f);
        MatStack.rotate(180.f, 0.f, 0.f, 1.f);
        MatStack.translate(xt, yt, 0.f);
    }

    MatStack.translate(1.0f*position.x, 1.0f*position.y, 0.0f);
    gluSphere(quadric, radius, 64, 32);

    LGL_ERROR;

    //set transformation to default
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    glPopAttrib();

    gluDeleteQuadric(quadric);
    LGL_ERROR;
}

void TouchTableOverlay::updateMenuCoordinates(vec2 buttonPos) {

    //set menu coordinates
    vec2 pos1, pos2;
    pos1 = buttonPos;

    int xLength = static_cast<int>(widgetsInMenu_.size() * (2 * widgetRadius_.get() + 10) + menuRadius_.get() + 10);
    xLength = std::max(xLength, (2 * widgetRadius_.get() + 10 + menuRadius_.get()));
    int yLength = 2 * widgetRadius_.get() + 10;

    pos2.x = (pos1.x - xLength <= 0) ? pos1.x + xLength : pos1.x - xLength;
    pos2.y = (pos1.y - yLength <= 0) ? pos1.y + yLength : pos1.y - yLength;

    overlayMenu_.setLL(tgt::min(pos1, pos2));
    overlayMenu_.setUR(tgt::max(pos1, pos2));

}

void TouchTableOverlay::onEvent(tgt::Event* e) {
    if(e->getEventType() == tgt::Event::TOUCHEVENT)
        handleTouchEvent((tgt::TouchEvent*) e);
}

void TouchTableOverlay::handleTouchEvent(tgt::TouchEvent* e) {

    const std::deque<TouchPoint>& tps = e->touchPoints();
    std::deque<tgt::TouchPoint> pointsToDistribute;
    std::deque<tgt::TouchPoint> pointsForMenuWidget;
    float radius = static_cast<float>(widgetRadius_.get());

    if (isInExclusiveMode_.get() && (currentMenuWidget_ != 0)){
        currentMenuWidget_->handleWidgetTouchPoints(tps);
        return;
    }

    //iterate over all TouchPoints
    std::deque<TouchPoint>::const_iterator i;
    for (i = tps.begin(); i != tps.end(); ++i) {
        tgt::TouchPoint::State state = i->state();

        //correct y corrdinate according to viewport size
        vec2 correctedPos = vec2(static_cast<float>(i->pos().x), static_cast<float>(outport_.getSize().y - i->pos().y));

        if (state == TouchPoint::TouchPointPressed) {

            //check if menu is open and if widgets in menu have been hit
            if(menuOpen_.get()){
                //is TouchPoint in Menu?
                if(isWidgetInMenu(correctedPos, 1)){
                    std::vector<TouchTableWidget*>::iterator widgetIter;
                    for(widgetIter = widgetsInMenu_.begin(); widgetIter != widgetsInMenu_.end(); ++widgetIter){

                        //get curren widget ,position and radius
                        TouchTableWidget* widget = *widgetIter;
                        vec2 pos = widget->getPosition();

                        if(tgt::length(correctedPos - pos) <= radius){
                            //if widget has been hit, make pair of corresponding touchpoint id and widget
                            touchPointMapping_.insert(std::make_pair(i->id(),widget));

                            //no other widget could have been hit ->break
                            break;
                        }
                    }//end for
                    continue;//done handling TouchPoint
                }
            }

            //check if TouchPoint hits menu of currentMenuWidget
            if((currentMenuWidget_ != 0) && currentMenuWidget_->isInMenu(correctedPos)){
                pointsForMenuWidget.push_back(*i);
                menuWidgetTouchPoints_.push_back(i->id());
                continue;
            }

            // run through all connected widget processors to check if one has been hit
            std::vector<TouchTableWidget*>::iterator iter;
            bool widgetHit = false;
            for (iter = widgetsOnTable_.begin(); iter != widgetsOnTable_.end(); ++iter) {

                // get current widget, position and radius
                TouchTableWidget* widget = *iter;
                vec2 pos = widget->getPosition();

                if (tgt::length(correctedPos - pos) <= radius) {
                    //if widget has been hit: insert corresponding "touch point id <-> widget" into map
                    touchPointMapping_.insert(std::make_pair(i->id(), widget));
                    //no other widget could have been hit -> break
                    widgetHit = true;
                    break;
                }

            }//end for
            if (!widgetHit)
                pointsToDistribute.push_back(*i);
        }
        else if (state == TouchPoint::TouchPointReleased) {

            //check if TouchPoint hits menu of currentMenuWidget
            std::vector<int>::iterator menuWidgetTouchPointIterator = std::find(menuWidgetTouchPoints_.begin(), menuWidgetTouchPoints_.end(), i->id());
            if((currentMenuWidget_) && (menuWidgetTouchPointIterator != menuWidgetTouchPoints_.end())) {
                pointsForMenuWidget.push_back(*i);
                menuWidgetTouchPoints_.erase(menuWidgetTouchPointIterator);
                continue;
            }

            bool menuHit = false;

            //check menu buttons for hit
            std::vector<ivec2>::iterator buttonIterator;
            for (buttonIterator = overlayMenu_.getButtons()->begin(); buttonIterator != overlayMenu_.getButtons()->end(); ++buttonIterator) {

                if (tgt::length(correctedPos - vec2(*buttonIterator)) <= menuRadius_.get()) {

                    menuHit = true;

                    if (menuOpen_.get() && (*buttonIterator == overlayMenu_.getButton())){

                        menuOpen_.set(false);

                    }else{
                        //update menu attributes
                        updateMenuCoordinates(*buttonIterator);
                        overlayMenu_.setButton(*buttonIterator);
                        menuOpen_.set(true);
                    }
                }
            }

                        // don't check widgets if menu button has been hit
            if (menuHit)
                continue;

            //find widgets in menu associated with the current touchpoint and release them
            std::map<int,TouchTableWidget*>::iterator iter;

            //find widgets associated with current touchpoint and release them
            iter = touchPointMapping_.find(i->id());

            //check if touch point is associated with menu
            if (menuOpen_.get() && (i->id() == overlayMenu_.getID())) {
                menuOpen_.set(false);
                touchPointMapping_.erase(iter);
                overlayMenu_.setID(-1);
                continue;
            }

            if (iter != touchPointMapping_.end()) {

                TouchTableWidget* widget = iter->second;

                if (!widget->isInMotion()) {
                    //widget has not been moved -> it has been selected
                    widget->pressed();
                }
                else {

                    //check if widget is in menu oder on table
                    if(menuOpen_.get()){

                        if(!isWidgetInMenu(widget->getPosition(), radius)){//release on Table
                            widget->setOnTable(true);

                        }else{ //release in Menu
                            widget->setOnTable(false);
                            //if widget was on Table before -> delete from widgetOnTable and insert in widgetInMenu
                            std::vector<TouchTableWidget*>::iterator widgetIter = std::find(widgetsOnTable_.begin(), widgetsOnTable_.end(), widget);
                            if(widgetIter != widgetsOnTable_.end()){
                                    widgetsOnTable_.erase(widgetIter);
                                    widgetsInMenu_.push_back(widget);
                                }

                        }


                        //refresh menu coordinates while keeping the same button position
                        updateMenuCoordinates(overlayMenu_.getButton());

                    }
                     //widget has been moved -> check if widget is in valid position, else set to last valid position
                    if  (widget->getPosition() != widget->getOldPosition())
                        widget->setPosition(widget->getOldPosition());
                }

                widget->setInMotion(false);
                widget->setPositionInvalid(false);
                touchPointMapping_.erase(iter);
            }
            else {
                pointsToDistribute.push_back(*i);
            }
        }
        else if(state == TouchPoint::TouchPointMoved || state == TouchPoint::TouchPointStationary) {

            //check if TouchPoint hits menu of currentMenuWidget
            std::vector<int>::iterator menuWidgetTouchPointIterator = std::find(menuWidgetTouchPoints_.begin(), menuWidgetTouchPoints_.end(), i->id());
            if((currentMenuWidget_) && (menuWidgetTouchPointIterator != menuWidgetTouchPoints_.end())) {
                pointsForMenuWidget.push_back(*i);
                continue;
            }

            std::map<int, TouchTableWidget*>::iterator iter;
            iter = touchPointMapping_.find(i->id());

            if (iter != touchPointMapping_.end()) {

                //get Position of the touchpoint and the widget
                vec2 tpPos= vec2(i->pos().x, outport_.getSize().y - i->pos().y);
                TouchTableWidget* widget = iter->second;

                // if no widget is associated: nothing to move
                if (!widget)
                    continue;

                //for widget not already moving: check if motion is above threshold, else return (do not move)
                if (!widget->isInMotion() && (tgt::length(tpPos - widget->getPosition()) < movementThreshold_.get()))
                   continue;


                //if widget in menu -> declare to be on Table
                if(!widget->isOnTable()){
                    widget->setOnTable(true);
                    std::vector<TouchTableWidget*>::iterator widgetIter = std::find(widgetsInMenu_.begin(), widgetsInMenu_.end(), widget);
                    if(widgetIter != widgetsInMenu_.end()){
                        widgetsInMenu_.erase(widgetIter);
                        widgetsOnTable_.push_back(widget);
                    }
                }

                widget->setInMotion(true);

                //create vector without currently moving widget, for collision handling
                std::vector<TouchTableWidget*> collisionWidgetVector = widgetsOnTable_;
                std::vector<TouchTableWidget*>::iterator widgetVecIter = std::find(collisionWidgetVector.begin(), collisionWidgetVector.end(), widget);
                if(widgetVecIter != collisionWidgetVector.end())
                    collisionWidgetVector.erase(widgetVecIter);

                //check for collisions with other widgets and the screen dimensions and check if the position is invalid
                handleWidgetCollision(tpPos, widget, collisionWidgetVector/*, &wallCollision*/);
                continue;

            }
            else{
                pointsToDistribute.push_back(*i);
            }
        }

    }//end for

    //let currentMenuWidget handle associated touchpoints
    if(currentMenuWidget_ != 0 && !pointsForMenuWidget.empty())
        currentMenuWidget_->handleWidgetTouchPoints(pointsForMenuWidget);

    //distribute touch points further in the network
    if (!pointsToDistribute.empty()) {

        if (pointsToDistribute.size() == 1) {
            //special case: only 1 touch point -> create mouse event for rotation / picking
            int posX = static_cast<int>(pointsToDistribute.at(0).pos().x);
            int posY = static_cast<int>(pointsToDistribute.at(0).pos().y);
            tgt::Event::Modifier modifier = tgt::Event::MODIFIER_NONE;
            if(moveView_.get())
                modifier = tgt::Event::SHIFT;
            tgt::MouseEvent::MouseAction action;
            if(pointsToDistribute.at(0).state() == tgt::TouchPoint::TouchPointMoved)
                action = tgt::MouseEvent::MOTION;
            else if(pointsToDistribute.at(0).state() == tgt::TouchPoint::TouchPointReleased)
                action = tgt::MouseEvent::RELEASED;
            else if(pointsToDistribute.at(0).state() == tgt::TouchPoint::TouchPointPressed)
                action = tgt::MouseEvent::PRESSED;

            tgt::MouseEvent* me = new tgt::MouseEvent(posX, posY, action, modifier, tgt::MouseEvent::MOUSE_BUTTON_LEFT, outport_.getSize());
            me->ignore();
            inport_.distributeEvent(me);
        }
        else if (pointsToDistribute.size() == 3) {
            //special case: 3 touch points - might be a 3-finger shift
            tgt::Bounds bounds(tgt::vec3(pointsToDistribute.at(0).pos(), 1.f), tgt::vec3(pointsToDistribute.at(1).pos(), 1.f));
            bounds.addPoint(tgt::vec3(pointsToDistribute.at(2).pos(), 1.f));

            if (tgt::length(bounds.diagonal()) < shiftDiameter_.get()) {
                // interpreted as 3-finger shift -> build mouse event

                //compute position and state
                int posX = static_cast<int>((pointsToDistribute.at(0).pos().x + pointsToDistribute.at(1).pos().x + pointsToDistribute.at(2).pos().x) / 3.f);
                int posY = static_cast<int>((pointsToDistribute.at(0).pos().y + pointsToDistribute.at(1).pos().y + pointsToDistribute.at(2).pos().y) / 3.f);

                tgt::Event::Modifier modifier = tgt::Event::SHIFT;
                tgt::MouseEvent::MouseAction action = tgt::MouseEvent::MOTION;

                if ((pointsToDistribute.at(0).state() == tgt::TouchPoint::TouchPointPressed)
                        || (pointsToDistribute.at(1).state() == tgt::TouchPoint::TouchPointPressed)
                        || (pointsToDistribute.at(2).state() == tgt::TouchPoint::TouchPointPressed))
                {
                    //create a new shift event
                    action = tgt::MouseEvent::PRESSED;
                    shiftTouchPoints_.clear();
                    shiftTouchPoints_.push_back(pointsToDistribute.at(0).id());
                    shiftTouchPoints_.push_back(pointsToDistribute.at(1).id());
                    shiftTouchPoints_.push_back(pointsToDistribute.at(2).id());
                }

                if ((pointsToDistribute.at(0).state() == tgt::TouchPoint::TouchPointReleased)
                        || (pointsToDistribute.at(1).state() == tgt::TouchPoint::TouchPointReleased)
                        || (pointsToDistribute.at(2).state() == tgt::TouchPoint::TouchPointReleased))
                {
                    //release the shift event
                    action = tgt::MouseEvent::RELEASED;
                    shiftTouchPoints_.clear();
                }

                if ((action == tgt::MouseEvent::MOTION) &&
                        ((std::find(shiftTouchPoints_.begin(), shiftTouchPoints_.end(), pointsToDistribute.at(0).id()) == shiftTouchPoints_.end()) ||
                         (std::find(shiftTouchPoints_.begin(), shiftTouchPoints_.end(), pointsToDistribute.at(1).id()) == shiftTouchPoints_.end()) ||
                         (std::find(shiftTouchPoints_.begin(), shiftTouchPoints_.end(), pointsToDistribute.at(2).id()) == shiftTouchPoints_.end())))
                {
                    //these touch points have not been registered as shift before -> distribute
                    tgt::TouchPoint::State states = pointsToDistribute.at(0).state();
                    tgt::TouchEvent* te = new tgt::TouchEvent(tgt::Event::MODIFIER_NONE, e->touchPointStates(),e->deviceType(), pointsToDistribute);
                    te->ignore();

                    // send touch event
                    inport_.distributeEvent(te);
                }else{

                    //distribute event
                    tgt::MouseEvent* me = new tgt::MouseEvent(posX, posY, action, modifier, tgt::MouseEvent::MOUSE_BUTTON_LEFT, outport_.getSize());
                    me->ignore();
                    inport_.distributeEvent(me);
                }
            }
            else {
                //not a 3-finger shift -> distribute events
                if ((pointsToDistribute.at(0).state() == tgt::TouchPoint::TouchPointReleased)
                        || (pointsToDistribute.at(1).state() == tgt::TouchPoint::TouchPointReleased)
                        || (pointsToDistribute.at(2).state() == tgt::TouchPoint::TouchPointReleased))
                    shiftTouchPoints_.clear();

                tgt::TouchPoint::State states = pointsToDistribute.at(0).state();
                tgt::TouchEvent* te = new tgt::TouchEvent(tgt::Event::MODIFIER_NONE, e->touchPointStates(),e->deviceType(), pointsToDistribute);
                te->ignore();

                // send touch event
                inport_.distributeEvent(te);
            }
        }
        else {

            tgt::TouchPoint::State states = pointsToDistribute.at(0).state();
            tgt::TouchEvent* te = new tgt::TouchEvent(tgt::Event::MODIFIER_NONE, e->touchPointStates(),e->deviceType(), pointsToDistribute);
            te->ignore();

            // send touch event
            inport_.distributeEvent(te);
        }
    }
}

void TouchTableOverlay::handleWidgetCollision(vec2 destPos, TouchTableWidget* widget, std::vector<TouchTableWidget*> collisionWidgetVector) {
    float radius = static_cast<float>(widgetRadius_.get());

    //clamp to screen position
    destPos = vec2(tgt::clamp( destPos.x, static_cast<float>(radius+1), static_cast<float>(outport_.getSize().x - radius - 1)),
                   tgt::clamp( destPos.y, static_cast<float>(radius+1), static_cast<float>(outport_.getSize().y - radius - 1)));


    std::vector<ivec2>::iterator iter;
    //iterate over menu buttons to check for collisions
    for(iter = overlayMenu_.getButtons()->begin(); iter != overlayMenu_.getButtons()->end(); ++iter){
        vec2 distVec = vec2(destPos.x - iter->x, destPos.y - iter->y);
        float absDist = tgt::length(distVec); //absolute distance between the two centers of the widgets

        if(absDist - (radius + menuRadius_.get()) < 0){
           destPos = vec2(iter->x + distVec.x / absDist * (menuRadius_.get()+radius + 1),
                          iter->y + distVec.y / absDist * (menuRadius_.get()+radius + 1));
        }
    }


    //iterate over widgets to check for collisions
    bool collision = false;
    std::vector<TouchTableWidget*>::iterator collisionWidgetIter;
    for(collisionWidgetIter=collisionWidgetVector.begin(); collisionWidgetIter!=collisionWidgetVector.end(); ++collisionWidgetIter){

        //calculate vector from center of the position destination of the widget to the center of the collision widget
        vec2 widgetDistVec = vec2((*collisionWidgetIter)->getPosition().x - destPos.x, (*collisionWidgetIter)->getPosition().y - destPos.y);

        float distBetweenWidgets = tgt::length(widgetDistVec); //absolute distance between the two centers of the widgets

        if(distBetweenWidgets - 2 * radius < 0){
            widget->setPositionInvalid(true);
            widget->setPosition(destPos);
            collision = true;
        }

    }//end for

    if (!collision) {
        //no collisions -> set position and old position
        widget->setPosition(destPos);
        widget->setOldPosition(destPos);
        widget->setPositionInvalid(false);
    }
}

bool TouchTableOverlay::isWidgetInMenu(vec2 pos, float radius) const {
    if((pos.x + radius) < overlayMenu_.getLL().x || (pos.y + radius) < overlayMenu_.getLL().y || (pos.x - radius) > overlayMenu_.getUR().x || (pos.y - radius) > overlayMenu_.getUR().y  ){//on table

        return false;

    }else{//in Menu
        return true;
    }

}

void TouchTableOverlay::setMenuWidget(TouchTableMenuWidget* widget){
    if(currentMenuWidget_)
        currentMenuWidget_->setActive(false);

    currentMenuWidget_ = widget;
    isInExclusiveMode_.set(false);

    menuWidgetTouchPoints_.clear();
}

void TouchTableOverlay::setExclusiveMode(bool mode){
    isInExclusiveMode_.set(mode);
}

tgt::ivec2 TouchTableOverlay::getOutportSize() const {
    return outport_.getSize();
}

int TouchTableOverlay::getWidgetRadius() const {
    return widgetRadius_.get();
}

bool TouchTableOverlay::getCameraMovement() const {
    return moveView_.get();
}

void TouchTableOverlay::setCameraMovement(bool move){

    moveView_.set(move);

}


} // namespace voreen
