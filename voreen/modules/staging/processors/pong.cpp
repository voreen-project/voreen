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

#include "pong.h"
#include "tgt/immediatemode/immediatemode.h"

#include "tgt/textureunit.h"
#include "tgt/event/timeevent.h"
#include "tgt/quadric.h"
#include "tgt/timer.h"

#include <deque>
#include <iterator>
#include <algorithm>
#include <time.h>

#include "voreen/core/voreenapplication.h"

using tgt::TextureUnit;
using tgt::vec2;
using tgt::TouchPoint;

namespace voreen {

Pong::Pong()
    : ImageProcessor("image/compositor")
    , inport_(Port::INPORT, "image.in", "Image Input", false, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_ORIGIN)
    , outport_(Port::OUTPORT, "image.out", "Image Output", true, Processor::INVALID_RESULT, RenderPort::RENDERSIZE_RECEIVER)
    , privatePort_(Port::OUTPORT, "image.tmp", "image.tmp", false)
    , color_("color", "Color", tgt::Color(1.f, 1.f, 1.f, 1.f))
    , opacity_("opacity", "Opacity", 1.f, 0.f, 1.f)
    , paddleHeight_("paddleHeight", "Paddle height", 0.2f, 0.05f, 1.f)
    , paddleWidth_("paddleWidth", "Paddle width", 0.04f, 0.01f, 0.1f)
    , ballSize_("ballSize", "Ball size", 0.05f, 0.01f, 0.1f)
    , speedFactor_("speedFac", "Ball Speed", 2.f, 1.f, 5.f)
    , paddleAutoSpeed_("paddleAutoSpeed", "Paddle Speed Factor", 5.f, 1.f, 10.f)
    , randomFactor_("randomFac", "Randomness", 2, 1, 5)
    , resetGame_("resetGame", "Reset Game")
    , crazyMode_("crazyMode", "Crazy Mode", false)
    , computerOpponent1_("computerOpponent1", "Use computer opponent left", false)
    , computerOpponent2_("computerOpponent2", "Use computer opponent right", false)
    , xmasMode_("xmasMode", "Happy X-mas!!", true)
    , copyShader_(0)
    , ballPos_(tgt::vec2(0.5f, 0.5f))
    , ballVel_(tgt::vec2(0.f))
    , paddle1Pos_(tgt::vec2(0.1f, 0.5f))
    , paddle2Pos_(tgt::vec2(0.9f, 0.5f))
    , paddle1Move_(false)
    , paddle2Move_(false)
    , upPressedL_(false)
    , downPressedL_(false)
    , leftPressedL_(false)
    , rightPressedL_(false)
    , upPressedR_(false)
    , downPressedR_(false)
    , leftPressedR_(false)
    , rightPressedR_(false)
    , leftTarget_(tgt::vec2(-1.f))
    , rightTarget_(tgt::vec2(-1.f))
    , angle_(0.f)
    , angleVel_(0.01f)
    , sleepTime_(50)
    , sleepCounter_(-1)
    , timer_(0)
    , renderScore_(false)
    , xmasBall_(0)
    , xmasFlake_(0)
{
    eventHandler_.addListenerToBack(this);

    addPort(inport_);
    addPort(outport_);
    addPrivateRenderPort(privatePort_);

    addProperty(color_);
    addProperty(opacity_);

    addProperty(resetGame_);
    addProperty(paddleHeight_);
    addProperty(paddleWidth_);
    addProperty(ballSize_);
    addProperty(speedFactor_);
    addProperty(paddleAutoSpeed_);
    addProperty(randomFactor_);
    addProperty(crazyMode_);
    addProperty(computerOpponent1_);
    addProperty(computerOpponent2_);
    addProperty(xmasMode_);

    resetGame_.onChange(MemberFunctionCallback<Pong>(this, &Pong::onReset));
    computerOpponent1_.onChange(MemberFunctionCallback<Pong>(this, &Pong::leftAIOnChange));
    computerOpponent2_.onChange(MemberFunctionCallback<Pong>(this, &Pong::rightAIOnChange));
    std::srand((unsigned)time(NULL));

    font_ = new tgt::Font(VoreenApplication::app()->getFontPath("Vera.ttf"), 36);

    for(size_t i = 0; i < 100; i++)
        snow_.push_back(std::pair<tgt::vec2, tgt::vec2>(tgt::vec2(std::abs(getRandomFloat()), std::abs(getRandomFloat())), tgt::vec2(getRandomFloat(), std::abs(getRandomFloat()) * speedFactor_.get() * 0.015f)));
}

Pong::~Pong() {
    delete font_;
}

bool Pong::isReady() const {
    return outport_.isReady();
}

void Pong::leftAIOnChange() {
    if(!computerOpponent1_.get())
        leftTarget_ = tgt::vec2(-1.f);
}

void Pong::rightAIOnChange() {
    if(!computerOpponent2_.get())
        rightTarget_ = tgt::vec2(-1.f);
}

void Pong::onReset() {
    score1_ = 0;
    score2_ = 0;
    ballPos_ = tgt::vec2(0.75f - paddleWidth_.get() * 0.5f - ballSize_.get() * 0.5f, 0.5f);
    ballVel_ = normalize(tgt::vec2(-2.f * std::abs(getRandomFloat()), getRandomFloat())) * speedFactor_.get() * 0.0015f;
    if(std::abs(ballVel_.x) < speedFactor_.get() * 0.001f)
        ballVel_.x = -speedFactor_.get() * 0.001f;

    angleVel_ = getRandomFloat() * 0.1f;

    sleepCounter_ = 0;
    renderScore_ = true;

    paddle1Pos_ = tgt::vec2(0.1f, 0.5f);
    paddle2Pos_ = tgt::vec2(0.9f, 0.5f);
    paddle1Move_ = false;
    paddle2Move_ = false;

    invalidate();
}

void Pong::initialize() {
    ImageProcessor::initialize();

    try {
        copyShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage.frag",
            ImageProcessor::generateHeader(), false);
        copyShader_->deactivate();
    } catch(tgt::Exception&) {
        processorState_ = PROCESSOR_STATE_NOT_INITIALIZED;
        throw VoreenException("Failed to load shader: passthrough.vert/copyimage.frag");
    }

    onReset();

    timer_ = VoreenApplication::app()->createTimer(&eventHandler_);
    if (timer_)
        timer_->start(30);

    xmasBall_ = TexMgr.load(VoreenApplication::app()->getModulePath("staging") + "/textures/xmas.png");
    xmasFlake_ = TexMgr.load(VoreenApplication::app()->getModulePath("staging") + "/textures/snowflake.png");
    if(!xmasBall_ || !xmasFlake_)
        throw VoreenException("Failed to load xmas textures");
}

void Pong::deinitialize() {
    TexMgr.dispose(xmasBall_);
    TexMgr.dispose(xmasFlake_);
    ShdrMgr.dispose(copyShader_);
    delete timer_;
    ImageProcessor::deinitialize();
}

void Pong::process() {

    if(inport_.getSize() != outport_.getSize())
        inport_.requestSize(outport_.getSize());
    if(privatePort_.getSize() != outport_.getSize())
        privatePort_.resize(outport_.getSize());

    privatePort_.activateTarget();
    privatePort_.clearTarget();

    MatStack.pushMatrix();
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.f, (float)outport_.getSize().x, (float)outport_.getSize().y, 0.f, 0.f, 1.f));

    glDepthFunc(GL_ALWAYS);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    tgt::vec2 paddle1PosScreen = paddle1Pos_ * tgt::vec2(outport_.getSize());
    tgt::vec2 paddle2PosScreen = paddle2Pos_ * tgt::vec2(outport_.getSize());
    float paddleWidthScreen = paddleWidth_.get() * outport_.getSize().x;
    float paddleHeightScreen = paddleHeight_.get() * outport_.getSize().y;

    IMode.color(color_.get().x, color_.get().y, color_.get().z, opacity_.get());
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(paddle1PosScreen.x - paddleWidthScreen * 0.5f, paddle1PosScreen.y - paddleHeightScreen * 0.5f);
        IMode.vertex(paddle1PosScreen.x + paddleWidthScreen * 0.5f, paddle1PosScreen.y - paddleHeightScreen * 0.5f);
        IMode.vertex(paddle1PosScreen.x + paddleWidthScreen * 0.5f, paddle1PosScreen.y + paddleHeightScreen * 0.5f);
        IMode.vertex(paddle1PosScreen.x - paddleWidthScreen * 0.5f, paddle1PosScreen.y + paddleHeightScreen * 0.5f);
    IMode.end();
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(paddle2PosScreen.x - paddleWidthScreen * 0.5f, paddle2PosScreen.y - paddleHeightScreen * 0.5f);
        IMode.vertex(paddle2PosScreen.x + paddleWidthScreen * 0.5f, paddle2PosScreen.y - paddleHeightScreen * 0.5f);
        IMode.vertex(paddle2PosScreen.x + paddleWidthScreen * 0.5f, paddle2PosScreen.y + paddleHeightScreen * 0.5f);
        IMode.vertex(paddle2PosScreen.x - paddleWidthScreen * 0.5f, paddle2PosScreen.y + paddleHeightScreen * 0.5f);
    IMode.end();

    if(xmasMode_.get()) {
        static std::vector<tgt::vec2> offsets1(10, tgt::vec2(paddle1PosScreen.x, paddle1PosScreen.y - paddleHeightScreen * 0.25f));
        static std::vector<tgt::vec2> offsets2(10, tgt::vec2(paddle2PosScreen.x, paddle2PosScreen.y - paddleHeightScreen * 0.25f));
        static int curOffsetIndex = 0;

        static char firstTen = 0;

        tgt::vec2 prevOffset1;
        tgt::vec2 prevOffset2;
        if(firstTen >= 10) {
            prevOffset1 = offsets1[curOffsetIndex];
            prevOffset2 = offsets2[curOffsetIndex];
        } else {
            firstTen++;
            prevOffset1 = tgt::vec2(paddle1PosScreen.x, paddle1PosScreen.y - paddleHeightScreen * 0.25f);// - paddleWidthScreen * 2.f);
            prevOffset2 = tgt::vec2(paddle2PosScreen.x, paddle2PosScreen.y - paddleHeightScreen * 0.25f);// - paddleWidthScreen * 2.f);
        }

        offsets1[curOffsetIndex] = tgt::vec2(paddle1PosScreen.x, paddle1PosScreen.y - paddleHeightScreen * 0.25f);// - paddleWidthScreen * 2.f);
        offsets2[curOffsetIndex] = tgt::vec2(paddle2PosScreen.x, paddle2PosScreen.y - paddleHeightScreen * 0.25f);// - paddleWidthScreen * 2.f);

        tgt::vec2 tip1Dist = prevOffset1 - tgt::vec2(paddle1PosScreen.x - paddleWidthScreen * 0.75f, paddle1PosScreen.y - paddleHeightScreen * 0.5f);
        if(length(tip1Dist) > 0.75f * paddleHeightScreen)
            prevOffset1 = tgt::vec2(paddle1PosScreen.x - paddleWidthScreen * 0.75f, paddle1PosScreen.y - paddleHeightScreen * 0.5f) + 0.75f * paddleHeightScreen * normalize(tip1Dist);
        tgt::vec2 tip2Dist = prevOffset2 - tgt::vec2(paddle2PosScreen.x - paddleWidthScreen * 0.75f, paddle2PosScreen.y - paddleHeightScreen * 0.5f);
        if(length(tip2Dist) > 0.75f * paddleHeightScreen)
            prevOffset2 = tgt::vec2(paddle2PosScreen.x - paddleWidthScreen * 0.75f, paddle2PosScreen.y - paddleHeightScreen * 0.5f) + 0.75f * paddleHeightScreen * normalize(tip2Dist);

        curOffsetIndex++;
        curOffsetIndex %= 10;

        //MatStack.pushMatrix();
        //MatStack.translate(ballPosScreen.x, ballPosScreen.y, 0.f);
        //MatStack.multMatrix(tgt::mat4::createRotationZ(10.f));

        // left hat
        IMode.color(1.f, 0.f, 0.f, opacity_.get());
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(paddle1PosScreen.x - paddleWidthScreen * 0.75f, paddle1PosScreen.y - paddleHeightScreen * 0.5f);
            IMode.vertex(paddle1PosScreen.x + paddleWidthScreen * 0.75f, paddle1PosScreen.y - paddleHeightScreen * 0.5f);
            IMode.vertex(paddle1PosScreen.x + paddleWidthScreen * 0.6f, paddle1PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen);
            IMode.vertex(paddle1PosScreen.x - paddleWidthScreen * 0.6f, paddle1PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen);
        IMode.end();
        IMode.color(0.8f, 0.8f, 0.8f, opacity_.get());
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(paddle1PosScreen.x - paddleWidthScreen * 0.8f, paddle1PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen * 0.5f);
            IMode.vertex(paddle1PosScreen.x + paddleWidthScreen * 0.8f, paddle1PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen * 0.5f);
            IMode.vertex(paddle1PosScreen.x + paddleWidthScreen * 0.8f, paddle1PosScreen.y - paddleHeightScreen * 0.5f + paddleWidthScreen * 0.5f);
            IMode.vertex(paddle1PosScreen.x - paddleWidthScreen * 0.8f, paddle1PosScreen.y - paddleHeightScreen * 0.5f + paddleWidthScreen * 0.5f);
        IMode.end();
        IMode.begin(tgt::ImmediateMode::TRIANGLES);
            IMode.color(1.f, 0.f, 0.f, opacity_.get());
            IMode.vertex(paddle1PosScreen.x - 0.6f* paddleWidthScreen, paddle1PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen);
            IMode.vertex(paddle1PosScreen.x + 0.6f* paddleWidthScreen, paddle1PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen);
            IMode.color(0.6f, 0.f, 0.f, opacity_.get());
            IMode.vertex(prevOffset1.x, prevOffset1.y);
        IMode.end();
        IMode.color(0.8f, 0.8f, 0.8f, opacity_.get());
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(prevOffset1.x - paddleWidthScreen * 0.1f, prevOffset1.y + 0.1f * paddleWidthScreen);
            IMode.vertex(prevOffset1.x + paddleWidthScreen * 0.1f, prevOffset1.y + 0.1f * paddleWidthScreen);
            IMode.vertex(prevOffset1.x + paddleWidthScreen * 0.1f, prevOffset1.y - 0.1f * paddleWidthScreen);
            IMode.vertex(prevOffset1.x - paddleWidthScreen * 0.1f, prevOffset1.y - 0.1f * paddleWidthScreen);
        IMode.end();

        // right hat
        IMode.color(1.f, 0.f, 0.f, opacity_.get());
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(paddle2PosScreen.x - paddleWidthScreen * 0.75f, paddle2PosScreen.y - paddleHeightScreen * 0.5f);
            IMode.vertex(paddle2PosScreen.x + paddleWidthScreen * 0.75f, paddle2PosScreen.y - paddleHeightScreen * 0.5f);
            IMode.vertex(paddle2PosScreen.x + paddleWidthScreen * 0.6f, paddle2PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen);
            IMode.vertex(paddle2PosScreen.x - paddleWidthScreen * 0.6f, paddle2PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen);
        IMode.end();
        IMode.color(0.8f, 0.8f, 0.8f, opacity_.get());
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(paddle2PosScreen.x - paddleWidthScreen * 0.8f, paddle2PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen * 0.5f);
            IMode.vertex(paddle2PosScreen.x + paddleWidthScreen * 0.8f, paddle2PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen * 0.5f);
            IMode.vertex(paddle2PosScreen.x + paddleWidthScreen * 0.8f, paddle2PosScreen.y - paddleHeightScreen * 0.5f + paddleWidthScreen * 0.5f);
            IMode.vertex(paddle2PosScreen.x - paddleWidthScreen * 0.8f, paddle2PosScreen.y - paddleHeightScreen * 0.5f + paddleWidthScreen * 0.5f);
        IMode.end();
        IMode.begin(tgt::ImmediateMode::TRIANGLES);
            IMode.color(1.f, 0.f, 0.f, opacity_.get());
            IMode.vertex(paddle2PosScreen.x - 0.6f* paddleWidthScreen, paddle2PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen);
            IMode.vertex(paddle2PosScreen.x + 0.6f* paddleWidthScreen, paddle2PosScreen.y - paddleHeightScreen * 0.5f - paddleWidthScreen);
            IMode.color(0.6f, 0.f, 0.f, opacity_.get());
            IMode.vertex(prevOffset2.x, prevOffset2.y);
        IMode.end();
        IMode.color(0.8f, 0.8f, 0.8f, opacity_.get());
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(prevOffset2.x - paddleWidthScreen * 0.1f, prevOffset2.y + 0.1f * paddleWidthScreen);
            IMode.vertex(prevOffset2.x + paddleWidthScreen * 0.1f, prevOffset2.y + 0.1f * paddleWidthScreen);
            IMode.vertex(prevOffset2.x + paddleWidthScreen * 0.1f, prevOffset2.y - 0.1f * paddleWidthScreen);
            IMode.vertex(prevOffset2.x - paddleWidthScreen * 0.1f, prevOffset2.y - 0.1f * paddleWidthScreen);
        IMode.end();
        //MatStack.popMatrix();

        IMode.color(1.f, 1.f, 1.f, opacity_.get());
        if(xmasFlake_) {
            xmasFlake_->enable();
            xmasFlake_->bind();
        }

        for(size_t i  = 0; i < snow_.size(); i++) {
            MatStack.pushMatrix();
            MatStack.translate(snow_[i].first.x * outport_.getSize().x, snow_[i].first.y * outport_.getSize().y, 0.f);

            IMode.begin(tgt::ImmediateMode::QUADS);
                IMode.texcoord(0.f, 0.f);
                IMode.vertex(-paddleWidthScreen * 0.2f,  0.2f * paddleWidthScreen);
                IMode.texcoord(1.f, 0.f);
                IMode.vertex( paddleWidthScreen * 0.2f,  0.2f * paddleWidthScreen);
                IMode.texcoord(1.f, 1.f);
                IMode.vertex( paddleWidthScreen * 0.2f, -0.2f * paddleWidthScreen);
                IMode.texcoord(0.f, 1.f);
                IMode.vertex(-paddleWidthScreen * 0.2f, -0.2f * paddleWidthScreen);
            IMode.end();
            MatStack.popMatrix();
        }
        if(xmasFlake_)
            xmasFlake_->disable();
    }

    IMode.color(color_.get().x, color_.get().y, color_.get().z, opacity_.get());

    tgt::vec2 ballPosScreen = ballPos_ * tgt::vec2(outport_.getSize());
    float ballRadiusScreen = ballSize_.get() * outport_.getSize().x * 0.5f;
    //tgt::Disk circle = tgt::Disk(0.f, 0.5f * ballSize_.get() * outport_.getSize().x, 4, 4);
    MatStack.pushMatrix();
    MatStack.translate(ballPosScreen.x, ballPosScreen.y, 0.f);
    MatStack.multMatrix(tgt::mat4::createRotationZ(angle_));

    if(xmasMode_.get() && xmasBall_) {
        xmasBall_->enable();
        xmasBall_->bind();
    }

    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.texcoord(0.f, 0.f);
        IMode.vertex(-ballRadiusScreen, -ballRadiusScreen);
        IMode.texcoord(1.f, 0.f);
        IMode.vertex(ballRadiusScreen, -ballRadiusScreen);
        IMode.texcoord(1.f, 1.f);
        IMode.vertex(ballRadiusScreen,  ballRadiusScreen);
        IMode.texcoord(0.f, 1.f);
        IMode.vertex(-ballRadiusScreen,  ballRadiusScreen);
    IMode.end();
    //circle.render();
    MatStack.popMatrix();
    if(xmasMode_.get() && xmasBall_)
        xmasBall_->disable();

    MatStack.popMatrix();

    if(renderScore_) {
        std::stringstream str1, str2;
        str1 << score1_;
        str2 << score2_;

        font_->setFontSize(48);
        font_->setPadding(7);
        tgt::ivec2 screensize = outport_.getSize();
        font_->setTextAlignment(tgt::Font::BottomLeft);
        font_->render(tgt::vec3(0.f, 0.8f *(float)outport_.getSize().y, 0), str1.str(), screensize);
        font_->setTextAlignment(tgt::Font::BottomRight);
        font_->render(tgt::vec3(outport_.getSize().x, 0.8f *(float)outport_.getSize().y, 0), str2.str(), screensize);

    }

    privatePort_.deactivateTarget();

    outport_.activateTarget();
    outport_.clearTarget();

    if(inport_.isReady()) {
        // blend pong image over incoming image
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
        program_->setUniform("depthTex0_", 0);
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

    // Restore GL state
    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDepthFunc(GL_LESS);
    IMode.color(1.0f, 1.0f, 1.0f, 1.0f);

    LGL_ERROR;
}

void Pong::touchEvent(tgt::TouchEvent* e) {

    e->accept();

    // average positions of left and right player's fingers
    tgt::vec2 leftPos = tgt::vec2(0.f);
    tgt::vec2 rightPos = tgt::vec2(0.f);

    unsigned char leftCounter = 0;
    unsigned char rightCounter = 0;
    bool leftReleased = true;
    bool rightReleased = true;

    for(std::deque<tgt::TouchPoint>::const_iterator it = e->touchPoints().begin(); it != e->touchPoints().end(); it++) {
        tgt::vec2 pos = it->pos() / tgt::vec2(e->getViewport());
        if(pos.x <= 0.5f) {
            leftCounter++;
            leftPos += pos;
            if(it->state() == tgt::TouchPoint::TouchPointPressed &&
               pos.x >= paddle1Pos_.x - 0.5f * paddleWidth_.get() && pos.x <= paddle1Pos_.x + 0.5f * paddleWidth_.get() &&
               pos.y >= paddle1Pos_.y - 0.5f * paddleHeight_.get() && pos.y <= paddle1Pos_.y + 0.5f * paddleHeight_.get()) {
                paddle1Move_ = true;
            }
            leftReleased &= (it->state() == tgt::TouchPoint::TouchPointReleased);
        }
        else {
            rightCounter++;
            rightPos += pos;
            if(it->state() == tgt::TouchPoint::TouchPointPressed &&
               pos.x >= paddle2Pos_.x - 0.5f * paddleWidth_.get() && pos.x <= paddle2Pos_.x + 0.5f * paddleWidth_.get() &&
               pos.y >= paddle2Pos_.y - 0.5f * paddleHeight_.get() && pos.y <= paddle2Pos_.y + 0.5f * paddleHeight_.get()) {
                paddle2Move_ = true;
            }
            rightReleased &= (it->state() == tgt::TouchPoint::TouchPointReleased);
        }
    }

    if(leftReleased)
        paddle1Move_ = false;
    if(rightReleased)
        paddle2Move_ = false;

    if(leftCounter > 0)
        leftPos /= (float)leftCounter;
    if(rightCounter > 0)
        rightPos /= (float)rightCounter;

    bool avoidDirectLeftMove = false;
    bool avoidDirectRightMove = false;
    if(!leftReleased && !paddle1Move_ && leftCounter > 0)
        leftTarget_ = leftPos;
    else
        leftTarget_ = tgt::vec2(-1.f);

    if(!rightReleased && !paddle2Move_ && rightCounter > 0)
        rightTarget_ = rightPos;
    else
        rightTarget_ = tgt::vec2(-1.f);

    if(leftTarget_.x >= 0.f)
        avoidDirectLeftMove = true;
    if(rightTarget_.x >= 0.f)
        avoidDirectRightMove = true;

    if(e->touchPointStates() & tgt::TouchPoint::TouchPointMoved) {
        if(!avoidDirectLeftMove && paddle1Move_ && leftCounter > 0) {
            if(leftPos.x < 0.25f && leftPos.x > paddleWidth_.get())
                paddle1Pos_.x = leftPos.x;
            if(leftPos.y + paddleHeight_.get() * 0.5f <= 1.f && leftPos.y - paddleHeight_.get() * 0.5f >= 0.f)
                paddle1Pos_.y = leftPos.y;
        }
        if(!avoidDirectRightMove && paddle2Move_ && rightCounter > 0) {
            if(rightPos.x > 0.75f && rightPos.x < 1.f - paddleWidth_.get())
                paddle2Pos_.x = rightPos.x;
            if(rightPos.y + paddleHeight_.get() * 0.5f <= 1.f && rightPos.y - paddleHeight_.get() * 0.5f >= 0.f)
                paddle2Pos_.y = rightPos.y;
        }
    }

    invalidate();
}

void Pong::keyEvent(tgt::KeyEvent* e) {

    if(e->keyCode() & tgt::KeyEvent::K_W ||
       e->keyCode() & tgt::KeyEvent::K_S ||
       e->keyCode() & tgt::KeyEvent::K_A ||
       e->keyCode() & tgt::KeyEvent::K_D)
        leftTarget_ = tgt::vec2(-1.f);

    if(e->keyCode() == tgt::KeyEvent::K_W)
        upPressedL_ = e->pressed();
    else if(e->keyCode() == tgt::KeyEvent::K_S)
        downPressedL_ = e->pressed();
    else if(e->keyCode() == tgt::KeyEvent::K_A)
        leftPressedL_ = e->pressed();
    else if(e->keyCode() == tgt::KeyEvent::K_D)
        rightPressedL_ = e->pressed();

    if(e->keyCode() & tgt::KeyEvent::K_UP ||
       e->keyCode() & tgt::KeyEvent::K_DOWN ||
       e->keyCode() & tgt::KeyEvent::K_LEFT||
       e->keyCode() & tgt::KeyEvent::K_RIGHT)
        rightTarget_ = tgt::vec2(-1.f);

    if(e->keyCode() == tgt::KeyEvent::K_UP)
        upPressedR_ = e->pressed();
    else if(e->keyCode() == tgt::KeyEvent::K_DOWN)
        downPressedR_ = e->pressed();
    else if(e->keyCode() == tgt::KeyEvent::K_LEFT)
        leftPressedR_ = e->pressed();
    else if(e->keyCode() == tgt::KeyEvent::K_RIGHT)
        rightPressedR_ = e->pressed();
}

void Pong::mouseDoubleClickEvent(tgt::MouseEvent* e) {
    onReset();
    e->accept();
}

void Pong::timerEvent(tgt::TimeEvent* e) {
    e->accept();

    for(size_t i  = 0; i < snow_.size(); i++) {
        snow_[i].first += snow_[i].second;
        snow_[i].second.x += getRandomFloat() * 0.001f * randomFactor_.get() * speedFactor_.get();
        snow_[i].second.y += getRandomFloat() * 0.001f * randomFactor_.get() * speedFactor_.get();
        snow_[i].second.x = tgt::clamp(snow_[i].second.x, -0.0015f * speedFactor_.get(), 0.0015f * speedFactor_.get());
        snow_[i].second.y = tgt::clamp(snow_[i].second.y, 0.0005f * speedFactor_.get(), 0.0015f * speedFactor_.get());

        if(snow_[i].first.y > 1.f || snow_[i].first.x < -0.3f || snow_[i].first.x > 1.3f)
            snow_[i] = std::pair<tgt::vec2, tgt::vec2>(tgt::vec2(std::abs(getRandomFloat()), 0.f), tgt::vec2(getRandomFloat(), std::abs(getRandomFloat())) * speedFactor_.get() * 0.015f);
    }

    if(computerOpponent1_.get()) {
        if(ballVel_.x < 0.f)
            leftTarget_.y = ballPos_.y;
        else
            leftTarget_.y = 0.5f;

        if(length(ballVel_) > speedFactor_.get() * 0.015f)
            leftTarget_.x = 0.05f;
        else
            leftTarget_.x = 0.2f;
    }

    if(computerOpponent2_.get()) {
        if(ballVel_.x > 0.f)
            rightTarget_.y = ballPos_.y;
        else
            rightTarget_.y = 0.5f;

        if(length(ballVel_) > speedFactor_.get() * 0.015f)
            rightTarget_.x = 0.95f;
        else
            rightTarget_.x = 0.8f;
    }

    //tgt::vec2 oldPaddle1Pos = paddle1Pos_;
    //tgt::vec2 oldPaddle2Pos = paddle2Pos_;

    // key press related paddle moving
    float keyInc = speedFactor_.get() * 0.0015f * paddleAutoSpeed_.get();
    if(upPressedL_)
        paddle1Pos_.y -= keyInc;
    if(downPressedL_)
        paddle1Pos_.y += keyInc;
    if(leftPressedL_)
        paddle1Pos_.x -= keyInc;
    if(rightPressedL_)
        paddle1Pos_.x += keyInc;
    if(upPressedR_)
        paddle2Pos_.y -= keyInc;
    if(downPressedR_)
        paddle2Pos_.y += keyInc;
    if(leftPressedR_)
        paddle2Pos_.x -= keyInc;
    if(rightPressedR_)
        paddle2Pos_.x += keyInc;

    // auto-move paddles if necessary
    if(leftTarget_.x >= 0.f) {
        tgt::vec2 dist = leftTarget_ - paddle1Pos_;
        if(length(dist) > speedFactor_.get() * 0.001f * paddleAutoSpeed_.get())
            paddle1Pos_ += normalize(dist) * std::pow(1.f + length(dist), 2.f) * speedFactor_.get() * 0.001f * paddleAutoSpeed_.get();
    }
    if(rightTarget_.x >= 0.f) {
        tgt::vec2 dist = rightTarget_ - paddle2Pos_;
        if(length(dist) > speedFactor_.get() * 0.001f * paddleAutoSpeed_.get())
            paddle2Pos_ += normalize(dist) * std::pow(1.f + length(dist), 2.f) * speedFactor_.get() * 0.001f * paddleAutoSpeed_.get();
    }

    tgt::vec2 bordersL = tgt::vec2(paddleWidth_.get() * 0.5f, paddleHeight_.get() * 0.5f);
    paddle1Pos_ = tgt::clamp(paddle1Pos_, bordersL, tgt::vec2(0.25f, 1.f - bordersL.y));
    tgt::vec2 bordersR = tgt::vec2(1.f - paddleWidth_.get() * 0.5f, 1.f - paddleHeight_.get() * 0.5f);
    paddle2Pos_ = tgt::clamp(paddle2Pos_, tgt::vec2(0.75f, 1.f - bordersR.y), bordersR);

    if(sleepCounter_ >= 0) {
        if(sleepCounter_ > sleepTime_) {
            sleepCounter_ = -1;
            renderScore_ = false;
        } else {
            sleepCounter_++;
            invalidate();
            return;
        }
    }

    // was a point scored?
    if(ballPos_.x - ballSize_.get() * 0.5f < 0.f) {
        score2_++;
        ballPos_ = tgt::vec2(0.75f - paddleWidth_.get() * 0.5f - ballSize_.get() * 0.5f, 0.5f);
        ballVel_ = normalize(tgt::vec2(-2.f * std::abs(getRandomFloat()), getRandomFloat())) * speedFactor_.get() * 0.0015f;
        if(std::abs(ballVel_.x) < speedFactor_.get() * 0.001f)
            ballVel_.x = -speedFactor_.get() * 0.001f;
        angleVel_ = getRandomFloat() * 0.1f;
        sleepCounter_ = 0;
        renderScore_ = true;
        invalidate();
        return;
    }

    if(ballPos_.x + ballSize_.get() * 0.5f > 1.f) {
        score1_++;
        ballPos_ = tgt::vec2(0.25f + paddleWidth_.get() * 0.5f + ballSize_.get() * 0.5f, 0.5f);
        ballVel_ = normalize(tgt::vec2(2.f * std::abs(getRandomFloat()), getRandomFloat())) * speedFactor_.get() * 0.0015f;
        if(ballVel_.x < speedFactor_.get() * 0.001f)
            ballVel_.x = speedFactor_.get() * 0.001f;
        angleVel_ = getRandomFloat() * 0.1f;
        sleepCounter_ = 0;
        renderScore_ = true;
        invalidate();
        return;
    }

    // check for ceiling collisions
    if(ballPos_.y + ballSize_.get() * 0.5f > 1.f) {
        ballPos_.y = 1.f - ballSize_.get() * 0.5f;
        ballVel_.y *= -1.f;
        variate(ballVel_);
        angleVel_ += ballVel_.x * 10.f;
        //ballVel_.x += angleVel_;
    }
    // check for floor collisions
    if(ballPos_.y - ballSize_.get() * 0.5f < 0.f) {
        ballPos_.y = ballSize_.get() * 0.5f;
        ballVel_.y *= -1.f;
        variate(ballVel_);
        angleVel_ += ballVel_.x * 10.f;
        //ballVel_.x += angleVel_;
    }

    //tgt::vec2 paddle1Vel = paddle1Pos_ - oldPaddle1Pos;
    //tgt::vec2 paddle2Vel = paddle2Pos_ - oldPaddle2Pos;

    // check for left paddle collision
    if(ballPos_.x - ballSize_.get() * 0.5f < paddle1Pos_.x + paddleWidth_.get() * 0.5f &&
       ballPos_.x + ballSize_.get() * 0.5f > paddle1Pos_.x - paddleWidth_.get() * 0.5f &&
       ballPos_.y < paddle1Pos_.y + paddleHeight_.get() * 0.5f &&
       ballPos_.y > paddle1Pos_.y - paddleHeight_.get() * 0.5f) {

        ballPos_.x = paddle1Pos_.x + paddleWidth_.get() * 0.5f + ballSize_.get() * 0.5f;

        ballVel_.x *= -1.f;

        //ballVel_ += tgt::vec2(std::max(0.f, paddle1Vel.x), paddle1Vel.y * angleVel_);

        variate(ballVel_);
        angleVel_ += ballVel_.y * 10.f;
    }

    // check for right paddle collision
    if(ballPos_.x + ballSize_.get() * 0.5f > paddle2Pos_.x - paddleWidth_.get() * 0.5f &&
       ballPos_.x - ballSize_.get() * 0.5f < paddle2Pos_.x + paddleWidth_.get() * 0.5f &&
       ballPos_.y < paddle2Pos_.y + paddleHeight_.get() * 0.5f &&
       ballPos_.y > paddle2Pos_.y - paddleHeight_.get() * 0.5f) {

        ballPos_.x = paddle2Pos_.x - paddleWidth_.get() * 0.5f - ballSize_.get() * 0.5f;

        //ballVel_ += tgt::vec2(std::min(0.f, paddle1Vel.x), paddle1Vel.y * angleVel_);

        ballVel_.x *= -1.f;
        variate(ballVel_);
        angleVel_ += ballVel_.y * 10.f;
    }

    if(length(ballVel_) < speedFactor_.get() * 0.01f)
        ballVel_ = normalize(ballVel_) * speedFactor_.get() * 0.01f;
    if(length(ballVel_) > speedFactor_.get() * 0.02f)
        ballVel_ = normalize(ballVel_) * speedFactor_.get() * 0.02f;

    // crazy mode
    if(crazyMode_.get())
        variate(ballVel_);

    ballPos_ += ballVel_;
    angleVel_ = tgt::clamp(angleVel_, -0.5f, 0.5f);
    angle_ += angleVel_;
    if(angle_ > 360.f)
        angle_ -= 360.f;

    invalidate();
}

void Pong::variate(tgt::vec2& v) {
    tgt::vec2 tmp = v;
    v.x += getRandomFloat() * speedFactor_.get() * 0.001f * randomFactor_.get();
    if(tgt::sign(v.x) != tgt::sign(tmp.x))
        v.x = tmp.x;
    v.y += getRandomFloat() * speedFactor_.get() * 0.001f * randomFactor_.get();
    if(tgt::sign(v.y) != tgt::sign(tmp.y))
        v.y = tmp.y;
}

std::string Pong::generateHeader(const tgt::GpuCapabilities::GlVersion* version) {
    std::string header = ImageProcessor::generateHeader(version);
    header += "#define MODE_ALPHA_BLENDING\n";
    return header;
}

} // namespace voreen
