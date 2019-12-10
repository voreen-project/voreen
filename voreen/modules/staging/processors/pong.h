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

#ifndef VRN_PONG_H
#define VRN_PONG_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/colorproperty.h"

#include "voreen/core/ports/renderport.h"

#include "tgt/event/touchevent.h"
#include "tgt/event/touchpoint.h"
#include "tgt/font.h"


namespace voreen {

class Pong : public ImageProcessor {

public:

    Pong();
    virtual ~Pong();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "Pong"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }

    virtual Processor* create() const { return new Pong(); }

    virtual void process();
    virtual bool isReady() const;
    virtual void initialize();
    virtual void deinitialize();

protected:
    virtual void setDescriptions() {
        setDescription("Play some Pong! Use your fingers to move the paddles (or, alternatively, WASD and arrow keys). Use mouth to generate beeping sounds.");
    }

    virtual void timerEvent(tgt::TimeEvent* e);
    virtual void touchEvent(tgt::TouchEvent* e);
    virtual void mouseDoubleClickEvent(tgt::MouseEvent* e);
    virtual void keyEvent(tgt::KeyEvent* e);
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version);

    void onReset();
    void leftAIOnChange();
    void rightAIOnChange();

    float getRandomFloat() const {
        float f = 0.f;
        while(f == 0.f) {
            f = (float)std::rand() / RAND_MAX;
        }
        return 2.f * f  - 1.f;
    }

    void variate(tgt::vec2& v);

    FloatProperty paddleHeight_;
    FloatProperty paddleWidth_;
    FloatProperty ballSize_;
    ButtonProperty resetGame_;
    BoolProperty crazyMode_;
    ColorProperty color_;
    FloatProperty opacity_;
    FloatProperty speedFactor_;
    FloatProperty paddleAutoSpeed_;
    IntProperty randomFactor_;
    BoolProperty computerOpponent1_;
    BoolProperty computerOpponent2_;
    BoolProperty xmasMode_;

    RenderPort inport_;
    RenderPort outport_;
    RenderPort privatePort_;

    tgt::vec2 ballPos_;
    tgt::vec2 ballVel_;
    tgt::vec2 paddle1Pos_;
    tgt::vec2 paddle2Pos_;
    bool paddle1Move_;
    bool paddle2Move_;
    tgt::vec2 leftTarget_;
    tgt::vec2 rightTarget_;
    float angle_;
    float angleVel_;

    tgt::EventHandler eventHandler_;    // A local eventhandler which is added to the timer.
    bool upPressedL_;
    bool downPressedL_;
    bool leftPressedL_;
    bool rightPressedL_;
    bool upPressedR_;
    bool downPressedR_;
    bool leftPressedR_;
    bool rightPressedR_;

    unsigned short score1_;
    unsigned short score2_;
    int sleepCounter_;
    int sleepTime_;
    tgt::Timer* timer_;
    bool renderScore_;
    tgt::Texture* xmasBall_;
    tgt::Texture* xmasFlake_;
    std::vector<std::pair<tgt::vec2, tgt::vec2> > snow_;

    tgt::Font* font_;

    tgt::Shader* copyShader_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
