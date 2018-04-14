/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "toucheventsimulator.h"

namespace voreen {

TouchEventSimulator::TouchEventSimulator()
    : ImageProcessor("copyimage")
    , inport_(Port::INPORT, "inport", "Inport")
    , outport_(Port::OUTPORT, "outport", "Outport")
{
    //register ports
    addPort(inport_);
    addPort(outport_);

    mouseEventProp_ = new EventProperty<TouchEventSimulator>(
        "toucheventsimulator.event", "Mouse Event", this,
        &TouchEventSimulator::onMouseEvent,
        tgt::MouseEvent::MOUSE_BUTTON_LEFT,
        tgt::MouseEvent::PRESSED | tgt::MouseEvent::RELEASED | tgt::MouseEvent::MOTION,
        tgt::Event::CTRL);

    addEventProperty(mouseEventProp_);

    shaderProp_.setVisibleFlag(false);
}

TouchEventSimulator::~TouchEventSimulator() {
    delete mouseEventProp_;
}

Processor* TouchEventSimulator::create() const {
    return new TouchEventSimulator();
}

std::string TouchEventSimulator::getClassName() const {
    return "TouchEventSimulator";
}

std::string TouchEventSimulator::getCategory() const {
    return "Touch Table";
}

void TouchEventSimulator::setDescriptions() {
    setDescription("Generates touch events from mouse events when STRG is pressed.");
}

void TouchEventSimulator::process() {

    // activate and clear output render target
    outport_.activateTarget();
    outport_.clearTarget();

    // bind input rendering to texture units
    tgt::TextureUnit colorUnit, depthUnit;
    inport_.bindTextures(colorUnit.getEnum(), depthUnit.getEnum());

    // activate shader and set uniforms
    tgtAssert(program_, "bypass shader not loaded");
    program_->activate();
    setGlobalShaderParameters(program_);
    inport_.setTextureParameters(program_, "texParams_");
    program_->setUniform("useTexcoords", false);
    program_->setUniform("colorTex_", colorUnit.getUnitNumber());
    program_->setUniform("depthTex_", depthUnit.getUnitNumber());

    // render screen aligned quad
    renderQuad();

    // cleanup
    program_->deactivate();
    outport_.deactivateTarget();
    tgt::TextureUnit::setZeroUnit();

    // check for OpenGL errors
    LGL_ERROR;
}

void TouchEventSimulator::onMouseEvent(tgt::MouseEvent* e) {

    // use continuous IDs for the touch events
    static int id = 0;

    // set attributes
    int deviceType = tgt::TouchEvent::TouchScreen;

    std::deque<tgt::TouchPoint> tps;

    tgt::vec2 pos = e->coord();

    tgt::TouchPoint::State state = static_cast<tgt::TouchPoint::State>(0);

    // translate mouse event state
    if (e->action() & tgt::MouseEvent::PRESSED) {
        state = tgt::TouchPoint::TouchPointPressed;
        id++;
    }
    else if(e->action() & tgt::MouseEvent::MOTION)
        state = tgt::TouchPoint::TouchPointMoved;
    else if (e->action() & tgt::MouseEvent::RELEASED)
        state = tgt::TouchPoint::TouchPointReleased;

    bool primary = true;

    tgt::TouchPoint::State states = state;

    // create touch point
    tgt::TouchPoint tp;
    tp.setId(id);
    tp.setPos(pos);
    tp.setState(state);
    tp.setPrimary(primary);

    tps.push_back(tp);

    // create event and add touch point
    tgt::TouchEvent* te = new tgt::TouchEvent(tgt::Event::MODIFIER_NONE, states, (tgt::TouchEvent::DeviceType)deviceType, tps, outport_.getSize());
    te->ignore();

    // send touch event
    inport_.distributeEvent(te);
}

} // namespace
