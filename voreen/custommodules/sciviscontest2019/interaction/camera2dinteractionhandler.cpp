/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2014 University of Muenster, Germany.                        *
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

#include "voreen/core/voreenapplication.h"
#include "camera2dinteractionhandler.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include "tgt/event/mouseevent.h"
#include "tgt/event/keyevent.h"
#include "tgt/event/timeevent.h"
#include "tgt/event/touchevent.h"
#include "tgt/event/touchpoint.h"
#include "tgt/timer.h"
#include "voreen/core/properties/eventproperty.h"

using tgt::Event;
using tgt::MouseEvent;
using tgt::KeyEvent;
using tgt::TimeEvent;
using tgt::TouchPoint;
using tgt::TouchEvent;

using tgt::vec2;
using tgt::vec3;
using tgt::ivec2;

namespace voreen {

Camera2DInteractionHandler::Camera2DInteractionHandler() :
    InteractionHandler("dummy", "dummy")
{}

Camera2DInteractionHandler::Camera2DInteractionHandler(const std::string& id, const std::string& guiName, FloatVec2Property* center, FloatProperty* width, const size_t* outPortWidth, bool sharing, bool enabled)
    : InteractionHandler(id, guiName)
    , center_(center)
    , width_(width)
    , outPortWidth_(outPortWidth)
{
    multiTouchEvent_ = new EventProperty<Camera2DInteractionHandler>(id + ".multitouch", guiName + "Multitouch", this,
        &Camera2DInteractionHandler::multiTouchEvent, sharing, enabled);

    addEventProperty(multiTouchEvent_);
    zoomEvent_ = new EventProperty<Camera2DInteractionHandler>(id + ".zoom", "Zoom", this,
        &Camera2DInteractionHandler::zoomEvent,
        MouseEvent::MOUSE_BUTTON_RIGHT,
        MouseEvent::ACTION_ALL,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    addEventProperty(zoomEvent_);

    shiftEvent_ = new EventProperty<Camera2DInteractionHandler>(id + ".shift", "Shift", this,
        &Camera2DInteractionHandler::shiftEvent,
        MouseEvent::MOUSE_BUTTON_LEFT,
        MouseEvent::ACTION_ALL,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    addEventProperty(shiftEvent_);

    wheelZoomEvent_ = new EventProperty<Camera2DInteractionHandler>(id + ".wheelZoom", "Wheel Zoom", this,
        &Camera2DInteractionHandler::zoomEvent,
        MouseEvent::MOUSE_WHEEL,
        MouseEvent::WHEEL,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    addEventProperty(wheelZoomEvent_);
}

Camera2DInteractionHandler::~Camera2DInteractionHandler() {
    // event properties are deleted by base class InteractionHandler
}

void Camera2DInteractionHandler::multiTouchEvent(tgt::TouchEvent* e) {
    if (e->touchPoints().size() == 2) {

        if (e->touchPointStates() & TouchPoint::TouchPointPressed) {
            vec2 touchPoint1 = e->touchPoints()[0].pos();
            vec2 touchPoint2 = e->touchPoints()[1].pos();

            lastDistance_ = length (touchPoint1 - touchPoint2);

            e->accept();
        }

        // Camera2DInteractionHandler only handles touch event to zoom
        else if (e->touchPointStates() & TouchPoint::TouchPointMoved) {

            vec2 pointPos1 = e->touchPoints()[0].pos();
            vec2 pointPos2 = e->touchPoints()[1].pos();
            vec2 middle = (pointPos1+pointPos2)*0.5f;
            vec2 middleOffset = center_->get() - middle;

            float newDistance = length(pointPos1 - pointPos2);
            float zoomFactor = lastDistance_ / newDistance;

            width_->set(width_->get() * zoomFactor);
            //Set center to zoom in on the middle of the two touch points
            center_->set(middle+middleOffset*zoomFactor);

            lastDistance_ = newDistance;
            e->accept();
        }
    }
}
void Camera2DInteractionHandler::zoomEvent(tgt::MouseEvent* e) {
    const float zoomSpeed = 0.05f;
    if (e->action() == MouseEvent::PRESSED) {
        QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
        lastMousePos_ = e->coord();
        e->accept();
    }
    else if (e->action() == MouseEvent::RELEASED) {
        QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
        lastMousePos_ = e->coord();
        e->accept();
    }
    else if (e->action() == MouseEvent::MOTION) {
        ivec2 mOffset = e->coord() - lastMousePos_;
        float mod = 1.0f + ((float)mOffset.y * zoomSpeed);

        width_->set(width_->get() * mod);

        lastMousePos_ = e->coord();
        e->accept();
    }
    else if (e->action() == MouseEvent::WHEEL) {
        float mod = 1.0f;
        if (e->button() == MouseEvent::MOUSE_WHEEL_UP)
            mod -= zoomSpeed;
        else
            mod += zoomSpeed;

        width_->set(width_->get() * mod);

        e->accept();
    }

    // invalidate processor and update camera prop widgets, if event has been accepted
    if (e->isAccepted()) {
        center_->invalidate();
    }
}

void Camera2DInteractionHandler::shiftEvent(tgt::MouseEvent* e) {

    if (e->action() == MouseEvent::PRESSED) {
        QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
        lastMousePos_ = e->coord();
        e->accept();
    }
    else if (e->action() == MouseEvent::RELEASED) {
        QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
        e->accept();
    }
    else if (e->action() == MouseEvent::MOTION) {
        ivec2 mOffset = e->coord() - lastMousePos_;
        tgt::vec2 newCenter(center_->get());
        float conversionFactor = width_->get() / *outPortWidth_;
        newCenter.x -= mOffset.x * conversionFactor;
        newCenter.y += mOffset.y * conversionFactor;
        center_->set(newCenter);

        lastMousePos_ = e->coord();
        e->accept();
    }
    else if (e->action() == MouseEvent::DOUBLECLICK) {
        e->accept();
    }

    // invalidate processor and update camera prop widgets, if event has been accepted
    if (e->isAccepted()) {
        center_->invalidate();
    }
}

void Camera2DInteractionHandler::onEvent(Event* eve) {

    // invalidate processor and update camera prop widgets, if event has been accepted
    if (eve->isAccepted()) {
        center_->invalidate();
    }
}

} // namespace
