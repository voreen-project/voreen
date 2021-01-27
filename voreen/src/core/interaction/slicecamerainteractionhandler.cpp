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

#include "voreen/core/voreenapplication.h"
#include "voreen/core/interaction/slicecamerainteractionhandler.h"
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

SliceCameraInteractionHandler::SliceCameraInteractionHandler() :
    InteractionHandler("dummy", "dummy")
{}

SliceCameraInteractionHandler::SliceCameraInteractionHandler(const std::string& id, const std::string& guiName,
                                                   CameraProperty* cameraProp, bool sharing, bool enabled)
    : InteractionHandler(id, guiName)
    , cameraProp_(cameraProp)
{
    tgtAssert(cameraProp, "No camera property");

    multiTouchEvent_ = new EventProperty<SliceCameraInteractionHandler>(id + ".multitouch", guiName + "Multitouch", this,
        &SliceCameraInteractionHandler::handleMultitouch, sharing, enabled);
    addEventProperty(multiTouchEvent_);

    zoomEvent_ = new EventProperty<SliceCameraInteractionHandler>(id + ".zoom", "Zoom", this,
        &SliceCameraInteractionHandler::zoomEvent,
        MouseEvent::MOUSE_BUTTON_RIGHT,
        MouseEvent::ACTION_ALL,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    addEventProperty(zoomEvent_);

    shiftEvent_ = new EventProperty<SliceCameraInteractionHandler>(id + ".shift", "Shift", this,
        &SliceCameraInteractionHandler::shiftEvent,
        MouseEvent::MOUSE_BUTTON_LEFT,
        MouseEvent::ACTION_ALL,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    addEventProperty(shiftEvent_);

    wheelZoomEvent_ = new EventProperty<SliceCameraInteractionHandler>(id + ".wheelZoom", "Wheel Zoom", this,
        &SliceCameraInteractionHandler::zoomEvent,
        MouseEvent::MOUSE_WHEEL,
        MouseEvent::WHEEL,
        tgt::Event::SHIFT, sharing, enabled);
    addEventProperty(wheelZoomEvent_);
}

SliceCameraInteractionHandler::~SliceCameraInteractionHandler() {
    // event properties are deleted by base class InteractionHandler
}

void SliceCameraInteractionHandler::handleMultitouch(tgt::TouchEvent* e) {
    tgtAssert(cameraProp_, "No camera property");

    if (e->touchPoints().size() == 2) {

        if (e->touchPointStates() & TouchPoint::TouchPointPressed) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            vec2 touchPoint1 = e->touchPoints()[0].pos();
            vec2 touchPoint2 = e->touchPoints()[1].pos();

            lastDistance_ = length (touchPoint1 - touchPoint2);

            e->accept();
            //cameraProp_->invalidate();
        }
        else if (e->touchPointStates() & TouchPoint::TouchPointReleased) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            e->accept();
            cameraProp_->invalidate();
        }

        // slicecamerainteractionhandler only handles touch event to zoom
        else if (e->touchPointStates() & TouchPoint::TouchPointMoved) {

            vec2 pointPos1 = e->touchPoints()[0].pos();
            vec2 pointPos2 = e->touchPoints()[1].pos();

            float newDistance = length(pointPos1 - pointPos2);
            float zoomFactor = lastDistance_ / newDistance;
            float mod = zoomFactor;

            tgt::Frustum f = cameraProp_->get().getFrustum();

            f.setLeft(f.getLeft()*mod);
            f.setRight(f.getRight()*mod);
            f.setBottom(f.getBottom()*mod);
            f.setTop(f.getTop()*mod);

            cameraProp_->setFrustum(f);

            lastDistance_ = newDistance;
            e->accept();

            cameraProp_->invalidate();
        }
    }
}

void SliceCameraInteractionHandler::zoomEvent(tgt::MouseEvent* e) {
    tgtAssert(cameraProp_, "No camera property");


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
        float mod = 1.0f + ((float)mOffset.y * 0.01f);

        tgt::Frustum f = cameraProp_->get().getFrustum();

        f.setLeft(f.getLeft()*mod);
        f.setRight(f.getRight()*mod);
        f.setBottom(f.getBottom()*mod);
        f.setTop(f.getTop()*mod);

        cameraProp_->setFrustum(f);

        lastMousePos_ = e->coord();
        e->accept();
    }
    else if (e->action() == MouseEvent::WHEEL) {
        float mod = 1.0f;
        if (e->button() == MouseEvent::MOUSE_WHEEL_UP)
            mod -= 0.01f;
        else
            mod += 0.01f;

        tgt::Frustum f = cameraProp_->get().getFrustum();

        f.setLeft(f.getLeft()*mod);
        f.setRight(f.getRight()*mod);
        f.setBottom(f.getBottom()*mod);
        f.setTop(f.getTop()*mod);

        cameraProp_->setFrustum(f);

        e->accept();
    }

    // invalidate processor and update camera prop widgets, if event has been accepted
    if (e->isAccepted()) {
        cameraProp_->invalidate();
    }
}

void SliceCameraInteractionHandler::shiftEvent(tgt::MouseEvent* e) {
    tgtAssert(cameraProp_, "No camera property");

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
        //perform panning motion:
        tgt::line3 lineCur = cameraProp_->get().getViewRay(e->viewport(), e->coord());
        tgt::line3 linePrev = cameraProp_->get().getViewRay(e->viewport(), lastMousePos_);
        vec3 offset = linePrev.getStart() - lineCur.getStart();

        cameraProp_->setPosition(cameraProp_->get().getPosition() + offset);
        cameraProp_->setFocus(cameraProp_->get().getFocus() + offset);

        lastMousePos_ = e->coord();
        e->accept();
    }
    else if (e->action() == MouseEvent::DOUBLECLICK) {
        tgt::line3 lineCur = cameraProp_->get().getViewRay(e->viewport(), e->coord());
        tgt::line3 lineCenter = cameraProp_->get().getViewRay(e->viewport(), e->viewport() / 2);
        vec3 offset = lineCur.getStart() - lineCenter.getStart();

        cameraProp_->setPosition(cameraProp_->get().getPosition() + offset);
        cameraProp_->setFocus(cameraProp_->get().getFocus() + offset);

        e->accept();
    }

    // invalidate processor and update camera prop widgets, if event has been accepted
    if (e->isAccepted()) {
        cameraProp_->invalidate();
    }
}

void SliceCameraInteractionHandler::onEvent(Event* eve) {
    tgtAssert(cameraProp_, "No camera property");

    // invalidate processor and update camera prop widgets, if event has been accepted
    if (eve->isAccepted()) {
        cameraProp_->invalidate();
    }
}

} // namespace
