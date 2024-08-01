/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/interaction/trackballnavigation.h"
#include "voreen/core/interaction/firstpersonnavigation.h"
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
using tgt::TouchEvent;
using tgt::TouchPoint;

namespace voreen {

CameraInteractionHandler::CameraInteractionHandler() :
    InteractionHandler("dummy", "dummy"),
    navigationMetaphor_("naviMetaphor", "Navigation Metaphor", Processor::VALID)
{}

CameraInteractionHandler::CameraInteractionHandler(const std::string& id, const std::string& guiName,
                                                   CameraProperty* cameraProp, bool sharing, bool enabled)
    : InteractionHandler(id, guiName)
    , cameraProp_(cameraProp)
    , navigationMetaphor_(id + ".interactionMetaphor", guiName + " Interaction", Processor::VALID)
{
    tgtAssert(cameraProp, "No camera property");
    cameraProp_ = cameraProp;

    navigationMetaphor_.addOption("trackball",     "Trackball",     TRACKBALL);
    navigationMetaphor_.addOption("first-person",  "First-Person",  FIRST_PERSON);
    navigationMetaphor_.select("trackball");
    navigationMetaphor_.onChange(MemberFunctionCallback<CameraInteractionHandler>(this, &CameraInteractionHandler::adjustWidgetStates));
    addProperty(navigationMetaphor_);

    // navigations
    tbNavi_ = new TrackballNavigation(cameraProp_, TrackballNavigation::ROTATE_MODE);
    fpNavi_ = new FirstPersonNavigation(cameraProp_);

    // event properties trackball
    multiTouchEvent_ = new EventProperty<CameraInteractionHandler>(id + ".multitouch", guiName + "Multitouch", this,
        &CameraInteractionHandler::handleMultitouch, sharing, enabled);
    addEventProperty(multiTouchEvent_);

    rotateEvent_ = new EventProperty<CameraInteractionHandler>(id + ".rotate", guiName + " Rotate", this,
        &CameraInteractionHandler::rotateEvent,
        MouseEvent::MOUSE_BUTTON_LEFT,
        MouseEvent::ACTION_ALL,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    addEventProperty(rotateEvent_);

    contRotateEvent_ = new EventProperty<CameraInteractionHandler>(id + ".rotatecont", guiName + " Rotate Continuously", this,
        &CameraInteractionHandler::rotateEvent,
        MouseEvent::MOUSE_BUTTON_LEFT,
        MouseEvent::ACTION_ALL,
        tgt::Event::CTRL, sharing, enabled);
    addEventProperty(contRotateEvent_);

    zoomEvent_ = new EventProperty<CameraInteractionHandler>(id + ".zoom", guiName + " Zoom", this,
        &CameraInteractionHandler::zoomEvent,
        MouseEvent::MOUSE_BUTTON_RIGHT,
        MouseEvent::ACTION_ALL,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    addEventProperty(zoomEvent_);

    shiftEvent_ = new EventProperty<CameraInteractionHandler>(id + ".shift", guiName + " Shift", this,
        &CameraInteractionHandler::shiftEvent,
        MouseEvent::MOUSE_BUTTON_LEFT,
        MouseEvent::ACTION_ALL,
        tgt::Event::SHIFT, sharing, enabled);
    addEventProperty(shiftEvent_);

    wheelZoomEvent_ = new EventProperty<CameraInteractionHandler>(id + ".wheelZoom", guiName + " Wheel Zoom", this,
        &CameraInteractionHandler::zoomEvent,
        MouseEvent::MOUSE_WHEEL,
        MouseEvent::WHEEL,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    addEventProperty(wheelZoomEvent_);

    // event properties first-person navigation
    moveForwardEvent_ = new EventProperty<CameraInteractionHandler>(id + ".moveForward", guiName + " Move Forward", this,
        &CameraInteractionHandler::keyEvent,
        tgt::KeyEvent::K_W,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    moveForwardEvent_->onChange(&CameraInteractionHandler::updateFpKeySettings);
    addEventProperty(moveForwardEvent_);

    moveLeftEvent_ = new EventProperty<CameraInteractionHandler>(id + ".moveLeft", guiName + " Move Left", this,
        &CameraInteractionHandler::keyEvent,
        tgt::KeyEvent::K_A,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    moveLeftEvent_->onChange(&CameraInteractionHandler::updateFpKeySettings);
    addEventProperty(moveLeftEvent_);

    moveBackwardEvent_ = new EventProperty<CameraInteractionHandler>(id + ".moveBackward", guiName + " Move Backward", this,
        &CameraInteractionHandler::keyEvent,
        tgt::KeyEvent::K_S,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    moveBackwardEvent_->onChange(&CameraInteractionHandler::updateFpKeySettings);
    addEventProperty(moveBackwardEvent_);

    moveRightEvent_ = new EventProperty<CameraInteractionHandler>(id + ".moveRight", guiName + " Move Right", this,
        &CameraInteractionHandler::keyEvent,
        tgt::KeyEvent::K_D,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    moveRightEvent_->onChange(&CameraInteractionHandler::updateFpKeySettings);
    addEventProperty(moveRightEvent_);

    moveUpEvent_ = new EventProperty<CameraInteractionHandler>(id + ".moveUp", guiName + " Move Up", this,
        &CameraInteractionHandler::keyEvent,
        tgt::KeyEvent::K_R,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    moveUpEvent_->onChange(&CameraInteractionHandler::updateFpKeySettings);
    addEventProperty(moveUpEvent_);

    moveDownEvent_ = new EventProperty<CameraInteractionHandler>(id + ".moveDown", guiName + " Move Down", this,
        &CameraInteractionHandler::keyEvent,
        tgt::KeyEvent::K_F,
        tgt::Event::MODIFIER_NONE, sharing, enabled);
    moveDownEvent_->onChange(&CameraInteractionHandler::updateFpKeySettings);
    addEventProperty(moveDownEvent_);

    // timer
    timerEventHandler_ = new tgt::EventHandler;
    motionTimer_ = VoreenApplication::app()->createTimer(timerEventHandler_);
    tgtAssert(timerEventHandler_, "No timerEventHandler");
    timerEventHandler_->addListenerToFront(this);

    adjustWidgetStates();
}

CameraInteractionHandler::~CameraInteractionHandler() {

    // event properties are deleted by base class InteractionHandler

    delete tbNavi_;
    delete fpNavi_;
    delete timerEventHandler_;
    delete motionTimer_;

    // Reset quality.
    QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
}

// handle multi touch event
// TODO: react on events with more than two touch points
void CameraInteractionHandler::handleMultitouch(tgt::TouchEvent* e){
    tgtAssert(cameraProp_, "No camera property");
    tgtAssert(fpNavi_, "No trackball navigation");
    if (e->touchPoints().size() == 1) {
        if (!navigationMetaphor_.isSelected("trackball"))
            return;

        tgtAssert(cameraProp_, "No camera property");
        tgtAssert(tbNavi_, "No trackball navigation");

        tgt::TouchPoint tp = e->touchPoints().front();

        // propagate event to trackball navigation
        if (tp.state() & TouchPoint::TouchPointPressed) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            tgt::MouseEvent event = tgt::MouseEvent((int)tp.pos().x, (int)tp.pos().y, tgt::MouseEvent::PRESSED, tgt::Event::MODIFIER_NONE, tgt::MouseEvent::MOUSE_BUTTON_LEFT, e->getViewport());
            tbNavi_->mousePressEvent(&event);
            e->accept();
            if(motionTimer_) motionTimer_->stop();
        }
        else if (tp.state() & TouchPoint::TouchPointReleased) {
            tgt::MouseEvent event = tgt::MouseEvent((int)tp.pos().x, (int)tp.pos().y, tgt::MouseEvent::RELEASED, tgt::Event::MODIFIER_NONE, tgt::MouseEvent::MOUSE_BUTTON_LEFT, e->getViewport());
            tbNavi_->mouseReleaseEvent(&event);
            e->accept();
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            cameraProp_->invalidate(); //< necessary to cause processor to re-render without interaction mode
        }
        else if (tp.state() & (TouchPoint::TouchPointMoved | TouchPoint::TouchPointStationary)) {
            tgt::MouseEvent event = tgt::MouseEvent((int)tp.pos().x, (int)tp.pos().y, tgt::MouseEvent::MOTION, tgt::Event::MODIFIER_NONE, tgt::MouseEvent::MOUSE_BUTTON_LEFT, e->getViewport());
            tbNavi_->mouseMoveEvent(&event);
            e->accept();
            cameraProp_->invalidate();
        }
    }

    if (e->touchPoints().size() == 2) {

        // propagate event to trackball navigation
        if (e->touchPointStates() & TouchPoint::TouchPointPressed) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            tbNavi_->touchPressEvent(e);
            e->accept();
            //cameraProp_->invalidate();
        }
        else if (e->touchPointStates() & TouchPoint::TouchPointReleased) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            tbNavi_->touchReleaseEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
        else if (e->touchPointStates() & (TouchPoint::TouchPointMoved | TouchPoint::TouchPointStationary)) {
            tbNavi_->touchMoveEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
    }
}

// TODO: left clicked mousedragging is catched by fpNavi, too. Maybe this method should get a proper name?!
void CameraInteractionHandler::rotateEvent(tgt::MouseEvent* e) {
    if (navigationMetaphor_.isSelected("trackball")){
        tgtAssert(cameraProp_, "No camera property");
        tgtAssert(tbNavi_, "No trackball navigation");

        // propagate event to trackball navigation
        if (e->action() == MouseEvent::PRESSED) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            tbNavi_->mousePressEvent(e);
            e->accept();
            if(motionTimer_) motionTimer_->stop();
        }
        else if (e->action() == MouseEvent::RELEASED) {
            tbNavi_->mouseReleaseEvent(e);
            e->accept();
            // if continuous modifier was used, start timer to coninuously spin the trackball
            if((e->modifiers() & contRotateEvent_->getModifier()) && motionTimer_)
                motionTimer_->start(10, 0);
            else
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            cameraProp_->invalidate(); //< necessary to cause processor to re-render without interaction mode
        }
        else if (e->action() == MouseEvent::MOTION) {
            tbNavi_->mouseMoveEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
        else if (e->action() == MouseEvent::DOUBLECLICK) {
            tbNavi_->mouseDoubleClickEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
        else if (e->action() == MouseEvent::WHEEL) {
            tbNavi_->wheelEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
    }
    else if (navigationMetaphor_.isSelected("first-person")){
        if (e->action() == MouseEvent::PRESSED){
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            fpNavi_->mousePressEvent(e);
            e->accept();
        }
        else if (e->action() == MouseEvent::RELEASED){
            if (!fpNavi_->isMoving()) {
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
                e->accept();
                cameraProp_->invalidate();
            }
        }
        else if (e->action() == MouseEvent::MOTION){
            if(motionTimer_) motionTimer_->stop();

            fpNavi_->mouseMoveEvent(e);
            cameraProp_->invalidate();

            // Restart motionTimer_, if necessary
            if (fpNavi_->isMoving() && motionTimer_)
                motionTimer_->start(40, 0);
        }
    }

}

// TODO: right clicked mousedragging is catched by fpNavi, too. Maybe this method should get a proper name?!
void CameraInteractionHandler::zoomEvent(tgt::MouseEvent* e) {
    if (navigationMetaphor_.isSelected("trackball")){
        tgtAssert(cameraProp_, "No camera property");
        tgtAssert(tbNavi_, "No trackball navigation");

        tbNavi_->setMode(TrackballNavigation::ZOOM_MODE);

        // propagate event to trackball navigation
        if (e->action() == MouseEvent::PRESSED) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            tbNavi_->mousePressEvent(e);
            e->accept();
            //cameraProp_->invalidate(); //< no need to issue an invalidation here, since camera has not changed yet
        }
        else if (e->action() == MouseEvent::RELEASED) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            tbNavi_->mouseReleaseEvent(e);
            e->accept();
            cameraProp_->invalidate(); //< necessary to cause processor to re-render without interaction mode
        }
        else if (e->action() == MouseEvent::MOTION) {
            tbNavi_->mouseMoveEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
        else if (e->action() == MouseEvent::DOUBLECLICK) {
            tbNavi_->mouseDoubleClickEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
        else if (e->action() == MouseEvent::WHEEL) {
            tbNavi_->wheelEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
        tbNavi_->setMode(TrackballNavigation::ROTATE_MODE);
    }
    else if (navigationMetaphor_.isSelected("first-person")){
        if (e->action() == MouseEvent::PRESSED){
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            fpNavi_->mousePressEvent(e);
            e->accept();
            //cameraProp_->invalidate();
        }
        else if (e->action() == MouseEvent::RELEASED){
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            e->accept();
            cameraProp_->invalidate();
        }
        else if (e->action() == MouseEvent::MOTION){
            if(motionTimer_) motionTimer_->stop();

            fpNavi_->mouseMoveEvent(e);
            e->accept();
            cameraProp_->invalidate();

            // Restart motionTimer_, if necessary
            if (fpNavi_->isMoving() && motionTimer_)
                motionTimer_->start(40, 0);
        }
    }
}

void CameraInteractionHandler::shiftEvent(tgt::MouseEvent* e) {
    if (navigationMetaphor_.isSelected("trackball")){
        tgtAssert(cameraProp_, "No camera property");
        tgtAssert(tbNavi_, "No trackball navigation");

        tbNavi_->setMode(TrackballNavigation::SHIFT_MODE);

        // propagate event to trackball navigation
        if (e->action() == MouseEvent::PRESSED) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
            tbNavi_->mousePressEvent(e);
            e->accept();
            //cameraProp_->invalidate(); //< no need to issue an invalidation here, since camera has not changed yet
        }
        else if (e->action() == MouseEvent::RELEASED) {
            QualityMode.requestQualityMode(VoreenQualityMode::RQ_DEFAULT, this);
            tbNavi_->mouseReleaseEvent(e);
            e->accept();
            cameraProp_->invalidate(); //< necessary to cause processor to re-render without interaction mode
        }
        else if (e->action() == MouseEvent::MOTION) {
            tbNavi_->mouseMoveEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
        else if (e->action() == MouseEvent::DOUBLECLICK) {
            tbNavi_->mouseDoubleClickEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }
        else if (e->action() == MouseEvent::WHEEL) {
            tbNavi_->wheelEvent(e);
            e->accept();
            cameraProp_->invalidate();
        }

        tbNavi_->setMode(TrackballNavigation::ROTATE_MODE);
    }
    else if (navigationMetaphor_.isSelected("first-person")){
        // Hack: prevents camerainteraction handler from executing shift events, if fpNavi is activated
        e->accept();
    }

}

void CameraInteractionHandler::keyEvent(tgt::KeyEvent* e){
    tgtAssert(cameraProp_, "No camera property");

    tgtAssert(fpNavi_, "No first person Navigation");

    if (navigationMetaphor_.isSelected("first-person")){
        fpNavi_->keyEvent(e);

        if (e->pressed()){
            if (fpNavi_->isMoving()){
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
                if(motionTimer_) motionTimer_->start(40, 0);
            }
            else
                if(motionTimer_) motionTimer_->stop();
        }
        else{
            if (!fpNavi_->isMoving()){
                QualityMode.requestQualityMode((pressedMouseButtons_.any() ? VoreenQualityMode::RQ_INTERACTIVE : VoreenQualityMode::RQ_DEFAULT), this);
                if(motionTimer_) motionTimer_->stop();
            }
        }

        e->accept();
    }
}

void CameraInteractionHandler::onEvent(Event* eve) {
    tgtAssert(cameraProp_, "No camera property");
    tgtAssert(fpNavi_, "No first person Navigation");

    if (tgt::TimeEvent* timeEve = dynamic_cast<tgt::TimeEvent*>(eve)){
        if (navigationMetaphor_.isSelected("first-person")){
            if (fpNavi_->isMoving()){
                QualityMode.requestQualityMode(VoreenQualityMode::RQ_INTERACTIVE, this);
                fpNavi_->timerEvent(timeEve);
                eve->accept();
            }
        } else if (navigationMetaphor_.isSelected("trackball")) {
            tgt::quat q = cameraProp_->getTrackball().getLastOrientationChange();
            cameraProp_->getTrackball().rotate(q);
            eve->accept();
        }
    }

    // invalidate processor and update camera prop widgets, if event has been accepted
    if (eve->isAccepted()) {
        cameraProp_->invalidate();
    }
}

void CameraInteractionHandler::setNavigationMode(NavigationMetaphor behaviour){
    navigationMetaphor_.selectByValue(behaviour);
}

CameraInteractionHandler::NavigationMetaphor CameraInteractionHandler::getNavigationMode() const {
    return navigationMetaphor_.getValue();
}

void CameraInteractionHandler::updateFpKeySettings() {
    tgtAssert(fpNavi_, "no first-person navigation");
    tgtAssert(moveForwardEvent_ && moveLeftEvent_ && moveBackwardEvent_ , "missing event property");
    tgtAssert(moveRightEvent_ && moveUpEvent_ && moveDownEvent_ , "missing event property");

    fpNavi_->setKeysetting(moveForwardEvent_->getKeyCode(), moveBackwardEvent_->getKeyCode(),
        moveLeftEvent_->getKeyCode(), moveRightEvent_->getKeyCode(),
        moveUpEvent_->getKeyCode(), moveDownEvent_->getKeyCode());
}

void CameraInteractionHandler::adjustWidgetStates() {
    tgtAssert(rotateEvent_ && zoomEvent_ && shiftEvent_ && wheelZoomEvent_, "missing event property");
    tgtAssert(moveForwardEvent_ && moveLeftEvent_ && moveBackwardEvent_ , "missing event property");
    tgtAssert(moveRightEvent_ && moveUpEvent_ && moveDownEvent_ , "missing event property");

    bool trackballSelected = navigationMetaphor_.isSelected("trackball");
    bool firstPersonSelected = navigationMetaphor_.isSelected("first-person");

    // trackball properties
    rotateEvent_->setVisibleFlag(trackballSelected);
    zoomEvent_->setVisibleFlag(trackballSelected);
    shiftEvent_->setVisibleFlag(trackballSelected);
    wheelZoomEvent_->setVisibleFlag(trackballSelected);

    // first-person properties
    moveForwardEvent_->setVisibleFlag(firstPersonSelected);
    moveLeftEvent_->setVisibleFlag(firstPersonSelected);
    moveBackwardEvent_->setVisibleFlag(firstPersonSelected);
    moveRightEvent_->setVisibleFlag(firstPersonSelected);
    moveUpEvent_->setVisibleFlag(firstPersonSelected);
    moveDownEvent_->setVisibleFlag(firstPersonSelected);
}

void CameraInteractionHandler::setVisible(bool state) {
    InteractionHandler::setVisible(state);
    navigationMetaphor_.setVisibleFlag(state);
    if (state)
        adjustWidgetStates();
}

} // namespace
