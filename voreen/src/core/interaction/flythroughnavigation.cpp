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

#include "voreen/core/interaction/flythroughnavigation.h"
#include "voreen/core/network/networkevaluator.h"

#include "tgt/spline.h"

#include <sstream>
#include <string>

namespace voreen {

using tgt::vec3;
using tgt::Navigation;
using tgt::Event;
using tgt::KeyEvent;
using tgt::MouseEvent;

const std::string FlythroughNavigation::loggerCat_ = "voreen.vis.FlythroughNavigation";

FlythroughNavigation::FlythroughNavigation(tgt::GLCanvas* canvas, tgt::Camera* camera)
    : Navigation(camera)
    , canvas_(canvas)
{

    wheelCounter_   = -1;
    moveCounter_    = -1;

    if (camera_) {
        initialPosition_ = camera_->getPosition();
        initialFocusPoint_ = camera_->getFocus();
        initialUpVector_ = camera_->getUpVector();
        initialCameraValid_ = true;
    }
    else {
        initialCameraValid_ = false;
    }

}

void FlythroughNavigation::keyEvent(tgt::KeyEvent* e) {

    bool cameraModified = false;
    if (e->pressed()) {

        float motionOffset = 0.01f;
        float motionOffsetBig = 0.06f;
        float angle = 0.01f;

        switch (e->keyCode()) {
        case KeyEvent::K_UP :
            if (e->modifiers() == 0)
                Navigation::moveCameraForward(motionOffset);
            else if (e->modifiers() == Event::SHIFT)
                Navigation::moveCameraUp(motionOffset);
            else if (e->modifiers() == Event::CTRL)
                Navigation::rollCameraVert(-angle);
            cameraModified = true;
            break;
        case KeyEvent::K_DOWN :
            if (e->modifiers() == 0)
                Navigation::moveCameraBackward(motionOffset);
            else if (e->modifiers() == Event::SHIFT)
                Navigation::moveCameraDown(motionOffset);
            else if (e->modifiers() == Event::CTRL)
                Navigation::rollCameraVert(angle);
            cameraModified = true;
            break;
        case KeyEvent::K_LEFT :
            if (e->modifiers() == 0 || e->modifiers() == Event::SHIFT)
                Navigation::moveCameraLeft(motionOffset);
            else if (e->modifiers() == Event::CTRL)
                Navigation::rotateViewHorz(angle);
            cameraModified = true;
            break;
        case KeyEvent::K_RIGHT :
            if (e->modifiers() == 0 || e->modifiers() == Event::SHIFT)
                Navigation::moveCameraRight(motionOffset);
            else if (e->modifiers() == Event::CTRL)
                Navigation::rotateViewHorz(-angle);
            cameraModified = true;
            break;
        case KeyEvent::K_PAGEUP :
            Navigation::moveCameraForward(motionOffsetBig);
            cameraModified = true;
            break;
        case KeyEvent::K_PAGEDOWN :
            Navigation::moveCameraBackward(motionOffsetBig);
            cameraModified = true;
            break;
        default:
            break;
        }
    }

    if (cameraModified) {
        //MsgDistr.postMessage( new BoolMsg(VoreenPainter::switchCoarseness_, true) );
        moveCounter_ = 0;
        canvas_->update();
    }

}


void FlythroughNavigation::mousePressEvent(tgt::MouseEvent* /*e*/) {

    /*
    if (!trackball_ || !trackball_->getCanvas() || !trackball_->getCanvas()->getPainter())
        return;

    if (trackballEnabled_) {
        postMessage( new BoolMsg(VoreenPainter::switchCoarseness_, true) );

        trackball_->mousePressEvent(e);
        if (e->button() == trackball_->getRotateButton() && trackball_->getContinuousSpin()) {
            spinit_ = true;
            moveCounter_ = 0;
        }
        // mid-button? what was supposed to happen in that case?
    }
    */

}

void FlythroughNavigation::mouseReleaseEvent(tgt::MouseEvent* /*e*/) {

    /*
    if (!trackball_ || !trackball_->getCanvas() || !trackball_->getCanvas()->getPainter())
        return;

    if (trackballEnabled_) {
            trackball_->mouseReleaseEvent(e);
            if (!trackball_->getContinuousSpin() || !spinit_)
                postMessage( new BoolMsg(VoreenPainter::switchCoarseness_, false) );

            moveCounter_ = -1;
            trackball_->getCanvas()->update();
    }
    */

}

void FlythroughNavigation::mouseMoveEvent(tgt::MouseEvent* /*e*/) {

    /*
    if (!trackball_ || !trackball_->getCanvas() || !trackball_->getCanvas()->getPainter())
        return;

    if (trackballEnabled_) {
        trackball_->mouseMoveEvent(e);

        // restrict distance within specified range
        if (trackball_->getCenterDistance() < minDistance_)
            trackball_->zoomAbsolute(minDistance_);
        if (trackball_->getCenterDistance() > maxDistance_)
            trackball_->zoomAbsolute(maxDistance_);

        postMessage( new CameraPtrMsg(VoreenPainter::cameraChanged_, trackball_->getCamera()) );
        trackball_->getCanvas()->update();

         if (trackball_->getContinuousSpin())
             moveCounter_ = 0;
    }
    */
}


void FlythroughNavigation::mouseDoubleClickEvent(tgt::MouseEvent* e) {

    // reset camera on right dbl-click
    if (e->button() == tgt::MouseEvent::MOUSE_BUTTON_RIGHT) {
        if (camera_ && initialCameraValid_) {
            camera_->setPosition(initialPosition_);
            camera_->setFocus(initialFocusPoint_);
            camera_->setUpVector(initialUpVector_);
            canvas_->repaint();
        }
    }
}

void FlythroughNavigation::wheelEvent(tgt::MouseEvent* e) {

    /*
    if (!trackball_ || !trackball_->getCanvas() || !trackball_->getCanvas()->getPainter())
        return;
    */

    float angle = 0.01f;

    bool up = (e->button() & MouseEvent::MOUSE_WHEEL_UP);
    bool cameraModified = false;
    if (e->modifiers() == 0) {
        if (up)
            Navigation::rollCameraHorz(angle);
        else
            Navigation::rollCameraHorz(-angle);
        cameraModified = true;
    }
    else if (e->modifiers() == Event::SHIFT) {
        if (up)
            Navigation::rollCameraVert(-angle);
        else
            Navigation::rollCameraVert(angle);
        cameraModified = true;
    }
    else if (e->modifiers() == Event::CTRL) {
        if (up)
            Navigation::rotateViewHorz(angle);
        else
            Navigation::rotateViewHorz(-angle);
        cameraModified = true;
    }

    /* not necessary, the same as rotateViewVert
    else if (e->modifiers() == Event::ALT) {
        if (up)
            Navigation::rollCameraVert(angle);
        else
            Navigation::rollCameraVert(-angle);
        cameraModified = true;
    } */

    if (cameraModified) {
        //MsgDistr.postMessage( new BoolMsg(VoreenPainter::switchCoarseness_, true) );
        wheelCounter_ = 0;
        canvas_->update();
    }

}

void FlythroughNavigation::timerEvent(tgt::TimeEvent* /*e*/) {

    /*
    if (!trackball_ || !trackball_->getCanvas() || !trackball_->getCanvas()->getPainter())
        return;
    */

    if (wheelCounter_ >= 0)
        wheelCounter_++;
    if (moveCounter_ >= 0)
        moveCounter_++;

    if (wheelCounter_ == 60) {
        //MsgDistr.postMessage( new BoolMsg(VoreenPainter::switchCoarseness_, false) );
        canvas_->update();
        wheelCounter_ = -1;
    }
    if (moveCounter_ == 60) {
        //MsgDistr.postMessage( new BoolMsg(VoreenPainter::switchCoarseness_, false) );
        canvas_->update();
        moveCounter_ = -1;
    }


}

} // namespace
