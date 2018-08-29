/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/navigation/trackball.h"

#include <cmath>

#include "tgt/assert.h"

namespace tgt {

Trackball::Trackball(GLCanvas* canvas, Camera* camera, bool defaultEventHandling, Timer* continuousSpinTimer)
    : Navigation(camera),
      canvas_(canvas),
      continuousSpin_(false),
      continuousSpinTimer_(continuousSpinTimer),
      continuousSpinStopwatch_(0),
      moveCenter_(true),
      size_(1.f),
      tracking_(false)
{

    if (camera_)
        center_ = camera_->getFocus();
    else
        center_ = vec3(0.f);

    saveCameraParameters();

    initializeEventHandling();

    if (defaultEventHandling) {
        setMouseRotate();
        setMouseMove();
        setMouseZoom();
        setMouseRoll();
        setMouseWheelZoom();
        setMouseWheelRoll();
    }

    if (continuousSpinTimer_)
        initializeContinuousSpin();
}

Trackball::~Trackball() {
    delete continuousSpinStopwatch_;
}

void Trackball::reset() {
    camera_->positionCamera(cameraPosition_, cameraFocus_, cameraUpVector_);
}

/*
 * Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
 * if we are away from the center of the sphere.
 */
vec3 Trackball::projectToSphere(const vec2 xy) const {

    static const float sqrt2 = sqrtf(2.f);
    vec3 result;

    float d = length(xy);
    if (d < size_ * sqrt2 / 2.f) {    /* Inside sphere */
        // The factor "sqrt2/2.f" make a smooth changeover from sphere to hyperbola. If we leave
        // factor 1/sqrt(2) away, the trackball would bounce at the changeover.
        result.z = sqrtf(size_ * size_ - d*d);
    } else {                         /* On hyperbola */
        float t = size_ / sqrt2;
        result.z = t*t / d;
    }

    result.x = xy.x;
    result.y = xy.y;
    return normalize(result);
}

vec3 Trackball::coordTransform(vec3 const axis) const {
    mat4 rotation;
    camera_->getRotateMatrix().invert(rotation);
    return rotation * axis;
}

void Trackball::startMouseDrag(MouseEvent* e) {
    mouse_ = scaleMouse( ivec2(e->x(), e->y()) );

    if (continuousSpin_ && continuousSpinTimer_) {
        if ( mouseRotateButton_ & e->button() &&
             (mouseRotateMod_ == e->modifiers() || mouseRotateMod_ & e->modifiers()) ) {
            // keep track how long last roation took
            continuousSpinTimer_->stop();
            continuousSpinStopwatch_->reset();
            continuousSpinLastOrientationChangeMSecs_ = 0;
            continuousSpinStopwatch_->start();
        }
    }

    tracking_ = true;
}

void Trackball::endMouseDrag(MouseEvent* e) {
    if (continuousSpin_ && continuousSpinStopwatch_) {
        // Make sure this endMouseDrag belongs to a rotation
        if ((mouseRotateButton_ & e->button()) &&
            (mouseRotateMod_ == e->modifiers() || mouseRotateMod_ & e->modifiers())) {
            // Make sure mouse has been moved within the last time.
            if (continuousSpinStopwatch_->getRuntime() < 2.*continuousSpinLastOrientationChangeMSecs_) {
                if (continuousSpinLastOrientationChangeMSecs_ < 1000.f / 25.f) {
                    // Last rotation took less than 40 msecs. As all rotations with more than
                    // 25 steps per second (less than 40 msecs) will look smooth, we do not
                    // want to rotate with a frequency that high. We will scale the last
                    // orientation change to be so large that rotation with 25 Hz will make
                    // the trackball rotate with the according speed.
                    float scale = (1000.f / 25.f) / static_cast<float>(continuousSpinLastOrientationChangeMSecs_);
//                    lastOrientationChange_.w = cosf( scale * acos( lastOrientationChange_.w ) );
                    lastOrientationChange_ = expQuat( scale * logQuat( lastOrientationChange_ ) );
                    lastOrientationChange_.normalize();
                    continuousSpinTimer_->start(1000 / 25, 0);
                } else {
                    continuousSpinTimer_->start(continuousSpinLastOrientationChangeMSecs_, 0);
                }
            }
            if (1 - lastOrientationChange_.w < 0.000002) {
                // Rotation is very small. Probably user did not mean to spin the trackball.
                lastOrientationChange_.x = 0.f;
                lastOrientationChange_.y = 0.f;
                lastOrientationChange_.z = 0.f;
                lastOrientationChange_.w = 1.f;
                continuousSpinTimer_->stop();
            }
        }
        continuousSpinStopwatch_->stop();
        continuousSpinStopwatch_->reset();
        continuousSpinLastOrientationChangeMSecs_ = 0;
    }

    tracking_ = false;
}

void Trackball::rotate(Quaternion<float> quat) {
    vec3 position = camera_->getPosition();
    position -= center_;
    position = quat::rotate(position, quat);
    position += center_;

    vec3 focus = camera_->getFocus();
    focus -= center_;
    // Usually focus - center == 0, so no need to rotate. But if we combine trackball
    // with some other navigations, this might be useful.
    focus = quat::rotate(focus, quat);
    focus += center_;

    vec3 upVector = camera_->getUpVector();
    upVector = quat::rotate(upVector, quat);

    camera_->positionCamera(position, focus, upVector);

    lastOrientationChange_ = quat;
    if (continuousSpin_ && continuousSpinStopwatch_) {
        // keep track how long last rotation took
        continuousSpinLastOrientationChangeMSecs_ = continuousSpinStopwatch_->getRuntime();
        continuousSpinStopwatch_->reset();
        continuousSpinStopwatch_->start();
    }

    getCanvas()->update();
}

void Trackball::rotate(vec3 axis, float phi) {
    // use coordTransform to get axis in world coordinates according to the axis given in camera coordinates
    rotate(quat::createQuat(phi, coordTransform(axis)));
}

void Trackball::rotate(vec2 newMouse) {

/* Project the points onto the virtual trackball,
 * then figure out the axis of rotation, which is the cross
 * product of P1-P2 and O-P1 (O is the center of the ball, 0,0,0)
 * Note:  This is a deformed trackball -- it is a trackball in the center,
 * but is deformed into a hyperbolic sheet of rotation away from the
 * center (projectToSphere does that job).  This particular function was
 * chosen after trying out several variations.
 *
 * It is assumed that the arguments to this routine are in the range
 * (-1.0 ... 1.0) */

    if (newMouse == mouse_) {
        /* Zero rotation -> do nothing */
        return;
    }

    /* First, figure out z-coordinates for projection of P1 and P2 to deformed sphere */
    vec3 p1 = projectToSphere(mouse_);
    vec3 p2 = projectToSphere(newMouse);

    /* Now, find the axis we are going to rotate about*/
    vec3 axis = cross(p2, p1);

    /*... and calculate the angle phi between the two vectors which is the
      angle we need to rotate about*/
    vec3 d = p1 - p2;
    float t = length(d) / (2.0f * getSize());

    /* avoid problems with out-of-control values... */
    if (t > 1.0f) t = 1.0f;
    if (t < -1.0f) t = -1.0f;

    float phi = 2.f * asinf(t);

    rotate(axis, phi);
}

void Trackball::move(float length, vec3 axis) {
    float frustFactor = camera_->getFocalLength() / camera_->getFrustum().getNearDist();
    float frustWidth  = ( camera_->getFrustum().getRight()
                          - camera_->getFrustum().getLeft() );
    float frustHeight = ( camera_->getFrustum().getTop()
                          - camera_->getFrustum().getBottom() );
    axis = normalize(axis) * length * frustFactor;
    axis.x *= frustWidth;
    axis.y *= frustHeight;

    // find axis in world coordinates according to the axis given in camera coordinates
    axis = coordTransform(axis);

    moveCamera(-axis);
    if (moveCenter_) { center_ -= axis; };
    getCanvas()->update();
}

void Trackball::move(vec2 newMouse) {
    vec2 mouseMotion = newMouse - mouse_;
    vec3 axis;
    axis.x = mouseMotion.x;
    axis.y = mouseMotion.y;
    axis.z = 0;

    move(length(mouseMotion), axis);
}

void Trackball::zoom(float factor) {
    // zoom factor is inverse proportional to scaling of the look vector, so invert:
    factor = 1.f / factor;
    camera_->setPosition( (1.f-factor) * camera_->getFocus()
                              + factor * camera_->getPosition());
    getCanvas()->update();
}

void Trackball::zoom(vec2 newMouse) {
    zoom( 1 + dot( (mouse_-newMouse), mouseZoomInDirection_) );
}

void Trackball::zoomAbsolute(float focallength) {
   camera_->setPosition( camera_->getFocus() - focallength * camera_->getLook());
   getCanvas()->update();
}

void Trackball::initializeContinuousSpin() {
    continuousSpinTimer_->stop();
    continuousSpinStopwatch_ = new Stopwatch();
    continuousSpinLastOrientationChangeMSecs_ = 0;
    setContinuousSpin(true);
}

float Trackball::getCenterDistance() {
    return dot( center_ - camera_->getPosition(), camera_->getLook() );
}

void Trackball::saveCameraParameters() {
    if (camera_) {
        cameraPosition_ = camera_->getPosition();
        cameraFocus_ = camera_->getFocus();
        cameraUpVector_ = camera_->getUpVector();
    }
    else {
        cameraPosition_ = vec3(0.f);
        cameraFocus_ = vec3(0.f, 0.f, -1.f);
        cameraUpVector_ = vec3(0.f, 1.f, 0.f);
    }
}

// Event Handling /////////////////////////////////////////////////////////////////////////////

void Trackball::setMouseRotate(MouseEvent::MouseButtons button, int mod) {
    mouseRotateButton_ = button;
    mouseRotateMod_ = mod;
}

void Trackball::setKeyRotate(float acuteness,
                             KeyEvent::KeyCode left,
                             KeyEvent::KeyCode right,
                             KeyEvent::KeyCode up,
                             KeyEvent::KeyCode down,
                             int mod, bool pressed) {
    if (acuteness != 0.f) {
        keyRotateAcuteness_ = acuteness;
        keyRotateLeft_ = left;
        keyRotateRight_ = right;
        keyRotateUp_ = up;
        keyRotateDown_ = down;
        keyRotateMod_ = mod;
        keyRotatePressed_ = pressed;
    } else {
        keyRotateLeft_ = keyRotateRight_ = keyRotateUp_ = keyRotateDown_ = KeyEvent::K_LAST;
    }
}

void Trackball::setMouseMove(MouseEvent::MouseButtons button, int mod) {
    mouseMoveButton_ = button;
    mouseMoveMod_ = mod;
}

void Trackball::setKeyMove(float acuteness,
                           KeyEvent::KeyCode left,
                           KeyEvent::KeyCode right,
                           KeyEvent::KeyCode up,
                           KeyEvent::KeyCode down,
                           int mod, bool pressed) {
    if (acuteness != 0.f) {
        keyMoveAcuteness_ = acuteness;
        keyMoveLeft_ = left;
        keyMoveRight_ = right;
        keyMoveUp_ = up;
        keyMoveDown_ = down;
        keyMoveMod_ = mod;
        keyMovePressed_ = pressed;
    } else {
        keyMoveLeft_ = keyMoveRight_ = keyMoveUp_ = keyMoveDown_ = KeyEvent::K_LAST;
    }
}

void Trackball::setMouseZoom(MouseEvent::MouseButtons button,
                             vec2 zoomInDirection,
                             int mod) {
    mouseZoomButton_ = button;
    mouseZoomInDirection_ = normalize(zoomInDirection);
    mouseZoomMod_ = mod;
}

void Trackball::setMouseWheelZoom(float acuteness, bool wheelUpZoomIn, int mod) {
    if (acuteness != 0.f) {
        mouseWheelZoom_ = true;
        mouseWheelZoomAcuteness_ = acuteness;
        mouseWheelUpZoomIn_ = wheelUpZoomIn;
        mouseWheelZoomMod_ = mod;
    } else {
        mouseWheelZoom_ = false;
    }
}

void Trackball::setKeyZoom(float acuteness, KeyEvent::KeyCode in, KeyEvent::KeyCode out,
                           int mod, bool pressed) {
    if (acuteness != 0.f) {
        keyZoomAcuteness_ = acuteness;
        keyZoomIn_ = in;
        keyZoomOut_ = out;
        keyZoomMod_ = mod;
        keyZoomPressed_ = pressed;
    } else {
        keyZoomIn_ = keyZoomOut_ = KeyEvent::K_LAST;
    }
}
void Trackball::setMouseRoll(MouseEvent::MouseButtons button, float acuteness, int mod) {
    mouseRollButton_ = button;
    mouseRollAcuteness_ = acuteness;
    mouseRollMod_ = mod;
}

void Trackball::setMouseWheelRoll(float acuteness, bool wheelUpRollLeft, int mod) {
    if (acuteness != 0.f) {
        mouseWheelRoll_ = true;
        mouseWheelRollAcuteness_ = acuteness;
        mouseWheelUpRollLeft_ = wheelUpRollLeft;
        mouseWheelRollMod_ = mod;
    } else {
        mouseWheelRoll_ = false;
    }
}

void Trackball::setKeyRoll(float acuteness, KeyEvent::KeyCode left, KeyEvent::KeyCode right,
                           int mod, bool pressed) {
    if (acuteness != 0.f) {
        keyRollAcuteness_ = acuteness;
        keyRollLeft_ = left;
        keyRollRight_ = right;
        keyRollMod_ = mod;
        keyRollPressed_ = pressed;
    } else {
        keyRollLeft_ = keyRollRight_ = KeyEvent::K_LAST;
    }
}

vec2 Trackball::scaleMouse(const ivec2& mouse) const {
    return vec2( static_cast<float>(mouse.x*2.f) / static_cast<float>(canvas_->getWidth()) - 1.f,
                 1.f - static_cast<float>(mouse.y*2.f) / static_cast<float>(canvas_->getHeight()) );
}

void Trackball::mousePressEvent(MouseEvent* e) {
    startMouseDrag( e );
    e->ignore();
}

void Trackball::mouseReleaseEvent(MouseEvent* e) {
    endMouseDrag( e );
    e->ignore();
}

void Trackball::mouseMoveEvent(MouseEvent* e) {
    if (!tracking_)
        return;

    vec2 newMouse = scaleMouse( ivec2(e->x(), e->y()) );

    if ( mouseRotateButton_ & e->button() &&
         (mouseRotateMod_ == e->modifiers() || mouseRotateMod_ & e->modifiers()) ) {
        rotate(newMouse);
    }
    if ( mouseMoveButton_ & e->button() &&
         (mouseMoveMod_ == e->modifiers() || mouseMoveMod_ & e->modifiers()) ) {
        move(newMouse);
    }
    if ( mouseZoomButton_ & e->button() &&
         (mouseZoomMod_ == e->modifiers() || mouseZoomMod_ & e->modifiers()) ) {
        zoom(newMouse);
    }
    if ( mouseRollButton_ & e->button() &&
         (mouseRollMod_ == e->modifiers() || mouseRollMod_ & e->modifiers()) ) {
        rollCameraHorz((newMouse.x-mouse_.x)/mouseRollAcuteness_);
        getCanvas()->update();
    }

    mouse_ = newMouse;

    e->ignore();
}

void Trackball::wheelEvent(MouseEvent* e) {
    if ( mouseWheelZoom_ &&
         (mouseWheelZoomMod_ == e->modifiers() || mouseWheelZoomMod_ & e->modifiers()) ) {
        zoom( getZoomFactor( mouseWheelZoomAcuteness_,
                             ( mouseWheelUpZoomIn_ && (e->button()&MouseEvent::MOUSE_WHEEL_UP)) ||
                             (!mouseWheelUpZoomIn_ && (e->button()&MouseEvent::MOUSE_WHEEL_DOWN))
                  ) );
    } else if ( mouseWheelRoll_ &&
                (mouseWheelRollMod_ == e->modifiers() || mouseWheelRollMod_ & e->modifiers()) ) {
        rollCameraHorz( getRollAngle( mouseWheelRollAcuteness_,
                             ( mouseWheelUpRollLeft_ && (e->button()&MouseEvent::MOUSE_WHEEL_UP)) ||
                             (!mouseWheelUpRollLeft_ && (e->button()&MouseEvent::MOUSE_WHEEL_DOWN))
                  ) );
        getCanvas()->update();
    }
    e->ignore();
}

void Trackball::keyEvent(KeyEvent* e) {
    if ((keyRotateMod_ == e->modifiers() || keyRotateMod_ & e->modifiers())
        && e->pressed() == keyRotatePressed_) {
        if (e->keyCode() == keyRotateLeft_) {
            rotate(vec3(0.f, 1.f, 0.f),  getRotationAngle(keyRotateAcuteness_));
        } else if (e->keyCode() == keyRotateRight_) {
            rotate(vec3(0.f, 1.f, 0.f), -getRotationAngle(keyRotateAcuteness_));
        } else if (e->keyCode() == keyRotateUp_) {
            rotate(vec3(1.f, 0.f, 0.f), -getRotationAngle(keyRotateAcuteness_));
        } else if (e->keyCode() == keyRotateDown_) {
            rotate(vec3(1.f, 0.f, 0.f),  getRotationAngle(keyRotateAcuteness_));
        }
    }
    if ((keyMoveMod_ == e->modifiers() || keyMoveMod_ & e->modifiers())
        && e->pressed() == keyMovePressed_) {
        if (e->keyCode() == keyMoveLeft_) {
            move(getMovementLength(keyMoveAcuteness_), vec3(-1.f, 0.f, 0.f));
        } else if (e->keyCode() == keyMoveRight_) {
            move(getMovementLength(keyMoveAcuteness_), vec3( 1.f, 0.f, 0.f));
        } else if (e->keyCode() == keyMoveUp_) {
            move(getMovementLength(keyMoveAcuteness_), vec3(0.f,  1.f, 0.f));
        } else if (e->keyCode() == keyMoveDown_) {
            move(getMovementLength(keyMoveAcuteness_), vec3(0.f, -1.f, 0.f));
        }
    }
    if ((keyZoomMod_ == e->modifiers() || keyZoomMod_ & e->modifiers())
        && e->pressed() == keyZoomPressed_) {
        if (e->keyCode() == keyZoomIn_) {
            zoom(getZoomFactor(keyZoomAcuteness_, true));
        } else if (e->keyCode() == keyZoomOut_) {
            zoom(getZoomFactor(keyZoomAcuteness_, false));
        }
    }
    if ((keyRollMod_ == e->modifiers() || keyRollMod_ & e->modifiers())
        && e->pressed() == keyRollPressed_) {
        if (e->keyCode() == keyRollLeft_) {
            rollCameraHorz(getRollAngle(keyRollAcuteness_, true));
            getCanvas()->update();
        } else if (e->keyCode() == keyRollRight_) {
            rollCameraHorz(getRollAngle(keyRollAcuteness_, false));
            getCanvas()->update();
        }
    }
    e->ignore();
}

float Trackball::getRotationAngle(const float& acuteness) const {
    return 1.f / acuteness;
}

float Trackball::getMovementLength(const float& acuteness) const {
    return getSize() / acuteness;
}

float Trackball::getZoomFactor(const float& acuteness, const bool& zoomIn) const {
    if (zoomIn) {
        return 1.f + 1.f/acuteness;
    } else {
        return 1.f - 1.f/acuteness;
    }
}

float Trackball::getRollAngle(const float& acuteness, const bool& left) const {
    if (left) {
        return - 1.f / acuteness;
    } else {
        return 1.f / acuteness;
    }
}

void Trackball::timerEvent(TimeEvent* e) {
    if (continuousSpin_ && continuousSpinTimer_ == e->getTimer() ) {
        rotate(getLastOrientationChange());
        canvas_->repaint();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////

} // namespace
