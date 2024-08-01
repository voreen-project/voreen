/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
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

#include "tgt/navigation/navigation.h"

//#include <cmath>

// #include "tgt/assert.h"
// #include "tgt/glmath.h"
// #include "tgt/quaternion.h"
// #include "tgt/spline.h"
// #include "tgt/tgt_gl.h"


namespace tgt {

// Rotate the Camera about an abitrary axis and angle given by the user
void Navigation::rotateView(float angle, const vec3& axis) {

    // for now, we use quaternion rotation
    vec3 look = normalize( quat::rotate(camera_->getLook(), angle, axis) );

    // we do not allow views that are too similar to the up-vector since this results
    // in "flickering"
    if ( fabs(dot( look, camera_->getUpVector() )) > .995f )
        return;

    // Set the new focus-point
    camera_->setFocus(camera_->getPosition() + camera_->getFocalLength()*look);
}

// Rotate the Camera about horizontal and vertical angles given by the user
void Navigation::rotateViewHV(float anglehorz, float anglevert) {

// FIXME: rotateViewHorz(anglehorz);rotateViewVert(anglevert) would be much more generic ... but maybe more expensive ...

    // for now, we use quaternion rotation... we'll have to benchmark this
    quat horz = quat::createQuat(anglehorz, camera_->getUpVector());
    quat vert = quat::createQuat(anglevert, camera_->getStrafe());

    // this way, we save calculations by multiplying the quaternions representing
    // the horizontal and vertical rotations instead of changing the view twice, once
    // for each angle and axis
    vec3 look = normalize( quat::rotate(camera_->getLook(), vert * horz) );

    // we do not allow views that are too similar to the up-vector since this results
    // in "flickering"
    if (fabs(dot(look, camera_->getUpVector())) > .995f)
        return;

    // Set the new focus-point
    camera_->setFocus(camera_->getPosition() + camera_->getFocalLength()*look);
}

// FIXME: the following few functions should be inline

void Navigation::rotateView(float angle, float x, float y, float z) {
    rotateView(angle, vec3(x, y, z));
}

// This is just for convenience since the rotation about the strafe-, the look-
// and the up-Vector are used in a lot of cases.
void Navigation::rotateViewVert(float angle) {
    rotateView(angle, camera_->getStrafe());
}

void Navigation::rotateViewHorz(float angle) {
    rotateView(angle, camera_->getUpVector());
}

// This can be used to rotate (or "roll") the camera left and right, angle in radian measure
void Navigation::rollCameraHorz(float angle) {
    vec3 up = normalize( quat::rotate(camera_->getUpVector(), angle, camera_->getLook()) );
    camera_->setUpVector(up);
}

// This can be used to rotate (or "roll") the camera forwards and backwards, angle in radian measure
void Navigation::rollCameraVert(float angle) {
    vec3 up = normalize( quat::rotate(camera_->getUpVector(), angle, camera_->getStrafe()) );
    camera_->setUpVector(up);
    vec3 look = cross(up, camera_->getStrafe());
    camera_->setFocus(camera_->getPosition() + look);
}

void Navigation::moveCameraForward(float length) {
    vec3 motionvector = length*camera_->getLook();
    camera_->setPosition( camera_->getPosition() + motionvector );
    camera_->setFocus( camera_->getFocus() + motionvector );
}

void Navigation::moveCameraBackward(float length) {
    moveCameraForward( -length );
}

void Navigation::moveCameraUp(float length) {
    vec3 motionvector = length*camera_->getUpVector();
    camera_->setPosition( camera_->getPosition() + motionvector );
    camera_->setFocus( camera_->getFocus() + motionvector );
}

void Navigation::moveCameraDown(float length) {
    moveCameraUp( -length );
}

void Navigation::moveCameraRight(float length) {
    vec3 motionvector = length*camera_->getStrafe();
    camera_->setPosition( camera_->getPosition() + motionvector );
    camera_->setFocus( camera_->getFocus() + motionvector );
}

void Navigation::moveCameraLeft(float length) {
    moveCameraRight( -length );
}

void Navigation::moveCamera(float length, float x, float y, float z) {
    moveCamera(length, vec3(x, y, z));
}

void Navigation::moveCamera(float length, const vec3& axis) {
    vec3 motionvector = length * normalize(axis);
    moveCamera(motionvector);
}

void Navigation::moveCamera(const vec3& motionvector) {
    camera_->setPosition( camera_->getPosition() + motionvector );
    camera_->setFocus( camera_->getFocus() + motionvector );
}

/***********************************************************************

    /// This can be used to look around the scene with the help of a mouse.  Middle is used
    /// as a reference value to calculate the angle by which the mouse has moved horizontally and vertically.
    void setViewByMouse(Mouse* mouse, const ivec2& middle, float speedfactor);


// Used to look around the scene with a mouse
void Navigation::setViewByMouse(Mouse* mouse, const ivec2& middle, float speedfactor) {

    float angleY;
    float angleZ;

    ivec2 current = mouse->getPosition();

    if ( (current.x == middle.x) && (current.y == middle.y) )
        return; //no change

    mouse->setPosition(middle);

    angleY = ((middle.x - current.x) * speedfactor) / (1000.0f * mouse->getRate());
    angleZ = ((middle.y - current.y) * speedfactor) / (1000.0f * mouse->getRate());

    rotateViewHV(angleY, angleZ);
}


    /// This method lets the camera follow a spline; the orientation and the upVector are also interpolated
    /// along given splines.
    void followSpline(const Spline& posSp, const Spline& focusSp, const Spline& upSp,
                      float t, float total = 1.f);


// This function takes the camera for a ride.  The camera's parameters are changed
// so the the view-vector is always looking in the direction of the focus-spline;
// the new up vector is also interpolated from the upSpline.

void Navigation::followSpline(const Spline& posSp, const Spline& focusSp, const Spline& upSp,
                      float t, float total) {
    // with less than 4 control-Points, we don't have engough information to form a B-Spline;
    // also, the running variable must not be smaller than zero or greater than the total-run-variable,
    // and all Splines that make up the Camera-run should have the same size.
    uint numberOfPoints = posSp.getControlPoints().size();

    if (numberOfPoints < 4 || (focusSp.getControlPoints().size() != numberOfPoints)
                           ||    (upSp.getControlPoints().size() != numberOfPoints)
                           || t < 0.f || t > total)
        return;

    vec3 pos = posSp.interpolate(t, total);

    vec3 focus_pos = focusSp.interpolate(t, total);

    vec3 upVector_pos = upSp.interpolate(t, total);

    setPosition(pos);
    setFocus(focus_pos);
    setUpVector(upVector_pos);
}

************************************************************************/

}; // namespace
