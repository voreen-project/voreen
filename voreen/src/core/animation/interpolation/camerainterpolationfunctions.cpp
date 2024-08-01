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

#include "voreen/core/animation/interpolation/camerainterpolationfunctions.h"
#include "voreen/core/animation/interpolation/basicfloatinterpolation.h"

#include "voreen/core/animation/interpolation/vec3interpolationfunctions.h"
#include "tgt/quaternion.h"
#include <cmath>
#include "voreen/core/animation/propertykeyvalue.h"

using tgt::vec3;
using tgt::Camera;

namespace voreen {

CameraLinearInterpolationFunction::CameraLinearInterpolationFunction() {}

std::string CameraLinearInterpolationFunction::getGuiName() const {
    return "linear interpolation";
}

std::string CameraLinearInterpolationFunction::getCategory() const {
    return "linear";
}

Camera CameraLinearInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    return interpolateInternal(startvalue, endvalue, time);
}

Camera CameraLinearInterpolationFunction::interpolateInternal(Camera startvalue, Camera endvalue, float time) {
    Vec3LinearInterpolationFunction* intfunc = new Vec3LinearInterpolationFunction();
    Vec3SphericalLinearInterpolationFunction* intfunc2 = new Vec3SphericalLinearInterpolationFunction();
    vec3 posvec = intfunc->interpolate(startvalue.getPosition(), endvalue.getPosition(), time);
    vec3 focvec = intfunc->interpolate(startvalue.getFocus(), endvalue.getFocus(), time);
    vec3 upvec = normalize(intfunc2->interpolate(startvalue.getUpVector(), endvalue.getUpVector(), time));
    //vec3 direction = normalize(endvalue.getPosition() - startvalue.getPosition());
    Camera node = Camera(startvalue);
    node.positionCamera(posvec, focvec, upvec);
/*    Camera node = new Camera(posvec, focvec, upvec, direction);
    if (startvalue->isTangential() && endvalue->isTangential())
        node->setTangential(true); */
    return node;
}

InterpolationFunction<Camera>* CameraLinearInterpolationFunction::create() const {
    return new CameraLinearInterpolationFunction();
}

//-----------------------------------------------------------------------------

CameraSmoothLinearInterpolationFunction::CameraSmoothLinearInterpolationFunction() {}

std::string CameraSmoothLinearInterpolationFunction::getGuiName() const {
    return "smooth linear interpolation";
}

std::string CameraSmoothLinearInterpolationFunction::getCategory() const {
    return "linear";
}

Camera CameraSmoothLinearInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    float smoothTime = BasicFloatInterpolation::inOutQuadInterpolation(0.0f, 1.0f, time);
    return CameraLinearInterpolationFunction::interpolateInternal(startvalue, endvalue, smoothTime);
}

InterpolationFunction<Camera>* CameraSmoothLinearInterpolationFunction::create() const {
    return new CameraSmoothLinearInterpolationFunction();
}

//-----------------------------------------------------------------------------

CameraRightRotationInterpolationFunction::CameraRightRotationInterpolationFunction() {}

std::string CameraRightRotationInterpolationFunction::getGuiName() const {
    return "rotation (right)";
}

std::string CameraRightRotationInterpolationFunction::getCategory() const {
    return "rotation";
}

Camera CameraRightRotationInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    return interpolateInternal(startvalue, endvalue, time);
}

Camera CameraRightRotationInterpolationFunction::interpolateInternal(Camera startvalue, Camera endvalue, float time) {
    vec3 focus = startvalue.getFocus();
    if(length(focus - endvalue.getFocus()) > 0.01f) {
        LERROR("focus mismatch!");
        return Camera();
    }

    vec3 up = startvalue.getUpVector();
    if(length(up - endvalue.getUpVector()) > 0.01f) {
        LERROR("upvector mismatch!");
        return Camera();
    }

    float r1 = length(focus - startvalue.getPosition());
    float r2 = length(focus - endvalue.getPosition());

    vec3 xVec = normalize(startvalue.getPosition() - focus);
    vec3 yVec = normalize(cross(xVec, up));

    float maxAngle = 2.0f * tgt::PIf;
    if(length(startvalue.getPosition() - endvalue.getPosition()) > 0.0001f) {
        bool smallAngle = dot(yVec, endvalue.getPosition() - focus) >= 0.0f;
        vec3 toEnd = normalize(endvalue.getPosition() - focus);

        if(smallAngle)
            maxAngle = acosf(dot(xVec, toEnd));
        else
            maxAngle = tgt::PIf + acosf(dot(-xVec, toEnd));
    }

    float angle = time * maxAngle;

    vec3 pos = focus + (sinf(angle) * yVec + cosf(angle) * xVec) * (time * r2 + (1.0f - time) * r1);

    Camera node = Camera(startvalue);
    node.positionCamera(pos, focus, up);
    return node;
}

InterpolationFunction<Camera>* CameraRightRotationInterpolationFunction::create() const {
        return new CameraRightRotationInterpolationFunction();
}

//-----------------------------------------------------------------------------

CameraLeftRotationInterpolationFunction::CameraLeftRotationInterpolationFunction() {}

std::string CameraLeftRotationInterpolationFunction::getGuiName() const {
    return "rotation (left)";
}

std::string CameraLeftRotationInterpolationFunction::getCategory() const {
    return "rotation";
}

Camera CameraLeftRotationInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    return interpolateInternal(startvalue, endvalue, time);
}

Camera CameraLeftRotationInterpolationFunction::interpolateInternal(Camera startvalue, Camera endvalue, float time) {
    vec3 focus = startvalue.getFocus();
    if(length(focus - endvalue.getFocus()) > 0.01f) {
        LERROR("focus mismatch!");
        return Camera();
    }

    vec3 up = startvalue.getUpVector();
    if(length(up - endvalue.getUpVector()) > 0.01f) {
        LERROR("upvector mismatch!");
        return Camera();
    }

    float r1 = length(focus - startvalue.getPosition());
    float r2 = length(focus - endvalue.getPosition());

    vec3 xVec = normalize(startvalue.getPosition() - focus);
    vec3 yVec = normalize(cross(xVec, up));

    float maxAngle = 2.0f * tgt::PIf;
    if(length(startvalue.getPosition() - endvalue.getPosition()) > 0.0001f) {
        bool smallAngle = dot(yVec, endvalue.getPosition() - focus) <= 0.0f;
        vec3 toEnd = normalize(endvalue.getPosition() - focus);

        if(smallAngle)
            maxAngle = acosf(dot(xVec, toEnd));
        else
            maxAngle = tgt::PIf + acosf(dot(-xVec, toEnd));
    }

    float angle = time * -maxAngle;

    vec3 pos = focus + (sinf(angle) * yVec + cosf(angle) * xVec) * (time * r2 + (1.0f - time) * r1);

    Camera node = Camera(startvalue);
    node.positionCamera(pos, focus, up);
    return node;
}

InterpolationFunction<Camera>* CameraLeftRotationInterpolationFunction::create() const {
        return new CameraLeftRotationInterpolationFunction();
}

//-----------------------------------------------------------------------------

CameraSmoothLeftRotationInterpolationFunction::CameraSmoothLeftRotationInterpolationFunction() {}

std::string CameraSmoothLeftRotationInterpolationFunction::getGuiName() const {
    return "smooth rotation (left)";
}

std::string CameraSmoothLeftRotationInterpolationFunction::getCategory() const {
    return "rotation";
}

Camera CameraSmoothLeftRotationInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    float smoothTime = BasicFloatInterpolation::inOutQuadInterpolation(0.0f, 1.0f, time);
    return CameraLeftRotationInterpolationFunction::interpolateInternal(startvalue, endvalue, smoothTime);
}

InterpolationFunction<Camera>* CameraSmoothLeftRotationInterpolationFunction::create() const {
    return new CameraSmoothLeftRotationInterpolationFunction();
}

//-----------------------------------------------------------------------------

CameraSmoothRightRotationInterpolationFunction::CameraSmoothRightRotationInterpolationFunction() {}

std::string CameraSmoothRightRotationInterpolationFunction::getGuiName() const {
    return "smooth rotation (right)";
}

std::string CameraSmoothRightRotationInterpolationFunction::getCategory() const {
    return "rotation";
}

Camera CameraSmoothRightRotationInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    float smoothTime = BasicFloatInterpolation::inOutQuadInterpolation(0.0f, 1.0f, time);
    return CameraRightRotationInterpolationFunction::interpolateInternal(startvalue, endvalue, smoothTime);
}

InterpolationFunction<Camera>* CameraSmoothRightRotationInterpolationFunction::create() const {
    return new CameraSmoothRightRotationInterpolationFunction();
}

//-----------------------------------------------------------------------------

CameraSphericalLinearInterpolationFunction::CameraSphericalLinearInterpolationFunction() {}

std::string CameraSphericalLinearInterpolationFunction::getGuiName() const {
    return "spherical linear interpolation";
}

std::string CameraSphericalLinearInterpolationFunction::getCategory() const {
    return "linear";
}

Camera CameraSphericalLinearInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    Vec3SphericalLinearInterpolationFunction* sphericalintfunc = new Vec3SphericalLinearInterpolationFunction();

    vec3 posvec = sphericalintfunc->interpolate(startvalue.getPosition(), endvalue.getPosition(), time);
    vec3 focvec = sphericalintfunc->interpolate(startvalue.getFocus(), endvalue.getFocus(), time);
    vec3 upvec  = normalize(sphericalintfunc->interpolate(startvalue.getUpVector(), endvalue.getUpVector(), time));

    vec3 direction;
//    float eps = 0.0001f;
    float dotVal = dot(normalize(startvalue.getPosition()), normalize(endvalue.getPosition()));
    if (dotVal > 0.9995){
        // small angle => linear interpolation
        direction = (endvalue.getPosition() - startvalue.getPosition());
    }
    Camera node(startvalue);
    node.positionCamera(posvec, focvec, upvec);
/*    float theta = acosf(dotVal);
    direction = (-1.0f * cosf((1.0f - time) * theta)/(sinf(theta) + eps)) * startvalue->getPosition()
                + (cosf(time * theta)/(sinf(theta) + eps)) * endvalue->getPosition();
    Camera node = new Camera(posvec, focvec, upvec, direction);
    if (startvalue->isTangential() && endvalue->isTangential())
        node->setTangential(true); */
    return node;
}

InterpolationFunction<Camera>* CameraSphericalLinearInterpolationFunction::create() const {
    return new CameraSphericalLinearInterpolationFunction();
}

CameraCubicSplineInterpolationFunction::CameraCubicSplineInterpolationFunction() {}

std::string CameraCubicSplineInterpolationFunction::getGuiName() const{
    return "cubic Bezier spline";
}

std::string CameraCubicSplineInterpolationFunction::getCategory() const{
    return "spline";
}

Camera CameraCubicSplineInterpolationFunction::interpolate(std::vector<PropertyKeyValue<Camera>*> controlpoints, float time) const{
    Camera startvalue;
    Camera endvalue;
    std::vector<PropertyKeyValue<Camera>*>::iterator it;
    it = controlpoints.begin();
    while ((*it)->getTime() < time)
        it++;

    float time2 = (*it)->getTime();
    tgt::vec3 f2 = (*it)->getValue().getFocus();
    endvalue = (*it)->getValue();
    it--;
    float time1 = (*it)->getTime();
    tgt::vec3 f1 = (*it)->getValue().getFocus();
    startvalue = (*it)->getValue();
    tgt::vec3 fm1;
    if (it!= controlpoints.begin()) {
        it--;
        //tgt::vec3 p0 = (*it)->getValue().getPosition();
        tgt::vec3 f0 = (*it)->getValue().getFocus();
        fm1 = 0.5f*(f2 - f0);
        it++;
    } else
        fm1 = f2 - f1;

    it++;
    it++;
    tgt::vec3 f3;//, p3;
    tgt::vec3 m2, fm2;
    if (it!= controlpoints.end()) {
        f3 = (*it)->getValue().getFocus();
        fm2 = 0.5f*(f3 - f1);
    } else
        fm2 = f2 - f1;

    float t = (time - time1) / (time2 - time1);

    vec3 p0 = startvalue.getPosition();
/*    vec3 p1 = startvalue->getPosition() + startvalue->getTangenteHead();
    vec3 p2 = endvalue->getPosition() - endvalue->getTangenteTail(); */
    vec3 p1 = startvalue.getPosition() + (endvalue.getPosition() - startvalue.getPosition())/3.f;
    vec3 p2 = startvalue.getPosition() + (endvalue.getPosition() - startvalue.getPosition())*2.f/3.f;
//    vec3 p2 = endvalue.getPosition() - endvalue.getTangenteTail();
    vec3 p3 = endvalue.getPosition();

    vec3 position;
    const float b03 = (1.0f - t) * (1.0f - t) * (1.0f - t);
    const float b13 = 3.0f * (1.0f - t) * (1.0f - t) * t;
    const float b23 = 3.0f * (1.0f - t) * t * t;
    const float b33 = t * t * t;
    position = b03 * p0 + b13 * p1 + b23 * p2 + b33 * p3;

    vec3 focus;
    const float h00 = (1.f + 2.f*t)*(1.f - t)*(1.f - t);
    const float h10 = t*(1.f- t)*(1.f - t);
    const float h01 = t*t*(3.f-2.f*t);
    const float h11 = t*t*(t - 1.f);
    focus = h00 * f1 + h10 * fm1 + h01 * f2 + h11 * fm2;

    vec3 upVector;
    Vec3SphericalLinearInterpolationFunction* intfunc2 = new Vec3SphericalLinearInterpolationFunction();
    upVector = normalize(intfunc2->interpolate(startvalue.getUpVector(), endvalue.getUpVector(), t));

    vec3 direction;
    const float db03 = -3.0f*(1.0f - t)*(1.0f - t);
    const float db13 = 1.0f - 4.0f*t + 3.0f*t*t;
    const float db23 = 2.0f*t - 3.0f*t*t;
    const float db33 = 3.0f*t*t;
    direction = db03 * p0 + db13 * p1 + db23 * p2 + db33 * p3;

    Camera node = Camera(startvalue);
    node.positionCamera(position, focus, upVector);
/*    Camera node = new Camera(position, focus, upVector, direction);
    if (startvalue->isTangential() && endvalue->isTangential())
        node->setTangential(true); */
    return node;
}
MultiPointInterpolationFunction<Camera>* CameraCubicSplineInterpolationFunction::create() const {
    return new CameraCubicSplineInterpolationFunction();
}

CameraStartInterpolationFunction::CameraStartInterpolationFunction() {}

std::string CameraStartInterpolationFunction::getGuiName() const {
    return "focus on startvalue";
}

std::string CameraStartInterpolationFunction::getCategory() const {
    return "boolean";
}

Camera CameraStartInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    if (time < 1.f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<Camera>* CameraStartInterpolationFunction::create() const {
    return new CameraStartInterpolationFunction();
}

CameraEndInterpolationFunction::CameraEndInterpolationFunction() {}

std::string CameraEndInterpolationFunction::getGuiName() const {
    return "focus on endvalue";
}
std::string CameraEndInterpolationFunction::getCategory() const {
    return "boolean";
}
Camera CameraEndInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    if (time > 0.f)
        return endvalue;
    else
        return startvalue;
}
InterpolationFunction<Camera>* CameraEndInterpolationFunction::create() const {
    return new CameraEndInterpolationFunction();
}

CameraStartEndInterpolationFunction::CameraStartEndInterpolationFunction() {}

std::string CameraStartEndInterpolationFunction::getGuiName() const {
    return "bisection";
}

std::string CameraStartEndInterpolationFunction::getCategory() const {
    return "boolean";
}

Camera CameraStartEndInterpolationFunction::interpolate(Camera startvalue, Camera endvalue, float time) const {
    if (time < 0.5f)
        return startvalue;
    else
        return endvalue;
}

InterpolationFunction<Camera>* CameraStartEndInterpolationFunction::create() const {
    return new CameraStartEndInterpolationFunction();
}

CameraCatmullRomInterpolationFunction::CameraCatmullRomInterpolationFunction() {}

std::string CameraCatmullRomInterpolationFunction::getGuiName() const {
    return "Catmull-Rom spline";
}

std::string CameraCatmullRomInterpolationFunction::getCategory() const {
    return "spline";
}

Camera CameraCatmullRomInterpolationFunction::interpolate(std::vector<PropertyKeyValue<Camera>*> controlpoints, float time) const {
    Camera camera1;
    Camera camera2;
    std::vector<PropertyKeyValue<Camera>*>::iterator it;
    it = controlpoints.begin();
    while ((*it)->getTime() < time)
        it++;

    float time2 = (*it)->getTime();
    tgt::vec3 p2 = (*it)->getValue().getPosition();
    tgt::vec3 f2 = (*it)->getValue().getFocus();
    camera2 = (*it)->getValue();
    it--;
    float time1 = (*it)->getTime();
    tgt::vec3 p1 = (*it)->getValue().getPosition();
    tgt::vec3 f1 = (*it)->getValue().getFocus();
    camera1 = (*it)->getValue();
    tgt::vec3 m1;
    tgt::vec3 fm1;
    if (it!= controlpoints.begin()) {
        it--;
        tgt::vec3 p0 = (*it)->getValue().getPosition();
        tgt::vec3 f0 = (*it)->getValue().getFocus();
        m1 = 0.5f*(p2 - p0);
        fm1 = 0.5f*(f2 - f0);
        it++;
    } else {
        m1 = p2 - p1;
        fm1 = f2 - f1;
    }
    it++;
    it++;
    tgt::vec3 p3, f3;
    tgt::vec3 m2, fm2;
    if (it!= controlpoints.end()) {
        p3 = (*it)->getValue().getPosition();
        f3 = (*it)->getValue().getFocus();
        m2 = 0.5f*(p3 - p1);
        fm2 = 0.5f*(f3 - f1);
    }
    else {
        m2 = p2 - p1;
        fm2 = f2 - f1;
    }

    float t = (time - time1)/(time2 - time1);

    tgt::vec3 position;
    const float h00 = (1.f + 2.f*t)*(1.f - t)*(1.f - t);
    const float h10 = t*(1.f- t)*(1.f - t);
    const float h01 = t*t*(3.f-2.f*t);
    const float h11 = t*t*(t - 1.f);
    position.x = h00 * p1.x + h10 * m1.x + h01 * p2.x + h11 * m2.x;
    position.y = h00 * p1.y + h10 * m1.y + h01 * p2.y + h11 * m2.y;
    position.z = h00 * p1.z + h10 * m1.z + h01 * p2.z + h11 * m2.z;

    vec3 focus;
    focus = h00 * f1 + h10 * fm1 + h01 * f2 + h11 * fm2;

    vec3 upVector;
    Vec3SphericalLinearInterpolationFunction* intfunc2 = new Vec3SphericalLinearInterpolationFunction();
    upVector = normalize(intfunc2->interpolate(camera1.getUpVector(), camera2.getUpVector(), t));

    vec3 direction;
    const float dh00 = 6.f*t*t - 6.f*t;
    const float dh10 = 3.f*t*t - 4.f*t + 1;
    const float dh01 = -6.f*t*t + 6.f*t;
    const float dh11 = 3.f*t*t - 2.f*t;
    direction.x = dh00 * p1.x + dh10 * m1.x + dh01 * p2.x + dh11 * m2.x;
    direction.y = dh00 * p1.y + dh10 * m1.y + dh01 * p2.y + dh11 * m2.y;
    direction.z = dh00 * p1.z + dh10 * m1.z + dh01 * p2.z + dh11 * m2.z;

    Camera node = Camera(camera1);
    node.positionCamera(position, focus, upVector);
/*    Camera node = new Camera(position, focus, upVector, direction);
    if (camera1->isTangential() && camera2->isTangential())
        node->setTangential(true); */
    return node;
}

MultiPointInterpolationFunction<Camera>* CameraCatmullRomInterpolationFunction::create() const {
    return new CameraCatmullRomInterpolationFunction();
}

CameraSquadInterpolationFunction::CameraSquadInterpolationFunction() {}

std::string CameraSquadInterpolationFunction::getGuiName() const {
    return "spherical cubic interpolation";
}

std::string CameraSquadInterpolationFunction::getCategory() const {
    return "spline";
}

Camera CameraSquadInterpolationFunction::interpolate(std::vector<PropertyKeyValue<Camera>*> controlpoints, float time) const {
    Camera camera1;
    Camera camera2;
    tgt::vec3 p0;
    tgt::vec3 p1;
    tgt::vec3 p2;
    tgt::vec3 p3;
    tgt::quat q0;
    tgt::quat q1;
    tgt::quat q2;
    tgt::quat q3;
    std::vector<PropertyKeyValue<Camera>*>::iterator it;
    it = controlpoints.begin();
    while ((*it)->getTime() < time)
        it++;

    float time2 = (*it)->getTime();
    p2 = (*it)->getValue().getPosition();
    q2 = tgt::quat(normalize(p2), 0.0f);
    camera2 = (*it)->getValue();
    it--;
    float time1 = (*it)->getTime();
    p1 = (*it)->getValue().getPosition();
    q1 = tgt::quat(normalize(p1), 0.0f);
    camera1 = (*it)->getValue();
    if (it!= controlpoints.begin()) {
        it--;
        p0 = (*it)->getValue().getPosition();
        q0 = tgt::quat(normalize(p0), 0.0f);
        it++;
    }
    else {
        q0 = q1;
    }
    it++;
    it++;
    tgt::vec3 m2;
    if (it!= controlpoints.end()) {
        p3 = (*it)->getValue().getPosition();
        q3 = tgt::quat(normalize(p3), 0.0f);
    }
    else
        q3 = q2;

    float t = (time - time1)/(time2 - time1);

    tgt::quat q1in = conjugate(q1);
    tgt::quat qh1 = q1 * tgt::expQuat((tgt::logQuat(q1in*q0) + tgt::logQuat(q1in*q2)) * -0.25f);//splineQuat(q0, q1, q2);
    tgt::quat q2in = conjugate(q2);
    tgt::quat qh2 = q2 * tgt::expQuat((tgt::logQuat(q2in*q1) + tgt::logQuat(q2in*q1)) * -0.25f);//splineQuat(q1, q2, q3);
    tgt::quat pos = slerpQuat(slerpQuat(q1, q2, t, false), slerpQuat(qh1, qh2, t, false), 2.0f*t*(1.0f-t), false);//squadQuat(q1, q2, qh1, qh2, t);

    tgt::vec3 position;
    float length1 = length(p1);
    float length2 = length(p2);
    float lengthRes = length1 * (1-t)+ length2 *t;
    position = pos.vec() * lengthRes;

    vec3 focus;
    Vec3SphericalLinearInterpolationFunction* intf = new Vec3SphericalLinearInterpolationFunction();
    focus = intf->interpolate(camera1.getFocus(), camera2.getFocus(), t);

    vec3 upVector;
    Vec3SphericalLinearInterpolationFunction* intfunc2 = new Vec3SphericalLinearInterpolationFunction();
    upVector = normalize(intfunc2->interpolate(camera1.getUpVector(), camera2.getUpVector(), t));

    Camera node = Camera(camera1);
    node.positionCamera(position, focus, upVector);
/*    vec3 direction;
    direction = intf->interpolate(camera1->getDirection(), camera2->getDirection(), t );

    Camera node = new Camera(position, focus, upVector, direction);
    if (camera1->isTangential() && camera2->isTangential())
        node->setTangential(true); */
    return node;
}
MultiPointInterpolationFunction<Camera>* CameraSquadInterpolationFunction::create() const {
    return new CameraSquadInterpolationFunction();
}

}
