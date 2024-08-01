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

#include "tgt/camera.h"

#include "tgt/assert.h"
#include "tgt/glmath.h"
#include "tgt/logmanager.h"
#include "tgt/quaternion.h"
#include "tgt/spline.h"
#include "tgt/tgt_gl.h"

#include <cmath>
#include <iostream>

namespace tgt {

// Constructor
Camera::Camera(const vec3& position, const vec3& focus, const vec3& up,
               float fovy, float ratio, float distn, float distf, ProjectionMode pm)
    : position_(position),
      focus_(focus),
      upVector_(normalize(up)),
      frust_(Frustum(fovy, ratio, distn, distf)),
      projectionMode_(pm),
      eyeSeparation_(65.f),
      stereoFocalLength_(600.f),
      stereoWidth_(450.f),
      eyeMode_(EYE_MIDDLE),
      stereoRelativeFocalLength_(0.05f),
      useRealWorldFrustum_(false),
      axisMode_(ON_AXIS),
      useOffset_(false),
      offset_(tgt::vec3(0.f)),
      useOrthoZoomFactor_(false),
      orthoZoomFactor_(0.01f)
{
    viewMatrix_ = mat4::createLookAt(position, focus, up);
}

Camera::Camera(const Camera& cam)
    : position_(cam.getPosition()),
      focus_(cam.getFocus()),
      upVector_(cam.getUpVector()),
      frust_(cam.getFrustum()),
      projectionMode_(cam.getProjectionMode()),
      eyeSeparation_(cam.getStereoEyeSeparation()),
      stereoFocalLength_(cam.getStereoFocalLength()),
      stereoWidth_(cam.getStereoWidth()),
      eyeMode_(cam.getStereoEyeMode()),
      stereoRelativeFocalLength_(cam.getStereoRelativeFocalLength()),
      useRealWorldFrustum_(cam.getUseRealWorldFrustum()),
      axisMode_(cam.getStereoAxisMode()),
      useOffset_(cam.isOffsetEnabled()),
      offset_(cam.getOffset()),
      useOrthoZoomFactor_(cam.getUseOrthoZoomFactorFlag()),
      orthoZoomFactor_(cam.getOrthoZoomFactor())
{
    viewMatrix_ = mat4::createLookAt(getPositionWithOffsets(), getFocusWithOffsets(), upVector_);
}

Camera::~Camera() {
}

Camera* Camera::clone() const {
    Camera* cam = new Camera(position_, focus_, upVector_, frust_.getFovy(),
        frust_.getRatio(), frust_.getNearDist(), frust_.getFarDist());
    cam->setStereoAxisMode(axisMode_);
    cam->setStereoEyeMode(eyeMode_, false);
    cam->setStereoEyeSeparation(eyeSeparation_, false);
    cam->setStereoFocalLength(stereoFocalLength_, false);
    cam->setStereoWidth(stereoWidth_, false);
    cam->setStereoRelativeFocalLength(stereoRelativeFocalLength_, false);
    cam->setUseRealWorldFrustum(useRealWorldFrustum_, false);
    cam->enableOffset(useOffset_);
    cam->setOffset(offset_);
    cam->setUseOrthoZoomFactorFlag(useOrthoZoomFactor_);
    cam->setOrthoZoomFactor(orthoZoomFactor_);
    return cam;
}

bool Camera::operator==(const Camera& rhs) const {
    return (rhs.position_ == position_) && (rhs.focus_ == focus_) && (rhs.upVector_ == upVector_) && (rhs.frust_ == frust_) &&
        (rhs.projectionMode_ == projectionMode_) && (rhs.useOrthoZoomFactor_ == useOrthoZoomFactor_) && (rhs.orthoZoomFactor_ == orthoZoomFactor_);
}

bool Camera::operator!=(const Camera& rhs) const {
    return !(*this == rhs);
}

// This is called to set up the Camera-View
void Camera::look(float windowRatio) {
    MatStack.matrixMode(MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    updateFrustum();
    MatStack.loadMatrix(getFrustumMatrix(windowRatio));
    //getProjectionMatrix();
    MatStack.matrixMode(MatrixStack::MODELVIEW);
    MatStack.loadIdentity();
    updateVM();
    MatStack.loadMatrix(viewMatrix_);
}

void Camera::look(ivec2 windowSize) {
    look(static_cast<float>(windowSize.x) / windowSize.y);
}

// Private method that updates the relevant frustum parameters
void Camera::updateFrustum() {
    frust_.update(this);
}

mat4 Camera::getViewMatrix() const {
    updateVM();
    return viewMatrix_;
}

void Camera::setViewMatrix(const mat4& mvMat) {
    mat4 inv;
    if (mvMat.invert(inv)) {
        // preserve the focallength
        float focallength = length(focus_ - position_);

        // calculate world-coordinates
        vec4 pos   = (inv * vec4(0.f, 0.f,  0.f, 1.f));
        vec4 look  = (inv * vec4(0.f, 0.f, -1.f, 0.f));
        vec4 focus = pos + focallength * look;
        vec4 up    = (inv * vec4(0.f, 1.f,  0.f, 0.f));

        positionCamera(pos.xyz(), focus.xyz(), up.xyz());

        updateVM();
    }
}

mat4 Camera::getRotateMatrix() const {
    updateVM();
    return viewMatrix_.getRotationalPart();
}

mat4 Camera::getViewMatrixInverse() const {
    updateVM();
    mat4 inv;
    if (viewMatrix_.invert(inv))
        return inv;
    else
        return mat4::identity;
}

tgt::Frustum Camera::getFrustumWithOffsets(float windowRatio) const {
    if(eyeMode_ != EYE_MIDDLE && axisMode_ == ON_AXIS)
        return stereoFrustumShift(windowRatio);
    else
        return frust_;
}

mat4 Camera::getFrustumMatrix(float windowRatio) const {
    Frustum f = getFrustumWithOffsets(windowRatio);
    return mat4::createFrustum(f.getLeft() * windowRatio, f.getRight() * windowRatio,
                               f.getBottom(), f.getTop(),
                               f.getNearDist(), f.getFarDist());
}

mat4 Camera::getFrustumMatrix(ivec2 windowSize) const {
    return getFrustumMatrix(static_cast<float>(windowSize.x) / windowSize.y);
}

mat4 Camera::getProjectionMatrix(ivec2 windowSize) const {
    return getProjectionMatrix(static_cast<float>(windowSize.x) / windowSize.y);
}

mat4 Camera::getProjectionMatrix(float windowRatio) const {
    if(projectionMode_ == ORTHOGRAPHIC) {
        tgt::mat4 scalingMatrix = useOrthoZoomFactor_ ? mat4::createScale(tgt::vec3(orthoZoomFactor_, orthoZoomFactor_, 1.f)) : mat4::identity;
        if(windowRatio > 1.0f)
            return scalingMatrix * mat4::createOrtho(frust_.getLeft() * windowRatio, frust_.getRight() * windowRatio,
                                     frust_.getBottom(), frust_.getTop(),
                                    -frust_.getNearDist(), frust_.getFarDist());
        else
            return scalingMatrix * mat4::createOrtho(frust_.getLeft(), frust_.getRight(),
                                     frust_.getBottom() * (1.0f/windowRatio), frust_.getTop() * (1.0f/windowRatio),
                                    -frust_.getNearDist(), frust_.getFarDist());
    } else if(projectionMode_ == PERSPECTIVE) {
        float fovy = frust_.getFovy();
        if(fovy < 6.f)
            fovy = 6.f;
        if(fovy > 175.f)
            fovy = 175.f;

        if(windowRatio >= 1.0f)
            return mat4::createPerspective(deg2rad(fovy), frust_.getRatio() * windowRatio , frust_.getNearDist(), frust_.getFarDist());
        else
            return mat4::createPerspective(atan(tan(deg2rad(fovy/2.f))/(windowRatio* frust_.getRatio()))*2, frust_.getRatio() * windowRatio , frust_.getNearDist(), frust_.getFarDist());
    }
    else
        return getFrustumMatrix(windowRatio);
}

float Camera::getFarDist(bool includeOffsets) const {
    if(!includeOffsets || eyeMode_ == EYE_MIDDLE || !useRealWorldFrustum_)
        return frust_.getFarDist();

    return std::max(frust_.getFarDist(), 2.f * stereoFocalLength_);
}

float Camera::getNearDist(bool includeOffsets) const {
    if(!includeOffsets || eyeMode_ == EYE_MIDDLE || !useRealWorldFrustum_)
        return frust_.getNearDist();

    return std::max(frust_.getNearDist(), 2.f * stereoFocalLength_ / 50000.f);
}

line3 Camera::getViewRay(ivec2 vp, ivec2 pixel) const {
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];

    tgt::mat4 projection_tgt = getProjectionMatrix(vp);
    tgt::mat4 modelview_tgt = getViewMatrix();
    for (int i = 0; i < 4; ++i) {
        modelview[i+0]   = modelview_tgt[i].x;
        modelview[i+4]   = modelview_tgt[i].y;
        modelview[i+8]   = modelview_tgt[i].z;
        modelview[i+12]  = modelview_tgt[i].w;
        projection[i+0]  = projection_tgt[i].x;
        projection[i+4]  = projection_tgt[i].y;
        projection[i+8]  = projection_tgt[i].z;
        projection[i+12] = projection_tgt[i].w;
    }
    viewport[0] = 0;
    viewport[1] = 0;
    viewport[2] = static_cast<GLint>(vp.x);
    viewport[3] = static_cast<GLint>(vp.y);

    GLfloat winX, winY, winZ;
    winX = (GLfloat)pixel.x;
    winY = (GLfloat)(vp.y - pixel.y);
    winZ = 1.0f;

    tgt::dvec3 posFar;
    tgt::gluUnProject(winX, winY, winZ, modelview, projection, viewport, &posFar.x, &posFar.y, &posFar.z);
    LGL_ERROR;

    winZ = 0.0f;
    tgt::dvec3 posNear;
    tgt::gluUnProject(winX, winY, winZ, modelview, projection, viewport, &posNear.x, &posNear.y, &posNear.z);
    LGL_ERROR;

    tgt::line3 l(posNear, posFar);
    return l;
}

vec3 Camera::project(ivec2 vp, vec3 point) const {
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];

    tgt::mat4 projection_tgt = getProjectionMatrix(vp);
    tgt::mat4 modelview_tgt = getViewMatrix();
    for (int i = 0; i < 4; ++i) {
        modelview[i+0]   = modelview_tgt[i].x;
        modelview[i+4]   = modelview_tgt[i].y;
        modelview[i+8]   = modelview_tgt[i].z;
        modelview[i+12]  = modelview_tgt[i].w;
        projection[i+0]  = projection_tgt[i].x;
        projection[i+4]  = projection_tgt[i].y;
        projection[i+8]  = projection_tgt[i].z;
        projection[i+12] = projection_tgt[i].w;
    }
    viewport[0] = 0;
    viewport[1] = 0;
    viewport[2] = vp.x;
    viewport[3] = vp.y;

    GLdouble pointProjectedGL[3];
    tgt::gluProject(point.x, point.y, point.z, modelview, projection, viewport,
               &pointProjectedGL[0], &pointProjectedGL[1], &pointProjectedGL[2]);

    return tgt::vec3(static_cast<float>(pointProjectedGL[0]),
                     static_cast<float>(pointProjectedGL[1]),
                     static_cast<float>(pointProjectedGL[2]));
}

bool Camera::setStereoEyeSeparation(float separation, bool updateCam) {
    eyeSeparation_ = separation;
    return updateCam;
}

bool Camera::setStereoFocalLength(float focallength, bool updateCam) {
    if(focallength <= 0.f || focallength == stereoFocalLength_)
        return false;

    stereoFocalLength_ = focallength;
    return updateCam;;
}

bool Camera::setStereoWidth(float width, bool updateCam) {
    if(width <= 0.f || width == stereoFocalLength_)
        return false;

    stereoWidth_ = width;
    return updateCam;;
}

bool Camera::setStereoEyeMode(StereoEyeMode mode, bool updateCam) {
    if(eyeMode_ == mode)
        return false;
    eyeMode_ = mode;
    return updateCam;
}

bool Camera::setStereoRelativeFocalLength(float stereoRelativeFocalLength, bool updateCam) {
    if(stereoRelativeFocalLength < 0.001f || stereoRelativeFocalLength > 1.f || stereoRelativeFocalLength == stereoRelativeFocalLength_)
        return false;
    stereoRelativeFocalLength_ = stereoRelativeFocalLength;
    return updateCam;
}

bool Camera::setUseRealWorldFrustum(bool useRealWorldFrustum, bool updateCam) {
    if(useRealWorldFrustum == useRealWorldFrustum_)
        return false;
    useRealWorldFrustum_ = useRealWorldFrustum;
    return updateCam;
}

tgt::vec3 Camera::getStereoShift() const {
    if(eyeMode_ == EYE_MIDDLE)
        return tgt::vec3(0.f);

    float eyeSep = eyeSeparation_;
    // heuristic: in non-real-world mode, adapt eye-separation to 1 / 1000.f of relative focal length times frustum size
    if(!useRealWorldFrustum_)
        eyeSep *= stereoRelativeFocalLength_ * (frust_.getFarDist() - frust_.getNearDist()) * 0.001f;

    if(eyeMode_ == EYE_LEFT)
        return getStrafeWithOffsets() * (-eyeSep * 0.5f);
    else
        return getStrafeWithOffsets() * ( eyeSep * 0.5f);

    //should not get here, but removes a warning
    return tgt::vec3(0.f);
}

tgt::vec3 Camera::getStrafeWithOffsets() const {
    tgt::vec3 pos = position_;
    if(useOffset_)
        pos += offset_;
    return normalize(cross(focus_ - pos, upVector_));
}

tgt::vec3 Camera::getFocusWithOffsets() const {
    tgt::vec3 foc = focus_;
    if(eyeMode_ != EYE_MIDDLE)
        foc += getStereoShift();
    // TODO necessary??
    //if(useOffset_);
        //foc += offset_;
    return foc;
}

tgt::vec3 Camera::getPositionWithOffsets() const {
    tgt::vec3 pos = position_;
    if(eyeMode_ != EYE_MIDDLE)
        pos += getStereoShift();
    if(useOffset_)
        pos += offset_;
    return pos;
}

tgt::Frustum Camera::stereoFrustumShift(float windowRatio) const {
    Frustum f = frust_;

    // for stereo, compute frustum based on near value, specified display distance and width
    float left, right;
    if(useRealWorldFrustum_) {
        right = f.getNearDist() * 0.5f * stereoWidth_ * (1.f / stereoFocalLength_);
        left  = -right;
        f.setFarDist(std::max(f.getFarDist(), 2.f * stereoFocalLength_));
        f.setNearDist(std::max(frust_.getNearDist(), 2.f * stereoFocalLength_ / 50000.f));
    } else {
        right = f.getRight();
        left  = f.getLeft();
    }

    float focalLength;
    if(useRealWorldFrustum_)
        focalLength = stereoFocalLength_;
    else
        focalLength = f.getNearDist() + (f.getFarDist() - f.getNearDist()) * stereoRelativeFocalLength_;

    tgt::vec3 shift = getStereoShift() / windowRatio;
    float moveXScalar = tgt::dot(shift, getStrafe());

    if(useRealWorldFrustum_) {
        f.setBottom(left);
        f.setTop(right);
    }
    f.setRight((right * focalLength - moveXScalar * f.getNearDist()) / focalLength);
    f.setLeft( ( left * focalLength - moveXScalar * f.getNearDist()) / focalLength);

    return f;
}

} // namespace tgt
