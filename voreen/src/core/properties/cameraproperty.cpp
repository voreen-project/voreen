/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/core/interaction/voreentrackball.h"
#include "voreen/core/properties/cameraproperty.h"

#include "tgt/camera.h"

namespace voreen {

using tgt::Camera;

CameraProperty::CameraProperty(const std::string& id, const std::string& guiText,
                       tgt::Camera const value, bool sceneAdjuster, bool adjustProjectionToViewport,
                       float maxValue,
                       int invalidationLevel, Property::LevelOfDetail lod)
        : TemplateProperty<Camera>(id, guiText, value, invalidationLevel, lod)
        , sceneAdjuster_(sceneAdjuster)
        , trackball_(this)
        , maxValue_(maxValue)
        , minValue_(maxValue / 50000.f)
        , currentSceneBounds_(tgt::Bounds(tgt::vec3(0.f), tgt::vec3(0.f)))
        , centerOption_(SCENE)
        , adaptOnChange_(true)
{}

CameraProperty::CameraProperty()
    : trackball_(this)
    , sceneAdjuster_(false)
    , maxValue_(100.f)
    , minValue_(100.f / 50000.f)
    , currentSceneBounds_(tgt::Bounds(tgt::vec3(0.f), tgt::vec3(0.f)))
    , adaptOnChange_(true)
{ }

CameraProperty::~CameraProperty() {
}

Property* CameraProperty::create() const {
    return new CameraProperty();
}

void CameraProperty::set(const tgt::Camera& camera) {
    value_ = Camera(camera);
    invalidate();
}

void CameraProperty::setPosition(const tgt::vec3& pos) {
    value_.setPosition(pos);
    invalidate();
}

void CameraProperty::setFocus(const tgt::vec3& focus) {
    value_.setFocus(focus);
    invalidate();
}

void CameraProperty::setUpVector(const tgt::vec3& up) {
    value_.setUpVector(up);
    invalidate();
}

void CameraProperty::setFrustum(const tgt::Frustum& frust) {
    value_.setFrustum(frust);
    invalidate();
}

void CameraProperty::setFarDist(float dist) {
    value_.setFarDist(dist);
    invalidate();
}

void CameraProperty::setNearDist(float dist) {
    value_.setNearDist(dist);
    invalidate();
}

void CameraProperty::setUseOrthoZoomFlag(bool b) {
    value_.setUseOrthoZoomFactorFlag(b);
    invalidate();
}

void CameraProperty::setOrthoZoomFactor(float zoomFactor) {
    value_.setOrthoZoomFactor(zoomFactor);
    invalidate();
}

bool CameraProperty::setStereoEyeSeparation(float separation) {
    bool result = value_.setStereoEyeSeparation(separation);
    if(result)
        invalidate();
    return result;
}

bool CameraProperty::setStereoRelativeFocalLength(float stereoRelativeFocalLength) {
    bool result = value_.setStereoRelativeFocalLength(stereoRelativeFocalLength);
    if(result)
        invalidate();
    return result;
}

bool CameraProperty::setUseRealWorldFrustum(bool useRealWorldFrustum) {
    bool result = value_.setUseRealWorldFrustum(useRealWorldFrustum);
    if(result)
        invalidate();
    return result;
}

bool CameraProperty::setStereoFocalLength(float focallength) {
    bool result = value_.setStereoFocalLength(focallength);
    if(result)
        invalidate();
    return result;
}

bool CameraProperty::setStereoWidth(float width) {
    bool result = value_.setStereoWidth(width);
    if(result)
        invalidate();
    return result;
}

bool CameraProperty::setStereoEyeMode(tgt::Camera::StereoEyeMode mode) {
    bool result = value_.setStereoEyeMode(mode);
    if(result)
        invalidate();
    return result;
}

//bool CameraProperty::stereoShift(tgt::vec3 shift) {
    //bool result = value_.stereoShift(shift);
    //if(result)
        //invalidate();
    //return result;
//}

bool CameraProperty::setStereoAxisMode(tgt::Camera::StereoAxisMode mode) {
    bool result = value_.setStereoAxisMode(mode);
    if(result)
        invalidate();
    return result;
}

void CameraProperty::setMinValue(float val) {
    minValue_ = val;
}

float CameraProperty::getMinValue() const {
    return minValue_;
}

void CameraProperty::setMaxValue(float val) {
    maxValue_ = val;
}

float CameraProperty::getMaxValue() const {
    return maxValue_;
}

void CameraProperty::serialize(Serializer& s) const {
    Property::serialize(s);

    s.serialize("projectionMode", (int)value_.getProjectionMode());

    s.serialize("position", value_.getPosition());
    s.serialize("focus", value_.getFocus());
    s.serialize("upVector", value_.getUpVector());

    s.serialize("maxValue", maxValue_);
    s.serialize("minValue", minValue_);

    s.serialize("frustLeft", value_.getFrustLeft());
    s.serialize("frustRight", value_.getFrustRight());
    s.serialize("frustBottom", value_.getFrustBottom());
    s.serialize("frustTop", value_.getFrustTop());
    s.serialize("frustNear", value_.getNearDist(false));
    s.serialize("frustFar", value_.getFarDist(false));

    s.serialize("fovy", value_.getFovy());

    s.serialize("eyeMode", (int)value_.getStereoEyeMode());
    s.serialize("eyeSeparation", value_.getStereoEyeSeparation());
    s.serialize("axisMode", (int)value_.getStereoAxisMode());
    s.serialize("stereoFocalLengh", value_.getStereoFocalLength());
    s.serialize("stereoWidth", value_.getStereoWidth());
    s.serialize("stereoRelativeFocalLength", value_.getStereoRelativeFocalLength());
    s.serialize("useRealWorldFrustum", value_.getUseRealWorldFrustum());

    s.serialize("trackball", trackball_);

    s.serialize("centerOption", (int)centerOption_);
    s.serialize("adaptOnChange", adaptOnChange_);
    s.serialize("sceneLLF", currentSceneBounds_.getLLF());
    s.serialize("sceneURB", currentSceneBounds_.getURB());

    s.serialize("useOrthoZoomFactor", value_.getUseOrthoZoomFactorFlag());
    s.serialize("orthoZoomFactor", value_.getOrthoZoomFactor());
}

void CameraProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);

    try {
        float left, right, bottom, top, nearP, farP;
        s.deserialize("frustLeft", left);
        s.deserialize("frustRight", right);
        s.deserialize("frustBottom", bottom);
        s.deserialize("frustTop", top);
        s.deserialize("frustNear", nearP);
        s.deserialize("frustFar", farP);
        value_.setFrustum(tgt::Frustum(left, right, bottom, top, nearP, farP));
    } catch(SerializationException&) {
        s.removeLastError();
    }

    int projMode;
    try {
        s.deserialize("projectionMode", projMode);
    } catch(SerializationException&) {
        s.removeLastError();
        projMode = tgt::Camera::PERSPECTIVE;
    }

    value_.setProjectionMode((tgt::Camera::ProjectionMode)projMode);

    tgt::vec3 vector;
    s.deserialize("position", vector);
    value_.setPosition(vector);
    s.deserialize("focus", vector);
    value_.setFocus(vector);
    s.deserialize("upVector", vector);
    value_.setUpVector(vector);

    try {
        float maxValue, minValue;
        s.deserialize("maxValue", maxValue);
        maxValue_ = maxValue;
        s.deserialize("minValue", minValue);
        minValue_ = minValue;
    }
    catch(SerializationException&) {
        s.removeLastError();
    }

    try {
        float fovy;
        s.deserialize("fovy", fovy);
        value_.setFovy(fovy);
    }
    catch(SerializationException&) {
        s.removeLastError();
    }

    try {
        int eyeMode, axisMode;
        float eyeSep, focalLength, focalWidth, stereoRelativeFocalLength;
        bool useRealWorldFrustum;
        s.deserialize("eyeMode", eyeMode);
        s.deserialize("eyeSeparation", eyeSep);
        s.deserialize("axisMode", axisMode);
        value_.setStereoEyeMode((tgt::Camera::StereoEyeMode)eyeMode, false);
        value_.setStereoEyeSeparation(eyeSep, false);
        value_.setStereoAxisMode((tgt::Camera::StereoAxisMode)axisMode, false);
        s.deserialize("stereoFocalLengh", focalLength);
        value_.setStereoFocalLength(focalLength, false);
        s.deserialize("stereoWidth", focalWidth);
        value_.setStereoWidth(focalWidth, false);
        s.deserialize("stereoRelativeFocalLength", stereoRelativeFocalLength);
        value_.setStereoRelativeFocalLength(stereoRelativeFocalLength, false);
        s.deserialize("useRealWorldFrustum", useRealWorldFrustum);
        value_.setUseRealWorldFrustum(useRealWorldFrustum, false);
    }
    catch(SerializationException&) {
        s.removeLastError();
    }

    try {
        int centerOption;
        s.deserialize("centerOption", centerOption);
        s.deserialize("adaptOnChange", adaptOnChange_);
        tgt::vec3 llf, urb;
        s.deserialize("sceneLLF", llf);
        s.deserialize("sceneURB", urb);
        currentSceneBounds_ = tgt::Bounds(llf, urb);
        setTrackballCenterBehaviour((TrackballCenter)centerOption);
    }
    catch(SerializationException&) {
        s.removeLastError();
    }

    // Deserialize trackball_ AFTER calling setTrackBallCenterBehavior
    // because it might overwrite trackball_.center_.
    try {
        s.deserialize("trackball", trackball_);
    }
    catch (SerializationException&) {
        s.removeLastError();
    }

    bool useOrthoZoomFactorFlag;
    float orthoZoomFactor;
    s.optionalDeserialize("useOrthoZoomFactor", useOrthoZoomFactorFlag, false);
    s.optionalDeserialize("orthoZoomFactor", orthoZoomFactor, 0.01f);
    value_.setUseOrthoZoomFactorFlag(useOrthoZoomFactorFlag);
    value_.setOrthoZoomFactor(orthoZoomFactor);
}

void CameraProperty::look(tgt::ivec2 viewportSize) {
    value_.look(viewportSize);
}

VoreenTrackball& CameraProperty::getTrackball() {
    return trackball_;
}

const VoreenTrackball& CameraProperty::getTrackball() const {
    return trackball_;
}

void CameraProperty::adaptInteractionToScene(const tgt::Bounds& bounds, float nearDist, bool forceUpdate) {
    if (!adaptOnChange_)
        return;

    if(!bounds.isDefined()) {
        LWARNING("Tried to adapt camera to undefined bounds, not adapting to scene");
        return;
    }

    tgt::Camera cam = get();

    bool cameraHasChanged = false;

    //center has been changed
    if ((forceUpdate || (bounds.center() != currentSceneBounds_.center())) && centerOption_ == SCENE) {
        tgt::vec3 distVec = currentSceneBounds_.center() - cam.getPosition();
        trackball_.setCenter(bounds.center());
        cam.setFocus(bounds.center());
        cam.setPosition(bounds.center() - distVec);
        cameraHasChanged = true;
    }

    float minDiag = std::min(tgt::length(bounds.diagonal()),tgt::length(currentSceneBounds_.diagonal()));
    float maxDiag = std::max(tgt::length(bounds.diagonal()),tgt::length(currentSceneBounds_.diagonal()));

    currentSceneBounds_ = bounds;

    //diagonal has been changed significantly
    if(forceUpdate || minDiag == 0.f || (maxDiag / minDiag > std::sqrt(2.1f))) { // sqrt(2)+epsilon
        float oldRelCamDist = (cam.getFocalLength() - getMinValue()) / (getMaxValue() - getMinValue());
        float maxSideLength = tgt::max(currentSceneBounds_.diagonal());

        // The factor 250 is derived from an earlier constant maxDist of 500 and a constant maximum cubeSize element of 2
        float newMaxDist = 250.f * maxSideLength;
        setMaxValue(newMaxDist);

        if(nearDist > 0.f) {
            if(cam.getProjectionMode() != tgt::Camera::FRUSTUM)
                cam.setNearDist(nearDist);
            setMinValue(nearDist);
        } else
            setMinValue(newMaxDist / 50000.f);

        float newAbsCamDist = getMinValue() + oldRelCamDist * (getMaxValue() - getMinValue());

        if(centerOption_ == CAMSHIFT) {
            tgt::vec3 newFocus = cam.getFocus() * (newAbsCamDist / cam.getFocalLength());
            tgt::vec3 newPos   = cam.getPosition() * (newAbsCamDist / cam.getFocalLength());
            cam.setFocus(newFocus);
            cam.setPosition(newPos);
        } else if(centerOption_ == SCENE) {
            tgt::vec3 newPos   = currentSceneBounds_.center() - cam.getLook() * newAbsCamDist;
            cam.setPosition(newPos);
        }
        cam.setFarDist(newMaxDist + maxSideLength);
        cameraHasChanged = true;
    }

    if(cameraHasChanged) {
        set(cam);
        adjustOrthogonalZoomToFocalLength();
        LDEBUGC("voreen.CameraProperty", "Adapting camera handling to new scene size...");
    }
}

void CameraProperty::adjustOrthogonalZoomToFocalLength() {
    tgt::Camera cam = get();
    cam.setOrthoZoomFactor(0.1f / cam.getFocalLength());
    set(cam);
}

void CameraProperty::resetCameraFocusToTrackballCenter() {
    tgt::Camera cam = get();
    tgt::vec3 center = trackball_.getCenter();
    cam.setPosition(center - cam.getFocalLength() * cam.getLook());
    cam.setFocus(center);
    set(cam);
}

void CameraProperty::setAdaptOnChange(bool b) {
    adaptOnChange_ = b;
}

bool CameraProperty::getAdaptOnChange() const {
    return adaptOnChange_;
}

void CameraProperty::setTrackballCenterBehaviour(TrackballCenter t) {
    switch(t) {
        case SCENE:
            trackball_.setMoveCenter(false);
            if (currentSceneBounds_.volume() != 0) {
                trackball_.setCenter(currentSceneBounds_.center());
            }
            break;
        case WORLD:
            trackball_.setMoveCenter(false);
            trackball_.setCenter(tgt::vec3(0.f));
            break;
        case CAMSHIFT:
            trackball_.setMoveCenter(true);
            // trackball center is updated by the trackball internally if setMoveCenter is set to true
            break;
        default:
            tgtAssert(false, "Unknown trackball center option");
    }
    centerOption_ = t;
}

CameraProperty::TrackballCenter CameraProperty::getTrackballCenterBehaviour() const {
    return centerOption_;
}

bool CameraProperty::isSceneAdjuster() const {
    return sceneAdjuster_;
}

const tgt::Bounds& CameraProperty::getSceneBounds() const {
    return currentSceneBounds_;
}

void CameraProperty::setSceneBounds(const tgt::Bounds& b) {
    currentSceneBounds_ = b;
}

} // namespace voreen
