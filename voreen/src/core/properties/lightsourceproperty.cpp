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

#include "voreen/core/properties/lightsourceproperty.h"

#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties//property.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

namespace voreen {

LightSourceProperty::LightSourceProperty()
    : FloatVec4Property()
    , camProp_(0)
{}

LightSourceProperty::LightSourceProperty(const std::string& id, const std::string& guiText,
                             const tgt::vec4& value, const tgt::vec4& minValue,
                             const tgt::vec4& maxValue, int invalidationLevel, Property::LevelOfDetail lod)
                             : FloatVec4Property(id, guiText, value, minValue, maxValue, invalidationLevel, NumericProperty<tgt::vec4>::STATIC,lod)
    , followCam_(false)
    , curCenter_(tgt::vec3(0.f))
    , maxDist_(10.f)
    , lightPos_(tgt::vec4(0.f, 0.f, 1.f, 1.f))
    , camProp_(0)
{}

// TODO force user to pass camprop pointer instead of blindly searching for it
CameraProperty* LightSourceProperty::getCamera() {
    if(camProp_)
        return camProp_;
    else {
        PropertyOwner* propOwner = getOwner();
        if(!propOwner)
            return 0;

        std::vector<Property*> props = propOwner->getProperties();
        std::vector<Property*>::iterator it = props.begin();
        while (it != props.end()) {
            if(CameraProperty* camProp = dynamic_cast<CameraProperty*>(*it)) {
                setCamera(camProp);
                break;
            }
            ++it;
        }
        return camProp_;
    }
}

bool LightSourceProperty::hasCamera() const {
    return camProp_ != 0;
}

void LightSourceProperty::setCamera(CameraProperty* cam) {
    camProp_ = cam;
    if(camProp_) {
        camProp_->onChange(MemberFunctionCallback<LightSourceProperty>(this, &LightSourceProperty::cameraUpdate));
        curCenter_ = camProp_->getTrackball().getCenter();
        maxDist_ = camProp_->getMaxValue() / 50.f;
        updateWidgets();
        cameraUpdate();
    }
}

void LightSourceProperty::setLightPos(const tgt::vec4& lightPos) {
    if(getCamera()) {
        const tgt::Camera& camera = getCamera()->get();
        set(tgt::vec4((camera.getViewMatrixInverse().getRotationalPart() * lightPos).xyz() + getCamera()->getTrackball().getCenter(), 1.f));
    } else {
        set(lightPos);
    }
    lightPos_ = lightPos;
    updateWidgets();
}

tgt::vec4 LightSourceProperty::getLightPos() const {
    return lightPos_;
}

void LightSourceProperty::setFollowCam(bool b) {
    followCam_ = b;
    updateWidgets();
}

bool LightSourceProperty::getFollowCam() const {
    return followCam_;
}

void LightSourceProperty::setMaxDist(float d) {
    maxDist_ = d;
    updateWidgets();
}

float LightSourceProperty::getMaxDist() const {
    return maxDist_;
}

void LightSourceProperty::cameraUpdate() {
    if(!getCamera())
        return;

    float maxVal = getCamera()->getMaxValue() / 50.f;
    tgt::vec3 newCenter = getCamera()->getTrackball().getCenter();

    if(fabs(maxVal - maxDist_) > 1.f || newCenter != curCenter_) {
        float oldRelDist = length(get().xyz() - curCenter_) / maxDist_;
        tgt::vec3 oldPos = get().xyz();
        FloatVec4Property::set(tgt::vec4(normalize(oldPos - curCenter_) * oldRelDist * maxVal + newCenter, 1.0));
        curCenter_ = newCenter;
        maxDist_ = maxVal;
        oldRelDist = length(get().xyz() - curCenter_) / maxDist_;
    }

    const tgt::Camera& camera = getCamera()->get();
    if (followCam_) {
        FloatVec4Property::set(tgt::vec4((camera.getViewMatrixInverse().getRotationalPart() * lightPos_).xyz() + curCenter_, 1.f));
    }
    else {
        lightPos_ = tgt::vec4((camera.getViewMatrix().getRotationalPart() * tgt::vec4(get().xyz() - curCenter_, 0.f)).xyz(), 1.f);
        updateWidgets();
    }
}

void LightSourceProperty::serialize(Serializer& s) const {
    FloatVec4Property::serialize(s);
    s.serialize("Center", curCenter_);
    s.serialize("LightPos", lightPos_);
    s.serialize("FollowCam", followCam_);
    s.serialize("MaxDist", maxDist_);
}

void LightSourceProperty::deserialize(Deserializer& s) {
    FloatVec4Property::deserialize(s);
    try {
        s.deserialize("Center", curCenter_);
        s.deserialize("LightPos", lightPos_);
        s.deserialize("FollowCam", followCam_);
        s.deserialize("MaxDist", maxDist_);
    }
    catch(SerializationException&) {
        // TODO necessary?
        lightPos_ = tgt::vec4(get().xyz(), 1.f);
        s.removeLastError();
    }
}

void LightSourceProperty::set(const tgt::vec4& value) {
    FloatVec4Property::set(value);
    if(!getCamera())
        return;
    lightPos_ = tgt::vec4((getCamera()->get().getViewMatrix().getRotationalPart() * tgt::vec4(get().xyz() - curCenter_, 0.f)).xyz(), 1.f);
    updateWidgets();
}

} // namespace voreen
