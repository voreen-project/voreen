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

#include "voreen/core/animation/propertykeyvalue.h"
#include "voreen/core/properties/volumeurllistproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/animation/interpolationfunction.h"
#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/transfunc2dprimitives.h"
#include "tgt/camera.h"

using tgt::Camera;

namespace voreen {

template <>
PropertyKeyValue<TransFunc1DKeys*>* PropertyKeyValue<TransFunc1DKeys*>::clone() const {
    PropertyKeyValue<TransFunc1DKeys*>* keyvalue = new PropertyKeyValue<TransFunc1DKeys*>();
    keyvalue->after_ = after_;
    keyvalue->before_ = before_;
    if (value_)
        keyvalue->value_ = static_cast<TransFunc1DKeys*>(value_->clone());
    else
        LERRORC("voreen.PropertyKeyValue<TransFunc1DKeys>", "No value");
    keyvalue->time_ = time_;
    keyvalue->smooth_ = smooth_;
    return keyvalue;
}

template <>
PropertyKeyValue<TransFunc2DPrimitives*>* PropertyKeyValue<TransFunc2DPrimitives*>::clone() const {
    PropertyKeyValue<TransFunc2DPrimitives*>* keyvalue = new PropertyKeyValue<TransFunc2DPrimitives*>();
    keyvalue->after_ = after_;
    keyvalue->before_ = before_;
    if (value_)
        keyvalue->value_ = static_cast<TransFunc2DPrimitives*>(value_->clone());
    else
        LERRORC("voreen.PropertyKeyValue<TransFunc2DPrimitives>", "No value");
    keyvalue->time_ = time_;
    keyvalue->smooth_ = smooth_;
    return keyvalue;
}

template <>
void PropertyKeyValue<tgt::Camera>::serialize(Serializer& s) const {
    // base class props
    s.serialize("time", time_);
    s.serialize("smooth", smooth_);
    s.serialize("before", before_);
    s.serialize("after", after_);

    // camera
    s.serialize("projectionMode", (int)value_.getProjectionMode());

    s.serialize("position", value_.getPosition());
    s.serialize("focus", value_.getFocus());
    s.serialize("upVector", value_.getUpVector());

    s.serialize("frustLeft", value_.getFrustLeft());
    s.serialize("frustRight", value_.getFrustRight());
    s.serialize("frustBottom", value_.getFrustBottom());
    s.serialize("frustTop", value_.getFrustTop());
    s.serialize("frustNear", value_.getNearDist());
    s.serialize("frustFar", value_.getFarDist());
    s.serialize("fovy", value_.getFovy());

}

template <>
void PropertyKeyValue<tgt::Camera>::deserialize(Deserializer& s) {
    // base class props
    s.deserialize("time", time_);
    s.deserialize("smooth", smooth_);
    s.deserialize("before", before_);
    s.deserialize("after", after_);

    // camera
    tgt::vec3 vector;
    s.deserialize("position", vector);
    value_.setPosition(vector);
    s.deserialize("focus", vector);
    value_.setFocus(vector);
    s.deserialize("upVector", vector);
    value_.setUpVector(vector);
    float winRatio = 1.f;
    try {
        s.deserialize("winRatio", winRatio);
    } catch(SerializationException&) {
        s.removeLastError();
    }

    // camera frustum
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
}

template class PropertyKeyValue<float>;
template class PropertyKeyValue<int>;
template class PropertyKeyValue<bool>;
template class PropertyKeyValue<tgt::ivec2>;
template class PropertyKeyValue<tgt::ivec3>;
template class PropertyKeyValue<tgt::ivec4>;
template class PropertyKeyValue<tgt::vec2>;
template class PropertyKeyValue<tgt::vec3>;
template class PropertyKeyValue<tgt::vec4>;
template class PropertyKeyValue<tgt::mat2>;
template class PropertyKeyValue<tgt::mat3>;
template class PropertyKeyValue<tgt::mat4>;
template class PropertyKeyValue<tgt::Camera>;
template class PropertyKeyValue<std::string>;
template class PropertyKeyValue<ShaderSource>;
template class PropertyKeyValue<TransFunc1DKeys*>;
template class PropertyKeyValue<TransFunc2DPrimitives*>;

} // namespace voreen
