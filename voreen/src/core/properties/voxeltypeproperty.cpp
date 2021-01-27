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

#include "voreen/core/properties/voxeltypeproperty.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "tgt/types.h"

namespace voreen {

const std::string VoxelTypeProperty::loggerCat_("voreen.VoxelTypeProperty");

VoxelTypeProperty::VoxelTypeProperty(const std::string& id, const std::string& guiText,
                    int invalidationLevel, Property::LevelOfDetail lod)
    : TemplateProperty<void*>(id, guiText, 0, invalidationLevel, lod)
    , volume_(0)
    , intVal_(0)
    , floatVal_(0.f)
    , minIntVal_(0)
    , maxIntVal_(1<<20)
    , minFloatVal_(tgt::vec4(0.f))
    , maxFloatVal_(tgt::vec4(1.f))
{}

VoxelTypeProperty::VoxelTypeProperty()
    : TemplateProperty<void*>("", "", 0, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , volume_(0)
    , intVal_(0)
    , floatVal_(0.f)
    , minIntVal_(0)
    , maxIntVal_(1<<20)
    , minFloatVal_(tgt::vec4(0.f))
    , maxFloatVal_(tgt::vec4(1.f))
{}

Property* VoxelTypeProperty::create() const {
    return new VoxelTypeProperty();
}

std::string VoxelTypeProperty::getGuiName() const {
    return Property::getGuiName();
}

std::string VoxelTypeProperty::getID() const {
    return Property::getID();
}

void VoxelTypeProperty::setGuiName(const std::string& guiname) {
    Property::setGuiName(guiname);
}

void VoxelTypeProperty::setVolume(const VolumeBase* volumeHandle, bool updateRange) {
    volume_ = volumeHandle;
    if (updateRange)
        updateMinMaxRange();

    updateWidgets();
}

const VolumeBase* VoxelTypeProperty::getVolume() const {
    return volume_;
}

void VoxelTypeProperty::updateMinMaxRange() {
    if (!volume_)
        return;

    const VolumeRAM* volume = volume_->getRepresentation<VolumeRAM>();
    tgtAssert(volume, "no CPU representation");

    // single-channel integer volume
    if (dynamic_cast<const VolumeRAM_UInt8*>(volume)  ||
            dynamic_cast<const VolumeRAM_UInt16*>(volume) ||
            dynamic_cast<const VolumeRAM_UInt32*>(volume) ||
            dynamic_cast<const VolumeRAM_Int8*>(volume)   ||
            dynamic_cast<const VolumeRAM_Int16*>(volume)  ||
            dynamic_cast<const VolumeRAM_Int32*>(volume)  ) {
        setMinValue<int>(tgt::iround(volume->elementRange().x));
        setMaxValue<int>(tgt::iround(volume->elementRange().y));
        setValue<int>(tgt::clamp(getValue<int>(), getMinValue<int>(), getMaxValue<int>()));
    }
    // single-channel float volume
    else if (dynamic_cast<const VolumeRAM_Float*>(volume)  ||
            dynamic_cast<const VolumeRAM_Double*>(volume) ) {
        setMinValue<float>(volume->elementRange().x);
        setMaxValue<float>(volume->elementRange().y);
        setValue<float>(tgt::clamp(getValue<float>(), getMinValue<float>(), getMaxValue<float>()));
    }
    // RGB integer volume
    else if (dynamic_cast<const VolumeRAM_3xUInt8*>(volume)  ||
            dynamic_cast<const VolumeRAM_3xUInt16*>(volume) ||
            dynamic_cast<const VolumeRAM_3xInt8*>(volume)   ||
            dynamic_cast<const VolumeRAM_3xInt16*>(volume)  ) {
        setMinValue<tgt::ivec3>(tgt::ivec3(tgt::iround(volume->elementRange().x)));
        setMaxValue<tgt::ivec3>(tgt::ivec3(tgt::iround(volume->elementRange().y)));
        setValue<tgt::ivec3>(tgt::clamp(getValue<tgt::ivec3>(), getMinValue<tgt::ivec3>(), getMaxValue<tgt::ivec3>()));
    }
    // RGB float volume
    else if (dynamic_cast<const VolumeRAM_3xFloat*>(volume)  ||
            dynamic_cast<const VolumeRAM_3xDouble*>(volume) ) {
        setMinValue<tgt::vec3>(tgt::vec3(volume->elementRange().x));
        setMaxValue<tgt::vec3>(tgt::vec3(volume->elementRange().y));
        setValue<tgt::vec3>(tgt::clamp(getValue<tgt::vec3>(), getMinValue<tgt::vec3>(), getMaxValue<tgt::vec3>()));
    }
    else {
        LWARNING("updateMinMaxRange() unknown volume type: " << typeid(volume).name());
    }

}

void VoxelTypeProperty::serialize(Serializer& s) const {
    Property::serialize(s);

    s.serialize("intVal", intVal_);
    s.serialize("minIntVal", minIntVal_);
    s.serialize("maxIntVal", maxIntVal_);

    s.serialize("floatVal", floatVal_);
    s.serialize("minFloatVal", minFloatVal_);
    s.serialize("maxFloatVal", maxFloatVal_);
}

void VoxelTypeProperty::deserialize(Deserializer& d) {
    Property::deserialize(d);

    d.deserialize("intVal", intVal_);
    d.deserialize("minIntVal", minIntVal_);
    d.deserialize("maxIntVal", maxIntVal_);

    d.deserialize("floatVal", floatVal_);
    d.deserialize("minFloatVal", minFloatVal_);
    d.deserialize("maxFloatVal", maxFloatVal_);
}


// template specialization definitions
// -----------------------------------

//
// get (scalar)
//

template<>
uint8_t VoxelTypeProperty::getValue<uint8_t>() const {
    return static_cast<uint8_t>(intVal_.x);
}

template<>
uint16_t VoxelTypeProperty::getValue<uint16_t>() const {
    return static_cast<uint16_t>(intVal_.x);
}

template<>
uint32_t VoxelTypeProperty::getValue<uint32_t>() const {
    return static_cast<uint32_t>(intVal_.x);
}

template<>
uint64_t VoxelTypeProperty::getValue<uint64_t>() const {
    return static_cast<uint64_t>(intVal_.x);
}

template<>
int8_t VoxelTypeProperty::getValue<int8_t>() const {
    return static_cast<int8_t>(intVal_.x);
}

template<>
int16_t VoxelTypeProperty::getValue<int16_t>() const {
    return static_cast<int16_t>(intVal_.x);
}

template<>
int32_t VoxelTypeProperty::getValue<int32_t>() const {
    return static_cast<int32_t>(intVal_.x);
}

template<>
int64_t VoxelTypeProperty::getValue<int64_t>() const {
    return static_cast<int64_t>(intVal_.x);
}

template<>
float VoxelTypeProperty::getValue<float>() const {
    return floatVal_.x;
}

template<>
double VoxelTypeProperty::getValue<double>() const {
    return static_cast<double>(floatVal_.x);
}

// get (vec2)

template<>
tgt::Vector2<uint8_t> VoxelTypeProperty::getValue<tgt::Vector2<uint8_t> >() const {
    return tgt::Vector2<uint8_t>(intVal_.xy());
}

template<>
tgt::Vector2<uint16_t> VoxelTypeProperty::getValue<tgt::Vector2<uint16_t> >() const {
    return tgt::Vector2<uint16_t>(intVal_.xy());
}

template<>
tgt::Vector2<uint32_t> VoxelTypeProperty::getValue<tgt::Vector2<uint32_t> >() const {
    return tgt::Vector2<uint32_t>(intVal_.xy());
}

template<>
tgt::Vector2<uint64_t> VoxelTypeProperty::getValue<tgt::Vector2<uint64_t> >() const {
    return tgt::Vector2<uint64_t>(intVal_.xy());
}

template<>
tgt::Vector2<int8_t> VoxelTypeProperty::getValue<tgt::Vector2<int8_t> >() const {
    return tgt::Vector2<int8_t>(intVal_.xy());
}

template<>
tgt::Vector2<int16_t> VoxelTypeProperty::getValue<tgt::Vector2<int16_t> >() const {
    return tgt::Vector2<int16_t>(intVal_.xy());
}

template<>
tgt::Vector2<int32_t> VoxelTypeProperty::getValue<tgt::Vector2<int32_t> >() const {
    return tgt::Vector2<int32_t>(intVal_.xy());
}

template<>
tgt::Vector2<int64_t> VoxelTypeProperty::getValue<tgt::Vector2<int64_t> >() const {
    return tgt::Vector2<int64_t>(intVal_.xy());
}

template<>
tgt::vec2 VoxelTypeProperty::getValue<tgt::vec2>() const {
    return tgt::vec2(floatVal_.xy());
}

template<>
tgt::dvec2 VoxelTypeProperty::getValue<tgt::dvec2>() const {
    return tgt::dvec2(floatVal_.xy());
}

// get (vec3)
template<>
tgt::Vector3<uint8_t> VoxelTypeProperty::getValue<tgt::Vector3<uint8_t> >() const {
    return tgt::Vector3<uint8_t>(intVal_.xyz());
}

template<>
tgt::Vector3<uint16_t> VoxelTypeProperty::getValue<tgt::Vector3<uint16_t> >() const {
    return tgt::Vector3<uint16_t>(intVal_.xyz());
}

template<>
tgt::Vector3<uint32_t> VoxelTypeProperty::getValue<tgt::Vector3<uint32_t> >() const {
    return tgt::Vector3<uint32_t>(intVal_.xyz());
}

template<>
tgt::Vector3<uint64_t> VoxelTypeProperty::getValue<tgt::Vector3<uint64_t> >() const {
    return tgt::Vector3<uint64_t>(intVal_.xyz());
}

template<>
tgt::Vector3<int8_t> VoxelTypeProperty::getValue<tgt::Vector3<int8_t> >() const {
    return tgt::Vector3<int8_t>(intVal_.xyz());
}

template<>
tgt::Vector3<int16_t> VoxelTypeProperty::getValue<tgt::Vector3<int16_t> >() const {
    return tgt::Vector3<int16_t>(intVal_.xyz());
}

template<>
tgt::Vector3<int32_t> VoxelTypeProperty::getValue<tgt::Vector3<int32_t> >() const {
    return tgt::Vector3<int32_t>(intVal_.xyz());
}

template<>
tgt::Vector3<int64_t> VoxelTypeProperty::getValue<tgt::Vector3<int64_t> >() const {
    return tgt::Vector3<int64_t>(intVal_.xyz());
}

template<>
tgt::vec3 VoxelTypeProperty::getValue<tgt::vec3>() const {
    return tgt::vec3(floatVal_.xyz());
}

template<>
tgt::dvec3 VoxelTypeProperty::getValue<tgt::dvec3>() const {
    return tgt::dvec3(floatVal_.xyz());
}

// get (vec4)
template<>
tgt::Vector4<uint8_t> VoxelTypeProperty::getValue<tgt::Vector4<uint8_t> >() const {
    return tgt::Vector4<uint8_t>(intVal_);
}

template<>
tgt::Vector4<uint16_t> VoxelTypeProperty::getValue<tgt::Vector4<uint16_t> >() const {
    return tgt::Vector4<uint16_t>(intVal_);
}

template<>
tgt::Vector4<uint32_t> VoxelTypeProperty::getValue<tgt::Vector4<uint32_t> >() const {
    return tgt::Vector4<uint32_t>(intVal_);
}

template<>
tgt::Vector4<uint64_t> VoxelTypeProperty::getValue<tgt::Vector4<uint64_t> >() const {
    return tgt::Vector4<uint64_t>(intVal_);
}

template<>
tgt::Vector4<int8_t> VoxelTypeProperty::getValue<tgt::Vector4<int8_t> >() const {
    return tgt::Vector4<int8_t>(intVal_);
}

template<>
tgt::Vector4<int16_t> VoxelTypeProperty::getValue<tgt::Vector4<int16_t> >() const {
    return tgt::Vector4<int16_t>(intVal_);
}

template<>
tgt::Vector4<int32_t> VoxelTypeProperty::getValue<tgt::Vector4<int32_t> >() const {
    return tgt::Vector4<int32_t>(intVal_);
}

template<>
tgt::Vector4<int64_t> VoxelTypeProperty::getValue<tgt::Vector4<int64_t> >() const {
    return tgt::Vector4<int64_t>(intVal_);
}

template<>
tgt::vec4 VoxelTypeProperty::getValue<tgt::vec4>() const {
    return floatVal_;
}

template<>
tgt::dvec4 VoxelTypeProperty::getValue<tgt::dvec4>() const {
    return tgt::dvec4(floatVal_);
}

//
// set (scalar)
//
template<>
void VoxelTypeProperty::setValue<uint8_t>(const uint8_t& value) {
    intVal_.x = static_cast<int>(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<uint16_t>(const uint16_t& value) {
    intVal_.x = static_cast<int>(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<uint32_t>(const uint32_t& value) {
    intVal_.x = static_cast<int>(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<uint64_t>(const uint64_t& value) {
    intVal_.x = static_cast<int>(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<int8_t>(const int8_t& value) {
    intVal_.x = static_cast<int>(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<int16_t>(const int16_t& value) {
    intVal_.x = static_cast<int>(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<int32_t>(const int32_t& value) {
    intVal_.x = static_cast<int>(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<int64_t>(const int64_t& value) {
    intVal_.x = static_cast<int>(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<float>(const float& value) {
    floatVal_.x = value;
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<double>(const double& value) {
    floatVal_.x = static_cast<float>(value);
    invalidate();
}

// set (vec2)
template<>
void VoxelTypeProperty::setValue<tgt::Vector2<uint8_t> >(const tgt::Vector2<uint8_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector2<uint16_t> >(const tgt::Vector2<uint16_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector2<uint32_t> >(const tgt::Vector2<uint32_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector2<uint64_t> >(const tgt::Vector2<uint64_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector2<int8_t> >(const tgt::Vector2<int8_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector2<int16_t> >(const tgt::Vector2<int16_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector2<int32_t> >(const tgt::Vector2<int32_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector2<int64_t> >(const tgt::Vector2<int64_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::vec2>(const tgt::vec2& value) {
    floatVal_.x = value.x;
    floatVal_.y = value.y;
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::dvec2>(const tgt::dvec2& value) {
    floatVal_.x = static_cast<float>(value.x);
    floatVal_.y = static_cast<float>(value.y);
    invalidate();
}

// set (vec3)
template<>
void VoxelTypeProperty::setValue<tgt::Vector3<uint8_t> >(const tgt::Vector3<uint8_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    intVal_.z = static_cast<int>(value.z);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector3<uint16_t> >(const tgt::Vector3<uint16_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    intVal_.z = static_cast<int>(value.z);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector3<uint32_t> >(const tgt::Vector3<uint32_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    intVal_.z = static_cast<int>(value.z);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector3<uint64_t> >(const tgt::Vector3<uint64_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    intVal_.z = static_cast<int>(value.z);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector3<int8_t> >(const tgt::Vector3<int8_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    intVal_.z = static_cast<int>(value.z);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector3<int16_t> >(const tgt::Vector3<int16_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    intVal_.z = static_cast<int>(value.z);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector3<int32_t> >(const tgt::Vector3<int32_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    intVal_.z = static_cast<int>(value.z);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector3<int64_t> >(const tgt::Vector3<int64_t>& value) {
    intVal_.x = static_cast<int>(value.x);
    intVal_.y = static_cast<int>(value.y);
    intVal_.z = static_cast<int>(value.z);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::vec3>(const tgt::vec3& value) {
    floatVal_.x = value.x;
    floatVal_.y = value.y;
    floatVal_.z = value.z;
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::dvec3>(const tgt::dvec3& value) {
    floatVal_.x = static_cast<float>(value.x);
    floatVal_.y = static_cast<float>(value.y);
    floatVal_.z = static_cast<float>(value.z);
    invalidate();
}

// set (vec4)
template<>
void VoxelTypeProperty::setValue<tgt::Vector4<uint8_t> >(const tgt::Vector4<uint8_t>& value) {
    intVal_ = tgt::ivec4(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector4<uint16_t> >(const tgt::Vector4<uint16_t>& value) {
    intVal_ = tgt::ivec4(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector4<uint32_t> >(const tgt::Vector4<uint32_t>& value) {
    intVal_ = tgt::ivec4(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector4<uint64_t> >(const tgt::Vector4<uint64_t>& value) {
    intVal_ = tgt::ivec4(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector4<int8_t> >(const tgt::Vector4<int8_t>& value) {
    intVal_ = tgt::ivec4(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector4<int16_t> >(const tgt::Vector4<int16_t>& value) {
    intVal_ = tgt::ivec4(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector4<int32_t> >(const tgt::Vector4<int32_t>& value) {
    intVal_ = tgt::ivec4(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::Vector4<int64_t> >(const tgt::Vector4<int64_t>& value) {
    intVal_ = tgt::ivec4(value);
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::vec4>(const tgt::vec4& value) {
    floatVal_ = value;
    invalidate();
}

template<>
void VoxelTypeProperty::setValue<tgt::dvec4>(const tgt::dvec4& value) {
    floatVal_ = tgt::vec4(value);
    invalidate();
}

// getMinValue
template<>
int VoxelTypeProperty::getMinValue<int>() const {
    return minIntVal_.x;
}

template<>
float VoxelTypeProperty::getMinValue<float>() const {
    return minFloatVal_.x;
}

template<>
tgt::ivec2 VoxelTypeProperty::getMinValue<tgt::ivec2>() const {
    return tgt::ivec2(minIntVal_.xy());
}

template<>
tgt::vec2 VoxelTypeProperty::getMinValue<tgt::vec2>() const {
    return tgt::vec2(minFloatVal_.xy());
}

template<>
tgt::ivec3 VoxelTypeProperty::getMinValue<tgt::ivec3>() const {
    return tgt::ivec3(minIntVal_.xyz());
}

template<>
tgt::vec3 VoxelTypeProperty::getMinValue<tgt::vec3>() const {
    return tgt::vec3(minFloatVal_.xyz());
}

template<>
tgt::ivec4 VoxelTypeProperty::getMinValue<tgt::ivec4>() const {
    return minIntVal_;
}

template<>
tgt::vec4 VoxelTypeProperty::getMinValue<tgt::vec4>() const {
    return minFloatVal_;
}

// setMinValue
template<>
void VoxelTypeProperty::setMinValue<uint8_t>(const uint8_t& value) {
    minIntVal_.x = static_cast<int>(value);
    if (minIntVal_.x > intVal_.x) {
        intVal_.x = minIntVal_.x;
    }
    invalidate();
}

template<>
void VoxelTypeProperty::setMinValue<uint16_t>(const uint16_t& value) {
    minIntVal_.x = static_cast<int>(value);
    if (minIntVal_.x > intVal_.x) {
        intVal_.x = minIntVal_.x;
    }
    invalidate();
}

template<>
void VoxelTypeProperty::setMinValue<uint32_t>(const uint32_t& value) {
    minIntVal_.x = static_cast<int>(value);
    if (minIntVal_.x > intVal_.x) {
        intVal_.x = minIntVal_.x;
    }
    invalidate();
}

template<>
void VoxelTypeProperty::setMinValue<uint64_t>(const uint64_t& value) {
    minIntVal_.x = static_cast<int>(value);
    if (minIntVal_.x > intVal_.x) {
        intVal_.x = minIntVal_.x;
    }
    invalidate();
}

template<>
void VoxelTypeProperty::setMinValue<int8_t>(const int8_t& value) {
    minIntVal_.x = static_cast<int>(value);
    if (minIntVal_.x > intVal_.x) {
        intVal_.x = minIntVal_.x;
    }
    invalidate();
}

template<>
void VoxelTypeProperty::setMinValue<int16_t>(const int16_t& value) {
    minIntVal_.x = static_cast<int>(value);
    if (minIntVal_.x > intVal_.x) {
        intVal_.x = minIntVal_.x;
    }
    invalidate();
}

template<>
void VoxelTypeProperty::setMinValue<int32_t>(const int32_t& value) {
    minIntVal_.x = static_cast<int>(value);
    if (minIntVal_.x > intVal_.x) {
        intVal_.x = minIntVal_.x;
    }
    invalidate();
}

template<>
void VoxelTypeProperty::setMinValue<int64_t>(const int64_t& value) {
    minIntVal_.x = static_cast<int>(value);
    if (minIntVal_.x > intVal_.x) {
        intVal_.x = minIntVal_.x;
    }
    invalidate();
}

template<>
void VoxelTypeProperty::setMinValue<float>(const float& value) {
    minFloatVal_.x = value;
    if (minFloatVal_.x > floatVal_.x) {
        floatVal_.x = minFloatVal_.x;
    }
    invalidate();
}

template<>
void VoxelTypeProperty::setMinValue<tgt::ivec2>(const tgt::ivec2& value) {
    minIntVal_.x = value.x;
    minIntVal_.y = value.y;
}

template<>
void VoxelTypeProperty::setMinValue<tgt::vec2>(const tgt::vec2& value) {
    minFloatVal_.x = value.x;
    minFloatVal_.y = value.y;
}

template<>
void VoxelTypeProperty::setMinValue<tgt::ivec3>(const tgt::ivec3& value) {
    minIntVal_.x = value.x;
    minIntVal_.y = value.y;
    minIntVal_.z = value.z;
}

template<>
void VoxelTypeProperty::setMinValue<tgt::vec3>(const tgt::vec3& value) {
    minFloatVal_.x = value.x;
    minFloatVal_.y = value.y;
    minFloatVal_.z = value.z;
}

template<>
void VoxelTypeProperty::setMinValue<tgt::ivec4>(const tgt::ivec4& value) {
    minIntVal_ = value;
}

template<>
void VoxelTypeProperty::setMinValue<tgt::vec4>(const tgt::vec4& value) {
    minFloatVal_ = value;
}

// getMaxValue
template<>
int VoxelTypeProperty::getMaxValue<int>() const {
    return maxIntVal_.x;
}

template<>
float VoxelTypeProperty::getMaxValue<float>() const {
    return maxFloatVal_.x;
}

template<>
tgt::ivec2 VoxelTypeProperty::getMaxValue<tgt::ivec2>() const {
    return tgt::ivec2(maxIntVal_.xy());
}

template<>
tgt::vec2 VoxelTypeProperty::getMaxValue<tgt::vec2>() const {
    return tgt::vec2(maxFloatVal_.xy());
}

template<>
tgt::ivec3 VoxelTypeProperty::getMaxValue<tgt::ivec3>() const {
    return tgt::ivec3(maxIntVal_.xyz());
}

template<>
tgt::vec3 VoxelTypeProperty::getMaxValue<tgt::vec3>() const {
    return tgt::vec3(maxFloatVal_.xyz());
}

template<>
tgt::ivec4 VoxelTypeProperty::getMaxValue<tgt::ivec4>() const {
    return maxIntVal_;
}

template<>
tgt::vec4 VoxelTypeProperty::getMaxValue<tgt::vec4>() const {
    return maxFloatVal_;
}

// setMaxValue
template<>
void VoxelTypeProperty::setMaxValue<uint8_t>(const uint8_t& value) {
    maxIntVal_.x = static_cast<int>(value);
}

template<>
void VoxelTypeProperty::setMaxValue<uint16_t>(const uint16_t& value) {
    maxIntVal_.x = static_cast<int>(value);
}

template<>
void VoxelTypeProperty::setMaxValue<uint32_t>(const uint32_t& value) {
    maxIntVal_.x = static_cast<int>(value);
}

template<>
void VoxelTypeProperty::setMaxValue<uint64_t>(const uint64_t& value) {
    maxIntVal_.x = static_cast<int>(value);
}

template<>
void VoxelTypeProperty::setMaxValue<int8_t>(const int8_t& value) {
    maxIntVal_.x = static_cast<int>(value);
}

template<>
void VoxelTypeProperty::setMaxValue<int16_t>(const int16_t& value) {
    maxIntVal_.x = static_cast<int>(value);
}

template<>
void VoxelTypeProperty::setMaxValue<int32_t>(const int32_t& value) {
    maxIntVal_.x = static_cast<int>(value);
}

template<>
void VoxelTypeProperty::setMaxValue<int64_t>(const int64_t& value) {
    maxIntVal_.x = static_cast<int>(value);
}

template<>
void VoxelTypeProperty::setMaxValue<float>(const float& value) {
    maxFloatVal_.x = value;
}

template<>
void VoxelTypeProperty::setMaxValue<tgt::ivec2>(const tgt::ivec2& value) {
    maxIntVal_.x = value.x;
    maxIntVal_.y = value.y;
}

template<>
void VoxelTypeProperty::setMaxValue<tgt::vec2>(const tgt::vec2& value) {
    maxFloatVal_.x = value.x;
    maxFloatVal_.y = value.y;
}

template<>
void VoxelTypeProperty::setMaxValue<tgt::ivec3>(const tgt::ivec3& value) {
    maxIntVal_.x = value.x;
    maxIntVal_.y = value.y;
    maxIntVal_.z = value.z;
}

template<>
void VoxelTypeProperty::setMaxValue<tgt::vec3>(const tgt::vec3& value) {
    maxFloatVal_.x = value.x;
    maxFloatVal_.y = value.y;
    maxFloatVal_.z = value.z;
}

template<>
void VoxelTypeProperty::setMaxValue<tgt::ivec4>(const tgt::ivec4& value) {
    maxIntVal_ = value;
}

template<>
void VoxelTypeProperty::setMaxValue<tgt::vec4>(const tgt::vec4& value) {
    maxFloatVal_ = value;
}

} // namespace voreen
