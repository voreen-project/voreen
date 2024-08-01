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

#ifndef VRN_VOXELTYPEPROPERTY_H
#define VRN_VOXELTYPEPROPERTY_H

#include "voreen/core/properties/property.h"
#include "voreen/core/properties/templateproperty.h"
#include "voreen/core/properties/propertyowner.h"
#include "tgt/types.h"

namespace voreen {

class VolumeBase;

/**
 * Property type to be used for parameters that depend on the voxel type (data type + channel count)
 * of a volume. The property adapts its used data type and its bounds to the voxel type of the
 * assigned volume.
 *
 * @note You have to call setValue<T>/getValue<T>() with the current volume's voxel type as template parameter
 *  to assign and retrieve the correct property value. Using a mis-matching template parameter yields undefined
 *  results. The inheritance from TemplateProperty is for technical reasons only, so do not use set()/get().
 *
 */
class VRN_CORE_API VoxelTypeProperty : public TemplateProperty<void*> {
public:
    VoxelTypeProperty(const std::string& id, const std::string& guiText,
        int invalidationLevel=Processor::INVALID_RESULT, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    VoxelTypeProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "VoxelTypeProperty"; }
    virtual std::string getTypeDescription() const { return "VoxelType"; }

    virtual std::string getGuiName() const;
    virtual std::string getID() const;
    virtual void setGuiName(const std::string& guiName);

    /**
     * Causes the property to adapt to the voxel type of the passed volume.
     *
     * @param volumeHandle the volume to adapt to
     * @param adjustMinMaxRange if true, the property's min/max range is adapted to
     *  the element range of the passed volume.
     */
    void setVolume(const VolumeBase* volumeHandle, bool adjustMinMaxRange = true);
    const VolumeBase* getVolume() const;

    template<typename T>
    void setValue(const T& value);

    template<typename T>
    T getValue() const;

    /**
     * Sets the minimum value the variable can take.
     *
     * @note min/max ranges are automatically updated
     *  on setVolume(), by default.
     */
    template<typename T>
    void setMinValue(const T& minValue);
    template<typename T>
    T getMinValue() const;

    /**
     * Sets the maximum value the variable can take.
     *
     * @note min/max ranges are automatically updated
     *  on setVolume(), by default.
     */
    template<typename T>
    void setMaxValue(const T& maxValue);
    template<typename T>
    T getMaxValue() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:
    /// Updates valid range according to the current voxel type.
    void updateMinMaxRange();

    /// The volume the property is currently adapted to.
    const VolumeBase* volume_;

    /*
     * The following members are used as replacement for dynamic typing:
     * - intVal_ is used for storing the property value, if the current voxel type
     *   is an integer type (uin8_t, uin16_t, 3xuint8_t, ...)
     * - floatVal_ is used for storing the property value, if the current voxel type
     *   is a float type (float, double, 3xfloat, ...)
     * If the current voxel type has less than four channels, the redundant channels are ignored.
     * E.g.: for voxel type uint8_t only the r-channel of intVal_ is used,
     * for voxel type 3xfloat the rgb-channels of floatVal_ are used.
     */
    tgt::ivec4 intVal_;
    tgt::vec4 floatVal_;

    /// min/max values (see above)
    tgt::ivec4 minIntVal_;
    tgt::ivec4 maxIntVal_;
    tgt::vec4 minFloatVal_;
    tgt::vec4 maxFloatVal_;

    static const std::string loggerCat_;
};


// template dummy definitions (see .cpp file for specializations)
// --------------------------------------------------------------
template<typename T>
T VoxelTypeProperty::getValue() const {
    LWARNING("get() unknown template parameter: " << typeid(T).name());
    return T();
}

template<typename T>
void VoxelTypeProperty::setValue(const T& /*value*/) {
    LWARNING("get() unknown template parameter: " << typeid(T).name());
}

template<typename T>
T VoxelTypeProperty::getMinValue() const {
    LWARNING("getMinValue() unknown template parameter: " << typeid(T).name());
    return T();
}

template<typename T>
void VoxelTypeProperty::setMinValue(const T& /*value*/) {
    LWARNING("setMinValue() unknown template parameter: " << typeid(T).name());
}

template<typename T>
T VoxelTypeProperty::getMaxValue() const {
    LWARNING("getMaxValue() unknown template parameter: " << typeid(T).name());
    return T();
}

template<typename T>
void VoxelTypeProperty::setMaxValue(const T& /*value*/) {
    LWARNING("setMaxValue() unknown template paramater: " << typeid(T).name());
}


// template specializations
// ------------------------

//
// get (scalar)
//
template<>
VRN_CORE_API uint8_t VoxelTypeProperty::getValue<uint8_t>() const;

template<>
VRN_CORE_API uint16_t VoxelTypeProperty::getValue<uint16_t>() const;

template<>
VRN_CORE_API uint32_t VoxelTypeProperty::getValue<uint32_t>() const;

template<>
VRN_CORE_API uint64_t VoxelTypeProperty::getValue<uint64_t>() const;

template<>
VRN_CORE_API int8_t VoxelTypeProperty::getValue<int8_t>() const;

template<>
VRN_CORE_API int16_t VoxelTypeProperty::getValue<int16_t>() const;

template<>
VRN_CORE_API int32_t VoxelTypeProperty::getValue<int32_t>() const;

template<>
VRN_CORE_API int64_t VoxelTypeProperty::getValue<int64_t>() const;

template<>
VRN_CORE_API float VoxelTypeProperty::getValue<float>() const;

template<>
VRN_CORE_API double VoxelTypeProperty::getValue<double>() const;

// get (vec2)
template<>
VRN_CORE_API tgt::Vector2<uint8_t> VoxelTypeProperty::getValue<tgt::Vector2<uint8_t> >() const;

template<>
VRN_CORE_API tgt::Vector2<uint16_t> VoxelTypeProperty::getValue<tgt::Vector2<uint16_t> >() const;

template<>
VRN_CORE_API tgt::Vector2<uint32_t> VoxelTypeProperty::getValue<tgt::Vector2<uint32_t> >() const;

template<>
VRN_CORE_API tgt::Vector2<uint64_t> VoxelTypeProperty::getValue<tgt::Vector2<uint64_t> >() const;

template<>
VRN_CORE_API tgt::Vector2<int8_t> VoxelTypeProperty::getValue<tgt::Vector2<int8_t> >() const;

template<>
VRN_CORE_API tgt::Vector2<int16_t> VoxelTypeProperty::getValue<tgt::Vector2<int16_t> >() const;

template<>
VRN_CORE_API tgt::Vector2<int32_t> VoxelTypeProperty::getValue<tgt::Vector2<int32_t> >() const;

template<>
VRN_CORE_API tgt::Vector2<int64_t> VoxelTypeProperty::getValue<tgt::Vector2<int64_t> >() const;

template<>
VRN_CORE_API tgt::vec2 VoxelTypeProperty::getValue<tgt::vec2>() const;

template<>
VRN_CORE_API tgt::dvec2 VoxelTypeProperty::getValue<tgt::dvec2>() const;

// get (vec3)
template<>
VRN_CORE_API tgt::Vector3<uint8_t> VoxelTypeProperty::getValue<tgt::Vector3<uint8_t> >() const;

template<>
VRN_CORE_API tgt::Vector3<uint16_t> VoxelTypeProperty::getValue<tgt::Vector3<uint16_t> >() const;

template<>
VRN_CORE_API tgt::Vector3<uint32_t> VoxelTypeProperty::getValue<tgt::Vector3<uint32_t> >() const;

template<>
VRN_CORE_API tgt::Vector3<uint64_t> VoxelTypeProperty::getValue<tgt::Vector3<uint64_t> >() const;

template<>
VRN_CORE_API tgt::Vector3<int8_t> VoxelTypeProperty::getValue<tgt::Vector3<int8_t> >() const;

template<>
VRN_CORE_API tgt::Vector3<int16_t> VoxelTypeProperty::getValue<tgt::Vector3<int16_t> >() const;

template<>
VRN_CORE_API tgt::Vector3<int32_t> VoxelTypeProperty::getValue<tgt::Vector3<int32_t> >() const;

template<>
VRN_CORE_API tgt::Vector3<int64_t> VoxelTypeProperty::getValue<tgt::Vector3<int64_t> >() const;

template<>
VRN_CORE_API tgt::vec3 VoxelTypeProperty::getValue<tgt::vec3>() const;

template<>
VRN_CORE_API tgt::dvec3 VoxelTypeProperty::getValue<tgt::dvec3>() const;

// get (vec4)
template<>
VRN_CORE_API tgt::Vector4<uint8_t> VoxelTypeProperty::getValue<tgt::Vector4<uint8_t> >() const;

template<>
VRN_CORE_API tgt::Vector4<uint16_t> VoxelTypeProperty::getValue<tgt::Vector4<uint16_t> >() const;

template<>
VRN_CORE_API tgt::Vector4<uint32_t> VoxelTypeProperty::getValue<tgt::Vector4<uint32_t> >() const;

template<>
VRN_CORE_API tgt::Vector4<uint64_t> VoxelTypeProperty::getValue<tgt::Vector4<uint64_t> >() const;

template<>
VRN_CORE_API tgt::Vector4<int8_t> VoxelTypeProperty::getValue<tgt::Vector4<int8_t> >() const;

template<>
VRN_CORE_API tgt::Vector4<int16_t> VoxelTypeProperty::getValue<tgt::Vector4<int16_t> >() const;

template<>
VRN_CORE_API tgt::Vector4<int32_t> VoxelTypeProperty::getValue<tgt::Vector4<int32_t> >() const;

template<>
VRN_CORE_API tgt::Vector4<int64_t> VoxelTypeProperty::getValue<tgt::Vector4<int64_t> >() const;

template<>
VRN_CORE_API tgt::vec4 VoxelTypeProperty::getValue<tgt::vec4>() const;

template<>
VRN_CORE_API tgt::dvec4 VoxelTypeProperty::getValue<tgt::dvec4>() const;

//
// set (scalar)
//
template<>
VRN_CORE_API void VoxelTypeProperty::setValue<uint8_t>(const uint8_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<uint16_t>(const uint16_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<uint32_t>(const uint32_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<uint64_t>(const uint64_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<int8_t>(const int8_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<int16_t>(const int16_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<int32_t>(const int32_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<int64_t>(const int64_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<float>(const float& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<double>(const double& value);

// set (vec2)
template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector2<uint8_t> >(const tgt::Vector2<uint8_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector2<uint16_t> >(const tgt::Vector2<uint16_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector2<uint32_t> >(const tgt::Vector2<uint32_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector2<uint64_t> >(const tgt::Vector2<uint64_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector2<int8_t> >(const tgt::Vector2<int8_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector2<int16_t> >(const tgt::Vector2<int16_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector2<int32_t> >(const tgt::Vector2<int32_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector2<int64_t> >(const tgt::Vector2<int64_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::vec2>(const tgt::vec2& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::dvec2>(const tgt::dvec2& value);

// set (vec3)
template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector3<uint8_t> >(const tgt::Vector3<uint8_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector3<uint16_t> >(const tgt::Vector3<uint16_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector3<uint32_t> >(const tgt::Vector3<uint32_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector3<uint64_t> >(const tgt::Vector3<uint64_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector3<int8_t> >(const tgt::Vector3<int8_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector3<int16_t> >(const tgt::Vector3<int16_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector3<int32_t> >(const tgt::Vector3<int32_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector3<int64_t> >(const tgt::Vector3<int64_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::vec3>(const tgt::vec3& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::dvec3>(const tgt::dvec3& value);

// set (vec4)
template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector4<uint8_t> >(const tgt::Vector4<uint8_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector4<uint16_t> >(const tgt::Vector4<uint16_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector4<uint32_t> >(const tgt::Vector4<uint32_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector4<uint64_t> >(const tgt::Vector4<uint64_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector4<int8_t> >(const tgt::Vector4<int8_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector4<int16_t> >(const tgt::Vector4<int16_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector4<int32_t> >(const tgt::Vector4<int32_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::Vector4<int64_t> >(const tgt::Vector4<int64_t>& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::vec4>(const tgt::vec4& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setValue<tgt::dvec4>(const tgt::dvec4& value);

// getMinValue
template<>
VRN_CORE_API int VoxelTypeProperty::getMinValue<int>() const;

template<>
VRN_CORE_API float VoxelTypeProperty::getMinValue<float>() const;

template<>
VRN_CORE_API tgt::ivec2 VoxelTypeProperty::getMinValue<tgt::ivec2>() const;

template<>
VRN_CORE_API tgt::vec2 VoxelTypeProperty::getMinValue<tgt::vec2>() const;

template<>
VRN_CORE_API tgt::ivec3 VoxelTypeProperty::getMinValue<tgt::ivec3>() const;

template<>
VRN_CORE_API tgt::vec3 VoxelTypeProperty::getMinValue<tgt::vec3>() const;

template<>
VRN_CORE_API tgt::ivec4 VoxelTypeProperty::getMinValue<tgt::ivec4>() const;

template<>
VRN_CORE_API tgt::vec4 VoxelTypeProperty::getMinValue<tgt::vec4>() const;

// setMinValue
template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<uint8_t>(const uint8_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<uint16_t>(const uint16_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<uint32_t>(const uint32_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<uint64_t>(const uint64_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<int8_t>(const int8_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<int16_t>(const int16_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<int32_t>(const int32_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<int64_t>(const int64_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<float>(const float& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<tgt::ivec2>(const tgt::ivec2& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<tgt::vec2>(const tgt::vec2& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<tgt::ivec3>(const tgt::ivec3& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<tgt::vec3>(const tgt::vec3& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<tgt::ivec4>(const tgt::ivec4& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMinValue<tgt::vec4>(const tgt::vec4& value);

// getMaxValue
template<>
VRN_CORE_API int VoxelTypeProperty::getMaxValue<int>() const;

template<>
VRN_CORE_API float VoxelTypeProperty::getMaxValue<float>() const;

template<>
VRN_CORE_API tgt::ivec2 VoxelTypeProperty::getMaxValue<tgt::ivec2>() const;

template<>
VRN_CORE_API tgt::vec2 VoxelTypeProperty::getMaxValue<tgt::vec2>() const;

template<>
VRN_CORE_API tgt::ivec3 VoxelTypeProperty::getMaxValue<tgt::ivec3>() const;

template<>
VRN_CORE_API tgt::vec3 VoxelTypeProperty::getMaxValue<tgt::vec3>() const;

template<>
VRN_CORE_API tgt::ivec4 VoxelTypeProperty::getMaxValue<tgt::ivec4>() const;

template<>
VRN_CORE_API tgt::vec4 VoxelTypeProperty::getMaxValue<tgt::vec4>() const;

// setMaxValue
template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<uint8_t>(const uint8_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<uint16_t>(const uint16_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<uint32_t>(const uint32_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<uint64_t>(const uint64_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<int8_t>(const int8_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<int16_t>(const int16_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<int32_t>(const int32_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<int64_t>(const int64_t& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<float>(const float& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<tgt::ivec2>(const tgt::ivec2& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<tgt::vec2>(const tgt::vec2& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<tgt::ivec3>(const tgt::ivec3& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<tgt::vec3>(const tgt::vec3& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<tgt::ivec4>(const tgt::ivec4& value);

template<>
VRN_CORE_API void VoxelTypeProperty::setMaxValue<tgt::vec4>(const tgt::vec4& value);

}  // namespace voreen

#endif // VRN_VOXELTYPEPROPERTY_H
