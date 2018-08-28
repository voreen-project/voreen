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

#ifndef VRN_PORTCONDITION_VOLUMETYPE_H
#define VRN_PORTCONDITION_VOLUMETYPE_H

#include "voreen/core/ports/conditions/portcondition.h"
#include "voreen/core/ports/volumeport.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

/**
 * Base class for port conditions that check for a specific volume type.
 */
class VRN_CORE_API PortConditionVolumeType : public PortCondition {
public:
    PortConditionVolumeType(const std::string& typeString, const std::string& volTypeName);

    virtual ~PortConditionVolumeType();

    /**
     * Returns true, if the inport has data and
     * the encapsulated volume is of type T.
     */
    virtual bool acceptsPortData() const;

protected:
    virtual void setCheckedPort(const Port* checkedPort);
    const std::string typeString_;
    const VolumePort* volumePort_;
};

// ------------------------------------------------------------------------------------------------

/**
 * Port condition that allows OpenGL-compatible volume types,
 * i.e., volumes that can be converted into an OpenGL texture.
 * These are all standard volume types except 64 integer and double:
 * - VolumeRAM_UInt8, VolumeRAM_Int8
 * - VolumeRAM_UInt16, VolumeRAM_Int16
 * - VolumeRAM_UInt32, VolumeRAM_Int32
 * - VolumeRAM_Float
 * - VolumeRAM_2xUInt8, VolumeRAM_2xInt8
 * - VolumeRAM_2xUInt16, VolumeRAM_2xInt16
 * - VolumeRAM_2xUInt32, VolumeRAM_2xInt32
 * - VolumeRAM_2xFloat
 * - VolumeRAM_3xUInt8, VolumeRAM_3xInt8
 * - VolumeRAM_3xUInt16, VolumeRAM_3xInt16
 * - VolumeRAM_3xUInt32, VolumeRAM_3xInt32
 * - VolumeRAM_3xFloat
 * - VolumeRAM_4xUInt8, VolumeRAM_4xInt8
 * - VolumeRAM_4xUInt16, VolumeRAM_4xInt16
 * - VolumeRAM_4xUInt32, VolumeRAM_4xInt32
 * - VolumeRAM_4xFloat
 */
class VRN_CORE_API PortConditionVolumeTypeGL : public PortCondition {
public:
    PortConditionVolumeTypeGL();

    virtual bool acceptsPortData() const;

protected:
    virtual void setCheckedPort(const Port* checkedPort);
    PortConditionLogicalOr orCondition_;
};

// ------------------------------------------------------------------------------------------------

/**
 * Port condition that allows float and double volume types:
 * - VolumeRAM_Float, VolumeRAM_Double
 * - VolumeRAM_2xFloat, VolumeRAM_2xDouble
 * - VolumeRAM_3xFloat, VolumeRAM_3xDouble
 * - VolumeRAM_4xFloat, VolumeRAM_4xDouble
 */
class VRN_CORE_API PortConditionVolumeTypeReal : public PortCondition {
public:
    PortConditionVolumeTypeReal();

    virtual bool acceptsPortData() const;

protected:
    virtual void setCheckedPort(const Port* checkedPort);
    PortConditionLogicalOr orCondition_;
};

// ------------------------------------------------------------------------------------------------

/**
 * Port condition that allows integer volume types:
 * - VolumeRAM_UInt8, VolumeRAM_Int8
 * - VolumeRAM_UInt16, VolumeRAM_Int16
 * - VolumeRAM_UInt32, VolumeRAM_Int32
 * - VolumeRAM_UInt64, VolumeRAM_Int64
 * - VolumeRAM_2xUInt8, VolumeRAM_2xInt8
 * - VolumeRAM_2xUInt16, VolumeRAM_2xInt16
 * - VolumeRAM_2xUInt32, VolumeRAM_2xInt32
 * - VolumeRAM_2xUInt64, VolumeRAM_2xInt64
 * - VolumeRAM_3xUInt8, VolumeRAM_3xInt8
 * - VolumeRAM_3xUInt16, VolumeRAM_3xInt16
 * - VolumeRAM_3xUInt32, VolumeRAM_3xInt32
 * - VolumeRAM_3xUInt64, VolumeRAM_3xInt64
 * - VolumeRAM_4xUInt8, VolumeRAM_4xInt8
 * - VolumeRAM_4xUInt16, VolumeRAM_4xInt16
 * - VolumeRAM_4xUInt32, VolumeRAM_4xInt32
 * - VolumeRAM_4xUInt64, VolumeRAM_4xInt64
 */
class VRN_CORE_API PortConditionVolumeTypeInteger : public PortCondition {
public:
    PortConditionVolumeTypeInteger();

    virtual bool acceptsPortData() const;

protected:
    virtual void setCheckedPort(const Port* checkedPort);
    PortConditionLogicalOr orCondition_;
};

// ------------------------------------------------------------------------------------------------

/**
 * Port condition that allows volumes with a specific number of channels.
 */
class VRN_CORE_API PortConditionVolumeChannelCount : public PortCondition {
public:
    PortConditionVolumeChannelCount(size_t numChannels);

    virtual bool acceptsPortData() const;

    virtual std::string getDescription() const;

protected:
    virtual void setCheckedPort(const Port* checkedPort);

    const VolumePort* volumePort_;
    size_t numChannels_;
};

// ------------------------------------------------------------------------------------------------

/**
 * Port condition that allows volumes with a specific value range.
 * Note: std::to_string() must have an implementation for type T.
 * Note: Needs a VolumeRAM representation to work properly.
 */
template<typename T>
class VRN_CORE_API PortConditionVolumeValueRange : public PortCondition {
public:
    PortConditionVolumeValueRange(const T& minValue, const T& maxValue);

    virtual bool acceptsPortData() const;

protected:
    virtual void setCheckedPort(const Port* checkedPort);

    const VolumePort* volumePort_;
    T minValue_, maxValue_;
};

// ------------------------------------------------------------------------------------------------
// implementation of PortConditionVolumeValueRange

template <typename T>
PortConditionVolumeValueRange<T>::PortConditionVolumeValueRange(const T& minValue, const T& maxValue)
    : PortCondition("Volume with values in range [" + std::to_string(minValue) + ", " + std::to_string(maxValue) + "] expected")
    , minValue_(minValue)
    , maxValue_(maxValue)
{
}

template <typename T>
bool PortConditionVolumeValueRange<T>::acceptsPortData() const
{
    if (!volumePort_ || !volumePort_->hasData())
        return false;

    const VolumeAtomic<T>* volume = dynamic_cast<const VolumeAtomic<T>*>(volumePort_->getData()->getRepresentation<VolumeRAM>());
    if (!volume)
        return false;

    return (volume->min() >= minValue_ && volume->max() <= maxValue_);
}

template <typename T>
void PortConditionVolumeValueRange<T>::setCheckedPort(const Port* checkedPort)
{
    if (!dynamic_cast<const VolumePort*>(checkedPort)) {
        LERRORC("voreen.PortConditionVolumeValueRange", "Assigned port is not a volume port");
    }
    else {
        volumePort_ = static_cast<const VolumePort*>(checkedPort);
    }
}

// ------------------------------------------------------------------------------------------------
// concrete subtypes of PortConditionVolumeType

// single channel

class VRN_CORE_API PortConditionVolumeTypeUInt8 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeUInt8() :
      PortConditionVolumeType("uint8","Volume_UInt8") {};
};

class VRN_CORE_API PortConditionVolumeTypeInt8 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeInt8() :
      PortConditionVolumeType("int8","Volume_Int8") {};
};

class VRN_CORE_API PortConditionVolumeTypeUInt12 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeUInt12() :
      PortConditionVolumeType("uint12","Volume_UInt12") {};
};

class VRN_CORE_API PortConditionVolumeTypeInt12 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeInt12() :
      PortConditionVolumeType("int12","Volume_Int12") {};
};

class VRN_CORE_API PortConditionVolumeTypeUInt16 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeUInt16() :
      PortConditionVolumeType("uint16","Volume_UInt16") {};
};

class VRN_CORE_API PortConditionVolumeTypeInt16 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeInt16() :
      PortConditionVolumeType("int16","Volume_Int16") {};
};

class VRN_CORE_API PortConditionVolumeTypeUInt32 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeUInt32() :
      PortConditionVolumeType("uint32","Volume_UInt32") {};
};

class VRN_CORE_API PortConditionVolumeTypeInt32 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeInt32() :
      PortConditionVolumeType("int32","Volume_Int32") {};
};

class VRN_CORE_API PortConditionVolumeTypeUInt64 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeUInt64() :
      PortConditionVolumeType("uint64","Volume_UInt64") {};
};

class VRN_CORE_API PortConditionVolumeTypeInt64 : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeInt64() :
      PortConditionVolumeType("int64","Volume_Int64") {};
};

class VRN_CORE_API PortConditionVolumeTypeFloat : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeFloat() :
      PortConditionVolumeType("float","Volume_Float") {};
};

class VRN_CORE_API PortConditionVolumeTypeDouble : public PortConditionVolumeType {
public:
    PortConditionVolumeTypeDouble() :
      PortConditionVolumeType("double","Volume_Double") {};
};

// two channels

class VRN_CORE_API PortConditionVolumeType2xUInt8 : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xUInt8() :
      PortConditionVolumeType("Vector2(uint8)","Volume_2xUInt8") {};
};

class VRN_CORE_API PortConditionVolumeType2xInt8 : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xInt8() :
      PortConditionVolumeType("Vector2(int8)","Volume_2xInt8") {};
};

class VRN_CORE_API PortConditionVolumeType2xUInt16 : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xUInt16() :
      PortConditionVolumeType("Vector2(uint16)","Volume_2xUInt16") {};
};

class VRN_CORE_API PortConditionVolumeType2xInt16 : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xInt16() :
      PortConditionVolumeType("Vector2(int16)","Volume_2xInt16") {};
};

class VRN_CORE_API PortConditionVolumeType2xUInt32 : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xUInt32() :
      PortConditionVolumeType("Vector2(uint32)","Volume_2xUInt32") {};
};

class VRN_CORE_API PortConditionVolumeType2xInt32 : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xInt32() :
      PortConditionVolumeType("Vector2(int32)","Volume_2xInt32") {};
};

class VRN_CORE_API PortConditionVolumeType2xUInt64 : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xUInt64() :
      PortConditionVolumeType("Vector2(uint64)","Volume_2xUInt64") {};
};

class VRN_CORE_API PortConditionVolumeType2xInt64 : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xInt64() :
      PortConditionVolumeType("Vector2(int64)","Volume_2xInt64") {};
};

class VRN_CORE_API PortConditionVolumeType2xFloat : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xFloat() :
      PortConditionVolumeType("Vector2(float)","Volume_2xFloat") {};
};

class VRN_CORE_API PortConditionVolumeType2xDouble : public PortConditionVolumeType {
public:
    PortConditionVolumeType2xDouble() :
      PortConditionVolumeType("Vector2(double)","Volume_2xDouble") {};
};

// three channels

class VRN_CORE_API PortConditionVolumeType3xUInt8 : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xUInt8() :
      PortConditionVolumeType("Vector3(uint8)","Volume_3xUInt8") {};
};

class VRN_CORE_API PortConditionVolumeType3xInt8 : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xInt8() :
      PortConditionVolumeType("Vector3(int8)","Volume_3xInt8") {};
};

class VRN_CORE_API PortConditionVolumeType3xUInt16 : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xUInt16() :
      PortConditionVolumeType("Vector3(uint16)","Volume_3xUInt16") {};
};

class VRN_CORE_API PortConditionVolumeType3xInt16 : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xInt16() :
      PortConditionVolumeType("Vector3(int16)","Volume_3xInt16") {};
};

class VRN_CORE_API PortConditionVolumeType3xUInt32 : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xUInt32() :
      PortConditionVolumeType("Vector3(uint32)","Volume_3xUInt32") {};
};

class VRN_CORE_API PortConditionVolumeType3xInt32 : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xInt32() :
      PortConditionVolumeType("Vector3(int32)","Volume_3xInt32") {};
};

class VRN_CORE_API PortConditionVolumeType3xUInt64 : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xUInt64() :
      PortConditionVolumeType("Vector3(uint64)","Volume_3xUInt64") {};
};

class VRN_CORE_API PortConditionVolumeType3xInt64 : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xInt64() :
      PortConditionVolumeType("Vector3(int64)","Volume_3xInt64") {};
};

class VRN_CORE_API PortConditionVolumeType3xFloat : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xFloat() :
      PortConditionVolumeType("Vector3(float)","Volume_3xFloat") {};
};

class VRN_CORE_API PortConditionVolumeType3xDouble : public PortConditionVolumeType {
public:
    PortConditionVolumeType3xDouble() :
      PortConditionVolumeType("Vector3(double)","Volume_3xDouble") {};
};

// four channels

class VRN_CORE_API PortConditionVolumeType4xUInt8 : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xUInt8() :
      PortConditionVolumeType("Vector4(uint8)","Volume_4xUInt8") {};
};

class VRN_CORE_API PortConditionVolumeType4xInt8 : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xInt8() :
      PortConditionVolumeType("Vector4(int8)","Volume_4xInt8") {};
};

class VRN_CORE_API PortConditionVolumeType4xUInt16 : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xUInt16() :
      PortConditionVolumeType("Vector4(uint16)","Volume_4xUInt16") {};
};

class VRN_CORE_API PortConditionVolumeType4xInt16 : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xInt16() :
      PortConditionVolumeType("Vector4(int16)","Volume_4xInt16") {};
};

class VRN_CORE_API PortConditionVolumeType4xUInt32 : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xUInt32() :
      PortConditionVolumeType("Vector4(uint32)","Volume_4xUInt32") {};
};

class VRN_CORE_API PortConditionVolumeType4xInt32 : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xInt32() :
      PortConditionVolumeType("Vector4(int32)","Volume_4xInt32") {};
};

class VRN_CORE_API PortConditionVolumeType4xUInt64 : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xUInt64() :
      PortConditionVolumeType("Vector4(uint64)","Volume_4xUInt64") {};
};

class VRN_CORE_API PortConditionVolumeType4xInt64 : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xInt64() :
      PortConditionVolumeType("Vector4(int64)","Volume_4xInt64") {};
};

class VRN_CORE_API PortConditionVolumeType4xFloat : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xFloat() :
      PortConditionVolumeType("Vector4(float)","Volume_4xFloat") {};
};

class VRN_CORE_API PortConditionVolumeType4xDouble : public PortConditionVolumeType {
public:
    PortConditionVolumeType4xDouble() :
      PortConditionVolumeType("Vector4(double)","Volume_4xDouble") {};
};

// ------------------------------------------------------------------------------------------------
// concrete subtypes of PortConditionVolumeType


} // namespace

#endif // VRN_PORTCONDITION_VOLUMETYPE_H
