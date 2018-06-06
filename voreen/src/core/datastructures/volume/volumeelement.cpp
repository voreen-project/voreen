/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/core/datastructures/volume/volumeelement.h"

namespace voreen {

template<>
float getTypeAsFloat(float value) {
    return value;
}

template<>
float getTypeAsFloat(double value) {
    return static_cast<float>(value);
}

template<>
float getTypeAsFloat(int8_t value) {
    if (value >= 0)
        return static_cast<float>(value) / VolumeElement<int8_t>::rangeMaxElement();
    else
        return static_cast<float>(value) / -VolumeElement<int8_t>::rangeMinElement();
}

template<>
float getTypeAsFloat(uint8_t value) {
    return static_cast<float>(value) / VolumeElement<uint8_t>::rangeMaxElement();
}

template<>
float getTypeAsFloat(int16_t value) {
    if (value >= 0)
        return static_cast<float>(value) / VolumeElement<int16_t>::rangeMaxElement();
    else
        return static_cast<float>(value) / -VolumeElement<int16_t>::rangeMinElement();
}

template<>
float getTypeAsFloat(uint16_t value) {
    return static_cast<float>(value) / VolumeElement<uint16_t>::rangeMaxElement();
}

template<>
float getTypeAsFloat(int32_t value) {
    if (value >= 0)
        return static_cast<float>(value) / VolumeElement<int32_t>::rangeMaxElement();
    else
        return static_cast<float>(value) / -VolumeElement<int32_t>::rangeMinElement();
}

template<>
float getTypeAsFloat(uint32_t value) {
    return static_cast<float>(value) / VolumeElement<uint32_t>::rangeMaxElement();
}

template<>
float getTypeAsFloat(int64_t value) {
    if (value >= 0)
        return static_cast<float>(value) / VolumeElement<int64_t>::rangeMaxElement();
    else
        return static_cast<float>(value) / -VolumeElement<int64_t>::rangeMinElement();
}

template<>
float getTypeAsFloat(uint64_t value) {
    return static_cast<float>(value) / VolumeElement<uint64_t>::rangeMaxElement();
}

//---------------------------------------------------------------------------------------

template<>
float getFloatAsType(float value) {
    return value;
}
template<>
double getFloatAsType(float value) {
    return static_cast<double>(value);
}

template<>
int8_t getFloatAsType(float value) {
    value = tgt::clamp(value, -1.0f, 1.0f);
    if(value >= 0.0f)
        return static_cast<int8_t>(value*VolumeElement<int8_t>::rangeMaxElement());
    else
        return static_cast<int8_t>(value*-VolumeElement<int8_t>::rangeMinElement());
}

template<>
uint8_t getFloatAsType(float value) {
    value = tgt::clamp(value, 0.0f, 1.0f);
    return static_cast<uint8_t>(value*VolumeElement<uint8_t>::rangeMaxElement());
}

template<>
int16_t getFloatAsType(float value) {
    value = tgt::clamp(value, -1.0f, 1.0f);
    if(value >= 0.0f)
        return static_cast<int16_t>(value*VolumeElement<int16_t>::rangeMaxElement());
    else
        return static_cast<int16_t>(value*-VolumeElement<int16_t>::rangeMinElement());
}
template<>
uint16_t getFloatAsType(float value) {
    value = tgt::clamp(value, 0.0f, 1.0f);
    return static_cast<uint16_t>(value*VolumeElement<uint16_t>::rangeMaxElement());
}

template<>
int32_t getFloatAsType(float value) {
    value = tgt::clamp(value, -1.0f, 1.0f);
    // int32_MAX does not fit into a float => We lose precision and an input value of 1.0f is mapped onto int32_MAX+1 == 0
    // Therefore we cast value to double first.
    if(value >= 0.0f)
        return static_cast<int32_t>(static_cast<double>(value)*std::numeric_limits<int32_t>::max());
    else
        return static_cast<int32_t>(static_cast<double>(-value)*std::numeric_limits<int32_t>::lowest());
}
template<>
uint32_t getFloatAsType(float value) {
    value = tgt::clamp(value, 0.0f, 1.0f);
    // uint32_MAX does not fit into a float => We lose precision and an input value of 1.0f is mapped onto uint32_MAX+1 == 0
    // Therefore we cast value to double first.
    return static_cast<uint32_t>(static_cast<double>(value)*std::numeric_limits<uint32_t>::max());
}

template<>
int64_t getFloatAsType(float value) {
    value = tgt::clamp(value, -1.0f, 1.0f);
    if(value >= 0.0f)
        return static_cast<int64_t>(value*VolumeElement<int64_t>::rangeMaxElement());
    else
        return static_cast<int64_t>(-value*VolumeElement<int64_t>::rangeMinElement());
}
template<>
uint64_t getFloatAsType(float value) {
    value = tgt::clamp(value, 0.0f, 1.0f);
    return static_cast<uint64_t>(value*VolumeElement<uint64_t>::rangeMaxElement());
}

} // namespace voreen
