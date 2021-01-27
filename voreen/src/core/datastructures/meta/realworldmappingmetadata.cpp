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

#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

namespace voreen {

RealWorldMapping::RealWorldMapping() : scale_(1.0f), offset_(0.0f), unit_("") {
}

RealWorldMapping::RealWorldMapping(float scale, float offset, std::string unit
)
    : scale_(scale), offset_(offset), unit_(unit) {
}
RealWorldMapping::RealWorldMapping(tgt::vec2 range, std::string unit)
    : scale_(range.y - range.x)
    , offset_(range.x)
    , unit_(unit)
{
    tgtAssert(range.x <= range.y, "Invalid range for RWM");
}

bool RealWorldMapping::operator==(const RealWorldMapping& other) const {
    return (scale_ == other.scale_ && offset_ == other.offset_ && unit_ == other.unit_);
}

bool RealWorldMapping::operator!=(const RealWorldMapping& other) const {
    return !(*this == other);
}

void RealWorldMapping::serialize(Serializer& s) const {
    s.serialize("scale", scale_);
    s.serialize("offset", offset_);
    s.serialize("unit", unit_);
}

void RealWorldMapping::deserialize(Deserializer& s) {
    s.deserialize("scale", scale_);
    s.deserialize("offset", offset_);
    s.deserialize("unit", unit_);
}

float RealWorldMapping::normalizedToRealWorld(float normalized) const {
    return (normalized * scale_) + offset_;
}

float RealWorldMapping::realWorldToNormalized(float realWorld) const {
    return (realWorld - offset_) / scale_;
}

std::string RealWorldMapping::getUnit() const {
    return unit_;
}

void RealWorldMapping::setUnit(std::string unit) {
    unit_ = unit;
}

float RealWorldMapping::getScale() const {
    return scale_;
}

void RealWorldMapping::setScale(float scale) {
    scale_ = scale;
}

float RealWorldMapping::getOffset() const {
    return offset_;
}

void RealWorldMapping::setOffset(float offset) {
    offset_ = offset;
}

RealWorldMapping RealWorldMapping::createDenormalizingMapping(const std::string& format) {
    if(format == "uint8") {
        return RealWorldMapping::createDenormalizingMapping<uint8_t>();
    } else if(format == "uint16") {
        return RealWorldMapping::createDenormalizingMapping<uint16_t>();
    } else if(format == "uint32") {
        return RealWorldMapping::createDenormalizingMapping<uint32_t>();
    } else if(format == "uint64") {
        return RealWorldMapping::createDenormalizingMapping<uint64_t>();
    } else if(format == "int8") {
        return RealWorldMapping::createDenormalizingMapping<int8_t>();
    } else if(format == "int16") {
        return RealWorldMapping::createDenormalizingMapping<int16_t>();
    } else if(format == "int32") {
        return RealWorldMapping::createDenormalizingMapping<int32_t>();
    } else if(format == "int64") {
        return RealWorldMapping::createDenormalizingMapping<int64_t>();
    } else if(format == "float") {
        return RealWorldMapping::createDenormalizingMapping<float>();
    } else if(format == "double") {
        return RealWorldMapping::createDenormalizingMapping<double>();
    } else {
        throw tgt::UnsupportedFormatException("Cannot generate real world mapping for volume format \"" + format + "\".");
    }
}

RealWorldMapping RealWorldMapping::getInverseMapping() const {
    return RealWorldMapping(1.f/scale_, -offset_/scale_, "1/" + unit_);
}

//--------------------------------------------------------------------------------------

RealWorldMappingMetaData::RealWorldMappingMetaData(float scale, float offset, std::string unit) : TemplateMetaData<RealWorldMapping>(RealWorldMapping(scale, offset, unit))
{}

std::string RealWorldMappingMetaData::toString() const {
    std::stringstream s;
    s << getValue();
    return s.str();
}

std::string RealWorldMappingMetaData::toString(const std::string& component) const {
    if (component == "offset")
        return ftos(getValue().getOffset());
    else if (component == "scale")
        return ftos(getValue().getScale());
    else if (component == "unit")
        return getValue().getUnit();
    else
        return toString();
}

} // namespace
