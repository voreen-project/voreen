/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

FloatVec2Property::FloatVec2Property(const std::string& id, const std::string& guiText,
                             const tgt::vec2& value, const tgt::vec2& minValue,
                             const tgt::vec2& maxValue, int invalidationLevel, NumericProperty<tgt::vec2>::BoundaryUpdatePolicy bup, Property::LevelOfDetail lod)
    : NumericProperty<tgt::vec2>(id, guiText, value, minValue, maxValue, tgt::vec2(0.01f), invalidationLevel, bup, lod)
{}

FloatVec2Property::FloatVec2Property()
    : NumericProperty<tgt::vec2>("", "", tgt::vec2(0.f), tgt::vec2(0.f), tgt::vec2(1.f), tgt::vec2(0.01f),
                                 Processor::INVALID_RESULT, NumericProperty<tgt::vec2>::STATIC, Property::LOD_DEFAULT)
{}

Property* FloatVec2Property::create() const {
    return new FloatVec2Property();
}

//---------------------------------------------------------------------------

FloatVec3Property::FloatVec3Property(const std::string& id, const std::string& guiText,
                             const tgt::vec3& value, const tgt::vec3& minValue,
                             const tgt::vec3& maxValue, int invalidationLevel, NumericProperty<tgt::vec3>::BoundaryUpdatePolicy bup, Property::LevelOfDetail lod)
    : NumericProperty<tgt::vec3>(id, guiText, value, minValue, maxValue, tgt::vec3(0.01f), invalidationLevel, bup, lod)
{}

FloatVec3Property::FloatVec3Property()
    : NumericProperty<tgt::vec3>("", "", tgt::vec3(0.f), tgt::vec3(0.f), tgt::vec3(1.f), tgt::vec3(0.01f),
        Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_DEFAULT)
{}

Property* FloatVec3Property::create() const {
    return new FloatVec3Property();
}

//---------------------------------------------------------------------------

FloatVec4Property::FloatVec4Property(const std::string& id, const std::string& guiText,
                             const tgt::vec4& value, const tgt::vec4& minValue,
                             const tgt::vec4& maxValue, int invalidationLevel, NumericProperty<tgt::vec4>::BoundaryUpdatePolicy bup, Property::LevelOfDetail lod)
    : NumericProperty<tgt::vec4>(id, guiText, value, minValue, maxValue, tgt::vec4(0.01f), invalidationLevel, bup, lod)
{}

FloatVec4Property::FloatVec4Property()
    : NumericProperty<tgt::vec4>("", "", tgt::vec4(0.f), tgt::vec4(0.f), tgt::vec4(1.f), tgt::vec4(0.01f),
        Processor::INVALID_RESULT, NumericProperty<tgt::vec4>::STATIC, Property::LOD_DEFAULT)
{}

Property* FloatVec4Property::create() const {
    return new FloatVec4Property();
}

//---------------------------------------------------------------------------

IntVec2Property::IntVec2Property(const std::string& id, const std::string& guiText,
                         const tgt::ivec2& value, const tgt::ivec2& minValue,
                         const tgt::ivec2& maxValue, int invalidationLevel, NumericProperty<tgt::ivec2>::BoundaryUpdatePolicy bup, Property::LevelOfDetail lod)
    : NumericProperty<tgt::ivec2>(id, guiText, value, minValue, maxValue, tgt::ivec2(1), invalidationLevel, bup, lod)
{}

IntVec2Property::IntVec2Property()
    : NumericProperty<tgt::ivec2>("", "", tgt::ivec2(0), tgt::ivec2(0), tgt::ivec2(100), tgt::vec2(1),
        Processor::INVALID_RESULT, NumericProperty<tgt::ivec2>::STATIC,Property::LOD_DEFAULT)
{}

Property* IntVec2Property::create() const {
    return new IntVec2Property();
}

//---------------------------------------------------------------------------

IntVec3Property::IntVec3Property(const std::string& id, const std::string& guiText,
                         const tgt::ivec3& value, const tgt::ivec3& minValue,
                         const tgt::ivec3& maxValue, int invalidationLevel, NumericProperty<tgt::ivec3>::BoundaryUpdatePolicy bup, Property::LevelOfDetail lod)
    : NumericProperty<tgt::ivec3>(id, guiText, value, minValue, maxValue, tgt::ivec3(1), invalidationLevel, bup, lod)
{}

IntVec3Property::IntVec3Property()
    : NumericProperty<tgt::ivec3>("", "", tgt::ivec3(0), tgt::ivec3(0), tgt::ivec3(100), tgt::vec3(1),
        Processor::INVALID_RESULT, NumericProperty<tgt::ivec3>::STATIC, Property::LOD_DEFAULT)
{}

Property* IntVec3Property::create() const {
    return new IntVec3Property();
}

//---------------------------------------------------------------------------

IntVec4Property::IntVec4Property(const std::string& id, const std::string& guiText,
                         const tgt::ivec4& value, const tgt::ivec4& minValue,
                         const tgt::ivec4& maxValue, int invalidationLevel, NumericProperty<tgt::ivec4>::BoundaryUpdatePolicy bup, Property::LevelOfDetail lod)
    : NumericProperty<tgt::ivec4>(id, guiText, value, minValue, maxValue, tgt::ivec4(1), invalidationLevel, bup, lod)
{}

IntVec4Property::IntVec4Property()
    : NumericProperty<tgt::ivec4>("", "", tgt::ivec4(0), tgt::ivec4(0), tgt::ivec4(100), tgt::vec4(1),
        Processor::INVALID_RESULT, NumericProperty<tgt::ivec4>::STATIC,Property::LOD_DEFAULT)
{}

Property* IntVec4Property::create() const {
    return new IntVec4Property();
}

} // namespace voreen
