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

#include "voreen/core/properties/matrixproperty.h"

namespace voreen {

FloatMat2Property::FloatMat2Property(const std::string& id, const std::string& guiText,
                                     const tgt::mat2& value, const tgt::mat2& minimum,
                                     const tgt::mat2& maximum, int invalidationLevel, NumericProperty<tgt::mat2>::BoundaryUpdatePolicy bup,
                                     Property::LevelOfDetail lod)
    : NumericProperty<tgt::mat2>(id, guiText, value, minimum, maximum, tgt::mat2(0.1f),
        invalidationLevel, bup, lod)
{}

FloatMat2Property::FloatMat2Property()
    : NumericProperty<tgt::mat2>("", "", tgt::mat2(0.f), tgt::mat2(-10.f), tgt::mat2(10.f), tgt::mat2(0.1f),
        Processor::INVALID_RESULT)
{}

Property* FloatMat2Property::create() const {
    return new FloatMat2Property();
}

//---------------------------------------------------------------------------

FloatMat3Property::FloatMat3Property( const std::string& id, const std::string& guiText,
                                     const tgt::mat3& value, const tgt::mat3& minimum,
                                     const tgt::mat3& maximum, int invalidationLevel, NumericProperty<tgt::mat3>::BoundaryUpdatePolicy bup,
                                     Property::LevelOfDetail lod)
    : NumericProperty<tgt::mat3>(id, guiText, value, minimum, maximum, tgt::mat3(0.1f),
        invalidationLevel, bup, lod)
{}

FloatMat3Property::FloatMat3Property()
    : NumericProperty<tgt::mat3>("", "", tgt::mat3(0.f), tgt::mat3(-10.f), tgt::mat3(10.f), tgt::mat3(0.1f),
        Processor::INVALID_RESULT)
{}

Property* FloatMat3Property::create() const {
    return new FloatMat3Property();
}

//---------------------------------------------------------------------------

FloatMat4Property::FloatMat4Property(const std::string& id, const std::string& guiText,
                                     const tgt::mat4& value, const tgt::mat4& minimum,
                                     const tgt::mat4& maximum, int invalidationLevel, NumericProperty<tgt::mat4>::BoundaryUpdatePolicy bup,
                                     Property::LevelOfDetail lod)
    : NumericProperty<tgt::mat4>(id, guiText, value, minimum, maximum, tgt::mat4(0.1f),
        invalidationLevel, bup, lod)
{}

FloatMat4Property::FloatMat4Property()
    : NumericProperty<tgt::mat4>("", "", tgt::mat4(0.f), tgt::mat4(-10.f), tgt::mat4(10.f), tgt::mat4(0.1f),
        Processor::INVALID_RESULT)
{}

Property* FloatMat4Property::create() const {
    return new FloatMat4Property();
}

} // namespace voreen
