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

#ifndef VRN_MATRIXROPERTY_H
#define VRN_MATRIXROPERTY_H

#include "voreen/core/properties/condition.h"
#include "voreen/core/properties/numericproperty.h"

#include <sstream>

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<tgt::mat2>;
#endif

class VRN_CORE_API FloatMat2Property : public NumericProperty<tgt::mat2> {
public:
    FloatMat2Property(const std::string& id, const std::string& guiText, const tgt::mat2& value,
        const tgt::mat2& minimum = tgt::mat2(-10.0f), const tgt::mat2& maximum = tgt::mat2(10.0f),
        int invalidationLevel=Processor::INVALID_RESULT, NumericProperty<tgt::mat2>::BoundaryUpdatePolicy bup = NumericProperty<tgt::mat2>::STATIC,
        Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    FloatMat2Property();
    virtual ~FloatMat2Property() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "FloatMat2Property"; }
    virtual std::string getTypeDescription() const { return "FloatMatrix2x2"; }
};

//---------------------------------------------------------------------------

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<tgt::mat3>;
#endif

class VRN_CORE_API FloatMat3Property : public NumericProperty<tgt::mat3> {
public:
    FloatMat3Property(const std::string& id, const std::string& guiText, const tgt::mat3& value,
        const tgt::mat3& minimum = tgt::mat3(-10.0f), const tgt::mat3& maximum = tgt::mat3(10.0f),
        int invalidationLevel=Processor::INVALID_RESULT, NumericProperty<tgt::mat3>::BoundaryUpdatePolicy bup = NumericProperty<tgt::mat3>::STATIC,
        Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    FloatMat3Property();
    virtual ~FloatMat3Property() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "FloatMat3Property"; }
    virtual std::string getTypeDescription() const { return "FloatMatrix3x3"; }
};

//---------------------------------------------------------------------------

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<tgt::mat4>;
#endif

class VRN_CORE_API FloatMat4Property : public NumericProperty<tgt::mat4> {
public:
    FloatMat4Property(const std::string& id, const std::string& guiText, const tgt::mat4& value,
        const tgt::mat4& minimum = tgt::mat4(-10.0f), const tgt::mat4& maximum = tgt::mat4(10.0f),
        int invalidationLevel=Processor::INVALID_RESULT, NumericProperty<tgt::mat4>::BoundaryUpdatePolicy bup = NumericProperty<tgt::mat4>::STATIC,
        Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    FloatMat4Property();
    virtual ~FloatMat4Property() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "FloatMat4Property"; }
    virtual std::string getTypeDescription() const { return "FloatMatrix4x4"; }
};

} // namespace voreen

#endif // VRN_MATRIXROPERTY_H
