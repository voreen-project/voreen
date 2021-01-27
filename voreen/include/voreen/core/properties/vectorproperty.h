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

#ifndef VRN_VECTORPROPERTY_H
#define VRN_VECTORPROPERTY_H

#include "voreen/core/properties/condition.h"
#include "voreen/core/properties/numericproperty.h"

#include <sstream>

namespace voreen {

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<tgt::vec2>;
#endif
class VRN_CORE_API FloatVec2Property : public NumericProperty<tgt::vec2> {
public:
    FloatVec2Property(const std::string& id, const std::string& guiText, const tgt::vec2& value,
                  const tgt::vec2& minimum = tgt::vec2(0.0f), const tgt::vec2& maximum = tgt::vec2(1.0f),
                  int invalidationLevel=Processor::INVALID_RESULT, NumericProperty<tgt::vec2>::BoundaryUpdatePolicy bup = NumericProperty<tgt::vec2>::STATIC,
                  Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    FloatVec2Property();
    virtual ~FloatVec2Property() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "FloatVec2Property"; }
    virtual std::string getTypeDescription() const { return "FloatVector2"; }
};

//---------------------------------------------------------------------------

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<tgt::vec3>;
#endif
class VRN_CORE_API FloatVec3Property : public NumericProperty<tgt::vec3> {
public:
    FloatVec3Property(const std::string& id, const std::string& guiText, const tgt::vec3& value,
                  const tgt::vec3& minimum = tgt::vec3(0.0f), const tgt::vec3& maximum = tgt::vec3(1.0f),
                  int invalidationLevel=Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::BoundaryUpdatePolicy bup = NumericProperty<tgt::vec3>::STATIC,
                  Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    FloatVec3Property();
    virtual ~FloatVec3Property() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "FloatVec3Property"; }
    virtual std::string getTypeDescription() const { return "FloatVector3"; }
};

//---------------------------------------------------------------------------

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<tgt::vec4>;
#endif
class VRN_CORE_API FloatVec4Property : public NumericProperty<tgt::vec4> {
public:
    FloatVec4Property();
    FloatVec4Property(const std::string& id, const std::string& guiText, const tgt::vec4& value,
                  const tgt::vec4& minimum = tgt::vec4(0.0f), const tgt::vec4& maximum = tgt::vec4(1.0f),
                  int invalidationLevel=Processor::INVALID_RESULT, NumericProperty<tgt::vec4>::BoundaryUpdatePolicy bup = NumericProperty<tgt::vec4>::STATIC,
                  Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    virtual ~FloatVec4Property() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "FloatVec4Property"; }
    virtual std::string getTypeDescription() const { return "FloatVector4"; }
};

//---------------------------------------------------------------------------

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<tgt::ivec2>;
#endif
class VRN_CORE_API IntVec2Property : public NumericProperty<tgt::ivec2> {
public:
    IntVec2Property();
    IntVec2Property(const std::string& id, const std::string& guiText, const tgt::ivec2& value,
                const tgt::ivec2& minimum = tgt::ivec2(0), const tgt::ivec2& maximum = tgt::ivec2(100),
                int invalidationLevel=Processor::INVALID_RESULT, NumericProperty<tgt::ivec2>::BoundaryUpdatePolicy bup = NumericProperty<tgt::ivec2>::STATIC,
                Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    virtual ~IntVec2Property() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "IntVec2Property"; }
    virtual std::string getTypeDescription() const { return "IntVector2"; }
};

//---------------------------------------------------------------------------

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<tgt::ivec3>;
#endif
class VRN_CORE_API IntVec3Property : public NumericProperty<tgt::ivec3> {
public:
    IntVec3Property(const std::string& id, const std::string& guiText, const tgt::ivec3& value,
                const tgt::ivec3& minimum = tgt::ivec3(0), const tgt::ivec3& maximum = tgt::ivec3(100),
                int invalidationLevel=Processor::INVALID_RESULT, NumericProperty<tgt::ivec3>::BoundaryUpdatePolicy bup = NumericProperty<tgt::ivec3>::STATIC,
                Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    IntVec3Property();
    virtual ~IntVec3Property() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "IntVec3Property"; }
    virtual std::string getTypeDescription() const { return "IntVector3"; }
};

//---------------------------------------------------------------------------

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API NumericProperty<tgt::ivec4>;
#endif
class VRN_CORE_API IntVec4Property : public NumericProperty<tgt::ivec4> {
public:
    IntVec4Property(const std::string& id, const std::string& guiText, const tgt::ivec4& value,
                const tgt::ivec4& minimum = tgt::ivec4(0), const tgt::ivec4& maximum = tgt::ivec4(100),
                int invalidationLevel=Processor::INVALID_RESULT, NumericProperty<tgt::ivec4>::BoundaryUpdatePolicy bup = NumericProperty<tgt::ivec4>::STATIC,
                Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    IntVec4Property();
    virtual ~IntVec4Property() {}

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "IntVec4Property"; }
    virtual std::string getTypeDescription() const { return "IntVector4"; }
};

} // namespace voreen

#endif // VRN_VECTORPROPERTY_H
