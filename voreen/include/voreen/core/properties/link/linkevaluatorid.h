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

#ifndef VRN_LINKEVALUATORID_H
#define VRN_LINKEVALUATORID_H

#include "voreen/core/properties/link/linkevaluatoridgeneric.h"

#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/lightsourceproperty.h"
#include "voreen/core/properties/volumeurlproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

namespace voreen {

///just a dummy class to enable auto conversion:
class VRN_CORE_API LinkEvaluatorId : public LinkEvaluatorBase {
public:
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorId(); }
    virtual std::string getClassName() const  { return "LinkEvaluatorId"; }
    virtual std::string getGuiName() const    { return "Identity"; }

    virtual void eval(Property*, Property*) {}
    virtual bool arePropertiesLinkable(const Property*, const Property*) const { return true; }
};

///just a dummy class to enable auto conversion:
class VRN_CORE_API LinkEvaluatorIdNormalized : public LinkEvaluatorBase {
public:
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIdNormalized(); }
    virtual std::string getClassName() const  { return "LinkEvaluatorIdNormalized"; }
    virtual std::string getGuiName() const     { return "Normalization"; }

    virtual void eval(Property*, Property*) {}
    virtual bool arePropertiesLinkable(const Property*, const Property*) const { return true; }
};

//-------------------------------------------------------------------------------------------------------
//Bool, Int, Float, Double:

class VRN_CORE_API LinkEvaluatorBoolId : public LinkEvaluatorIdGeneric<bool> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorBoolId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorBoolId(); }
};

class VRN_CORE_API LinkEvaluatorIntId : public LinkEvaluatorIdGeneric<int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntId(); }
};

class VRN_CORE_API LinkEvaluatorFloatId : public LinkEvaluatorIdGeneric<float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatId(); }
};

class VRN_CORE_API LinkEvaluatorDoubleId : public LinkEvaluatorIdGeneric<double> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDoubleId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDoubleId(); }
};

class VRN_CORE_API LinkEvaluatorIntIdBounds : public LinkEvaluatorIdNumeric<int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntIdBounds"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntIdBounds(); }
};

class VRN_CORE_API LinkEvaluatorFloatIdBounds : public LinkEvaluatorIdNumeric<float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatIdBounds"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatIdBounds(); }
};

class VRN_CORE_API LinkEvaluatorDoubleIdBounds : public LinkEvaluatorIdNumeric<double> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDoubleIdBounds"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDoubleIdBounds(); }
};

// id conversion
class VRN_CORE_API LinkEvaluatorDoubleFloatId : public LinkEvaluatorIdGenericConversion<double, float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDoubleFloatId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDoubleFloatId(); }
};

class VRN_CORE_API LinkEvaluatorDoubleIntId : public LinkEvaluatorIdGenericConversion<double, int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDoubleIntId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDoubleIntId(); }
};

class VRN_CORE_API LinkEvaluatorDoubleBoolId : public LinkEvaluatorIdGenericConversion<double, bool> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDoubleBoolId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDoubleBoolId(); }
};

class VRN_CORE_API LinkEvaluatorFloatIntId : public LinkEvaluatorIdGenericConversion<float, int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatIntId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatIntId(); }
};

class VRN_CORE_API LinkEvaluatorFloatBoolId : public LinkEvaluatorIdGenericConversion<float, bool> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatBoolId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoolId(); }
};

class VRN_CORE_API LinkEvaluatorIntBoolId : public LinkEvaluatorIdGenericConversion<int, bool> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntBoolId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoolId(); }
};

// id normalized
class VRN_CORE_API LinkEvaluatorIntIdNormalized : public LinkEvaluatorIdNormalizedGeneric<int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntIdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntIdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorFloatIdNormalized : public LinkEvaluatorIdNormalizedGeneric<float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatIdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatIdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDoubleIdNormalized : public LinkEvaluatorIdNormalizedGeneric<double> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDoubleIdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDoubleIdNormalized(); }
};

// id normalized conversion
class VRN_CORE_API LinkEvaluatorDoubleFloatIdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<double, float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDoubleFloatIdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDoubleFloatIdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDoubleIntIdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<double, int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDoubleIntIdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDoubleIntIdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorFloatIntIdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<float, int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatIntIdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatIntIdNormalized(); }
};

//-------------------------------------------------------------------------------------------------------
// Vector properties:

// id
class VRN_CORE_API LinkEvaluatorIVec2Id : public LinkEvaluatorIdGeneric<tgt::ivec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIVec2Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIVec2Id(); }
};

class VRN_CORE_API LinkEvaluatorIVec3Id : public LinkEvaluatorIdGeneric<tgt::ivec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIVec3Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIVec3Id(); }
};

class VRN_CORE_API LinkEvaluatorIVec4Id : public LinkEvaluatorIdGeneric<tgt::ivec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIVec4Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIVec4Id(); }
};

class VRN_CORE_API LinkEvaluatorVec2Id : public LinkEvaluatorIdGeneric<tgt::vec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec2Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec2Id(); }
};

class VRN_CORE_API LinkEvaluatorVec3Id : public LinkEvaluatorIdGeneric<tgt::vec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec3Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec3Id(); }
};

class VRN_CORE_API LinkEvaluatorVec4Id : public LinkEvaluatorIdGeneric<tgt::vec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec4Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec4Id(); }
};

class VRN_CORE_API LinkEvaluatorDVec2Id : public LinkEvaluatorIdGeneric<tgt::dvec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec2Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec2Id(); }
};

class VRN_CORE_API LinkEvaluatorDVec3Id : public LinkEvaluatorIdGeneric<tgt::dvec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec3Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec3Id(); }
};

class VRN_CORE_API LinkEvaluatorDVec4Id : public LinkEvaluatorIdGeneric<tgt::dvec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec4Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec4Id(); }
};

// id conversion
class VRN_CORE_API LinkEvaluatorDVec2Vec2Id : public LinkEvaluatorIdGenericConversion<tgt::dvec2, tgt::vec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec2Vec2Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec2Vec2Id(); }
};

class VRN_CORE_API LinkEvaluatorDVec3Vec3Id : public LinkEvaluatorIdGenericConversion<tgt::dvec3, tgt::vec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec3Vec3Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec3Vec3Id(); }
};

class VRN_CORE_API LinkEvaluatorDVec4Vec4Id : public LinkEvaluatorIdGenericConversion<tgt::dvec4, tgt::vec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec4Vec4Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec4Vec4Id(); }
};

class VRN_CORE_API LinkEvaluatorDVec2IVec2Id : public LinkEvaluatorIdGenericConversion<tgt::dvec2, tgt::ivec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec2IVec2Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec2IVec2Id(); }
};

class VRN_CORE_API LinkEvaluatorDVec3IVec3Id : public LinkEvaluatorIdGenericConversion<tgt::dvec3, tgt::ivec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec3IVec3Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec3IVec3Id(); }
};

class VRN_CORE_API LinkEvaluatorDVec4IVec4Id : public LinkEvaluatorIdGenericConversion<tgt::dvec4, tgt::ivec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec4IVec4Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec4IVec4Id(); }
};

class VRN_CORE_API LinkEvaluatorVec2IVec2Id : public LinkEvaluatorIdGenericConversion<tgt::vec2, tgt::ivec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec2IVec2Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec2IVec2Id(); }
};

class VRN_CORE_API LinkEvaluatorVec3IVec3Id : public LinkEvaluatorIdGenericConversion<tgt::vec3, tgt::ivec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec3IVec3Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec3IVec3Id(); }
};

class VRN_CORE_API LinkEvaluatorVec4IVec4Id : public LinkEvaluatorIdGenericConversion<tgt::vec4, tgt::ivec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec4IVec4Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec4IVec4Id(); }
};

// id normalized
class VRN_CORE_API LinkEvaluatorIVec2IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::ivec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIVec2IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIVec2IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorIVec3IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::ivec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIVec3IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIVec3IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorIVec4IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::ivec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIVec4IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIVec4IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorVec2IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::vec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec2IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec2IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorVec3IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::vec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec3IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec3IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorVec4IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::vec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec4IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec4IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDVec2IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::dvec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec2IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec2IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDVec3IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::dvec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec3IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec3IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDVec4IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::dvec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec4IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec4IdNormalized(); }
};

// id normalized conversion
class VRN_CORE_API LinkEvaluatorDVec2Vec2IdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<tgt::dvec2, tgt::vec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec2Vec2IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec2Vec2IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDVec3Vec3IdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<tgt::dvec3, tgt::vec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec3Vec3IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec3Vec3IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDVec4Vec4IdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<tgt::dvec4, tgt::vec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec4Vec4IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec4Vec4IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDVec2IVec2IdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<tgt::dvec2, tgt::ivec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec2IVec2IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec2IVec2IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDVec3IVec3IdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<tgt::dvec3, tgt::ivec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec3IVec3IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec3IVec3IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorDVec4IVec4IdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<tgt::dvec4, tgt::ivec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDVec4IVec4IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDVec4IVec4IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorVec2IVec2IdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<tgt::vec2, tgt::ivec2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec2IVec2IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec2IVec2IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorVec3IVec3IdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<tgt::vec3, tgt::ivec3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec3IVec3IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec3IVec3IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorVec4IVec4IdNormalized : public LinkEvaluatorIdNormalizedGenericConversion<tgt::vec4, tgt::ivec4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorVec4IVec4IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorVec4IVec4IdNormalized(); }
};

//-------------------------------------------------------------------------------------------------------
//Links between scalars and vectors

// Programmatically construct the class name of vec to scalar link classes (value)
#define VEC_TO_SCALAR_LINK_CLASSNAME(scalarType, vecType, vecDim, vecComponent) \
    LinkEvaluator ## scalarType ## Vec ## vecType ## vecDim ## Component ## vecComponent ## Id

// Programmatically construct the class name of vec to scalar link classes (normalized)
#define VEC_TO_SCALAR_LINK_CLASSNAME_NORMALIZED(scalarType, vecType, vecDim, vecComponent) \
    LinkEvaluator ## scalarType ## Vec ## vecType ## vecDim ## Component ## vecComponent ## IdNormalized


// Define a macro that executes another macro for all types of vecs
#define OP_VEC_TO_SCALAR_ALL(op, scalarType, vecType) \
    op(scalarType,vecType,2,0) \
    op(scalarType,vecType,2,1) \
    op(scalarType,vecType,3,0) \
    op(scalarType,vecType,3,1) \
    op(scalarType,vecType,3,2) \
    op(scalarType,vecType,4,0) \
    op(scalarType,vecType,4,1) \
    op(scalarType,vecType,4,2) \
    op(scalarType,vecType,4,3)


// Helper: Convert a token to a string
#define TOKEN_TO_STR(s) #s

// Define a macro to declare a link class (for values) for vec and scalar types, a given dimension (-> vec2/3/4) and a component (0 <= vecComponent < vecDim)
// This is intended to be used in OP_VEC_TO_SCALAR_ALL
#define DECLARE_VEC_TO_SCALAR_LINK(scalarType,vecType,vecDim,vecComponent) \
class VRN_CORE_API VEC_TO_SCALAR_LINK_CLASSNAME(scalarType, vecType, vecDim, vecComponent) : public LinkEvaluatorIdGenericVecComponentConversion<tgt::Vector##vecDim < vecType > , scalarType , vecComponent> { \
public: \
    virtual std::string getClassName() const { return TOKEN_TO_STR(VEC_TO_SCALAR_LINK_CLASSNAME(scalarType, vecType, vecDim, vecComponent)); } \
    virtual LinkEvaluatorBase* create() const { return new VEC_TO_SCALAR_LINK_CLASSNAME(scalarType, vecType, vecDim, vecComponent) (); } \
};

// Finally use the macros to declare 4x4x(2+3+4) = 144 link classes
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,bool,bool)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,bool,int)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,bool,float)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,bool,double)

OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,int,bool)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,int,int)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,int,float)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,int,double)

OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,float,bool)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,float,int)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,float,float)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,float,double)

OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,double,bool)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,double,int)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,double,float)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK,double,double)

// We are done using the macro to declare the classes
#undef DECLARE_VEC_TO_SCALAR


// Define a macro to declare a link class (normalized) for vec and scalar types, a given dimension (-> vec2/3/4) and a component (0 <= vecComponent < vecDim)
// This is intended to be used in OP_VEC_TO_SCALAR_ALL
#define DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED(scalarType,vecType,vecDim,vecComponent) \
class VRN_CORE_API VEC_TO_SCALAR_LINK_CLASSNAME_NORMALIZED(scalarType, vecType, vecDim, vecComponent) : public LinkEvaluatorIdNormalizedGenericVecComponentConversion<tgt::Vector##vecDim < vecType > , scalarType , vecComponent> { \
public: \
    virtual std::string getClassName() const { return TOKEN_TO_STR(VEC_TO_SCALAR_LINK_CLASSNAME_NORMALIZED(scalarType, vecType, vecDim, vecComponent)); } \
    virtual LinkEvaluatorBase* create() const { return new VEC_TO_SCALAR_LINK_CLASSNAME_NORMALIZED(scalarType, vecType, vecDim, vecComponent) (); } \
};

// Finally use the macros to declare 3x3x(2+3+4) = 81 link classes
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED,int,int)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED,int,float)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED,int,double)

OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED,float,int)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED,float,float)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED,float,double)

OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED,double,int)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED,double,float)
OP_VEC_TO_SCALAR_ALL(DECLARE_VEC_TO_SCALAR_LINK_NORMALIZED,double,double)

// We are done using the macro to declare the classes
#undef DECLARE_VEC_TO_SCALAR_NORMALIZED

//-------------------------------------------------------------------------------------------------------
//Matrix properties:

class VRN_CORE_API LinkEvaluatorMat2Id : public LinkEvaluatorIdGeneric<tgt::mat2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorMat2Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorMat2Id(); }
};

class VRN_CORE_API LinkEvaluatorMat3Id : public LinkEvaluatorIdGeneric<tgt::mat3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorMat3Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorMat3Id(); }
};

class VRN_CORE_API LinkEvaluatorMat4Id : public LinkEvaluatorIdGeneric<tgt::mat4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorMat4Id"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorMat4Id(); }
};

// id normalized
class VRN_CORE_API LinkEvaluatorMat2IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::mat2> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorMat2IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorMat2IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorMat3IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::mat3> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorMat3IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorMat3IdNormalized(); }
};

class VRN_CORE_API LinkEvaluatorMat4IdNormalized : public LinkEvaluatorIdNormalizedGeneric<tgt::mat4> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorMat4IdNormalized"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorMat4IdNormalized(); }
};

//-------------------------------------------------------------------------------------------------------
//String (+FileDialog)

class VRN_CORE_API LinkEvaluatorStringId : public LinkEvaluatorIdGeneric<std::string> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorStringId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorStringId(); }
};


class VRN_CORE_API LinkEvaluatorIntStringId : public LinkEvaluatorGenericStringId<int> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntStringId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntStringId(); }
};

class VRN_CORE_API LinkEvaluatorFloatStringId : public LinkEvaluatorGenericStringId<float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatStringId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatStringId(); }
};

class VRN_CORE_API LinkEvaluatorDoubleStringId : public LinkEvaluatorGenericStringId<double> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorDoubleStringId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorDoubleStringId(); }
};


//-------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorShaderId : public LinkEvaluatorIdGeneric<ShaderSource> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorShaderId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorShaderId(); }
};

//-------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorCameraId : public LinkEvaluatorBase {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorCameraId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorCameraId(); }
    virtual std::string getGuiName() const { return "Camera Identity"; }

    ///Special implementation to only link position focus and up vector
    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

class VRN_CORE_API LinkEvaluatorColorSwitchId : public LinkEvaluatorBase {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorColorSwitchId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorColorSwitchId(); }
    virtual std::string getGuiName() const { return "ColorSwitch Identity"; }

    ///Special implementation to only link position focus and up vector
    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//-------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorCameraOrientationId : public LinkEvaluatorBase {
public:
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorCameraOrientationId(); }
    virtual std::string getClassName() const { return "LinkEvaluatorCameraOrientationId"; }
    virtual std::string getGuiName() const { return "Camera Orientation Identity"; }

    ///Special implementation to only link position focus and up vector
    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//------------------------------------------------------------------------------------------------------
//
class VRN_CORE_API LinkEvaluatorCameraPosId : public LinkEvaluatorBase {
public:
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorCameraPosId(); }
    virtual std::string getClassName() const { return "LinkEvaluatorCameraPosId"; }
    virtual std::string getGuiName() const { return "Camera Position Identity"; }

    ///Special implementation to only link position with FloatVec3Property
    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//------------------------------------------------------------------------------------------------------
//
class VRN_CORE_API LinkEvaluatorCameraLookId : public LinkEvaluatorBase {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorCameraLookId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorCameraLookId(); }
    virtual std::string getGuiName() const { return "Camera Look Identity"; }

    ///Special implementation to only link look vector with FloatVec3Property
    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorCameraFocusId : public LinkEvaluatorBase {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorCameraFocusId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorCameraFocusId(); }
    virtual std::string getGuiName() const { return "Camera Focus Identity"; }

    ///Special implementation to only link focus vector with FloatVec3Property
    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorCameraFrustumId : public LinkEvaluatorBase {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorCameraFrustumId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorCameraFrustumId(); }
    virtual std::string getGuiName() const { return "Camera Frustum Identity"; }

    ///Special implementation to only link focus vector with FloatVec3Property
    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//-------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorTransFunc1DKeysId : public LinkEvaluatorBase {
public:
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorTransFunc1DKeysId(); }
    virtual std::string getClassName() const  { return "LinkEvaluatorTransFunc1DKeysId"; }
    virtual std::string getGuiName() const    { return "Identity"; }

    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//-------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorTransFunc1DGaussianId : public LinkEvaluatorBase {
public:
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorTransFunc1DGaussianId(); }
    virtual std::string getClassName() const  { return "LinkEvaluatorTransFunc1DGaussianId"; }
    virtual std::string getGuiName() const    { return "Identity"; }

    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//-------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorTransFunc2DPrimitivesId : public LinkEvaluatorBase {
public:
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorTransFunc2DPrimitivesId(); }
    virtual std::string getClassName() const  { return "LinkEvaluatorTransFunc2DPrimitivesId"; }
    virtual std::string getGuiName() const    { return "Identity"; }

    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//-------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorButtonId : public LinkEvaluatorBase {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorButtonId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorButtonId(); }
    virtual std::string getGuiName() const { return "Button"; }

    ///Special implementation to only link position focus and up vector
    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//-------------------------------------------------------------------------------------------------------

class VRN_CORE_API LinkEvaluatorLightSourceId : public LinkEvaluatorBase {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorLightSourceId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorLightSourceId(); }
    virtual std::string getGuiName() const { return "LightSource"; }

    ///Special implementation to make sure light source property widgets are updated
    virtual void eval(Property* src, Property* dst);

    bool arePropertiesLinkable(const Property* p1, const Property* p2) const;
};

//---------------------------------------------------------------------------------------------------------
// Interval and Boundingbox propertys
template<typename T>
class TemplateLinkEvaluatorIntervalId : public LinkEvaluatorBase {
    virtual bool arePropertiesLinkable( const Property* p1, const Property* p2 ) const{
        return dynamic_cast<const IntervalProperty<T>*> (p1)
            && dynamic_cast<const IntervalProperty<T>*> (p2);
    }

    virtual void eval( Property* src, Property* dst ) {
        IntervalProperty<T>* srcP = dynamic_cast<IntervalProperty<T>*>(src);
        IntervalProperty<T>* dstP = dynamic_cast<IntervalProperty<T>*>(dst);
        dstP->set(srcP->get());
    }
};

class VRN_CORE_API LinkEvaluatorFloatIntervalId : public TemplateLinkEvaluatorIntervalId<float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatIntervalId"; }
    virtual std::string getGuiName() const { return "Float Interval"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatIntervalId(); }
};

class VRN_CORE_API LinkEvaluatorIntIntervalId : public TemplateLinkEvaluatorIntervalId<int>{
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntIntervalBoxId"; }
    virtual std::string getGuiName() const { return "Int Interval"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntIntervalId(); }
};

template<typename T>
class TemplateLinkEvaluatorBoundingBoxId : public LinkEvaluatorBase {
    virtual bool arePropertiesLinkable( const Property* p1, const Property* p2 ) const{
        return dynamic_cast<const TemplateBoundingBoxProperty<T>*> (p1)
            && dynamic_cast<const TemplateBoundingBoxProperty<T>*> (p2);
    }

    virtual void eval( Property* src, Property* dst ) {
        TemplateBoundingBoxProperty<T>* srcP = dynamic_cast<TemplateBoundingBoxProperty<T>*>(src);
        TemplateBoundingBoxProperty<T>* dstP = dynamic_cast<TemplateBoundingBoxProperty<T>*>(dst);
        dstP->set(srcP->get());
    }
};

class VRN_CORE_API LinkEvaluatorFloatBoundingBoxId : public TemplateLinkEvaluatorBoundingBoxId<float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatBoundingBoxId"; }
    virtual std::string getGuiName() const { return "Float Boundingbox"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxId(); }
};

class VRN_CORE_API LinkEvaluatorIntBoundingBoxId : public TemplateLinkEvaluatorBoundingBoxId<int>{
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntBoundingBoxId"; }
    virtual std::string getGuiName() const { return "Int Boundingbox"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxId(); }
};


// with limits
template<typename T>
class TemplateLinkEvaluatorIntervalIdWithLimits : public LinkEvaluatorBase {
    virtual bool arePropertiesLinkable( const Property* p1, const Property* p2 ) const{
        return dynamic_cast<const IntervalProperty<T>*> (p1)
            && dynamic_cast<const IntervalProperty<T>*> (p2);
    }

    virtual void eval( Property* src, Property* dst ) {
        IntervalProperty<T>* srcP = dynamic_cast<IntervalProperty<T>*>(src);
        IntervalProperty<T>* dstP = dynamic_cast<IntervalProperty<T>*>(dst);
        dstP->setMinValue(srcP->getMinValue());
        dstP->setMaxValue(srcP->getMaxValue());
        dstP->setMinRange(srcP->getMinRange());
        dstP->setMaxRange(srcP->getMaxRange());
        dstP->set(srcP->get());
    }
};

class VRN_CORE_API LinkEvaluatorFloatIntervalIdWithLimits : public TemplateLinkEvaluatorIntervalIdWithLimits<float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatIntervalIdWithLimits"; }
    virtual std::string getGuiName() const { return "Float Interval with limits"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatIntervalIdWithLimits(); }
};

class VRN_CORE_API LinkEvaluatorIntIntervalIdWithLimits : public TemplateLinkEvaluatorIntervalIdWithLimits<int>{
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntIntervalIdWithLimits"; }
    virtual std::string getGuiName() const { return "Int Interval with limits"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntIntervalIdWithLimits(); }
};

template<typename T>
class TemplateLinkEvaluatorBoundingBoxIdWithLimits : public LinkEvaluatorBase {
    virtual bool arePropertiesLinkable( const Property* p1, const Property* p2 ) const{
        return dynamic_cast<const TemplateBoundingBoxProperty<T>*> (p1)
            && dynamic_cast<const TemplateBoundingBoxProperty<T>*> (p2);
    }

    virtual void eval( Property* src, Property* dst ) {
        TemplateBoundingBoxProperty<T>* srcP = dynamic_cast<TemplateBoundingBoxProperty<T>*>(src);
        TemplateBoundingBoxProperty<T>* dstP = dynamic_cast<TemplateBoundingBoxProperty<T>*>(dst);
        dstP->setMinValue(srcP->getMinValue());
        dstP->setMaxValue(srcP->getMaxValue());
        dstP->setMinRange(srcP->getMinRange());
        dstP->setMaxRange(srcP->getMaxRange());
        dstP->set(srcP->get());
    }
};

class VRN_CORE_API LinkEvaluatorFloatBoundingBoxIdWithLimits : public TemplateLinkEvaluatorBoundingBoxIdWithLimits<float> {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorFloatBoundingBoxIdWithLimits"; }
    virtual std::string getGuiName() const { return "Float Boundingbox with limits"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorFloatBoundingBoxIdWithLimits(); }
};

class VRN_CORE_API LinkEvaluatorIntBoundingBoxIdWithLimits : public TemplateLinkEvaluatorBoundingBoxIdWithLimits<int>{
public:
    virtual std::string getClassName() const { return "LinkEvaluatorIntBoundingBoxIdWithLimits"; }
    virtual std::string getGuiName() const { return "Int Boundingbox with limits"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorIntBoundingBoxIdWithLimits(); }
};



} // namespace

#endif // VRN_LINKEVALUATORID_H
