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

#ifndef VRN_GLSLVALUE_H
#define VRN_GLSLVALUE_H

#include <string>

namespace voreen {

namespace glslparser {

/**
 * Abstract base class for any value which can be the result from evaluating a
 * GLSL expression.
 */
class GLSLValue
{
public:
    enum GLSLValueType {
        VALUE_INT, VALUE_UINT, VALUE_FLOAT, VALUE_BOOL, VALUE_ADDRESS
    };

public:
    GLSLValue(const GLSLValueType type)
        : type_(type)
    {
    }

    virtual ~GLSLValue() = 0;

    GLSLValueType getType() const { return type_; }

private:
    const GLSLValueType type_;
};  // class GLSLValue

inline GLSLValue::~GLSLValue() {}

// ============================================================================

class GLSLValueInt : public GLSLValue {
public:
    GLSLValueInt(const int val)
        : GLSLValue(VALUE_INT),
        value_(val)
    {
    }

    int value_;
};  // class GLSLValueInt

// ============================================================================

class GLSLValueUint : public GLSLValue {
public:
    GLSLValueUint(const unsigned int val)
        : GLSLValue(VALUE_UINT),
        value_(val)
    {
    }

    unsigned int value_;
};  // class GLSLValueUint

// ============================================================================

class GLSLValueFloat : public GLSLValue {
public:
    GLSLValueFloat(const float val)
        : GLSLValue(VALUE_FLOAT),
        value_(val)
    {
    }

    float value_;
};  // class GLSLValueFloat

// ============================================================================

class GLSLValueBool : public GLSLValue {
public:
    GLSLValueBool(const bool val)
        : GLSLValue(VALUE_BOOL),
        value_(val)
    {
    }

    bool value_;
};  // class GLSLValueBool

// ============================================================================


class GLSLValueAddress : public GLSLValue
{
public:
    GLSLValueAddress(const std::string& symbolName)
        : GLSLValue(VALUE_ADDRESS),
        symbolName_(symbolName),
        arrayIndex_(-1)
    {
    }

    GLSLValueAddress(const std::string& symbolName, const int arrayIndex)
        : GLSLValue(VALUE_ADDRESS),
        symbolName_(symbolName),
        arrayIndex_(arrayIndex)
    {
    }

    bool isArray() const { return (arrayIndex_ >= 0); }

    const std::string& symbolName_;
    int arrayIndex_;
};  // class GLSLValueAddress

// ============================================================================

}   // namespace glslparser

}   // namespace voreen

#endif  // VRN_GLSLVALUE_H
