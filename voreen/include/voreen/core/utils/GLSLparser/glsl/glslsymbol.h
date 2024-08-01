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

#ifndef VRN_GLSLSYMBOL_H
#define VRN_GLSLSYMBOL_H

#include "voreen/core/utils/GLSLparser/symbol.h"

#include "voreen/core/utils/GLSLparser/glslannotation.h"

namespace voreen {

namespace glslparser {

class GLSLSymbol : public Symbol {
public:
    enum InternalType {
        INTERNAL_VOID,
        INTERNAL_BOOL,
        INTERNAL_FLOAT,
        INTERNAL_INT,
        INTERNAL_UINT,
        INTERNAL_SAMPLER,
        INTERNAL_STRUCT
    };

    enum PrecisionQualifier {
        PRECQ_NONE,
        PRECQ_LOW,
        PRECQ_MEDIUM,
        PRECQ_HIGH
    };

public:
    GLSLSymbol(const std::string& identifier, const InternalType type, const int elementCount,
        const bool isDeclared = true);

    virtual ~GLSLSymbol() = 0;

    const std::vector<GLSLAnnotation*>& getAnnotations() const { return annotations_; }
    const GLSLAnnotation* getAnnotation(std::string name) const;
    const std::string getAnnotationValueString(std::string name) const;

    InternalType getInternalType() const { return internalType_; }
    int getNumInternalElements() const { return elementCount_; }
    PrecisionQualifier getPrecision() const { return precision_; }

    bool isDeclared() const { return isDeclared_; }

    void setAnnotations(const std::vector<GLSLAnnotation*>& annotations) {
        annotations_ = annotations;
    }
    void setInternalType(const InternalType internalType) { internalType_ = internalType; }
    void setPrecision(const PrecisionQualifier precision) { precision_ = precision; }

protected:
    InternalType internalType_;
    PrecisionQualifier precision_;
    int elementCount_;  // indicates number of elements of internal type (3 for vec3, 4 for vec4, etc.)
    std::vector<GLSLAnnotation*> annotations_;
    bool isDeclared_;
};

// ============================================================================

class GLSLVariableSymbol : public GLSLSymbol {
public:
    enum StorageQualifier {
        SQ_NONE,
        SQ_CONST,
        SQ_ATTRIBUTE,
        SQ_VARYING,
        SQ_CENTROID_VARYING,
        SQ_IN,
        SQ_OUT,
        SQ_CENTROID_IN,
        SQ_CENTROID_OUT,
        SQ_UNIFORM
    };

    enum InterpolationQualifier {
        IQ_NONE,
        IQ_FLAT,
        IQ_SMOOTH,
        IQ_NOPERSPECTIVE
    };

public:
    GLSLVariableSymbol(const std::string& identifier, const InternalType type,
        const int elementCount, const bool isDeclared = true);

    InterpolationQualifier getInterpolationQualifier() const { return interpolation_; }
    int getNumArrayElements() const { return numArrayElements_; }
    StorageQualifier getStorageQualifier() const { return storage_; }
    bool isArray() const { return isArray_; }
    bool isInvariant() const { return hasInvariantQualifier_; }

    void setInterpolation(const InterpolationQualifier interpolation) { interpolation_ = interpolation; }
    void setInvariant(const bool invariant) { hasInvariantQualifier_ = invariant; }
    void setIsArray(const bool isArray) { isArray_ = isArray; }
    void setNumArrayElements(const int numArrayElements) { numArrayElements_ = numArrayElements; }
    void setStorage(const StorageQualifier storage) { storage_ = storage; }

protected:
    StorageQualifier storage_;
    InterpolationQualifier interpolation_;

    bool hasInvariantQualifier_;
    bool isArray_;
    int numArrayElements_;
    // initializer?
};

/*
// ============================================================================

class GLSLParameterSymbol : public GLSLSymbol {
public:
    enum ParameterQualifier { PQ_NONE, PQ_IN, PQ_OUT, PQ_INOUT };

public:
    GLSLParameterSymbol(const std::string& identifier);
};
*/

}   // namespace glslparser

}   // namespace voreen

#endif
