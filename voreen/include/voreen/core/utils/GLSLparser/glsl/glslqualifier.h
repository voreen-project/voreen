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

#ifndef VRN_GLSLQUALIFIER_H
#define VRN_GLSLQUALIFIER_H

#include "voreen/core/utils/GLSLparser/glsl/glslparsetreenode.h"

namespace voreen {

namespace glslparser {

class GLSLQualifier : public GLSLNode {
public:
    GLSLQualifier(const Token& token)
        : GLSLNode(token.getTokenID()),
        token_(token.getCopy())
    {
    }

    virtual ~GLSLQualifier() {
        delete token_;
    }

    Token* getToken() const { return token_; }

protected:
    Token* const token_;
};

// ============================================================================

class GLSLInterpolationQualifier : public GLSLQualifier {
public:
    GLSLInterpolationQualifier(Token* const token)
        : GLSLQualifier(*token)
    {
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_INTERPOLATION_QUALIFIER; }
};

// ============================================================================

class GLSLLayoutQualifier : public GLSLQualifier {
public:
    GLSLLayoutQualifier(IdentifierToken* const identifier, ConstantToken* const value)
        : GLSLQualifier(Token(GLSLTerminals::ID_LAYOUT))
    {
        addElement(identifier, value);
    }

    virtual ~GLSLLayoutQualifier() {
        for (size_t i = 0; i < elements_.size(); ++i) {
            delete elements_[i].first;
            delete elements_[i].second;
        }
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_LAYOUT_QUALIFIER; }

    void addElement(IdentifierToken* const identifier, ConstantToken* const value) {
        if (identifier != 0)
            elements_.push_back(std::make_pair(identifier, value));
    }

    const std::vector<std::pair<IdentifierToken*, ConstantToken*> >& getElements() const {
        return elements_;
    }

protected:
    std::vector<std::pair<IdentifierToken*, ConstantToken*> > elements_;
};

// ============================================================================

class GLSLParameterQualifier : public GLSLQualifier {
public:
    GLSLParameterQualifier(Token* const token)
        : GLSLQualifier(*token)
    {
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_PARAMETER_QUALIFIER; }
};

// ============================================================================

class GLSLPrecisionQualifier : public GLSLQualifier {
public:
    GLSLPrecisionQualifier(Token* const token)
        : GLSLQualifier(*token)
    {
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_PRECISION_QUALIFIER; }
};

// ============================================================================

class GLSLStorageQualifier : public GLSLQualifier {
public:
    GLSLStorageQualifier(Token* const token, const bool centroid)
        : GLSLQualifier(*token),
        hasCentroid_(centroid)
    {
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_STORAGE_QUALIFIER; }

    bool hasCentroid() const { return hasCentroid_; }

protected:
    bool hasCentroid_;
};

// ============================================================================

class GLSLTypeQualifier : public GLSLQualifier {
public:
    explicit GLSLTypeQualifier(GLSLStorageQualifier* const storage,
                      GLSLInterpolationQualifier* const interpolation,
                      GLSLLayoutQualifier* const layout,
                      const bool hasInvariant)
        : GLSLQualifier(Token(GLSLTerminals::ID_UNKNOWN)),
        interpolation_(interpolation),
        storage_(storage),
        layout_(layout),
        hasInvariant_(hasInvariant)
    {
    }

    virtual ~GLSLTypeQualifier() {
        delete interpolation_;
        delete storage_;
        delete layout_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_TYPE_QUALIFIER; }

    GLSLInterpolationQualifier* getInterpolationQualifier() const {
        return interpolation_;
    }
    GLSLStorageQualifier* getStorageQualifier() const { return storage_; }
    GLSLLayoutQualifier* getLayoutQualifier() const { return layout_; }
    bool hasInvariantQualifier() const { return hasInvariant_; }

protected:
    GLSLInterpolationQualifier* const interpolation_;
    GLSLStorageQualifier* const storage_;
    GLSLLayoutQualifier* const layout_;
    const bool hasInvariant_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
