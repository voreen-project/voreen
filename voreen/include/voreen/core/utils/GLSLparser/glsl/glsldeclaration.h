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

#ifndef VRN_GLSLDECLARATION_H
#define VRN_GLSLDECLARATION_H

#include "voreen/core/utils/GLSLparser/glsl/glslexternaldeclaration.h"
#include "voreen/core/utils/GLSLparser/glsl/glslfunctionprototype.h"

namespace voreen {

namespace glslparser {

class GLSLDeclaration : public GLSLExternalDeclaration {
public:
    GLSLDeclaration(const int symbolID)
        //: GLSLNode(symbolID),
        : GLSLExternalDeclaration(symbolID),
        leadingAnnotation_(0),
        trailingAnnotation_(0)
    {
    }

    virtual ~GLSLDeclaration() {
        delete leadingAnnotation_;
        delete trailingAnnotation_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_DECLARATION; }

    AnnotationToken* getLeadingAnnotation() const { return leadingAnnotation_; }
    AnnotationToken* getTrailingAnnotation() const { return trailingAnnotation_; }

    void setLeadingAnnotation(AnnotationToken* const leading) {
        leadingAnnotation_ = dynamic_cast<AnnotationToken*>(leading->getCopy());
    }

    void setTrailingAnnotation(AnnotationToken* const trailing) {
        trailingAnnotation_ = dynamic_cast<AnnotationToken*>(trailing->getCopy());
    }

protected:
    AnnotationToken* leadingAnnotation_;
    AnnotationToken* trailingAnnotation_;
};

// ============================================================================

class GLSLDeclarationList : public GLSLDeclaration {
public:
    GLSLDeclarationList(GLSLTypeSpecifier* const typeSpecifier, GLSLVariable* const var)
        : GLSLDeclaration(var->getSymbolID()),
        typeSpecifier_(typeSpecifier)
    {
        addVariable(var);
    }

    virtual ~GLSLDeclarationList() {
        delete typeSpecifier_;
        for (size_t i = 0; i < vars_.size(); ++i)
            delete vars_[i];
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_DECLARATION_LIST; }

    void addVariable(GLSLVariable* const var) {
        if (var != 0)
            vars_.push_back(var);
    }

    GLSLTypeSpecifier* getTypeSpecifier() const { return typeSpecifier_; }
    const std::vector<GLSLVariable*>& getVariables() const { return vars_; }

    void setTypeSpecifier(GLSLTypeSpecifier* const typeSpec) {
        typeSpecifier_ = typeSpec;
    }

protected:
    GLSLTypeSpecifier* typeSpecifier_;
    std::vector<GLSLVariable*> vars_;
};

// ============================================================================

class GLSLStructDeclaratorList : public GLSLDeclaration {
public:
    GLSLStructDeclaratorList(GLSLDeclarationList* const declarator)
        : GLSLDeclaration(declarator->getSymbolID())
    {
        addDeclarator(declarator);
    }

    virtual ~GLSLStructDeclaratorList() {
        for (size_t i = 0; i < declarators_.size(); ++i)
            delete declarators_[i];
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_STRUCT_DECLARATION; }

    /**
     * The declaration withint the list may not contain any initializers.
     */
    void addDeclarator(GLSLDeclarationList* const declarator) {
        if (declarator != 0)
            declarators_.push_back(declarator);
    }

    const std::vector<GLSLDeclarationList*>& getDeclators() const {
        return declarators_;
    }

protected:
    std::vector<GLSLDeclarationList*> declarators_;
};

// ============================================================================

class GLSLFieldDeclaration : public GLSLDeclaration {
public:
    GLSLFieldDeclaration(GLSLTypeQualifier* const typeQualifier, IdentifierToken* const name,
        GLSLStructDeclaratorList* const structDecls)
        : GLSLDeclaration(name->getTokenID()),
        typeQualifier_(typeQualifier),
        name_(dynamic_cast<IdentifierToken*>(name->getCopy())),
        structDecls_(structDecls),
        var_(0)
    {
    }

    virtual ~GLSLFieldDeclaration() {
        delete typeQualifier_;
        delete name_;
        delete structDecls_;
        delete var_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_FIELD_DECLARATION; }

    void setVariable(GLSLVariable* const var) { var_ = var; }

protected:
    GLSLTypeQualifier* const typeQualifier_;
    IdentifierToken* const name_;
    GLSLStructDeclaratorList* const structDecls_;
    GLSLVariable* var_;

};

// ============================================================================

class GLSLFunctionDeclaration : public GLSLDeclaration {
public:
    /**
     * @param   funcProto   Must not be NULL!
     */
    GLSLFunctionDeclaration(GLSLFunctionPrototype* const funcProto)
        : GLSLDeclaration(funcProto->getName()->getTokenID()),
        funcProto_(funcProto)
    {
    }

    virtual ~GLSLFunctionDeclaration() {
        delete funcProto_;
    }

    virtual int getNodeType() const { return GLSLNodeTypes::NODE_FUNCTION_DECLARATION; }

    GLSLFunctionPrototype* getFunctionPrototype() { return funcProto_; }
protected:
    GLSLFunctionPrototype* const funcProto_;
};  // class GLSLFunctionDeclaration

}   // namespace glslparser

}   // namespace voreen

#endif
