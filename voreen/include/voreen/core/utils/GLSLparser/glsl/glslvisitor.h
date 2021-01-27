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

#ifndef VRN_GLSLVISITOR_H
#define VRN_GLSLVISITOR_H

#include "voreen/core/utils/GLSLparser/symboltable.h"
#include "voreen/core/utils/GLSLparser/glsl/glslfunctioncall.h"
#include "voreen/core/utils/GLSLparser/glsl/glslfunctiondefinition.h"
#include "voreen/core/utils/GLSLparser/glsl/glsltoken.h"
#include "voreen/core/utils/GLSLparser/glsl/glslstructspecifier.h"
#include "voreen/core/utils/GLSLparser/glsl/glsltranslation.h"
#include "voreen/core/utils/GLSLparser/glsl/glslsymbol.h"

#include <set>
#include <sstream>

namespace voreen {

namespace glslparser {

class GLSLVisitor : public ParseTreeVisitor {
public:
    GLSLVisitor();
    virtual ~GLSLVisitor();

    virtual bool visit(ParseTreeNode* const node);

    void printGlobalSymbolTable();

    void printSymbolTable(const std::string& tableName, const std::map<std::string, GLSLSymbol*>& table,
        const unsigned int depth) const;

    unsigned int getNumWarnings() const { return numWarnings_; }

    /**
     * Returns all declarations of uniform symbols within this GLSL program.
     *
     * The returned symbols are no long kept within the symbol table if the
     * parameter keepInTable is set to true and the caller must delete the
     * returned pointers using C++ delete operator in order to prevent memory
     * leaks.
     * If keepInTable is false, the pointers are kept within the table, but
     * become deleted invalid as soon as the dtor of GLSLVistitor is called,
     * for it deletes its symbol table and all contained symbols.
     */
    std::vector<GLSLVariableSymbol*> getUniforms(const bool keepInTable);

    /**
     * Returns all declarations of out symbols within this GLSL program.
     *
     * The returned symbols are no long kept within the symbol table if the
     * parameter keepInTable is set to true and the caller must delete the
     * returned pointers using C++ delete operator in order to prevent memory
     * leaks.
     * If keepInTable is false, the pointers are kept within the table, but
     * become deleted invalid as soon as the dtor of GLSLVistitor is called,
     * for it deletes its symbol table and all contained symbols.
     */
    std::vector<GLSLVariableSymbol*> getOuts(const bool keepInTable);

    /**
     * Returns a set with the indices of the special array gl_FragData which have
     * been used during assignments as LHS operands.
     */
    const std::set<unsigned int>& getReferencedGlFragData() const { return glFragDataElements_; }

protected:
    std::string getPrefixLower(const std::string& input, const unsigned int len) const;

    std::vector<GLSLVariableSymbol*> getStorageQualifiedVars(const bool keepInTable,
        const GLSLVariableSymbol::StorageQualifier storageQualifier);

    void pushSymbolTable(const std::string& newName);

    void popSymbolTable();

private:
    typedef std::map<std::string, GLSLSymbol*> SymbolMap;

private:
    // Declarations
    //
    void visitNode(GLSLDeclaration* const decl);
    void visitNode(GLSLDeclarationList* const decls);
    //void visitNode(GLSLStructDeclaratorList* const decls);
    //void visitNode(GLSLFieldDeclaration* const decl);
    //void visitNode(GLSLFunctionDeclaration* const decl);

    // Expressions
    //
    GLSLValue* visitNode(GLSLExpression* const expr);
    GLSLValue* visitNode(GLSLExpressionList* const exprLst);
    GLSLValue* visitNode(GLSLAssignmentExpression* const assign);
    GLSLValue* visitNode(GLSLBinaryExpression* const bin);
    GLSLValue* visitNode(GLSLConditionalExpression* const cond);
    GLSLValue* visitNode(GLSLFunctionCall* const funCall);
    GLSLValue* visitNode(GLSLUnaryExpression* const unary);
    GLSLValue* visitNode(GLSLPrimaryExpression* const primaryExpr);
    GLSLValue* visitNode(GLSLPostfixExpression* const postfix);

    // Qualifiers
    //
    GLSLVariableSymbol::InterpolationQualifier visitNode(GLSLInterpolationQualifier* const interpolation);
    GLSLSymbol::PrecisionQualifier visitNode(GLSLPrecisionQualifier* const precision);
    GLSLVariableSymbol::StorageQualifier visitNode(GLSLStorageQualifier* const storage);

    // Specifiers
    //
    GLSLVariableSymbol visitNode(GLSLTypeSpecifier* const typeSpec);

    // Statements
    //
    void visitNode(GLSLStatement* const statement);
    void visitNode(GLSLCompoundStatement* const statements);
    void visitNode(GLSLSimpleStatement* const statement);
    void visitNode(GLSLCaseLabel* const lbl);
    void visitNode(GLSLDeclarationStatement* const decl);
    void visitNode(GLSLDoWhileStatement* const dwhl);
    void visitNode(GLSLExpressionStatement* const expr);
    void visitNode(GLSLForStatement* const fr);
    void visitNode(GLSLJumpStatement* const jmp);
    void visitNode(GLSLSelectionStatement* const sel);
    void visitNode(GLSLSwitchStatement* const swtch);
    void visitNode(GLSLWhileStatement* const whl);

    // Misc
    //
    GLSLValue* visitNode(GLSLCondition* const cond);
    void visitNode(GLSLExternalDeclaration* const extDecl);
    void visitNode(GLSLFunctionDefinition* const funcDef);
    void visitNode(GLSLFunctionPrototype* const funcProto);
    void visitNode(GLSLParameter* const param);
    void visitNode(GLSLTranslation* const trans);

    std::vector<GLSLAnnotation*> processAnnotation(AnnotationToken* const annotation);

protected:
    typedef SymbolTable<GLSLSymbol> GLSLSymbolMap;
    GLSLSymbolMap globalSymbols_;

    GLSLSymbolMap* activeSymbols_;

private:
    GLSLTerminals terminals_;
    unsigned int numWarnings_;
    std::set<unsigned int> glFragDataElements_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
