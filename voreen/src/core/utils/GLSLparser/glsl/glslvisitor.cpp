/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/utils/GLSLparser/glsl/glslvisitor.h"

#include "voreen/core/utils/GLSLparser/annotations/annotationparser.h"
#include "voreen/core/utils/GLSLparser/annotations/annotationvisitor.h"

#include <algorithm>

namespace voreen {

namespace glslparser {

GLSLVisitor::GLSLVisitor()
    : ParseTreeVisitor(),
    globalSymbols_("globals"),
    activeSymbols_(&globalSymbols_),
    terminals_(),
    numWarnings_(0)
{
}

GLSLVisitor::~GLSLVisitor() {
}

void GLSLVisitor::printGlobalSymbolTable() {
    unsigned int depth = 0;

    std::list<GLSLSymbolMap*> children;
    GLSLSymbolMap* child = &globalSymbols_;
    while (child != 0)
    {
        const SymbolMap& symbols = child->getSymbolsMap();
        printSymbolTable(child->getName(), symbols, depth);

        const std::list<GLSLSymbolMap*>& c = child->getChildTables();
        if (! c.empty())
            children.insert(children.end(), c.begin(), c.end());

        if (! children.empty())
        {
            child = children.front();
            children.pop_front();
        }
        else
            child = 0;
    }
}

void GLSLVisitor::printSymbolTable(const std::string& tableName, const std::map<std::string, GLSLSymbol*>& table,
                                   const unsigned int depth) const
{
    std::string indent(depth, ' ');

    std::cout << indent << "Table '" << tableName << "':\n------------------------\n";
    for (SymbolMap::const_iterator it = table.begin(); it != table.end(); ++it) {
        GLSLSymbol* const sym = it->second;

        std::cout << indent << "  name = '" << it->first << "', type = '" << sym->getInternalType() << "'";
        std::cout << ", #elements = " << sym->getNumInternalElements() << "\n";
    }
    std::cout << "\n";
}

std::vector<GLSLVariableSymbol*> GLSLVisitor::getUniforms(const bool keepInTable) {
    return getStorageQualifiedVars(keepInTable, GLSLVariableSymbol::SQ_UNIFORM);
}

std::vector<GLSLVariableSymbol*> GLSLVisitor::getOuts(const bool keepInTable) {
    return getStorageQualifiedVars(keepInTable, GLSLVariableSymbol::SQ_OUT);
}

bool GLSLVisitor::visit(ParseTreeNode* const node) {
    if (node == 0)
        return false;

    bool result = true;
    switch (node->getNodeType())
    {
        case GLSLNodeTypes::NODE_TRANSLATION:
            visitNode(static_cast<GLSLTranslation* const>(node));
            break;

        default:
            log_ << "GLSLVisitor::visit(): error! Initial parse tree node is not of type NODE_TRANSLATION!" << std::endl;
            log_ << "  node type: " << node->getNodeType() << std::endl;
            break;
    }

    return result;
}

// protected methods
//

std::string GLSLVisitor::getPrefixLower(const std::string& input, const unsigned int len) const
{
    if (input.empty())
        return "";

    std::string sub = input.substr(0, len);
    std::transform(sub.begin(), sub.end(), sub.begin(), ::tolower);
    return sub;
}

std::vector<GLSLVariableSymbol*> GLSLVisitor::getStorageQualifiedVars(const bool keepInTable,
    const GLSLVariableSymbol::StorageQualifier storageQualifier)
{
    const SymbolMap& symbols = globalSymbols_.getSymbolsMap();

    std::vector<GLSLVariableSymbol*> vars;
    for (SymbolMap::const_iterator it = symbols.begin(); it != symbols.end(); )
    {
        GLSLVariableSymbol* const varSymbol = dynamic_cast<GLSLVariableSymbol* const>(it->second);
        if (varSymbol == 0)
            continue;

        if (varSymbol->getStorageQualifier() == storageQualifier) {
            vars.push_back(varSymbol);
            if (keepInTable)
                ++it;
            else
                globalSymbols_.removeSymbol((it++)->second);

        } else
            ++it;
    }

    return vars;
}

void GLSLVisitor::pushSymbolTable(const std::string& newName)
{
    GLSLSymbolMap* const active = new GLSLSymbolMap(newName, activeSymbols_);
    activeSymbols_ = active;
}

void GLSLVisitor::popSymbolTable()
{
    GLSLSymbolMap* const parentTable = activeSymbols_->getParentTable();
    if (parentTable)
        activeSymbols_ = parentTable;
}

// private
//

// Expressions
//

GLSLValue* GLSLVisitor::visitNode(GLSLExpression* const expr)
{
    if (! expr)
        return 0;

    switch (expr->getNodeType())
    {
        case GLSLNodeTypes::NODE_ASSIGNMENT_EXPRESSION:
            return visitNode(static_cast<GLSLAssignmentExpression* const>(expr));

        case GLSLNodeTypes::NODE_BINARY_EXPRESSION:
            return visitNode(static_cast<GLSLBinaryExpression* const>(expr));

        case GLSLNodeTypes::NODE_CONDITIONAL_EXPRESSION:
            return visitNode(static_cast<GLSLConditionalExpression* const>(expr));

        case GLSLNodeTypes::NODE_EXPRESSION_LIST:
            return visitNode(static_cast<GLSLExpressionList* const>(expr));

        case GLSLNodeTypes::NODE_FUNCTION_CALL:
            return visitNode(static_cast<GLSLFunctionCall* const>(expr));

        case GLSLNodeTypes::NODE_POSTFIX_EXPRESSION:
            return visitNode(static_cast<GLSLPostfixExpression* const>(expr));

        case GLSLNodeTypes::NODE_PRIMARY_EXPRESSION:
            return visitNode(static_cast<GLSLPrimaryExpression* const>(expr));

        case GLSLNodeTypes::NODE_UNARY_EXPRESSION:
            return visitNode(static_cast<GLSLUnaryExpression* const>(expr));

        default:
            break;
    }   // switch(expr->getNodeType()

    return 0;
}

GLSLValue* GLSLVisitor::visitNode(GLSLExpressionList* const exprLst)
{
    if (! exprLst)
        return 0;

    GLSLValue* res = 0;
    std::vector<GLSLExpression*> lst = exprLst->getExpressions();
    for (size_t i = 0; i < lst.size(); ++i)
        res = visitNode(lst[i]);
    return res;
}

GLSLValue* GLSLVisitor::visitNode(GLSLAssignmentExpression* const assign)
{
    if (! assign)
        return 0;

    // Perform semantic action for assignment expressions:
    // The action in this case is currently only to evaluate check the lefthand side value
    // (a unary expression) and check if its token is an IDENTIFIER. If so, check if
    // the the token's name starts with "gl_" and insert it into the symbol table.
    // In other cases check if a corresponding symbol is containing within the symbol table
    // and emit a warning that the variable has not been declared, if no symbol is yet present.
    //
    // NOTE: currently any undeclared variable will be inserted into the symbol table. In
    // addition, a warning will be issued (not an error!), unless the variable name starts
    // with "gl_" (in lower or upper case).
    //
    GLSLValue* rhs = visitNode(assign->getRValue());
    GLSLValue* lhs = visitNode(assign->getLValue());
    if ((! rhs) || (! lhs))
        return 0;

    if (lhs->getType() == GLSLValue::VALUE_ADDRESS)
    {
        GLSLValueAddress* const addr = static_cast<GLSLValueAddress* const>(lhs);
        const std::string& varName = addr->symbolName_;

        if (! activeSymbols_->findSymbol(varName))
        {
            GLSLVariableSymbol* const sym = new GLSLVariableSymbol(varName, GLSLSymbol::INTERNAL_VOID, 1, false);

            if (getPrefixLower(varName, 3) != "gl_")
            {
                activeSymbols_->insertSymbol(sym);

                ++numWarnings_;
                log_ << "warning: Variable '" << sym->getID() << "' has been assigned without previous declaration!" << std::endl;
            }
            else
                // add gl_*** symbols ALWAYS to global symbol table.
                globalSymbols_.insertSymbol(sym);
        }

        if ((varName == "gl_FragData") && (addr->isArray()))
            glFragDataElements_.insert(addr->arrayIndex_);
    }
    else
    {
        ++numWarnings_;

        log_ << "warning: lhs operand must be an lvalue !" << std::endl;
    }

    return lhs;
}

GLSLValue* GLSLVisitor::visitNode(GLSLBinaryExpression* const /*bin*/)
{
    // TODO: implement
    return 0;
}

GLSLValue* GLSLVisitor::visitNode(GLSLConditionalExpression* const /*cond*/)
{
    // TODO: implement
    return 0;
}

GLSLValue* GLSLVisitor::visitNode(GLSLUnaryExpression* const unary)
{
    if (! unary)
        return 0;

    GLSLExpression* const operand = unary->getExpression();
    if (! operand)
        return 0;

    GLSLValue* res = 0;
    switch (operand->getNodeType())
    {
        case GLSLNodeTypes::NODE_POSTFIX_EXPRESSION:
            res = visitNode(static_cast<GLSLPostfixExpression* const>(operand));
            break;

        case GLSLNodeTypes::NODE_UNARY_EXPRESSION:
            res = visitNode(static_cast<GLSLUnaryExpression*>(operand));
            break;
    }   // switch

    // TODO: to "really" evaluate this expression, check opToken, value and
    // apply that operator.
    //Token* const opToken = unary->getToken();

    return res;
}

GLSLValue* GLSLVisitor::visitNode(GLSLPostfixExpression* const postfix)
{
    if (! postfix)
        return 0;

    GLSLExpression* const operand = postfix->getOperand();
    if (! operand)
        return 0;

    // 1. Evaluate "nested" postfix expression
    //
    GLSLValue* res = 0;
    switch (operand->getNodeType())
    {
        case GLSLNodeTypes::NODE_EXPRESSION:
            res = visitNode(operand);
            break;

        case GLSLNodeTypes::NODE_PRIMARY_EXPRESSION:
            res = visitNode(static_cast<GLSLPrimaryExpression* const>(operand));
            break;

        case GLSLNodeTypes::NODE_FUNCTION_CALL:
            res = visitNode(static_cast<GLSLFunctionCall* const>(operand));
            break;

        case GLSLNodeTypes::NODE_POSTFIX_EXPRESSION:    // recursive call
            res = visitNode(static_cast<GLSLPostfixExpression* const>(operand));
            break;
    }   // switch

    if (! res)
        return 0;

    // 2. Investigate token
    //
    Token* const token = postfix->getToken();
    if (token) {
        switch (token->getTokenID())
        {
            case GLSLTerminals::ID_DEC_OP:
            case GLSLTerminals::ID_INC_OP:
                // TODO: handle sementics for cases
                // ID  62:  [postfix-expression] ::= [postfix-expression] ++
                // ID  63:  [postfix-expression] ::= [postfix-expression] --
                break;

            case GLSLTerminals::ID_DOT:
                // TODO: handle sementics for case
                //ID  61:   [postfix-expression] ::= [postfix-expression] . FIELD-SELECTION
                break;

            // Subscript operator
            case GLSLTerminals::ID_LBRACKET:
            {
//std::cout << "visitNode(GLSLPostfixExpression): evaluating \"[]\" expression...\n";

                // Evaluate integer expression
                //
                GLSLValue* const idx = visitNode(postfix->getIntExpression());
                if ((idx != 0) && (idx->getType() != GLSLValue::VALUE_ADDRESS)
                    && (res->getType() == GLSLValue::VALUE_ADDRESS))
                {
                    GLSLValueAddress* const addr = static_cast<GLSLValueAddress* const>(res);
                    GLSLValueInt* const i = static_cast<GLSLValueInt* const>(idx);
                    addr->arrayIndex_ = i->value_;
                } else
                    std::cout << "no address result returned from operand evaluation!\n";
            }
                break;

            default:
                break;
        }   // switch
    }   // if (token

    return res;
}

GLSLValue* GLSLVisitor::visitNode(GLSLFunctionCall* const /*funCall*/)
{
    // TODO: add semantics. does currently do nothing.
    return 0;
}

GLSLValue* GLSLVisitor::visitNode(GLSLPrimaryExpression* const primaryExpr)
{
    GLSLExpression* const nested = primaryExpr->getExpression();
    if (nested)
        return visitNode(nested);

    Token* const token = primaryExpr->getToken();
    if (! token)
        return 0;

    switch (token->getTokenID())
    {
        case GLSLTerminals::ID_FALSE:
            return new GLSLValueBool(false);

        case GLSLTerminals::ID_TRUE:
            return new GLSLValueBool(true);

        case GLSLTerminals::ID_IDENTIFIER:
            if (IdentifierToken* const id = dynamic_cast<IdentifierToken* const>(token))
                return new GLSLValueAddress(id->getValue());

        case GLSLTerminals::ID_INTCONST:
            if (ConstantToken* const ct = dynamic_cast<ConstantToken* const>(token))
                return new GLSLValueInt(ct->convert<int>());
            break;

        case GLSLTerminals::ID_UINTCONST:
            if (ConstantToken* const ct = dynamic_cast<ConstantToken* const>(token))
                return new GLSLValueUint(ct->convert<unsigned int>());

        case GLSLTerminals::ID_FLOATCONST:
            if (ConstantToken* const ct = dynamic_cast<ConstantToken* const>(token))
                return new GLSLValueFloat(ct->convert<float>());

        default:
            break;
    }

    return 0;
}

// Delcarations
//

void GLSLVisitor::visitNode(GLSLDeclaration* const decl) {
    if (decl == 0)
        return;

    switch (decl->getNodeType())
    {
        case GLSLNodeTypes::NODE_DECLARATION_LIST:
            visitNode(static_cast<GLSLDeclarationList* const>(decl));
            break;

        case GLSLNodeTypes::NODE_STRUCT_DECLARATION:
            //visitNode(static_cast<GLSLStructDeclaration* const>(decl));
            break;

        case GLSLNodeTypes::NODE_FIELD_DECLARATION:
            //visitNode(static_cast<GLSLFieldDeclaration* const>(decl));
            break;

        case GLSLNodeTypes::NODE_FUNCTION_DECLARATION:
            //visitNode(static_cast<GLSLFunctionDeclaration* const>(decl));
            break;

        default:
            break;
    }   // switch (decl->getNodeType()
}

void GLSLVisitor::visitNode(GLSLDeclarationList* const decls) {
    if (decls == 0)
        return;

    GLSLTypeSpecifier* const typeSpec = decls->getTypeSpecifier();
    if (typeSpec == 0)
        throw std::runtime_error("visitNode(GLSLDeclarationList): given declarations are not specified!");

    try {
        GLSLVariableSymbol metaSymbol = visitNode(typeSpec);
        const std::vector<GLSLVariable*>& vars = decls->getVariables();

        for (size_t i = 0; i < vars.size(); ++i) {
            GLSLVariableSymbol* const symbol = new GLSLVariableSymbol(metaSymbol);

            symbol->setID(vars[i]->getName());
            symbol->setIsArray(vars[i]->isArray());
            int numArrayElements = 0;   // TODO: evaluate expression from vars[i] to determine array size
            symbol->setNumArrayElements(numArrayElements);

            std::vector<GLSLAnnotation*> leading = processAnnotation(decls->getLeadingAnnotation());
            std::vector<GLSLAnnotation*> trailing = processAnnotation(decls->getTrailingAnnotation());

            leading.insert(leading.end(), trailing.begin(), trailing.end());
            symbol->setAnnotations(leading);

            activeSymbols_->insertSymbol(symbol);
        }
    } catch (std::runtime_error& e) {
        log_ << "  Exeption: " << e.what() << "\n";
    }
}

// Qualifiers
//

GLSLVariableSymbol::InterpolationQualifier GLSLVisitor::visitNode(
    GLSLInterpolationQualifier* const interpolation)
{
    if (interpolation == 0)
        return GLSLVariableSymbol::IQ_NONE;

    Token* const token = interpolation->getToken();
    if (token != 0) {
        switch (token->getTokenID()) {
            case GLSLTerminals::ID_FLAT:
                return GLSLVariableSymbol::IQ_FLAT;

            case GLSLTerminals::ID_SMOOTH:
                return GLSLVariableSymbol::IQ_SMOOTH;

            case GLSLTerminals::ID_NONPERSPECTIVE:
                return GLSLVariableSymbol::IQ_NOPERSPECTIVE;
        }
    }

    return GLSLVariableSymbol::IQ_NONE;
}

GLSLSymbol::PrecisionQualifier GLSLVisitor::visitNode(GLSLPrecisionQualifier* const precision) {
    if (precision == 0)
        return GLSLSymbol::PRECQ_NONE;

    Token* const token = precision->getToken();
    if (token != 0) {
        switch (token->getTokenID()) {
            case GLSLTerminals::ID_HIGHP:
                return GLSLSymbol::PRECQ_HIGH;

            case GLSLTerminals::ID_MEDIUMP:
                return GLSLSymbol::PRECQ_MEDIUM;

            case GLSLTerminals::ID_LOWP:
                return GLSLSymbol::PRECQ_LOW;
        }
    }

    return GLSLSymbol::PRECQ_NONE;
}

GLSLVariableSymbol::StorageQualifier GLSLVisitor::visitNode(GLSLStorageQualifier* const storage) {
    if (storage == 0)
        return GLSLVariableSymbol::SQ_NONE;

    Token* const token = storage->getToken();
    if (token != 0) {
        switch (token->getTokenID()) {
            case GLSLTerminals::ID_ATTRIBUTE:
                return GLSLVariableSymbol::SQ_ATTRIBUTE;

            case GLSLTerminals::ID_CONST:
                return GLSLVariableSymbol::SQ_CONST;

            case GLSLTerminals::ID_VARYING:
                return ((storage->hasCentroid()) ? GLSLVariableSymbol::SQ_CENTROID_VARYING
                    : GLSLVariableSymbol::SQ_VARYING);

            case GLSLTerminals::ID_IN:
                return ((storage->hasCentroid()) ? GLSLVariableSymbol::SQ_CENTROID_IN
                    : GLSLVariableSymbol::SQ_IN);

            case GLSLTerminals::ID_OUT:
                return ((storage->hasCentroid()) ? GLSLVariableSymbol::SQ_CENTROID_OUT
                    : GLSLVariableSymbol::SQ_OUT);

            case GLSLTerminals::ID_UNIFORM:
                return GLSLVariableSymbol::SQ_UNIFORM;
        }
    }

    return GLSLVariableSymbol::SQ_NONE;
}

// Specifiers
//

GLSLVariableSymbol GLSLVisitor::visitNode(GLSLTypeSpecifier* const typeSpec) {
    if (typeSpec == 0)
        throw std::runtime_error("visitNode(GLSLTypeSpecifier): given specifier is NULL!");

    Token* const token = typeSpec->getTypeToken();
    if (token == 0)
        throw std::runtime_error("visitNode(GLSLTypeSpecifier): specifier does not specify a type!");

    // Determine internal type of the declaration
    //
    GLSLSymbol::InternalType type = GLSLSymbol::INTERNAL_VOID;
    switch (token->getTokenID()) {
        case GLSLTerminals::ID_VOID:
            break;

        case GLSLTerminals::ID_STRUCT:
            type = GLSLSymbol::INTERNAL_STRUCT;
            break;

        case GLSLTerminals::ID_BOOL:
        case GLSLTerminals::ID_BVEC2:
        case GLSLTerminals::ID_BVEC3:
        case GLSLTerminals::ID_BVEC4:
            type = GLSLSymbol::INTERNAL_BOOL;
            break;

        case GLSLTerminals::ID_FLOAT:
        case GLSLTerminals::ID_VEC2:
        case GLSLTerminals::ID_VEC3:
        case GLSLTerminals::ID_VEC4:
        case GLSLTerminals::ID_MAT2:
        case GLSLTerminals::ID_MAT3:
        case GLSLTerminals::ID_MAT4:
        case GLSLTerminals::ID_MAT2X2:
        case GLSLTerminals::ID_MAT2X3:
        case GLSLTerminals::ID_MAT2X4:
        case GLSLTerminals::ID_MAT3X2:
        case GLSLTerminals::ID_MAT3X3:
        case GLSLTerminals::ID_MAT3X4:
        case GLSLTerminals::ID_MAT4X2:
        case GLSLTerminals::ID_MAT4X3:
        case GLSLTerminals::ID_MAT4X4:
            type = GLSLSymbol::INTERNAL_FLOAT;
            break;

        case GLSLTerminals::ID_INT:
        case GLSLTerminals::ID_IVEC2:
        case GLSLTerminals::ID_IVEC3:
        case GLSLTerminals::ID_IVEC4:
            type = GLSLSymbol::INTERNAL_INT;
            break;

        case GLSLTerminals::ID_UINT:
        case GLSLTerminals::ID_UVEC2:
        case GLSLTerminals::ID_UVEC3:
        case GLSLTerminals::ID_UVEC4:
            type = GLSLSymbol::INTERNAL_UINT;
            break;

        default:
            if ((token->getTokenID() >= GLSLTerminals::ID_SAMPLER1D)
                && (token->getTokenID() <= GLSLTerminals::ID_USAMPLER2DMSARRAY))
            {
                type = GLSLSymbol::INTERNAL_SAMPLER;
            }
            break;
    }

    // Determine number of elements of internal type
    //
    int numElements = 1;
    switch (token->getTokenID()) {
        case GLSLTerminals::ID_STRUCT:
            numElements = 0;
            break;

        case GLSLTerminals::ID_BVEC2:
        case GLSLTerminals::ID_VEC2:
        case GLSLTerminals::ID_IVEC2:
        case GLSLTerminals::ID_UVEC2:
            numElements = 2;
            break;

        case GLSLTerminals::ID_BVEC3:
        case GLSLTerminals::ID_VEC3:
        case GLSLTerminals::ID_IVEC3:
        case GLSLTerminals::ID_UVEC3:
            numElements = 3;
            break;

        case GLSLTerminals::ID_BVEC4:
        case GLSLTerminals::ID_VEC4:
        case GLSLTerminals::ID_IVEC4:
        case GLSLTerminals::ID_UVEC4:
        case GLSLTerminals::ID_MAT2:
        case GLSLTerminals::ID_MAT2X2:
            numElements = 4;
            break;

        case GLSLTerminals::ID_MAT2X3:
        case GLSLTerminals::ID_MAT3X2:
            numElements = 6;
            break;

        case GLSLTerminals::ID_MAT2X4:
        case GLSLTerminals::ID_MAT4X2:
            numElements = 8;
            break;

        case GLSLTerminals::ID_MAT3:
        case GLSLTerminals::ID_MAT3X3:
            numElements = 9;
            break;

        case GLSLTerminals::ID_MAT3X4:
        case GLSLTerminals::ID_MAT4X3:
            numElements = 12;
            break;

        case GLSLTerminals::ID_MAT4:
        case GLSLTerminals::ID_MAT4X4:
            numElements = 16;
            break;

        case GLSLTerminals::ID_SAMPLER1D:
        case GLSLTerminals::ID_SAMPLER1DSHADOW:
        case GLSLTerminals::ID_SAMPLER1DARRAY:
        case GLSLTerminals::ID_SAMPLER1DARRAYSHADOW:
        case GLSLTerminals::ID_ISAMPLER1D:
        case GLSLTerminals::ID_ISAMPLER1DARRAY:
        case GLSLTerminals::ID_USAMPLER1D:
        case GLSLTerminals::ID_USAMPLER1DARRAY:
            numElements = 1;
            break;

        case GLSLTerminals::ID_SAMPLER2D:
        case GLSLTerminals::ID_SAMPLER2DSHADOW:
        case GLSLTerminals::ID_SAMPLER2DARRAY:
        case GLSLTerminals::ID_SAMPLER2DARRAYSHADOW:
        case GLSLTerminals::ID_ISAMPLER2D:
        case GLSLTerminals::ID_ISAMPLER2DARRAY:
        case GLSLTerminals::ID_USAMPLER2D:
        case GLSLTerminals::ID_USAMPLER2DARRAY:
        case GLSLTerminals::ID_SAMPLER2DRECT:
        case GLSLTerminals::ID_SAMPLER2DRECTSHADOW:
        case GLSLTerminals::ID_ISAMPLER2DRECT:
        case GLSLTerminals::ID_USAMPLER2DRECT:
        case GLSLTerminals::ID_SAMPLER2DMS:
        case GLSLTerminals::ID_ISAMPLER2DMS:
        case GLSLTerminals::ID_USAMPLER2DMS:
        case GLSLTerminals::ID_SAMPLER2DMSARRY:
        case GLSLTerminals::ID_ISAMPLER2DMSARRAY:
        case GLSLTerminals::ID_USAMPLER2DMSARRAY:
            numElements = 2;
            break;

        case GLSLTerminals::ID_SAMPLER3D:
        case GLSLTerminals::ID_ISAMPLER3D:
        case GLSLTerminals::ID_USAMPLER3D:
            numElements = 3;
            break;

        case GLSLTerminals::ID_SAMPLERCUBE:
        case GLSLTerminals::ID_ISAMPLERCUBE:
        case GLSLTerminals::ID_USAMPLERCUBE:
            numElements = 6;
            break;

        case GLSLTerminals::ID_SAMPLERBUFFER:
        case GLSLTerminals::ID_ISAMPLERBUFFER:
        case GLSLTerminals::ID_USAMPLEBUFFER:
            numElements = 0;
            break;
    }

    GLSLVariableSymbol metaSymbol("", type, numElements);

    GLSLTypeQualifier* const typeQuali = typeSpec->getTypeQualifier();
    if (typeQuali != 0) {
        GLSLVariableSymbol::InterpolationQualifier iq = visitNode(typeQuali->getInterpolationQualifier());
        GLSLVariableSymbol::StorageQualifier sq = visitNode(typeQuali->getStorageQualifier());

        metaSymbol.setInterpolation(iq);
        metaSymbol.setStorage(sq);
    }

    GLSLSymbol::PrecisionQualifier pq = visitNode(typeSpec->getPrecisionQualifier());
    metaSymbol.setPrecision(pq);

    return metaSymbol;
}

// Statement nodes
//

void GLSLVisitor::visitNode(GLSLStatement* const statement)
{
    if (! statement)
        return;

    switch (statement->getNodeType())
    {
        case GLSLNodeTypes::NODE_COMPOUND_STATEMENT:
            visitNode(static_cast<GLSLCompoundStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_SIMPLE_STATEMENT:
            visitNode(static_cast<GLSLSimpleStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_CASE_LABEL:
            visitNode(static_cast<GLSLCaseLabel* const>(statement));
            break;

        case GLSLNodeTypes::NODE_DECLARATION_STATEMENT:
            visitNode(static_cast<GLSLDeclarationStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_DO_WHILE_STATEMENT:
            visitNode(static_cast<GLSLDoWhileStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_FOR_STATEMENT:
            visitNode(static_cast<GLSLForStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_EXPRESSION_STATEMENT:
            visitNode(static_cast<GLSLExpressionStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_JUMP_STATEMENT:
            visitNode(static_cast<GLSLJumpStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_SELECTION_STATEMENT:
            visitNode(static_cast<GLSLSelectionStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_SWITCH_STATEMENT:
            visitNode(static_cast<GLSLSwitchStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_WHILE_STATEMENT:
            visitNode(static_cast<GLSLWhileStatement* const>(statement));
            break;

        default:
            break;
    }   // switch
}

void GLSLVisitor::visitNode(GLSLCompoundStatement* const statements)
{
    if (! statements)
        return;

    // TODO: add handling for creation of new scopes
    GLSLStatementList* const list = statements->getStatementList();
    if (list)
    {
        const std::vector<GLSLStatement*>& stmts = list->getStatements();
        for (size_t i = 0; i < stmts.size(); ++i)
            visitNode(stmts[i]);
    }
}

void GLSLVisitor::visitNode(GLSLSimpleStatement* const statement)
{
    if (! statement)
        return;

    switch (statement->getNodeType())
    {
        case GLSLNodeTypes::NODE_CASE_LABEL:
            visitNode(static_cast<GLSLCaseLabel* const>(statement));
            break;

        case GLSLNodeTypes::NODE_DECLARATION_STATEMENT:
            visitNode(static_cast<GLSLDeclarationStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_DO_WHILE_STATEMENT:
            visitNode(static_cast<GLSLDoWhileStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_FOR_STATEMENT:
            visitNode(static_cast<GLSLForStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_EXPRESSION_STATEMENT:
            visitNode(static_cast<GLSLExpressionStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_JUMP_STATEMENT:
            visitNode(static_cast<GLSLJumpStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_SELECTION_STATEMENT:
            visitNode(static_cast<GLSLSelectionStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_SWITCH_STATEMENT:
            visitNode(static_cast<GLSLSwitchStatement* const>(statement));
            break;

        case GLSLNodeTypes::NODE_WHILE_STATEMENT:
            visitNode(static_cast<GLSLWhileStatement* const>(statement));
            break;

        default:
            break;
    }   // switch
}

void GLSLVisitor::visitNode(GLSLCaseLabel* const lbl)
{
    if (! lbl)
        return;
    // TODO: add semantics
}

void GLSLVisitor::visitNode(GLSLDeclarationStatement* const decl)
{
    if (decl)
        visitNode(decl->getDeclaration());
}

void GLSLVisitor::visitNode(GLSLDoWhileStatement* const dwhl)
{
    if (dwhl)
    {
        visitNode(dwhl->getBody());
        delete visitNode(dwhl->getCondition());
    }
}

void GLSLVisitor::visitNode(GLSLExpressionStatement* const expr)
{
    if (expr)
        delete visitNode(expr->getExpression());
}

void GLSLVisitor::visitNode(GLSLForStatement* const fr)
{
    if (fr)
    {
        visitNode(fr->getInit());
        delete visitNode(fr->getCondition());
        delete visitNode(fr->getIterationExpr());
        visitNode(fr->getBody());
    }
}

void GLSLVisitor::visitNode(GLSLJumpStatement* const jmp)
{
    if (jmp)
        delete visitNode(jmp->getExpression());
}

void GLSLVisitor::visitNode(GLSLSelectionStatement* const sel)
{
    if (sel)
    {
        //GLSLExpression* const res =
            delete visitNode(sel->getCondition());
        // TODO: decide whether both case need to be evaluated / visited
        // depending on the evaluation of the condition. This can enable
        // short-circuit evaluation, which may be unwanted.
        visitNode(sel->getTrueStatement());
        visitNode(sel->getFalseStatement());
    }
}

void GLSLVisitor::visitNode(GLSLSwitchStatement* const swtch)
{
    if (swtch)
    {
        delete visitNode(swtch->getExpression());
        GLSLStatementList* const lst = swtch->getStatements();
        if (! lst)
            return;
        const std::vector<GLSLStatement*>& stmts = lst->getStatements();
        for (size_t i = 0; i < stmts.size(); ++i)
            visitNode(stmts[i]);
    }
}

void GLSLVisitor::visitNode(GLSLWhileStatement* const whl)
{
    if (whl)
    {
        delete visitNode(whl->getCondition());
        visitNode(whl->getBody());
    }
}

// Misc nodes
//

GLSLValue* GLSLVisitor::visitNode(GLSLCondition* const cond)
{
    if (! cond)
        return 0;

    GLSLTypeSpecifier* const spec = cond->getTypeSpecifier();
    IdentifierToken* const id = cond->getIdentifier();
    if ((spec) && (id))
    {
        try {
            GLSLVariableSymbol symbol = visitNode(spec);
            symbol.setID(id->getValue());
            symbol.setIsArray(false);
            symbol.setNumArrayElements(0);

            activeSymbols_->insertSymbol(new GLSLVariableSymbol(symbol));
        } catch (std::runtime_error& e) {
            log_ << "  Exeption in GLSLVisitor::visitNode(GLSLCondition*): " << e.what() << "\n";
        }
    }

    // If spec and id where both not NULL, the expression is the initializer
    //
    // NOTE: initializer is however only an assignemnt expression, not
    // an arbitrary expression to be more precise.
    return visitNode(cond->getExpression());
}

void GLSLVisitor::visitNode(GLSLExternalDeclaration* const extDecl)
{
    if (! extDecl)
        return;

    switch (extDecl->getNodeType())
    {
        case GLSLNodeTypes::NODE_DECLARATION:
            visitNode(static_cast<GLSLDeclaration* const>(extDecl));
            break;

        case GLSLNodeTypes::NODE_DECLARATION_LIST:
            visitNode(static_cast<GLSLDeclarationList* const>(extDecl));
            break;

        case GLSLNodeTypes::NODE_FUNCTION_DEFINITION:
            visitNode(static_cast<GLSLFunctionDefinition* const>(extDecl));
            break;

        default:
            std::cout << "GLSLVisitor::visitNode(): Error! Unknown external declaration type!\n";
            std::cout << "  node type = " << extDecl->getNodeType() << std::endl;
            break;
    }   // switch
}

void GLSLVisitor::visitNode(GLSLFunctionDefinition* const funcDef)
{
    if (! funcDef)
        return;

    GLSLFunctionPrototype* const proto = funcDef->getFunctionPrototype();
    if (proto)
    {
        visitNode(proto);

        if (proto->getName())
            pushSymbolTable(proto->getName()->getValue());
        else
            pushSymbolTable("anonymous-function");

        const std::vector<GLSLParameter*>& params = proto->getParameters();
        for (size_t i = 0; i < params.size(); ++i)
            visitNode(params[i]);

        GLSLCompoundStatement* const statements = funcDef->getStatements();
        visitNode(statements);
        popSymbolTable();
    }
}

void GLSLVisitor::visitNode(GLSLFunctionPrototype* const funcProto)
{
    if (! funcProto)
        return;

    // TODO: extract function name, return value, parameter count etc.
    // and insert function symbol into symboltable.
    //IdentifierToken* const funcName = funcProto->getName();
}

void GLSLVisitor::visitNode(GLSLParameter* const param)
{
    if (! param)
        return;

    GLSLTypeSpecifier* const spec = param->getSpecifier();
    GLSLVariable* const var = param->getName();

    if ((spec) && (var))
    {
        // TODO: handle qualifier and variable

        GLSLVariableSymbol meta = visitNode(param->getSpecifier());
        int arrCount = 0; // TODO: evaluate expression var->getNumArrayElements();
        GLSLVariableSymbol* sym = new GLSLVariableSymbol(var->getName(),
            meta.getInternalType(), arrCount , true);

        sym->setIsArray(var->isArray());

        activeSymbols_->insertSymbol(sym);
    }
}

void GLSLVisitor::visitNode(GLSLTranslation* const trans)
{
    if (! trans)
        return;

    const std::vector<GLSLExternalDeclaration*>& ext = trans->getExternalDeclarations();
    for (size_t i = 0; i < ext.size(); ++i)
        visitNode(ext[i]);
}

// ----------------------------------------------------------------------------

std::vector<GLSLAnnotation*> GLSLVisitor::processAnnotation(AnnotationToken* const annotation)
{
    std::vector<GLSLAnnotation*> annotations;

    if (annotation != 0) {
        AnnotationParser annotationParser(annotation->getParts());
        //annotationParser.setDebugging(true);

        ParseTreeNode* const root = annotationParser.parse();
        if (root != 0) {
            AnnotationVisitor av;
            annotations = av.getAnnotations(root);
            delete root;
        } else {
            log_ << "Annotation parser log:\n------------------------\n";
            log_ << annotationParser.getLog().str() << "\n";
        }
    }

    return annotations;
}

}   // namespace glslparser

}   // namespace voreen
