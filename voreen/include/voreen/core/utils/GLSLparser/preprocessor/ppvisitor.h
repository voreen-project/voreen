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

#ifndef VRN_PPVISITOR_H
#define VRN_PPVISITOR_H

#include "voreen/core/utils/GLSLparser/symboltable.h"
#include "voreen/core/utils/GLSLparser/preprocessor/ppmacro.h"
#include "voreen/core/utils/GLSLparser/preprocessor/ppstatement.h"
#include "voreen/core/utils/GLSLparser/preprocessor/ppsymbol.h"

#include <sstream>

namespace voreen {

namespace glslparser {

class PreprocessorVisitor : public ParseTreeVisitor {
public:
    PreprocessorVisitor();
    virtual ~PreprocessorVisitor();

    virtual bool visit(ParseTreeNode* const node);

    const std::string& getIncludePath() const { return includePath_; }
    const std::string& getShaderHeader() const { return shaderHeader_; }

    void setIncludePath(const std::string& includePath);
    void setShaderHeader(const std::string& shaderHeader) { shaderHeader_ = shaderHeader; }

    std::ostringstream& translate(std::istream* const is, const std::string& path = "");
    std::ostringstream& translate(ParseTreeNode* const node);

private:
    int visitNode(Expression* const n);
    int visitNode(ExpressionList* const n);
    int visitNode(ParenthesisExpression* const n);
    int visitNode(ArithmeticExpression* const n);
    int visitNode(IntConstant* const n);
    int visitNode(UnaryExpression* const n);
    int visitNode(BinaryExpression* const n);
    int visitNode(LogicalExpression* const n);
    int visitNode(LogicalBinaryExpression* const n);
    int visitNode(DefinedOperator* const n);

    int visitNode(TextNode* const n);
    int visitNode(TokenList* const n);
    int visitNode(Macro* const n);

    int visitNode(Statement* const n);
    int visitNode(IdentifierList* const n);
    bool visitNode(DefineDirective* const n);
    void visitNode(ErrorDirective* const n);
    void visitNode(ExtensionDirective* const n);
    bool visitNode(IncludeDirective* const n);
    int visitNode(LineDirective* const n);
    void visitNode(PragmaDirective* const n);
    bool visitNode(UndefineDirective* const n);
    int visitNode(VersionDirective* const n);
    void visitNode(ConditionalDirective* const n);

    /**
     * Adjusts the given file name to match the requirements of the OS and
     * returns an OS-dependent absolute file name.
     *
     * @param   fileName    String containg a (relative) path to the file with
     *                      slashes ('/') as path separators.
     * @return String containing the absolute path and the file name with OS
     *          specific path separators.
     */
    std::string adjustIncludeFile(const std::string& fileName);

    /**
     * Expands all occurences of the macro given by symbol within thestring
     * passed.
     *
     * @param   input   The input string containing macro-invocations which shall
     *                  be expanded.
     * @param   symbol  The symbol from the symbol table representing the macro
     *                  definition.
     * @return  true if at least one expansion has been perfomed, false otherwise
     */
    bool expandMacro(std::string& input, PreprocessorSymbol* const symbol) const;

    /**
     * Extracts the parameters from macro-calls. The parameters are are fromatted
     * in the following fashion:
     *
     * <paramter-list> ::= <whitespace> ( <whitespace> <paramters> <whitespace> )
     * <parameters> ::= EMPTY | PARAMETER | <parameters> <whitespace> , <whitespace> PARAMETER
     *
     * @param   input   String to scane for parameters
     * @param   offset  Offset indicating where to start within the input
     *
     * @return  A vector containing the names of the parameters in the order of there
     *          appearance within the input string.
     */
    std::vector<std::string> readParameters(const std::string& input,
        size_t& offset) const;

    /**
     * Replaces all occurrences of string 'search' within string 'in' by the string
     * 'replacement'.
     *
     * @return  A string derived from 'in' where all occurrences of string 'search'
     * are replaced by string 'replacement'. If 'search' does not appear within 'in',
     * the retrun value is an exact copy of 'in'.
     */
    std::string replace(const std::string& in, const std::string& search,
        const std::string& replacement) const;

private:
    std::ostringstream translation_;

    typedef SymbolTable<PreprocessorSymbol> PPSymbolMap;
    PPSymbolMap symbols_;

    PreprocessorTerminals terminals_;

    std::string includePath_;
    std::string shaderHeader_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
