/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/core/utils/GLSLparser/generator/glslgrammar.h"

namespace voreen {

namespace glslparser {

GLSLGrammar::GLSLGrammar()
    : setup_(false)
{
    setup(terminals_);
}

void GLSLGrammar::setup(const TerminalsMap& terminals) {
    if (setup_)
        return;

    setup_ = true;

    // Do not call inherited symbol for we do not need to have all the
    // reserved keyword in the grammar and we want to use a surrogate terminal
    // NATIVE-TYPE instead of all terminals for int, float, vec2, etc...
    //
    const std::map<std::string, int>& terminalSymbols = terminals.getLexemeMap();
    for (std::map<std::string, int>::const_iterator it = terminalSymbols.begin();
        it != terminalSymbols.end(); ++it)
    {
        if (it->second == GLSLTerminals::ID_RESERVED_KEYWORD)
            continue;
        else if ((it->second >= GLSLTerminals::ID_BOOL) && (it->second <= GLSLTerminals::ID_USAMPLER2DMSARRAY))
            continue;
        else
            addTerminal(it->first, it->second);
    }

    // add abstract symbol to make the parser report an error if a reserved keyword is used.
    //
    addTerminal("RESERVED-KEYWORD", GLSLTerminals::ID_RESERVED_KEYWORD);

    addTerminal("$START$", GLSLTerminals::ID_START, GrammarSymbol::SYMBOL_START);

    addProduction("$START$", "program");

    addProduction("program", "translation-unit");

    addProduction("translation-unit", "external-declaration");
    addProduction("translation-unit", "translation-unit external-declaration");

    addProduction("external-declaration", "function-definition");
    addProduction("external-declaration", "declaration");
    addProduction("external-declaration", "enhanced-declaration");

    addProduction("function-definition", "function-prototype compound-statement-no-new-scope");

    addProduction("statement", "compound-statement");
    addProduction("statement", "simple-statement");

    // Grammar Note: labeled statements for SWITCH only; 'goto' is not supported.
    addProduction("simple-statement", "declaration-statement");
    addProduction("simple-statement", "expression-statement");
    addProduction("simple-statement", "selection-statement");
    addProduction("simple-statement", "switch-statement");
    addProduction("simple-statement", "case-label");
    addProduction("simple-statement", "iteration-statement");
    addProduction("simple-statement", "jump-statement");

    addProduction("compound-statement", "{ }");
    addProduction("compound-statement", "{ statement-list }");

    addProduction("declaration-statement", "declaration");
    addProduction("declaration-statement", "enhanced-declaration");

    addProduction("statement-no-new-scope", "compound-statement-no-new-scope");
    addProduction("statement-no-new-scope", "simple-statement");

    addProduction("compound-statement-no-new-scope", "{ }");
    addProduction("compound-statement-no-new-scope", "{ statement-list }");

    addProduction("statement-list", "statement");
    addProduction("statement-list", "statement-list statement");

    addProduction("expression-statement", ";");
    addProduction("expression-statement", "expression ;");

    addProduction("selection-statement", "if ( expression ) statement");
    addProduction("selection-statement", "if ( expression ) statement else statement");

    addProduction("condition", "expression");
    addProduction("condition", "fully-specified-type IDENTIFIER = initializer");

    addProduction("switch-statement", "switch ( expression ) { statement-list }");
    addProduction("switch-statement", "switch ( expression ) { }");

    addProduction("case-label", "case expression :");
    addProduction("case-label", "default :");

    addProduction("iteration-statement", "while ( condition ) statement-no-new-scope");
    addProduction("iteration-statement", "do statement while ( expression ) ;");
    addProduction("iteration-statement",
        "for ( for-init-statement for-rest-statement ) statement-no-new-scope");

    addProduction("for-init-statement", "expression-statement");
    addProduction("for-init-statement", "declaration-statement");

    addProduction("for-rest-statement", ";");
    addProduction("for-rest-statement", "condition ;");
    addProduction("for-rest-statement", "; expression");
    addProduction("for-rest-statement", "condition ; expression");

    addProduction("jump-statement", "continue ;");
    addProduction("jump-statement", "break ;");
    addProduction("jump-statement", "return ;");
    addProduction("jump-statement", "return expression ;");
    addProduction("jump-statement", "discard ;"); // Fragment shader only

    addProduction("primary-expression", "IDENTIFIER");
    addProduction("primary-expression", "INTCONSTANT");
    addProduction("primary-expression", "UINTCONSTANT");
    addProduction("primary-expression", "FLOATCONSTANT");
    addProduction("primary-expression", "false");
    addProduction("primary-expression", "true");
    addProduction("primary-expression", "( expression )");

    addProduction("postfix-expression", "primary-expression");
    addProduction("postfix-expression", "postfix-expression [ integer-expression ]");
    addProduction("postfix-expression", "function-call-or-method");
    addProduction("postfix-expression", "postfix-expression . FIELD-SELECTION");
    addProduction("postfix-expression", "postfix-expression ++");
    addProduction("postfix-expression", "postfix-expression --");

    addProduction("integer-expression", "expression");

    addProduction("function-call-or-method", "function-call-generic");
    addProduction("function-call-or-method", "postfix-expression . function-call-generic");

    addProduction("function-call-generic", "function-call-header void )");
    addProduction("function-call-generic", "function-call-header )");
    addProduction("function-call-generic", "function-call-header-with-parameters )");

    addProduction("function-call-header-with-parameters", "function-call-header assignment-expression");
    addProduction("function-call-header-with-parameters",
        "function-call-header-with-parameters , assignment-expression");

    addProduction("function-call-header", "type-specifier (");
    addProduction("function-call-header", "IDENTIFIER (");
    addProduction("function-call-header", "FIELD-SELECTION (");

    addProduction("unary-expression", "postfix-expression");
    addProduction("unary-expression", "++ unary-expression");
    addProduction("unary-expression", "-- unary-expression");
    addProduction("unary-expression", "+ unary-expression");
    addProduction("unary-expression", "- unary-expression");
    addProduction("unary-expression", "! unary-expression");
    addProduction("unary-expression", "~ unary-expression");

    addProduction("multiplicative-expression", "unary-expression");
    addProduction("multiplicative-expression", "multiplicative-expression * unary-expression");
    addProduction("multiplicative-expression", "multiplicative-expression / unary-expression");
    addProduction("multiplicative-expression", "multiplicative-expression % unary-expression");

    addProduction("additive-expression", "multiplicative-expression");
    addProduction("additive-expression", "additive-expression + multiplicative-expression");
    addProduction("additive-expression", "additive-expression - multiplicative-expression");

    addProduction("shift-expression", "additive-expression");
    addProduction("shift-expression", "shift-expression << additive-expression");
    addProduction("shift-expression", "shift-expression >> additive-expression");

    addProduction("relational-expression", "shift-expression");
    addProduction("relational-expression", "relational-expression < shift-expression");
    addProduction("relational-expression", "relational-expression > shift-expression");
    addProduction("relational-expression", "relational-expression <= shift-expression");
    addProduction("relational-expression", "relational-expression >= shift-expression");

    addProduction("equality-expression", "relational-expression");
    addProduction("equality-expression", "equality-expression == relational-expression");
    addProduction("equality-expression", "equality-expression != relational-expression");

    addProduction("and-expression", "equality-expression");
    addProduction("and-expression", "and-expression & equality-expression");

    addProduction("exclusive-or-expression", "and-expression");
    addProduction("exclusive-or-expression", "exclusive-or-expression ^ and-expression");

    addProduction("inclusive-or-expression", "exclusive-or-expression");
    addProduction("inclusive-or-expression", "inclusive-or-expression | exclusive-or-expression");

    addProduction("logical-and-expression", "inclusive-or-expression");
    addProduction("logical-and-expression", "logical-and-expression && inclusive-or-expression");

    addProduction("logical-xor-expression", "logical-and-expression");
    addProduction("logical-xor-expression", "logical-xor-expression ^^ logical-and-expression");

    addProduction("logical-or-expression", "logical-xor-expression");
    addProduction("logical-or-expression", "logical-or-expression || logical-xor-expression");

    addProduction("conditional-expression", "logical-or-expression");
    addProduction("conditional-expression", "logical-or-expression ? expression : assignment-expression");

    addProduction("assignment-expression", "conditional-expression");
    addProduction("assignment-expression", "unary-expression = assignment-expression");
    addProduction("assignment-expression", "unary-expression *= assignment-expression");
    addProduction("assignment-expression", "unary-expression /= assignment-expression");
    addProduction("assignment-expression", "unary-expression %= assignment-expression");
    addProduction("assignment-expression", "unary-expression += assignment-expression");
    addProduction("assignment-expression", "unary-expression -= assignment-expression");
    addProduction("assignment-expression", "unary-expression <<= assignment-expression");
    addProduction("assignment-expression", "unary-expression >>= assignment-expression");
    addProduction("assignment-expression", "unary-expression &= assignment-expression");
    addProduction("assignment-expression", "unary-expression ^= assignment-expression");
    addProduction("assignment-expression", "unary-expression |= assignment-expression");

    addProduction("expression", "assignment-expression");
    addProduction("expression", "expression , assignment-expression");

    addProduction("constant-expression", "conditional-expression");

    addProduction("enhanced-declaration", "LEADING-ANNOTATION declaration");
    addProduction("enhanced-declaration", "declaration TRAILING-ANNOTATION");
    addProduction("enhanced-declaration", "LEADING-ANNOTATION declaration TRAILING-ANNOTATION");

    addProduction("declaration", "function-prototype ;");
    addProduction("declaration", "init-declarator-list ;");
    addProduction("declaration", "precision precision-qualifier type-specifier-no-prec ;");
    addProduction("declaration", "type-qualifier IDENTIFIER { struct-declaration-list } ;");
    addProduction("declaration", "type-qualifier IDENTIFIER { struct-declaration-list } IDENTIFIER ;");
    addProduction("declaration", "type-qualifier IDENTIFIER { struct-declaration-list } IDENTIFIER [ ] ;");
    addProduction("declaration",
        "type-qualifier IDENTIFIER { struct-declaration-list } IDENTIFIER [ constant-expression ] ;");
    addProduction("declaration", "type-qualifier ;");

    addProduction("function-prototype", "function-header )");
    addProduction("function-prototype", "function-header-with-parameters )");

    addProduction("function-header-with-parameters", "function-header parameter-declaration");
    addProduction("function-header-with-parameters", "function-header-with-parameters , parameter-declaration");

    addProduction("function-header", "fully-specified-type IDENTIFIER (");

    addProduction("parameter-declarator", "type-specifier IDENTIFIER");
    addProduction("parameter-declarator", "type-specifier IDENTIFIER [ constant-expression ]");

    addProduction("parameter-declaration", "const parameter-declarator");
    addProduction("parameter-declaration", "const parameter-qualifier parameter-declarator");
    addProduction("parameter-declaration", "parameter-declarator");
    addProduction("parameter-declaration", "parameter-qualifier parameter-declarator");
    addProduction("parameter-declaration", "const parameter-type-specifier");
    addProduction("parameter-declaration", "const parameter-qualifier parameter-type-specifier");
    addProduction("parameter-declaration", "parameter-type-specifier");
    addProduction("parameter-declaration", "parameter-qualifier parameter-type-specifier");

    addProduction("parameter-qualifier", "in");
    addProduction("parameter-qualifier", "out");
    addProduction("parameter-qualifier", "inout");

    addProduction("parameter-type-specifier", "type-specifier");

    addProduction("init-declarator-list", "single-declaration");
    addProduction("init-declarator-list", "single-declaration");
    addProduction("init-declarator-list", "init-declarator-list , IDENTIFIER");
    addProduction("init-declarator-list", "init-declarator-list , IDENTIFIER [ ]");
    addProduction("init-declarator-list", "init-declarator-list , IDENTIFIER [ constant-expression ]");
    addProduction("init-declarator-list", "init-declarator-list , IDENTIFIER [ ] = initializer");
    addProduction("init-declarator-list",
        "init-declarator-list , IDENTIFIER [ constant-expression ] = initializer");
    addProduction("init-declarator-list", "init-declarator-list , IDENTIFIER = initializer");

    addProduction("single-declaration", "fully-specified-type");
    addProduction("single-declaration", "fully-specified-type IDENTIFIER");
    addProduction("single-declaration", "fully-specified-type IDENTIFIER [ ]");
    addProduction("single-declaration", "fully-specified-type IDENTIFIER [ constant-expression ]");
    addProduction("single-declaration", "fully-specified-type IDENTIFIER [ ] = initializer");
    addProduction("single-declaration",
        "fully-specified-type IDENTIFIER [ constant-expression ] = initializer");
    addProduction("single-declaration", "fully-specified-type IDENTIFIER = initializer");

    addProduction("fully-specified-type", "type-specifier");
    addProduction("fully-specified-type", "type-qualifier type-specifier");

    addProduction("interpolation-qualifier", "smooth");
    addProduction("interpolation-qualifier", "flat");
    addProduction("interpolation-qualifier", "noperspective");

    addProduction("layout-qualifier", "layout ( layout-qualifier-id-list )");

    addProduction("layout-qualifier-id-list", "IDENTIFIER");
    addProduction("layout-qualifier-id-list", "IDENTIFIER = INTCONSTANT");
    addProduction("layout-qualifier-id-list", "layout-qualifier-id-list , IDENTIFIER");
    addProduction("layout-qualifier-id-list", "layout-qualifier-id-list , IDENTIFIER = INTCONSTANT");

    addProduction("type-qualifier", "storage-qualifier");
    addProduction("type-qualifier", "layout-qualifier");
    addProduction("type-qualifier", "layout-qualifier storage-qualifier");
    addProduction("type-qualifier", "interpolation-qualifier storage-qualifier");
    addProduction("type-qualifier", "interpolation-qualifier");
    addProduction("type-qualifier", "invariant storage-qualifier");
    addProduction("type-qualifier", "invariant interpolation-qualifier storage-qualifier");
    addProduction("type-qualifier", "invariant");

    addProduction("storage-qualifier", "const");
    addProduction("storage-qualifier", "attribute"); // Vertex only.
    addProduction("storage-qualifier", "varying");
    addProduction("storage-qualifier", "centroid varying");
    addProduction("storage-qualifier", "in");
    addProduction("storage-qualifier", "out");
    addProduction("storage-qualifier", "centroid in");
    addProduction("storage-qualifier", "centroid out");
    addProduction("storage-qualifier", "uniform");

    addProduction("type-specifier", "type-specifier-no-prec");
    addProduction("type-specifier", "precision-qualifier type-specifier-no-prec");

    addProduction("precision-qualifier", "highp");
    addProduction("precision-qualifier", "mediump");
    addProduction("precision-qualifier", "lowp");

    addProduction("type-specifier-no-prec", "type-specifier-nonarray");
    addProduction("type-specifier-no-prec", "type-specifier-nonarray [ ]");
    addProduction("type-specifier-no-prec", "type-specifier-nonarray [ constant-expression ]");

    addProduction("type-specifier-nonarray", "void");
    addProduction("type-specifier-nonarray", "NATIVE-TYPE");
    addProduction("type-specifier-nonarray", "struct-specifier");
    addProduction("type-specifier-nonarray", "TYPE-NAME");

    addProduction("struct-specifier", "struct IDENTIFIER { struct-declaration-list }");
    addProduction("struct-specifier", "struct { struct-declaration-list }");

    addProduction("struct-declaration-list", "struct-declaration");
    addProduction("struct-declaration-list", "struct-declaration-list struct-declaration");

    addProduction("struct-declaration", "type-specifier struct-declarator-list ;");
    addProduction("struct-declaration", "type-qualifier type-specifier struct-declarator-list ;");

    addProduction("struct-declarator-list", "struct-declarator");
    addProduction("struct-declarator-list", "struct-declarator-list , struct-declarator");

    addProduction("struct-declarator", "IDENTIFIER");
    addProduction("struct-declarator", "IDENTIFIER [ ]");
    addProduction("struct-declarator", "IDENTIFIER [ constant-expression ]");

    addProduction("initializer", "assignment-expression");
}

}   // namespace glslparser

}   // namespace voreen
