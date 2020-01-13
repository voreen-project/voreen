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

#include "voreen/core/utils/GLSLparser/generator/ppgrammar.h"

#include "voreen/core/utils/GLSLparser/preprocessor/ppterminals.h"

namespace voreen {

namespace glslparser {

PreprocessorGrammar::PreprocessorGrammar() : setup_(false) {
    setup(PreprocessorTerminals());
}

void PreprocessorGrammar::setup(const TerminalsMap& terminals) {
    if (setup_)
        return;

    setup_ = true;
    Grammar::setup(terminals);
    addTerminal("$START$", PreprocessorTerminals::ID_START, GrammarSymbol::SYMBOL_START);

    addProduction("$START$", "program");

    addProduction("program", "particle-list");
    addProduction("program", "expression-parsing");

    addProduction("expression-parsing", "pirate constant-expression");

    addProduction("particle-list", "particle");
    addProduction("particle-list", "particle-list particle");

    addProduction("particle", "directive");
    addProduction("particle", "text");
    addProduction("particle", "newline");

    addProduction("directive", "def-directive");
    addProduction("directive", "undef-directive");
    addProduction("directive", "error-directive");
    addProduction("directive", "extension-directive");
    addProduction("directive", "include-directive");
    addProduction("directive", "line-directive");
    addProduction("directive", "null-directive");
    addProduction("directive", "pragma-directive");
    addProduction("directive", "version-directive");
    addProduction("directive", "conditional-directive");

    addProduction("def-directive", "#define identifier newline");
    addProduction("def-directive", "#define identifier token-list newline");
    addProduction("def-directive", "#define identifier ( ) newline");
    addProduction("def-directive", "#define identifier ( ) token-list newline");
    addProduction("def-directive", "#define identifier ( formals-list ) newline");
    addProduction("def-directive", "#define identifier ( formals-list ) token-list newline");

    addProduction("token-list", "token");
    addProduction("token-list", "token-list token");

    addProduction("formals-list", "identifier");
    addProduction("formals-list", "formals-list , identifier");

    addProduction("error-directive", "#error newline");
    addProduction("error-directive", "#error token-list newline");

    addProduction("extension-directive", "#extension newline");
    addProduction("extension-directive", "#extension token-list newline");

    addProduction("include-directive", "#include \" string \" newline");

    addProduction("line-directive", "#line int-constant newline");
    addProduction("line-directive", "#line int-constant \" string \" newline");

    addProduction("null-directive", "# newline");

    addProduction("pragma-directive", "#pragma token-list newline");

    addProduction("undef-directive", "#undef identifier newline");

    addProduction("version-directive", "#version int-constant newline");
    addProduction("version-directive", "#version int-constant identifier newline");

    addProduction("conditional-directive", "if-part endif-directive");
    addProduction("conditional-directive", "if-part else-part endif-directive");

    addProduction("if-part", "if-directive particle-list");
    addProduction("if-part", "ifdef-directive particle-list");
    addProduction("if-part", "ifndef-directive particle-list");

    addProduction("elif-parts", "elif-part");
    addProduction("elif-parts", "elif-parts elif-part");

    addProduction("elif-part", "elif-directive particle-list");

    addProduction("else-part", "elif-parts");
    addProduction("else-part", "elif-parts else-directive particle-list");
    addProduction("else-part", "else-directive particle-list");

    addProduction("if-directive", "#if constant-expression newline");
    addProduction("ifdef-directive", "#ifdef identifier newline");
    addProduction("ifndef-directive", "#ifndef identifier newline");
    addProduction("else-directive", "#else newline");
    addProduction("elif-directive", "#elif constant-expression newline");
    addProduction("endif-directive", "#endif newline");

    addProduction("constant-expression", "logical-or-expression");

    addProduction("macro-evaluation", "identifier");
    addProduction("macro-evaluation", "identifier ( )");
    addProduction("macro-evaluation", "identifier ( parameter-list )");

    addProduction("parameter-list", "constant-expression");
    addProduction("parameter-list", "parameter-list , constant-expression");

    addProduction("defined-operator", "defined ( identifier )");
    addProduction("defined-operator", "defined identifier");

    addProduction("logical-or-expression", "logical-and-expression");
    addProduction("logical-or-expression", "logical-or-expression || logical-and-expression");

    addProduction("logical-and-expression", "inclusive-or-expression");
    addProduction("logical-and-expression", "logical-and-expression && inclusive-or-expression");

    addProduction("inclusive-or-expression", "exclusive-or-expression");
    addProduction("inclusive-or-expression", "inclusive-or-expression | exclusive-or-expression");

    addProduction("exclusive-or-expression", "and-expression");
    addProduction("exclusive-or-expression", "exclusive-or-expression ^ and-expression");

    addProduction("and-expression", "equality-expression");
    addProduction("and-expression", "and-expression & equality-expression");

    addProduction("equality-expression", "relational-expression");
    addProduction("equality-expression", "equality-expression == relational-expression");
    addProduction("equality-expression", "equality-expression != relational-expression");

    addProduction("relational-expression", "shift-expression");
    addProduction("relational-expression", "relational-expression < shift-expression");
    addProduction("relational-expression", "relational-expression > shift-expression");
    addProduction("relational-expression", "relational-expression <= shift-expression");
    addProduction("relational-expression", "relational-expression >= shift-expression");

    addProduction("shift-expression", "additive-expression");
    addProduction("shift-expression", "shift-expression << additive-expression");
    addProduction("shift-expression", "shift-expression >> additive-expression");

    addProduction("additive-expression", "multiplicative-expression");
    addProduction("additive-expression", "additive-expression + multiplicative-expression");
    addProduction("additive-expression", "additive-expression - multiplicative-expression");

    addProduction("multiplicative-expression", "unary-expression");
    addProduction("multiplicative-expression", "multiplicative-expression * unary-expression");
    addProduction("multiplicative-expression", "multiplicative-expression / unary-expression");
    addProduction("multiplicative-expression", "multiplicative-expression % unary-expression");

    addProduction("unary-expression", "primary-expression");
    addProduction("unary-expression", "+ unary-expression");
    addProduction("unary-expression", "- unary-expression");
    addProduction("unary-expression", "! unary-expression");
    addProduction("unary-expression", "~ unary-expression");

    addProduction("primary-expression", "defined-operator");
    addProduction("primary-expression", "macro-evaluation");
    addProduction("primary-expression", "int-constant");
    addProduction("primary-expression", "( constant-expression )");
}

}   // namespace glslparser

}   // namespace voreen
