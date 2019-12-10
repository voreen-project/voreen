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

#include "plotfunctiongrammar.h"

#include "plotfunctionterminals.h"

namespace voreen {

namespace glslparser {

PlotFunctionGrammar::PlotFunctionGrammar()
    : setup_(false)
{
    setup(PlotFunctionTerminals());
}

void PlotFunctionGrammar::setup(const TerminalsMap& terminals) {
    if (setup_)
        return;

    setup_ = true;
    Grammar::setup(terminals);
    addTerminal("$START$" , PlotFunctionTerminals::ID_START , GrammarSymbol::SYMBOL_START);
    //addTerminal("$END$"   , PlotFunctionTerminals::ID_END   , GrammarSymbol::SYMBOL_END);
    //addTerminal("$EMPTY$" , PlotFunctionTerminals::ID_EMPTY , GrammarSymbol::SYMBOL_EMPTY);

    addProduction("$START$", "expression");

    addProduction("expression", "term");
    addProduction("expression", "term-list");
    addProduction("expression", "term-list ; term");
    addProduction("expression", "EMPTY");

    addProduction("term-list", "termwithinterval");
    addProduction("term-list", "term-list ; termwithinterval ");

    addProduction("termwithinterval", "term : interval");

    addProduction("interval", "range");
    addProduction("interval", "interval | range");

    addProduction("range", "- range");
    addProduction("range" , "leftrange , rightrange");

    addProduction("leftrange","( interval-term");
    addProduction("leftrange","[ interval-term");

    addProduction("rightrange","interval-term )");
    addProduction("rightrange","interval-term ]");

    addProduction("term", "term + addend");
    addProduction("term", "term - addend");
    addProduction("term", "+ addend");
    addProduction("term", "- addend");
    addProduction("term", "- - addend");
    addProduction("term", "term + - addend");
    addProduction("term", "term - + addend");
    addProduction("term", "term - - addend");
    addProduction("term", "addend");

    addProduction("addend", "factor");
    addProduction("addend", "addend * factor");
    addProduction("addend", "addend / factor");

    addProduction("factor", "variable");
    addProduction("factor", "variable ^ factor");

    //addProduction("variable", "+ variable");
    //addProduction("variable", "- variable");
    addProduction("variable", "atom");
    addProduction("variable", "bracket-term");

    addProduction("bracket-term", "( term )");

    addProduction("atom","function bracket-term");
    addProduction("atom", "single-value");
    //addProduction("atom", "- atom");
    //addProduction("atom", "+ atom");
    addProduction("single-value", "INTCONST");
    addProduction("single-value", "FLOATCONST");
    addProduction("single-value", "VARIABLE");

    addProduction("function", "FUNCTION");

    addProduction("interval-term", "interval-term + interval-addend");
    addProduction("interval-term", "interval-term - interval-addend");
    addProduction("interval-term", "+ interval-addend");
    addProduction("interval-term", "- interval-addend");
    addProduction("interval-term", "- - interval-addend");
    addProduction("interval-term", "interval-term + - interval-addend");
    addProduction("interval-term", "interval-term - + interval-addend");
    addProduction("interval-term", "interval-term - - interval-addend");
    addProduction("interval-term", "interval-addend");

    addProduction("interval-addend", "interval-factor");
    addProduction("interval-addend", "interval-addend * interval-factor");
    addProduction("interval-addend", "interval-addend / interval-factor");

    addProduction("interval-factor", "interval-variable");
    addProduction("interval-factor", "interval-variable ^ interval-factor");

    //addProduction("interval-variable", "+ interval-variable");
    //addProduction("interval-variable", "- interval-variable");
    addProduction("interval-variable", "interval-atom");
    addProduction("interval-variable", "interval-bracket-term");

    addProduction("interval-bracket-term", "( interval-term )");

    addProduction("interval-atom","function interval-bracket-term");
    addProduction("interval-atom", "interval-single-value");
    //addProduction("interval-atom", "+ interval-atom");
    //addProduction("interval-atom", "- interval-atom");
    addProduction("interval-single-value", "INTCONST");
    addProduction("interval-single-value", "FLOATCONST");
}

}   // namespace glslparser

}   // namespace voreen
