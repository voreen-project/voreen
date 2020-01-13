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

#include "voreen/core/utils/GLSLparser/generator/annotationgrammar.h"

#include "voreen/core/utils/GLSLparser/annotations/annotationterminals.h"

namespace voreen {

namespace glslparser {

AnnotationGrammar::AnnotationGrammar()
    : setup_(false)
{
    setup(AnnotationTerminals());
}

void AnnotationGrammar::setup(const TerminalsMap& terminals) {
    if (setup_)
        return;

    setup_ = true;
    Grammar::setup(terminals);
    addTerminal("$START$", AnnotationTerminals::ID_START, GrammarSymbol::SYMBOL_START);

    addProduction("$START$", "annotation");

    addProduction("annotation", "tag-list");

    addProduction("tag-list", "tag");
    addProduction("tag-list", "tag-list tag");
    addProduction("tag-list", "tag-list , tag");

    addProduction("tag", "tag-assignment [ value-list ]");
    addProduction("tag", "tag-assignment single-value");

    addProduction("tag-assignment", "tag-name =");
    addProduction("tag-assignment", "tag-assignment tag-name = ");

    addProduction("tag-name", "@ IDENTIFIER");

    addProduction("value-list", "single-value");
    addProduction("value-list", "value-list , single-value");

    addProduction("single-value", "true");
    addProduction("single-value", "false");
    addProduction("single-value", "INTCONST");
    addProduction("single-value", "FLOATCONST");
    addProduction("single-value", "\" STRING \"");
    addProduction("single-value", "' CHARACTER '");
}

}   // namespace glslparser

}   // namespace voreen
