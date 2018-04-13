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

#include "voreen/core/utils/GLSLparser/preprocessor/ppmacro.h"

#include "voreen/core/utils/GLSLparser/lexer.h"
#include "voreen/core/utils/GLSLparser/preprocessor/ppsymbol.h"

namespace voreen {

namespace glslparser {

void Macro::tokenize(TokenList* const list) const {
    Expression::tokenize(list);
    if (parameters_ != 0) {
        for (size_t i = 0; i < parameters_->size(); ++i)
            (*parameters_)[i]->tokenize(list);
    }
}

std::list<Token*> Macro::expandMacro(PreprocessorSymbol* const symbol) {
    typedef std::list<Token*> TList;

    if ((symbol == 0) || (parameters_ == 0) || (symbol->getNumFormals() != parameters_->size()))
        return TList();

    TList body = symbol->getBody()->getCopy();
    const std::vector<std::string>& formals = symbol->getFormals();

    for (size_t i = 0; i < parameters_->size(); ++i) {
        TokenList* tokenList = new TokenList();
        (*parameters_)[i]->tokenize(tokenList);

        for (TList::iterator it = body.begin(); it != body.end(); ) {
            if (IdentifierToken* id = dynamic_cast<IdentifierToken*>(*it)) {
                if (id->getValue() == formals[i]) {
                    TList expansion = tokenList->getCopy();
                    body.insert(it, expansion.begin(), expansion.end());
                    body.erase(it++);
                } else
                    ++it;
            } else
                ++it;
        }

        delete tokenList;
    }

    return body;
}

}   // namespace glslparser

}   // namespace voreen
