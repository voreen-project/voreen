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

#include "voreen/core/utils/GLSLparser/terminalsmap.h"

namespace voreen {

namespace glslparser {

TerminalsMap::TerminalsMap()
   : lexMap_()
{
}

TerminalsMap::~TerminalsMap() {
}

int TerminalsMap::findSymbolID(const std::string& lexeme) const {
    LexemeMap::const_iterator it = lexMap_.find(lexeme);
    if (it != lexMap_.end())
        return it->second;

    return -1;
}

std::string TerminalsMap::findLexeme(const int symbolID) const {
    IDMap::const_iterator it = idMap_.find(symbolID);
    if (it != idMap_.end())
        return it->second;

    return "";
}

// protected methods
//

bool TerminalsMap::insertSymbol(const std::string& lexeme, const int symbolID, const bool isKeyword)
{
    std::pair<LexemeMap::iterator, bool> res = lexMap_.insert(std::make_pair(lexeme, symbolID));
    if (isKeyword)
        keywords_.insert(std::make_pair(lexeme, Word(symbolID, lexeme)));

    std::pair<IDMap::iterator, bool> res2 = idMap_.insert(std::make_pair(symbolID, lexeme));
    return (res.second && res2.second);
}

}   // namespace glslparser

}   // namespace voreen
