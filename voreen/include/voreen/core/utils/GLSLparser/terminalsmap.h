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

#ifndef VRN_TERMINALSMAP_H
#define VRN_TERMINALSMAP_H

//#define GETSYMBOLIC(symbol) std::make_pair<int, std::string>(symbol, #symbol)

#include "voreen/core/utils/GLSLparser/token.h"

#include <map>
#include <string>

namespace voreen {

namespace glslparser {

/**
 * Class to setup a grammar and finally the tables for a parser for
 * that grammar.
 *
 * It is intended to be derived and hold a public enum for all token
 * IDs. Those tokens are considered as terminals in the grammar and
 * should therefore have the same IDs as the terminal symbols.
 */
class TerminalsMap {
public:
    TerminalsMap();
    virtual ~TerminalsMap() = 0;    // ensures that the class is abstract

    /**
     * Returns the ID of the symbol corresponding to the given lexeme if
     * the lexeme is known or -1 otherwise.
     */
    int findSymbolID(const std::string& lexeme) const;

    std::string findLexeme(const int symbolID) const;

    const std::map<std::string, Word>& getKeywords() const { return keywords_; }

    /**
     * Returns the map mapping the lexemes to their IDs.
     */
    const std::map<std::string, int>& getLexemeMap() const { return lexMap_; }

    int getNumTerminals() const { return static_cast<int>(lexMap_.size()); }

protected:
    /**
     * Adds the given lexeme with the given ID into the map if the lexeme is
     * is not already contained.
     *
     * NOTE: this method does not check whether the given ID is already associated
     * with any other lexeme within the map!
     */
    bool insertSymbol(const std::string& lexeme, const int symbolID, const bool isKeyword = false);

private:
    typedef std::map<std::string, int> LexemeMap;
    LexemeMap lexMap_;

    typedef std::map<int, std::string> IDMap;
    IDMap idMap_;

    typedef std::map<std::string, Word> KeywordMap;
    KeywordMap keywords_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
