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

#ifndef VRN_PRODUCTION_H
#define VRN_PRODUCTION_H

#include <string>
#include <vector>

namespace voreen {

namespace glslparser {

class GrammarSymbol;

/**
 * A -> B C D ...
 */
class Production {
public:
    Production(const int id, GrammarSymbol* const head, const std::vector<GrammarSymbol*>& body)
        : productionID_(id),
        head_(head),
        body_(body)
    {
    }

    bool operator==(const Production& rhs) const;
    bool operator!=(const Production& rhs) const;

    bool appendSymbol(GrammarSymbol* const symbol);
    bool removeSymbol(const int position);

    const std::vector<GrammarSymbol*>& getBody() const { return body_; }
    std::vector<int> getBodyIDs() const;
    const GrammarSymbol* getHead() const { return head_; }
    int getHeadID() const;
    int getProductionID() const { return productionID_; }
    bool isEmpty() const { return body_.empty(); }

    Production getLeftPart(const int position) const;
    Production getRightPart(const int position) const;

    std::string toString() const;

private:
    int productionID_;
    GrammarSymbol* head_;
    std::vector<GrammarSymbol*> body_;
};

}   // namespace glslparser

}   // namespace voreen

#endif
