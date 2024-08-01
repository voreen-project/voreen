/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "voreen/core/utils/GLSLparser/production.h"

#include "voreen/core/utils/GLSLparser/grammarsymbol.h"
#include <sstream>

namespace voreen {

namespace glslparser {

bool Production::operator!=(const Production& rhs) const {
    if (body_.size() != rhs.body_.size())
        return true;

    bool result = false;
    for (size_t i = 0; ((result == false) && (i < body_.size())); ++i)
        result = (*(body_[i]) != *(rhs.body_[i]));
    return result;
}

bool Production::operator==(const Production& rhs) const {
    if (body_.size() != rhs.body_.size())
        return false;

    bool result = true;
    for (size_t i = 0; ((result == true) && (i < body_.size())); ++i)
        result = (*(body_[i]) == *(rhs.body_[i]));
    return result;
}

std::vector<int> Production::getBodyIDs() const {
    std::vector<int> bodyIDs;
    for (size_t i = 0; i < body_.size(); ++i)
        bodyIDs.push_back(body_[i]->getSymbolID());
    return bodyIDs;
}

int Production::getHeadID() const {
    return head_->getSymbolID();
}

bool Production::appendSymbol(GrammarSymbol* const symbol) {
    if (symbol == 0)
        return false;

    body_.push_back(symbol);
    return true;
}

bool Production::removeSymbol(const int position) {
    if ( (position < 0) || (position >= static_cast<int>(body_.size())) )
        return false;

    body_.erase(body_.begin() + position);
    return true;
}

Production Production::getLeftPart(const int position) const {
    std::vector<GrammarSymbol*> split;

    if ((position < static_cast<int>(body_.size())) && (position < 0))
        split.insert(split.end(), body_.begin(), (body_.end() - position));

    return Production(productionID_, head_, split);
}

Production Production::getRightPart(const int position) const {
    std::vector<GrammarSymbol*> split;

    if ((position < static_cast<int>(body_.size())) && (position >= 0))
        split.insert(split.end(), (body_.begin() + position), body_.end());

    return Production(productionID_, head_, split);
}

std::string Production::toString() const {
    std::ostringstream oss;

    oss << head_->toString();
    for (size_t i = 0; i < body_.size(); ++i) {
        if (i > 0)
            oss << " ";
        else
            oss << " ::= ";

        oss << body_[i]->toString();
    }

    return oss.str();
}

}   // namespace glslparser

}   // namespace voreen
