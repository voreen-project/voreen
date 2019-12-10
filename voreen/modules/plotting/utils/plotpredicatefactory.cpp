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

#include "plotpredicatefactory.h"

#include "../datastructures/plotpredicate.h"

namespace voreen {

PlotPredicateFactory* PlotPredicateFactory::instance_ = 0;

PlotPredicateFactory::PlotPredicateFactory() {
}

PlotPredicateFactory* PlotPredicateFactory::getInstance() {
    if (!instance_)
        instance_ = new PlotPredicateFactory();

    return instance_;
}

std::vector<std::string> PlotPredicateFactory::getAllTypeStrings(PredicateClasses predicateClasses) {
    std::vector<std::string> toReturn;

    toReturn.push_back("PlotPredicateLess");
    toReturn.push_back("PlotPredicateEqual");
    toReturn.push_back("PlotPredicateGreater");
    if (predicateClasses == ALL_PREDICATES ||
        predicateClasses == ALL_VALUE_PREDICATES_ONLY ||
        predicateClasses == ALL_TAG_PREDICATES_ONLY)
        toReturn.push_back("PlotPredicateNotBetween");
    toReturn.push_back("PlotPredicateBetween");
    if (predicateClasses == ALL_PREDICATES ||
        predicateClasses == ALL_VALUE_PREDICATES_ONLY ||
        predicateClasses == ALL_TAG_PREDICATES_ONLY)
        toReturn.push_back("PlotPredicateNotBetweenOrEqual");
    toReturn.push_back("PlotPredicateBetweenOrEqual");
    toReturn.push_back("PlotPredicateEmpty");
    toReturn.push_back("PlotPredicateNotEmpty");

    if (predicateClasses == ALL_PREDICATES ||
        predicateClasses == ALL_TAG_PREDICATES_ONLY ||
        predicateClasses == SIMPLE_TAG_PREDICATES_ONLY) {
        toReturn.push_back("PlotPredicateIsSubStr");
        toReturn.push_back("PlotPredicateAlphaNumeric");
        toReturn.push_back("PlotPredicateNotAlphaNumeric");
    }

    return toReturn;
}

std::string PlotPredicateFactory::getSerializableTypeString(const std::type_info& type) const {
    if (type == typeid(PlotPredicateLess))
        return "PlotPredicateLess";
    else if (type == typeid(PlotPredicateEqual))
        return "PlotPredicateEqual";
    else if (type == typeid(PlotPredicateGreater))
        return "PlotPredicateGreater";
    else if (type == typeid(PlotPredicateBetween))
        return "PlotPredicateBetween";
    else if (type == typeid(PlotPredicateNotBetween))
        return "PlotPredicateNotBetween";
    else if (type == typeid(PlotPredicateBetweenOrEqual))
        return "PlotPredicateBetweenOrEqual";
    else if (type == typeid(PlotPredicateNotBetweenOrEqual))
        return "PlotPredicateNotBetweenOrEqual";
    else if (type == typeid(PlotPredicateEmpty))
        return "PlotPredicateEmpty";
    else if (type == typeid(PlotPredicateNotEmpty))
        return "PlotPredicateNotEmpty";
    else if (type == typeid(PlotPredicateIsSubStr))
        return "PlotPredicateIsSubStr";
    else if (type == typeid(PlotPredicateAlphaNumeric))
        return "PlotPredicateAlphaNumeric";
    else if (type == typeid(PlotPredicateNotAlphaNumeric))
        return "PlotPredicateNotAlphaNumeric";
    else
        return "";
}

Serializable* PlotPredicateFactory::createSerializableType(const std::string& typeString) const {
    if (typeString == "PlotPredicateLess")
        return new PlotPredicateLess();
    else if (typeString == "PlotPredicateEqual")
        return new PlotPredicateEqual();
    else if (typeString == "PlotPredicateGreater")
        return new PlotPredicateGreater();
    else if (typeString == "PlotPredicateBetween")
        return new PlotPredicateBetween();
    else if (typeString == "PlotPredicateNotBetween")
        return new PlotPredicateNotBetween();
    else if (typeString == "PlotPredicateBetweenOrEqual")
        return new PlotPredicateBetweenOrEqual();
    else if (typeString == "PlotPredicateNotBetweenOrEqual")
        return new PlotPredicateNotBetweenOrEqual();
    else if (typeString == "PlotPredicateEmpty")
        return new PlotPredicateEmpty();
    else if (typeString == "PlotPredicateNotEmpty")
        return new PlotPredicateNotEmpty();
    else if (typeString == "PlotPredicateIsSubStr")
        return new PlotPredicateIsSubStr();
    else if (typeString == "PlotPredicateAlphaNumeric")
        return new PlotPredicateAlphaNumeric();
    else if (typeString == "PlotPredicateNotAlphaNumeric")
        return new PlotPredicateNotAlphaNumeric();
    else
        return 0;
}

} // namespace voreen
