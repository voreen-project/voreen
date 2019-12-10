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

#include "voreen/core/properties/link/linkevaluatorhelper.h"
#include "voreen/core/voreenapplication.h"

namespace voreen {

std::vector<std::string> LinkEvaluatorHelper::listFunctionNames() {
    if (!VoreenApplication::app()) {
        LERRORC("voreen.LinkEvaluatorHelper", "VoreenApplication not instantiated");
        return std::vector<std::string>();
    }

    std::vector<const LinkEvaluatorBase*> linkEvaluators = VoreenApplication::app()->getSerializableTypes<LinkEvaluatorBase>();
    std::vector<std::string> result;
    for (size_t i=0; i<linkEvaluators.size(); i++)
        result.push_back(linkEvaluators[i]->getClassName());
    sort(result.begin(), result.end());
    return result;
}

LinkEvaluatorBase* LinkEvaluatorHelper::createEvaluator(const std::string& typeString) {
    if (!VoreenApplication::app()) {
        LERRORC("voreen.LinkEvaluatorHelper", "VoreenApplication not instantiated");
        return 0;
    }

    Serializable* serializable = VoreenApplication::app()->createSerializableType(typeString);
    if (!serializable)
        return 0;
    else if (LinkEvaluatorBase* evaluator = dynamic_cast<LinkEvaluatorBase*>(serializable))
        return evaluator;
    else {
        LERRORC("voreen.LinkEvaluatorHelper", "Type '" << typeString << "' is not a LinkEvaluatorBase");
        delete serializable;
        return 0;
    }
}

bool LinkEvaluatorHelper::arePropertiesLinkable(const Property* p1, const Property* p2, bool bidirectional) {
    tgtAssert(p1, "null pointer passed");
    tgtAssert(p2, "null pointer passed");
    if (!VoreenApplication::app()) {
        LERRORC("voreen.LinkEvaluatorHelper", "VoreenApplication not instantiated");
        return false;
    }

    std::vector<const LinkEvaluatorBase*> linkEvaluators = VoreenApplication::app()->getSerializableTypes<LinkEvaluatorBase>();
    for (size_t i=0; i<linkEvaluators.size(); i++) {
        if (linkEvaluators[i]->arePropertiesLinkable(p1, p2)) {
            if (!bidirectional)
                return true;
            else
                return linkEvaluators[i]->arePropertiesLinkable(p2, p1);
        }
    }
    return false;
}

std::vector<std::pair<std::string, std::string> > LinkEvaluatorHelper::getCompatibleLinkEvaluators(const Property* p1, const Property* p2) {
    tgtAssert(p1, "null pointer passed");
    tgtAssert(p2, "null pointer passed");
    if (!VoreenApplication::app()) {
        LERRORC("voreen.LinkEvaluatorHelper", "VoreenApplication not instantiated");
        return std::vector<std::pair<std::string, std::string> >();
    }

    std::vector<std::pair<std::string, std::string> > compatible;
    std::vector<const LinkEvaluatorBase*> linkEvaluators = VoreenApplication::app()->getSerializableTypes<LinkEvaluatorBase>();
    for (size_t i=0; i<linkEvaluators.size(); i++) {
        if (linkEvaluators[i]->arePropertiesLinkable(p1, p2))
            compatible.push_back(std::pair<std::string, std::string>(linkEvaluators[i]->getClassName(), linkEvaluators[i]->getGuiName()));
    }
    return compatible;
}

} // namespace
