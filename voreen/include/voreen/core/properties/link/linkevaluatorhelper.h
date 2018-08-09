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

#ifndef VRN_LINKEVALUATORHELPER_H
#define VRN_LINKEVALUATORHELPER_H

#include "voreen/core/properties/link/linkevaluatorbase.h"

#include <vector>
#include <string>

namespace voreen {

/**
 * Provides helper functions that internally retrieve the registered LinkEvaluators
 * from the VoreenApplication.
 */
class VRN_CORE_API LinkEvaluatorHelper {
public:
    /// Returns a new instance of the class corresponding to the given typeString.
    static LinkEvaluatorBase* createEvaluator(const std::string& typeString);

    ///Returns the complete list of functions registered with the factory.
    static std::vector<std::string> listFunctionNames();

    /**
     * Checks if the properties p1 and p2 are linkable.
     * A property is linkable if at least one registered LinkEvaluator can link these properties.
     *
     * \return true, if p1 and p2 are compatible, false otherwise
     */
    static bool arePropertiesLinkable(const Property* p1, const Property* p2, bool bidirectional = false);

    /**
     * Get all link evaluators that can link p1 to p2.
     *
     * @return Vector of compatible link evaluators, in the form of <Classname, name> pairs.
     */
    static std::vector<std::pair<std::string, std::string> > getCompatibleLinkEvaluators(const Property* p1, const Property* p2);
};

} // namespace

#endif // VRN_LINKEVALUATORHELPER_H
