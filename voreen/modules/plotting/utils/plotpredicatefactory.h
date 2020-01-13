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

#ifndef VRN_PLOTPREDICATEFACTORY_H
#define VRN_PLOTPREDICATEFACTORY_H

#include "voreen/core/io/serialization/serialization.h"

#include <string>
#include <vector>

namespace voreen {

/**
 * The @c PlotPredicateFactory class is responsible for PlotPredicate creation
 * during serialization process.
 *
 * @see SerializableFactory
 */
class VRN_CORE_API PlotPredicateFactory : public SerializableFactory {
public:
    enum PredicateClasses {
        ALL_PREDICATES,         ///< all predicates
        SIMPLE_VALUE_PREDICATES_ONLY,  ///< only predicates suitable for numeric values (this predicates are representable by one interval)
        ALL_VALUE_PREDICATES_ONLY,  ///< only predicates suitable for numeric values
        SIMPLE_TAG_PREDICATES_ONLY,     ///< only predicates suitable for tags (strings) (this predicates are representable by one interval)
        ALL_TAG_PREDICATES_ONLY     ///< only predicates suitable for tags (strings)
    };

    /**
     * Returns the singleton instance of the tranfer function factory.
     *
     * @returns the instance.
     */
    static PlotPredicateFactory* getInstance();

    /**
     * Returns a vector of type strings of all predicates suitable for \a predicateClasses.
     *
     * \param   predicateClasses    class of predicates to filter
     **/
    static std::vector<std::string> getAllTypeStrings(PredicateClasses predicateClasses);

    /**
     * @see SerializableFactory::getTypeString
     */
    virtual std::string getSerializableTypeString(const std::type_info& type) const;

    /**
     * @see SerializableFactory::createType
     */
    virtual Serializable* createSerializableType(const std::string& typeString) const;

private:
    /**
     * Default constructor.
     */
    PlotPredicateFactory();

    /**
     * Singleton instance of the transfer function factory.
     */
    static PlotPredicateFactory* instance_;
};

} // namespace voreen

#endif // VRN_TRANSFUNCFACTORY_H
