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

#ifndef VRN_PROPERTYLINK_H
#define VRN_PROPERTYLINK_H

#include "voreen/core/io/serialization/serialization.h"
#include <vector>
#include <map>

namespace voreen {

class Property;
class LinkEvaluatorBase;

/**
 * This class establishes a link between two properties.
 *
 * When the value of the source property has changed the onChange method is called.
 * OnChange publishes the new value to the destination property by executing the
 * link evaluator assigned to the link.
 *
 * @see LinkEvaluatorBase
 */
class VRN_CORE_API PropertyLink : public Serializable {
    /** To get access to visitedProperties_ for loop detection. */
    friend class Property;
public:
    /**
     * Registers itself as an listener of the source property.
     * If no linkEvaluator is passed, an id-link is created.
     *
     * Source and destination properties must not be null.
     * Takes ownership of the LinkEvaluator.
     */
    PropertyLink(Property* src, Property* dest, LinkEvaluatorBase* linkEvaluator = 0);

    /**
     * Removes this link from the source property.
     * Deletes the LinkEvaluator;
     */
    virtual ~PropertyLink();

    /**
     * Callback-function when source property value has changed. It may modify
     * the new value and publish it to the destination property. Moreover, it
     * detects circular linking and stops publishing on detection.
     *
     * @throws VoreenException, if link execution failed. This usually
     *      because of incompatible data types.
     */
    void executeLink();

    /**
     * Forces the propagation with the current value thus testing functionality.
     *
     * @warn This function is deliberately not const, since it actually executes
     *      the link and might therefore affect the destination property.
     *
     * @return true, if the propagation succeeded. If false is returned,
     *      linking the respective properties is not possible.
     */
    bool testPropertyLink();

    /**
     * Assigns a new evaluator to enable transparent transition to a new
     * evaluator without recreating the link.
     * Takes ownership of the LinkEvaluator.
     */
    void setLinkEvaluator(LinkEvaluatorBase* evaluator);

    /**
     * Returns the link's evaluator.
     */
    LinkEvaluatorBase* getLinkEvaluator() const;

    /**
     * Returns the link's source property.
     */
    Property* getSourceProperty() const;

    /**
     * Returns the link's destination property.
     */
    Property* getDestinationProperty() const;

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

private:
    /// Default constructor needed for serialization.
    PropertyLink();
    friend class XmlDeserializer;

    Property* src_;                 ///< Source property of the link.
    Property* dest_;                ///< Destination property of the link.
    LinkEvaluatorBase* evaluator_;  ///< Evaluator responsible for actually executing the link.

    /// Used for cycle prevention during link execution. Is cleared in Property::executeLinks
    static std::vector<Property*> visitedProperties_;
    /// Category used in logging
    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_PROPERTYLINK_H
