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

#ifndef VRN_PORTCONDITION_H
#define VRN_PORTCONDITION_H

#include "voreen/core/voreencoreapi.h"

#include <string>
#include <vector>

namespace voreen {

class Port;

/**
 * Port conditions are used to restrict the allowed set of port data types.
 *
 * @see Port::addCondition
 */
class VRN_CORE_API PortCondition {

public:
    PortCondition(const std::string& description);

    virtual ~PortCondition();

    /**
     * Supposed to return whether the port data currently stored
     * in the checked port is accepted by this condition.
     */
    virtual bool acceptsPortData() const = 0;

    /**
     * Returns a textual description of the condition.
     */
    virtual std::string getDescription() const;

    /**
     * Assigns the port to be checked by this condition.
     *
     * @note Should usually not be called directly, since
     *  it is already called by Port::addCondition().
     */
    virtual void setCheckedPort(const Port* checkedPort);

protected:
    const Port* checkedPort_;
    std::string description_;
};

// ---------------------------------------------------------------------

class VRN_CORE_API PortConditionLogicalOr : public PortCondition {

public:
    PortConditionLogicalOr();

    virtual ~PortConditionLogicalOr();

    /**
     * Adds a condition to be linked by the or-operator.
     *
     * @note Takes ownership of the passed object.
     */
    void addLinkedCondition(PortCondition* condition);

    virtual bool acceptsPortData() const;

    virtual std::string getDescription() const;

    /// Propagates the assigned port to the linked conditions.
    virtual void setCheckedPort(const Port* checkedPort);

protected:
    std::vector<PortCondition*> linkedConditions_;
};


} // namespace voreen

#endif // VRN_PORTCONDITION_H
