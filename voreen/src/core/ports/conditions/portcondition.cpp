/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "voreen/core/ports/conditions/portcondition.h"

#include "tgt/logmanager.h"

namespace voreen {

PortCondition::PortCondition(const std::string& description)
    : description_(description)
{}

PortCondition::~PortCondition() {
}

std::string PortCondition::getDescription() const {
    return description_;
}

void PortCondition::setCheckedPort(const Port* checkedPort) {
    checkedPort_ = checkedPort;
}

//----------------------------------------------------------------------------------

PortConditionLogicalOr::PortConditionLogicalOr() :
    PortCondition("")
{}


PortConditionLogicalOr::~PortConditionLogicalOr() {
    for (size_t i=0; i<linkedConditions_.size(); i++)
        delete linkedConditions_.at(i);
    linkedConditions_.clear();
}

void PortConditionLogicalOr::addLinkedCondition(PortCondition* condition) {
    linkedConditions_.push_back(condition);
}

bool PortConditionLogicalOr::acceptsPortData() const {
    if (linkedConditions_.empty()) {
        LWARNINGC("voreen.PortConditionLogicalOr", "acceptsPortData(): No linked conditions");
        return false;
    }

    for (size_t i=0; i<linkedConditions_.size(); i++) {
        if (linkedConditions_.at(i)->acceptsPortData())
            return true;
    }
    return false;
}

std::string PortConditionLogicalOr::getDescription() const {
    std::string description = "PortConditionLogicalOr: None of the linked conditions is met: ";
    for (size_t i=0; i<linkedConditions_.size(); i++) {
        description += linkedConditions_.at(i)->getDescription();
        if (i < linkedConditions_.size()-1)
            description += ", ";
    }
    return description;
}

void PortConditionLogicalOr::setCheckedPort(const Port* checkedPort) {
    if (linkedConditions_.empty()) {
        LWARNINGC("voreen.PortConditionLogicalOr", "setCheckedPort(): No linked conditions");
        return;
    }

    for (size_t i=0; i<linkedConditions_.size(); i++)
        linkedConditions_.at(i)->setCheckedPort(checkedPort);
}

} // namespace
