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

#include "portconditionensemble.h"

#include "voreen/core/utils/stringutils.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/utils/stringutils.h"

namespace voreen {

PortConditionEnsemble::PortConditionEnsemble(const std::string& message)
    : PortCondition(message)
{}

PortConditionEnsemble::~PortConditionEnsemble()
{}

void PortConditionEnsemble::setCheckedPort(const Port* checkedPort) {
    if (!dynamic_cast<const EnsembleDatasetPort*>(checkedPort)) {
        LERRORC("voreen.PortConditionEnsemble", "Assigned port is not an ensemble port");
    }
    else {
        ensemblePort_ = static_cast<const EnsembleDatasetPort*>(checkedPort);
    }
}

PortConditionEnsembleFieldName::PortConditionEnsembleFieldName(const std::string& fieldName)
    : PortConditionEnsemble("Ensemble with field " + fieldName + " expected")
    , fieldNames_(1, fieldName)
{}

PortConditionEnsembleFieldName::PortConditionEnsembleFieldName(const std::vector<std::string>& fieldNames)
    : PortConditionEnsemble("Ensemble with fields: " + strJoin(fieldNames, ", ") + " expected")
    , fieldNames_(fieldNames)
{}

PortConditionEnsembleFieldName::~PortConditionEnsembleFieldName()
{}

bool PortConditionEnsembleFieldName::acceptsPortData() const  {
    if (!ensemblePort_ || !ensemblePort_->hasData())
        return false;

    if(ensemblePort_->getData()->getCommonFieldNames().empty())
        return false;

    for(const std::string& field1 : ensemblePort_->getData()->getCommonFieldNames()) {
        bool found = false;
        for(const std::string& field2 : fieldNames_) {
            if(field1 == field2) {
                found = true;
                break;
            }
        }
        if(!found)
            return false;
    }

    return true;
}

PortConditionEnsembleFieldCount::PortConditionEnsembleFieldCount(size_t fieldCount)
    : PortConditionEnsemble("Ensemble with " + std::to_string(fieldCount) + " field(s) expected")
    , fieldCount_(fieldCount)
{}

PortConditionEnsembleFieldCount::~PortConditionEnsembleFieldCount()
{}

bool PortConditionEnsembleFieldCount::acceptsPortData() const  {
    if (!ensemblePort_ || !ensemblePort_->hasData())
        return false;

    if(ensemblePort_->getData()->getCommonFieldNames().size() != fieldCount_)
        return false;

    return true;
}

PortConditionEnsembleSingleTimeStep::PortConditionEnsembleSingleTimeStep()
    : PortConditionEnsemble("Only one time step expected!")
{}

PortConditionEnsembleSingleTimeStep::~PortConditionEnsembleSingleTimeStep()
{}

bool PortConditionEnsembleSingleTimeStep::acceptsPortData() const  {
    if (!ensemblePort_ || !ensemblePort_->hasData())
        return false;

    if(ensemblePort_->getData()->getMinNumTimeSteps() != 1)
        return false;

    return true;
}

} // namespace
