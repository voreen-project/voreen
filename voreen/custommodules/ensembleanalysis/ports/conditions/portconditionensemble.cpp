/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

PortConditionEnsembleChannelName::PortConditionEnsembleChannelName(const std::string& channelName)
    : PortConditionEnsemble("Ensemble with channel " + channelName + " expected")
    , channelNames_(1, channelName)
{}

PortConditionEnsembleChannelName::PortConditionEnsembleChannelName(const std::vector<std::string>& channelNames)
    : PortConditionEnsemble("Ensemble with channels: " + strJoin(channelNames, ", ") + " expected")
    , channelNames_(channelNames)
{}

PortConditionEnsembleChannelName::~PortConditionEnsembleChannelName()
{}

bool PortConditionEnsembleChannelName::acceptsPortData() const  {
    if (!ensemblePort_ || !ensemblePort_->hasData())
        return false;

    if(ensemblePort_->getData()->getCommonChannels().empty())
        return false;

    for(const std::string& channel1 : ensemblePort_->getData()->getCommonChannels()) {
        bool found = false;
        for(const std::string& channel2 : channelNames_) {
            if(channel1 == channel2) {
                found = true;
                break;
            }
        }
        if(!found)
            return false;
    }

    return true;
}

PortConditionEnsembleChannelCount::PortConditionEnsembleChannelCount(size_t channelCount)
    : PortConditionEnsemble("Ensemble with " + std::to_string(channelCount) + " channel(s) expected")
    , channelCount_(channelCount)
{}

PortConditionEnsembleChannelCount::~PortConditionEnsembleChannelCount()
{}

bool PortConditionEnsembleChannelCount::acceptsPortData() const  {
    if (!ensemblePort_ || !ensemblePort_->hasData())
        return false;

    if(ensemblePort_->getData()->getCommonChannels().size() != channelCount_)
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
