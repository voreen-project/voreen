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

#include "utils.h"

namespace voreen {

TimeStepMapper::TimeStepMapper(const EnsembleDataset* dataset)
    : dataset_(dataset)
    , minDuration_(std::numeric_limits<float>::max())
    , maxTime_(std::numeric_limits<float>::lowest())
    , maxNumTimeSteps_(0)
{
    init();
}

TimeStepMapper::~TimeStepMapper() {
}

void TimeStepMapper::init() {
    for(EnsembleDataset::Run run : dataset_->getRuns()) {
        for(EnsembleDataset::TimeStep timeStep : run.timeSteps_) {
            minDuration_ = std::min(minDuration_, timeStep.duration_);
            maxTime_ = std::max(maxTime_, timeStep.time_);
        }
    }
    float globalDuration = 0.0f;
    while(globalDuration < maxTime_) {
        for(EnsembleDataset::Run run : dataset_->getRuns()) {
            for(EnsembleDataset::TimeStep timeStep : run.timeSteps_) {
                if(timeStep.time_ > globalDuration) {
                    std::pair<std::string, EnsembleDataset::TimeStep> timeStepPair;
                    timeStepPair = std::make_pair(run.name_, timeStep);
                    timeStepMap_[globalDuration].push_back(timeStepPair);
                    break;
                }
            }
        }
        globalDuration += minDuration_;
    }
    maxNumTimeSteps_ = timeStepMap_.size();
}

std::map<float, std::vector<std::pair<std::string, EnsembleDataset::TimeStep>>> TimeStepMapper::getTimeSteps() {
    return timeStepMap_;
}


void TimeStepMapper::iterate(const std::function <void(std::string, EnsembleDataset::TimeStep, float)>& callback) {
    for(std::pair<float, std::vector<std::pair<std::string, EnsembleDataset::TimeStep>>> ensembleTimeSteps : getTimeSteps()) {
        float time = ensembleTimeSteps.first;

        for(std::pair<std::string, EnsembleDataset::TimeStep> timeStepPair : ensembleTimeSteps.second) {
            std::string runName = timeStepPair.first;
            EnsembleDataset::TimeStep timeStep = timeStepPair.second;

            callback(runName, timeStep, time);
        }
     }
}

}


