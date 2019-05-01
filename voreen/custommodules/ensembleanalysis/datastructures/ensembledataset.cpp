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

#include "ensembledataset.h"

#include "tgt/assert.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

namespace voreen {

EnsembleDataset::EnsembleDataset()
    : minNumTimeSteps_(std::numeric_limits<size_t>::max())
    , maxNumTimeSteps_(0)
    , totalNumTimeSteps_(0)
    , maxTimeStepDuration_(0.0f)
    , minTimeStepDuration_(std::numeric_limits<float>::max())
    , startTime_(std::numeric_limits<float>::max())
    , endTime_(0.0f)
    , commonTimeInterval_(endTime_, startTime_)
    , bounds_()
    , commonBounds_()
    , roi_()
{
}

EnsembleDataset::EnsembleDataset(const EnsembleDataset& origin)
    : EnsembleDataset()
{
    // Adding runs set's attributes accordingly.
    for(const Run& run : origin.runs_)
        addRun(run);

    // Set Roi first since it might has been modified.
    setRoi(origin.getRoi());
}

void EnsembleDataset::addRun(const Run& run) {

    // Skip empty runs.
    if (run.timeSteps_.empty())
        return;

    // Notify Observers
    notifyPendingDataInvalidation();

    minNumTimeSteps_ = std::min(run.timeSteps_.size(), minNumTimeSteps_);
    maxNumTimeSteps_ = std::max(run.timeSteps_.size(), maxNumTimeSteps_);
    totalNumTimeSteps_ += run.timeSteps_.size();
    startTime_ = std::min(startTime_, run.timeSteps_.front().time_);
    endTime_   = std::max(endTime_,   run.timeSteps_.back().time_);

    RunMetaData metaData;

    for (size_t t = 0; t < run.timeSteps_.size(); t++) {
        std::vector<std::string> channels;
        for (const auto& channel : run.timeSteps_[t].channels_) {

            const std::string& channelName = channel.first;
            channels.push_back(channelName);

            const VolumeBase* volume = channel.second;
            // Bounds are stored in physical space, so don't transform to world space.
            tgt::Bounds bounds = volume->getBoundingBox(false).getBoundingBox();
            tgt::vec2 minMax;
            if(volume->getNumChannels() == 1) {
                VolumeMinMax* derivedData = volume->getDerivedData<VolumeMinMax>();
                minMax = tgt::vec2(derivedData->getMin(), derivedData->getMax());
            }
            else { // bigger than 1
                VolumeMinMaxMagnitude* derivedData = volume->getDerivedData<VolumeMinMaxMagnitude>();
                minMax = tgt::vec2(derivedData->getMinMagnitude(), derivedData->getMaxMagnitude());
            }

            bool firstChannelElement = channelMetaData_.find(channelName) == channelMetaData_.end();
            ChannelMetaData& channelMetaData = channelMetaData_[channelName];
            if(firstChannelElement) {
                channelMetaData.valueRange_ = minMax;
                channelMetaData.numChannels_ = volume->getNumChannels();
            }
            else {
                channelMetaData.valueRange_.x = std::min(channelMetaData.valueRange_.x, minMax.x);
                channelMetaData.valueRange_.y = std::max(channelMetaData.valueRange_.y, minMax.y);
                if(channelMetaData.numChannels_ != volume->getNumChannels()) {
                    LERRORC("voreen.EnsembleDataSet", "Number of channels inside channel differs, taking min.");
                    channelMetaData.numChannels_ = std::min(channelMetaData.numChannels_, volume->getNumChannels());
                }
            }

            metaData.timeStepDurationStats_.addSample(run.timeSteps_[t].duration_);

            if (!bounds_.isDefined()) {
                if(!commonBounds_.isDefined()) {
                    commonBounds_.addVolume(bounds);
                }
            }
            else if(commonBounds_.isDefined()) {
                commonBounds_.intersectVolume(bounds);
                if(!commonBounds_.isDefined()) {
                    LWARNINGC("voreen.EnsembeDataSet", "There is no overlap between the bounds of Run " << run.name_ << " and the previously defined bounds");
                }
            }
            bounds_.addVolume(bounds);
        }

        // Calculate common channels.
        if (!commonChannels_.empty()) {
            std::vector<std::string> intersection;
            std::set_intersection(
                commonChannels_.begin(),
                commonChannels_.end(),
                channels.begin(),
                channels.end(),
                std::back_inserter(intersection)
            );

            if (commonChannels_.size() != intersection.size() && !runs_.empty()) {
                LWARNINGC("voreen.EnsembeDataSet", "Time Step " << t << " of Run " << run.name_ << " has less channels than the previously added Run " << runs_.back().name_);
            }

            commonChannels_ = intersection;
        }
        else if (runs_.empty()) {
            commonChannels_ = channels;
        }

        // Calculate times and durations.
        if (t < run.timeSteps_.size() - 1) {
            maxTimeStepDuration_ = std::max(maxTimeStepDuration_, run.timeSteps_[t].duration_);
            minTimeStepDuration_ = std::min(minTimeStepDuration_, run.timeSteps_[t].duration_);
        }
    }

    // Reset roi to current common bounds.
    roi_ = commonBounds_;

    commonTimeInterval_.x = std::max(commonTimeInterval_.x, run.timeSteps_.front().time_);
    commonTimeInterval_.y = std::min(commonTimeInterval_.y, run.timeSteps_.back().time_+run.timeSteps_.back().duration_);

    if(commonTimeInterval_.x > commonTimeInterval_.y) {
        LWARNINGC("voreen.EnsembleDataSet", "The time interval of the currently added Run " << run.name_ << " does not overlap with the previous interval");
        commonTimeInterval_ = tgt::vec2::zero;
    }

    runMetaData_.push_back(metaData);
    runs_.push_back(run);
}

const std::vector<EnsembleDataset::Run>& EnsembleDataset::getRuns() const {
    return runs_;
}

size_t EnsembleDataset::getMinNumTimeSteps() const {
    return minNumTimeSteps_;
}

size_t EnsembleDataset::getMaxNumTimeSteps() const {
    return maxNumTimeSteps_;
}

size_t EnsembleDataset::getTotalNumTimeSteps() const {
    return totalNumTimeSteps_;
}

const Statistics& EnsembleDataset::getTimeStepDurationStats(size_t runIdx) const {
    tgtAssert(runIdx < runs_.size(), "Run not available");
    return runMetaData_[runIdx].timeStepDurationStats_;
}

const tgt::vec3& EnsembleDataset::getColor(size_t runIdx) const {
    tgtAssert(runIdx < runs_.size(), "Run not available");
    return runs_[runIdx].color_;
}

float EnsembleDataset::getMinTimeStepDuration() const {
    return minTimeStepDuration_;
}

float EnsembleDataset::getMaxTimeStepDuration() const {
    return maxTimeStepDuration_;
}

float EnsembleDataset::getStartTime() const {
    return startTime_;
}

float EnsembleDataset::getEndTime() const {
    return endTime_;
}

float EnsembleDataset::getMaxTotalDuration() const {
    return getEndTime() - getStartTime();
}

const tgt::vec2& EnsembleDataset::getCommonTimeInterval() const {
    return commonTimeInterval_;
}

const tgt::Bounds& EnsembleDataset::getBounds() const {
    return bounds_;
}

const tgt::Bounds& EnsembleDataset::getCommonBounds() const {
    return commonBounds_;
}

const tgt::Bounds& EnsembleDataset::getRoi() const {
    return roi_;
}

void EnsembleDataset::setRoi(tgt::Bounds roi) {
    roi.intersectVolume(commonBounds_);
    if(!roi.isDefined()) {
        notifyPendingDataInvalidation();
        roi_ = roi;
    }
    else {
        LWARNINGC("voreen.EnsembleDataSet", "Roi must overlap with common domain bounds");
    }
}

const tgt::vec2& EnsembleDataset::getValueRange(const std::string& channel) const {
    tgtAssert(channelMetaData_.find(channel) != channelMetaData_.end(), "Channel not available");
    return channelMetaData_.at(channel).valueRange_;
}

size_t EnsembleDataset::getNumChannels(const std::string& channel) const {
    tgtAssert(channelMetaData_.find(channel) != channelMetaData_.end(), "Channel not available");
    return channelMetaData_.at(channel).numChannels_;
}

const std::vector<std::string>& EnsembleDataset::getCommonChannels() const {
    return commonChannels_;
}

std::vector<const VolumeBase*> EnsembleDataset::getVolumes() const {
    std::vector<const VolumeBase*> result;
    for(const Run& run : runs_) {
        for(const TimeStep& timeStep : run.timeSteps_) {
            for(const auto& channel : timeStep.channels_) {
                result.push_back(channel.second);
            }
        }
    }
    return result;
}

size_t EnsembleDataset::pickTimeStep(size_t runIdx, float time) const {
    tgtAssert(runIdx < runs_.size(), "run index too big");

    if (runs_[runIdx].timeSteps_.empty())
        return -1;

    size_t t = 0;
    while (t < runs_[runIdx].timeSteps_.size()-1 && runs_[runIdx].timeSteps_[t].time_ < time) t++;
    return t;

}

}   // namespace
