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
    , dimensions_(tgt::svec3::zero)
    , spacing_(tgt::vec3::zero)
    , roi_()
{
}

EnsembleDataset::EnsembleDataset(const EnsembleDataset& origin)
    : EnsembleDataset(&origin)
{
}

EnsembleDataset::EnsembleDataset(const EnsembleDataset* const origin)
    : EnsembleDataset()
{
    tgtAssert(origin, "Origin was null");

    // Set Roi first since it might has been modified.
    setRoi(origin->getRoi());

    // Adding runs set's attributes accordingly.
    for(const Run& run : origin->runs_)
        addRun(run);
}

EnsembleDataset::~EnsembleDataset()
{
}
    //----------------
    //  Access
    //----------------

void EnsembleDataset::addRun(const Run& run) {

    // Skip empty runs.
    if (run.timeSteps_.empty())
        return;

    // Notify Observers
    notifyPendingDataInvalidation();

    minNumTimeSteps_ = std::min(run.timeSteps_.size(), minNumTimeSteps_);
    maxNumTimeSteps_ = std::max(run.timeSteps_.size(), maxNumTimeSteps_);
    totalNumTimeSteps_ += run.timeSteps_.size();

    for (size_t t = 0; t < run.timeSteps_.size(); t++) {
        std::vector<std::string> channels;
        for (const auto& channel : run.timeSteps_[t].channels_) {

            const std::string& channelName = channel.first;
            channels.push_back(channelName);

            tgtAssert(channel.second->hasDerivedData<VolumeMinMax>(), "Derived data min max not available");
            VolumeMinMax* minMax = channel.second->getDerivedData<VolumeMinMax>();

            if (valueRange_.find(channelName) != valueRange_.end())
                valueRange_[channelName] = tgt::vec2(std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest());

            valueRange_[channelName].x = std::min(valueRange_[channelName].x, minMax->getMin());
            valueRange_[channelName].y = std::max(valueRange_[channelName].y, minMax->getMax());

            if (runs_.empty()) {
                dimensions_ = channel.second->getDimensions();
                spacing_ = channel.second->getSpacing();
            }

            /*
            //TODO: handle shrinking and expanding dimensions at the same time.
            if(dimensions_ != channel.second->getDimensions())
                //LWARNINGC("voreen.ensembledataset", "Dimensions of run '" << run.name_ << "' do not match dimension of previously set runs");
            }
            */
            tgtAssert(dimensions_ == channel.second->getDimensions(), "Dimensions do not match!");

            // Assuming dimensions don't change, we can init the roi here.
            if (!roi_.isDefined()) {
                roi_.addPoint(tgt::ivec3::zero);
                roi_.addPoint(tgt::ivec3(dimensions_) - tgt::ivec3::one);
            }

            const tgt::vec3& spacing = channel.second->getSpacing();
            if (spacing_ != spacing) {
                //LWARNINGC("voreen.ensembledataset", "Spacing does not match for run '" << run.name_ << "'");
                spacing_ = tgt::min(spacing_, spacing);
            }
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

            if (commonChannels_.size() != intersection.size() && !runs_.empty())
                LWARNINGC("voreen.EnsembeDataSet", "TimeStep " << t << " of Run " << run.name_ << " has less channels than the previously added Run " << runs_.back().name_);

            commonChannels_ = intersection;
        }
        else if (runs_.empty())
            commonChannels_ = channels;

        // Calculate times and durations.
        if (t < run.timeSteps_.size() - 1) {
            maxTimeStepDuration_ = std::max(maxTimeStepDuration_, run.timeSteps_[t].duration_);
            minTimeStepDuration_ = std::min(minTimeStepDuration_, run.timeSteps_[t].duration_);
        }

        startTime_ = std::min(startTime_, run.timeSteps_[t].time_);
        endTime_   = std::max(endTime_,   run.timeSteps_[t].time_+run.timeSteps_[t].duration_);
        commonTimeInterval_.x = std::max(commonTimeInterval_.x, run.timeSteps_[t].time_);
        commonTimeInterval_.y = std::min(commonTimeInterval_.y, run.timeSteps_[t].time_+run.timeSteps_[t].duration_);
    }

    if(commonTimeInterval_.x > commonTimeInterval_.y) {
        LWARNINGC("voreen.EnsembleDataSet", "The time interval of the currently added Run " << run.name_ << " does not overlap with the prior interval");
        commonTimeInterval_ = tgt::vec2::zero;
    }

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

const tgt::svec3& EnsembleDataset::getDimensions() const {
    return dimensions_;
}

const tgt::vec3& EnsembleDataset::getSpacing() const {
    return spacing_;
}

const tgt::IntBounds& EnsembleDataset::getRoi() const {
    return roi_;
}

void EnsembleDataset::setRoi(const tgt::IntBounds& roi) {
    notifyPendingDataInvalidation();
    roi_ = roi;
}

const tgt::vec2& EnsembleDataset::getValueRange(const std::string& channel) const {
    tgtAssert(valueRange_.find(channel) != valueRange_.end(), "Channel not available");
    return valueRange_.at(channel);
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

float EnsembleDataset::pickSample(const VolumeRAM_Float* volume, const tgt::vec3& spacing, tgt::vec3 sample, VolumeRAM::Filter filter) const {
    tgtAssert(volume, "volume null");
    tgtAssert(volume->getDimensions() == getDimensions(), "Dimensions do not match");

    tgt::vec3 offset = (spacing - getSpacing()) * tgt::vec3(getDimensions()) / (2.0f * spacing);
    sample = sample * (getSpacing() / spacing) + offset;

    // switch between filtering options
    switch (filter) {
    case VolumeRAM::NEAREST: {
        tgt::svec3 index = tgt::clamp(tgt::svec3(tgt::iround(sample)), tgt::svec3::zero, volume->getDimensions() - tgt::svec3::one);
        return volume->voxel(index); // round and do the lookup
    }
    case VolumeRAM::LINEAR: {
        // clamp to volume dimensions
        sample = tgt::clamp(sample, tgt::vec3::zero, tgt::vec3(volume->getDimensions() - tgt::svec3::one));

        // get decimal part and lower / upper voxel
        tgt::vec3 p = sample - tgt::floor(sample);
        tgt::ivec3 llb(sample);
        tgt::ivec3 urf(tgt::ceil(sample));

        // clamp again for safety so the lookups do not exceed the dimensions
        llb = tgt::max(llb, tgt::ivec3::zero);
        urf = tgt::min(urf, tgt::ivec3(volume->getDimensions()) - 1);

        //interpolate linearly
        return
                (volume->voxel(llb.x, llb.y, llb.z)) * (1.f-p.x)*(1.f-p.y)*(1.f-p.z)  // llB
              + (volume->voxel(urf.x, llb.y, llb.z)) * (    p.x)*(1.f-p.y)*(1.f-p.z)  // lrB
              + (volume->voxel(urf.x, urf.y, llb.z)) * (    p.x)*(    p.y)*(1.f-p.z)  // urB
              + (volume->voxel(llb.x, urf.y, llb.z)) * (1.f-p.x)*(    p.y)*(1.f-p.z)  // ulB
              + (volume->voxel(llb.x, llb.y, urf.z)) * (1.f-p.x)*(1.f-p.y)*(    p.z)  // llF
              + (volume->voxel(urf.x, llb.y, urf.z)) * (    p.x)*(1.f-p.y)*(    p.z)  // lrF
              + (volume->voxel(urf.x, urf.y, urf.z)) * (    p.x)*(    p.y)*(    p.z)  // urF
              + (volume->voxel(llb.x, urf.y, urf.z)) * (1.f-p.x)*(    p.y)*(    p.z); // ulF
    }
    case VolumeRAM::CUBIC: {
        // clamp to volume dimensions
        sample = tgt::clamp(sample, tgt::vec3::zero, tgt::vec3(volume->getDimensions() - tgt::svec3::one));
        return volume->getVoxelNormalizedCubic(sample);
    }
    default:
        return 0.0f;
    }
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
