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
    // Adding runs sets attributes accordingly.
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
        std::vector<std::string> fields;
        for (const auto& field : run.timeSteps_[t].fieldNames_) {

            const std::string& fieldName = field.first;
            fields.push_back(fieldName);

            // Retrieve volume.
            const VolumeBase* volume = field.second;

            // Gather parameters (take first time step as representative).
            if(t==0) {
                for (const std::string& key : volume->getMetaDataKeys()) {
                    if (key.find("Parameter") != std::string::npos) {
                        allParameters_.insert(key);
                    }
                }
            }

            // Bounds are stored in physical space, so don't transform to world space.
            tgt::Bounds bounds = volume->getBoundingBox(false).getBoundingBox();
            VolumeMinMax* vmm = volume->getDerivedData<VolumeMinMax>();
            tgt::vec2 minMax(std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest());
            for(size_t c = 0; c < vmm->getNumChannels(); c++) {
                minMax.x = std::min(minMax.x, vmm->getMin(c));
                minMax.y = std::max(minMax.y, vmm->getMax(c));
            }

            /*
            // In further applications it might be useful to force derived data calculations
            // to improve responsivity after loading the data.
            if(volume->getNumChannels() > 1) {
                volume->getDerivedData<VolumeMinMaxMagnitude>();
            }
            */

            bool firstFieldElement = fieldMetaData_.find(fieldName) == fieldMetaData_.end();
            FieldMetaData& fieldMetaData = fieldMetaData_[fieldName];
            if(firstFieldElement) {
                fieldMetaData.valueRange_ = minMax;
                fieldMetaData.numChannels_ = volume->getNumChannels();
            }
            else {
                fieldMetaData.valueRange_.x = std::min(fieldMetaData.valueRange_.x, minMax.x);
                fieldMetaData.valueRange_.y = std::max(fieldMetaData.valueRange_.y, minMax.y);
                if(fieldMetaData.numChannels_ != volume->getNumChannels()) {
                    LERRORC("voreen.EnsembleDataSet", "Number of channels inside channel differs, taking min.");
                    fieldMetaData.numChannels_ = std::min(fieldMetaData.numChannels_, volume->getNumChannels());
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

        // Calculate common fields.
        if (!commonFieldNames_.empty()) {
            std::vector<std::string> intersection;
            std::set_intersection(
                commonFieldNames_.begin(),
                commonFieldNames_.end(),
                fields.begin(),
                fields.end(),
                std::back_inserter(intersection)
            );

            if (commonFieldNames_.size() != intersection.size() && !runs_.empty()) {
                LWARNINGC("voreen.EnsembeDataSet", "Time Step " << t << " of Run " << run.name_ << " has less fields than the previously added Run " << runs_.back().name_);
            }

            commonFieldNames_ = intersection;
        }
        else if (runs_.empty()) {
            commonFieldNames_ = fields;
        }

        // Update all fields.
        std::vector<std::string> fieldUnion;
        std::set_union(uniqueFieldNames_.begin(), uniqueFieldNames_.end(), fields.begin(), fields.end(), std::back_inserter(fieldUnion));
        uniqueFieldNames_ = fieldUnion;

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
    if(roi.isDefined()) {
        notifyPendingDataInvalidation();
        roi_ = roi;
    }
    else {
        LWARNINGC("voreen.EnsembleDataSet", "Roi must overlap with common domain bounds, ignoring.");
    }
}

const tgt::vec2& EnsembleDataset::getValueRange(const std::string& field) const {
    tgtAssert(fieldMetaData_.find(field) != fieldMetaData_.end(), "Field not available");
    return fieldMetaData_.at(field).valueRange_;
}

size_t EnsembleDataset::getNumChannels(const std::string& field) const {
    tgtAssert(fieldMetaData_.find(field) != fieldMetaData_.end(), "Field not available");
    return fieldMetaData_.at(field).numChannels_;
}

const std::vector<std::string>& EnsembleDataset::getUniqueFieldNames() const {
    return uniqueFieldNames_;
}

const std::vector<std::string>& EnsembleDataset::getCommonFieldNames() const {
    return commonFieldNames_;
}

std::vector<const VolumeBase*> EnsembleDataset::getVolumes() const {
    std::vector<const VolumeBase*> result;
    for(const Run& run : runs_) {
        for(const TimeStep& timeStep : run.timeSteps_) {
            for(const auto& field : timeStep.fieldNames_) {
                result.push_back(field.second);
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

std::string EnsembleDataset::toHTML() const {

    std::stringstream stream;

    stream << "<html><head>"
              "<meta content=\"text/html;charset=utf-8\" http-equiv=\"Content-Type\">\n"
              "<meta content=\"utf-8\" http-equiv=\"encoding\">\n"
              "<link src=\"https://cdn.datatables.net/1.10.20/css/jquery.dataTables.min.css\" rel=\"stylesheet\">\n"
              "<script src=\"https://code.jquery.com/jquery-3.4.1.min.js\"></script>\n"
              "<script src=\"https://cdn.datatables.net/1.10.20/js/jquery.dataTables.min.js\"></script>\n"
              "<style>table,th,td {border: 1px solid black;}</style></head>"
              "<body><table id=\"ensemble\"><thead>";
    // Parameter names
    stream << "  <tr>\n";
    // Run Name and Color are mandatory.
    stream << "    <th>Name</th>\n";
    stream << "    <th>Color</th>\n";
    stream << "    <th>Num. Time Steps</th>\n";
    stream << "    <th>Start Time</th>\n";
    stream << "    <th>End Time</th>\n";
    for(const std::string& parameter : allParameters_) {
        stream << "    <th>" << parameter << "</th>\n";

    }
    stream << "  </tr></thead><tbody>\n";

    // Runs and their parameters.
    for(size_t i=0; i<runs_.size(); i++) {
        stream << "  <tr>\n";
        stream << "    <th>" << runs_[i].name_ << "</th>\n";
        tgt::ivec3 color(runs_[i].color_ * 255.0f);
        stream << "    <th style=\"background-color: rgb(" << color.r << ", " << color.g << ", " << color.b << ")\"></th>\n";
        stream << "    <th>" << runs_[i].timeSteps_.size() << "</th>\n";
        stream << "    <th>" << runs_[i].timeSteps_.front().time_ << "</th>\n";
        stream << "    <th>" << runs_[i].timeSteps_.back().time_ << "</th>\n";

        for(const std::string& parameter : allParameters_) {
            // TODO: assumes that all fields contain the same parameters.
            const VolumeBase* reference = runs_[i].timeSteps_.front().fieldNames_.begin()->second;
            stream << "    <th>";
            if(reference->hasMetaData(parameter)) {
                stream << reference->getMetaData(parameter)->toString();
            }
            stream << "</th>\n";
        }
        stream << "  </tr>\n";
    }

    stream <<"</tbody></table>"
             "<script>$(document).ready( function () {$('#ensemble').DataTable({paging: false});} );</script>"
             "</body></html>";
    return stream.str();
}

}   // namespace
