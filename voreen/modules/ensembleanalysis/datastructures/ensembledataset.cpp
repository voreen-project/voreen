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

#include "ensembledataset.h"

#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"

namespace voreen {

//////////// Time Step


EnsembleDataset::TimeStep::VolumeCache::VolumeCache()
{
}

EnsembleDataset::TimeStep::VolumeCache::VolumeCache(const std::map<std::string, const VolumeBase*>& volumeData) {
    for(const auto& vol : volumeData) {
        const VolumeBase* volume = vol.second;
        volume->Observable<VolumeObserver>::addObserver(static_cast<const VolumeObserver*>(this));
        cacheEntries_.insert(std::make_pair(volume->getOrigin().getURL(), VolumeCacheEntry{volume, false}));
    }
}

EnsembleDataset::TimeStep::VolumeCache::~VolumeCache() {
    for(const auto& entry : cacheEntries_) {
        if(entry.second.owned_) {
            delete entry.second.volume_;
        }
    }
}

void EnsembleDataset::TimeStep::VolumeCache::volumeDelete(const VolumeBase* source) {
    // If the volume gets deleted by the owner, we remove it from the cache.
    // Otherwise, we would have a dangling pointer.
    for(auto iter = cacheEntries_.begin(); iter != cacheEntries_.end(); iter++) {
        if(iter->second.volume_ == source) {
            cacheEntries_.erase(iter);
            return;
        }
    }
}

void EnsembleDataset::TimeStep::VolumeCache::volumeChange(const VolumeBase* source) {
}

const VolumeBase* EnsembleDataset::TimeStep::VolumeCache::requestVolume(const VolumeURL& url) {
    boost::lock_guard<boost::mutex> lockGuard(volumeDataMutex_);

    // First query volume data.
    std::string urlString = url.getURL();
    auto volIter = cacheEntries_.find(urlString);
    if(volIter != cacheEntries_.end()) {
        return volIter->second.volume_;
    }

    // If not available, load it using the stored url.
    const VolumeBase* volume = VolumeSerializerPopulator().getVolumeSerializer()->read(url);
    if(!volume) {
        LERROR("Could not load volume: " << urlString);
        return nullptr;
    }

    // Cache result.
    cacheEntries_.insert(std::make_pair(urlString, VolumeCacheEntry{volume, true}));

    return volume;
}


EnsembleDataset::TimeStep::TimeStep()
    : TimeStep(std::map<std::string, const VolumeBase*>(), 0.0f, 0.0f)
{}

EnsembleDataset::TimeStep::TimeStep(const std::map<std::string, const VolumeBase*>& volumeData,
                                    float time, float duration)
    : time_(time)
    , duration_(duration)
    , volumeCache_(new VolumeCache(volumeData))
{
    for(const auto& vol : volumeData) {
        urls_.insert(std::make_pair(vol.first, vol.second->getOrigin()));
    }
}

float EnsembleDataset::TimeStep::getTime() const {
    return time_;
}

float EnsembleDataset::TimeStep::getDuration() const {
    return duration_;
}

std::vector<std::string> EnsembleDataset::TimeStep::getFieldNames() const {
    std::vector<std::string> fieldNames;
    for(const auto& volumeData : urls_) {
        fieldNames.push_back(volumeData.first);
    }
    return fieldNames;
}

const VolumeBase* EnsembleDataset::TimeStep::getVolume(const std::string& fieldName) const {
    auto urlIter = urls_.find(fieldName);
    if(urlIter != urls_.end()) {
        const VolumeBase* volume = volumeCache_->requestVolume(urlIter->second);
        if(!volume) {
            LERROR("Could not load volume" << urlIter->second.getURL());
            return nullptr;
        }

        return volume;
    }

    return nullptr;
}

VolumeURL EnsembleDataset::TimeStep::getURL(const std::string& fieldName) const {
    auto iter = urls_.find(fieldName);
    if(iter != urls_.end()) {
        return iter->first;
    }

    return VolumeURL();
}

void EnsembleDataset::TimeStep::serialize(Serializer& s) const {
    s.serialize("time", time_);
    s.serialize("duration", duration_);
    s.serialize("urls", urls_);
}

void EnsembleDataset::TimeStep::deserialize(Deserializer& s) {
    s.deserialize("time", time_);
    s.deserialize("duration", duration_);
    s.deserialize("urls", urls_);

    volumeCache_.reset(new VolumeCache());
}


//////////// Run

EnsembleDataset::Run::Run()
    : Run("", tgt::vec3::zero, std::vector<TimeStep>())
{}

EnsembleDataset::Run::Run(const std::string& name, const tgt::vec3& color, const std::vector<TimeStep>& timeSteps)
    : name_(name)
    , color_(color)
    , timeSteps_(timeSteps)
    , timeStepDurationStats_(false)
{
    for(const TimeStep& timeStep : timeSteps_) {
        timeStepDurationStats_.addSample(timeStep.getDuration());
    }
}

const std::string& EnsembleDataset::Run::getName() const {
    return name_;
}

const tgt::vec3& EnsembleDataset::Run::getColor() const {
    return color_;
}

const std::vector<EnsembleDataset::TimeStep>& EnsembleDataset::Run::getTimeSteps() const {
    return timeSteps_;
}

size_t EnsembleDataset::Run::getTimeStep(float time) const {
    if(timeSteps_.empty())
        return -1;

    size_t t = 0;
    while (t < getTimeSteps().size()-1 && getTimeSteps()[t].getTime() < time) t++;
    return t;
}

const Statistics& EnsembleDataset::Run::getTimeStepDurationStats() const {
    return timeStepDurationStats_;
}

void EnsembleDataset::Run::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("color", color_);
    s.serialize("timeSteps", timeSteps_);
}

void EnsembleDataset::Run::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.deserialize("color", color_);
    s.deserialize("timeSteps", timeSteps_);

    timeStepDurationStats_.reset();
    for(const TimeStep& timeStep : timeSteps_) {
        timeStepDurationStats_.addSample(timeStep.getDuration());
    }
}

//////////// Field

EnsembleDataset::FieldMetaData::FieldMetaData()
    : valueRange_(tgt::vec2::zero)
    , magnitudeRange_(tgt::vec2::zero)
    , numChannels_(0)
    , baseType_("")
{}

void EnsembleDataset::FieldMetaData::serialize(Serializer& s) const {
    s.serialize("valueRange", valueRange_);
    s.serialize("magnitudeRange", magnitudeRange_);
    s.serialize("numChannels", numChannels_);
    s.serialize("baseType", baseType_);
}

void EnsembleDataset::FieldMetaData::deserialize(Deserializer& s) {
    s.deserialize("valueRange", valueRange_);
    s.deserialize("magnitudeRange", magnitudeRange_);
    s.deserialize("numChannels", numChannels_);
    s.deserialize("baseType", baseType_);
}

//////////// EnsembleDataset

const std::string EnsembleDataset::loggerCat_ = "voreen.ensembleanalysis.EnsembleDataSet";

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
{
}

EnsembleDataset::EnsembleDataset(const EnsembleDataset& origin)
    : runs_(origin.runs_)
    , uniqueFieldNames_(origin.uniqueFieldNames_)
    , commonFieldNames_(origin.commonFieldNames_)
    , fieldMetaData_(origin.fieldMetaData_)
    , allParameters_(origin.allParameters_)
    , minNumTimeSteps_(origin.minNumTimeSteps_)
    , maxNumTimeSteps_(origin.maxNumTimeSteps_)
    , totalNumTimeSteps_(origin.totalNumTimeSteps_)
    , maxTimeStepDuration_(origin.maxTimeStepDuration_)
    , minTimeStepDuration_(origin.minTimeStepDuration_)
    , startTime_(origin.startTime_)
    , endTime_(origin.endTime_)
    , commonTimeInterval_(origin.commonTimeInterval_)
    , bounds_(origin.bounds_)
    , commonBounds_(origin.commonBounds_)
{
}

void EnsembleDataset::addRun(const Run& run) {

    // Skip empty runs.
    if (run.getTimeSteps().empty()) {
        LERROR("Can't add empty run");
        return;
    }

    // Notify Observers
    notifyPendingDataInvalidation();

    minNumTimeSteps_ = std::min(run.getTimeSteps().size(), minNumTimeSteps_);
    maxNumTimeSteps_ = std::max(run.getTimeSteps().size(), maxNumTimeSteps_);
    totalNumTimeSteps_ += run.getTimeSteps().size();
    startTime_ = std::min(startTime_, run.getTimeSteps().front().getTime());
    endTime_   = std::max(endTime_,   run.getTimeSteps().back().getTime());

    for (size_t t = 0; t < run.getTimeSteps().size(); t++) {
        std::vector<std::string> fields;
        for (const std::string& fieldName : run.getTimeSteps()[t].getFieldNames()) {

            fields.push_back(fieldName);

            // Retrieve volume.
            const VolumeBase* volume = run.getTimeSteps()[t].getVolume(fieldName);

            // Gather parameters (take first time step as representative).
            if(t==0) {
                for (const std::string& key : volume->getMetaDataKeys()) {
                    if (key.find("Parameter") != std::string::npos) {
                        allParameters_.insert(key);
                    }
                }
            }

            // In further applications it might be useful to force derived data calculations
            // to improve responsivity after loading the data.
            // TODO: Think about calculating on demand.

            VolumeMinMax* vmm = volume->getDerivedData<VolumeMinMax>();
            tgt::vec2 minMax(std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest());
            for(size_t c = 0; c < vmm->getNumChannels(); c++) {
                minMax.x = std::min(minMax.x, vmm->getMin(c));
                minMax.y = std::max(minMax.y, vmm->getMax(c));
            }

            tgt::vec2 minMaxMagnitude = minMax;
            if(volume->getNumChannels() > 1) {
                VolumeMinMaxMagnitude* vmmm = volume->getDerivedData<VolumeMinMaxMagnitude>();
                minMaxMagnitude.x = vmmm->getMinMagnitude();
                minMaxMagnitude.y = vmmm->getMaxMagnitude();
            }

            bool firstFieldElement = fieldMetaData_.find(fieldName) == fieldMetaData_.end();
            FieldMetaData& fieldMetaData = fieldMetaData_[fieldName];
            if(firstFieldElement) {
                fieldMetaData.valueRange_ = minMax;
                fieldMetaData.magnitudeRange_ = minMaxMagnitude;
                fieldMetaData.numChannels_ = volume->getNumChannels();
                fieldMetaData.baseType_ = volume->getBaseType();
            }
            else {
                fieldMetaData.valueRange_.x = std::min(fieldMetaData.valueRange_.x, minMax.x);
                fieldMetaData.valueRange_.y = std::max(fieldMetaData.valueRange_.y, minMax.y);
                fieldMetaData.magnitudeRange_.x = std::min(fieldMetaData.magnitudeRange_.x, minMaxMagnitude.x);
                fieldMetaData.magnitudeRange_.y = std::max(fieldMetaData.magnitudeRange_.y, minMaxMagnitude.y);
                if(fieldMetaData.numChannels_ != volume->getNumChannels()) {
                    LERRORC("voreen.EnsembleDataSet", "Number of channels differs per field, taking min.");
                    fieldMetaData.numChannels_ = std::min(fieldMetaData.numChannels_, volume->getNumChannels());
                }
                if(fieldMetaData.baseType_ != volume->getBaseType()) {
                    LERRORC("voreen.EnsembleDataSet", "Base type differs per field, taking first.");
                }
            }

            tgt::Bounds bounds = volume->getBoundingBox().getBoundingBox();
            if (!bounds_.isDefined()) {
                if(!commonBounds_.isDefined()) {
                    commonBounds_.addVolume(bounds);
                }
            }
            else if(commonBounds_.isDefined()) {
                commonBounds_.intersectVolume(bounds);
                if(!commonBounds_.isDefined()) {
                    LWARNINGC("voreen.EnsembeDataSet", "There is no overlap between the bounds of Run " << run.getName() << " and the previously defined bounds");
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
                LWARNINGC("voreen.EnsembeDataSet", "Time Step " << t << " of Run " << run.getName() << " has less fields than the previously added Run " << runs_.back().getName());
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
        if (t < run.getTimeSteps().size() - 1) {
            maxTimeStepDuration_ = std::max(maxTimeStepDuration_, run.getTimeSteps()[t].getDuration());
            minTimeStepDuration_ = std::min(minTimeStepDuration_, run.getTimeSteps()[t].getDuration());
        }
    }

    commonTimeInterval_.x = std::max(commonTimeInterval_.x, run.getTimeSteps().front().getTime());
    commonTimeInterval_.y = std::min(commonTimeInterval_.y, run.getTimeSteps().back().getTime()+run.getTimeSteps().back().getDuration());

    if(commonTimeInterval_ != tgt::vec2::zero && commonTimeInterval_.x > commonTimeInterval_.y) {
        LWARNINGC("voreen.EnsembleDataSet", "The time interval of the currently added Run " << run.getName() << " does not overlap with the previous interval");
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

const tgt::Bounds& EnsembleDataset::getBounds() const {
    return bounds_;
}

const tgt::Bounds& EnsembleDataset::getCommonBounds() const {
    return commonBounds_;
}

const tgt::vec2& EnsembleDataset::getValueRange(const std::string& field) const {
    tgtAssert(fieldMetaData_.find(field) != fieldMetaData_.end(), "Field not available");
    return fieldMetaData_.at(field).valueRange_;
}

const tgt::vec2& EnsembleDataset::getMagnitudeRange(const std::string& field) const {
    tgtAssert(fieldMetaData_.find(field) != fieldMetaData_.end(), "Field not available");
    return fieldMetaData_.at(field).magnitudeRange_;
}

size_t EnsembleDataset::getNumChannels(const std::string& field) const {
    tgtAssert(fieldMetaData_.find(field) != fieldMetaData_.end(), "Field not available");
    return fieldMetaData_.at(field).numChannels_;
}

const std::string& EnsembleDataset::getBaseType(const std::string& field) const {
    tgtAssert(fieldMetaData_.find(field) != fieldMetaData_.end(), "Field not available");
    return fieldMetaData_.at(field).baseType_;
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
        for(const TimeStep& timeStep : run.getTimeSteps()) {
            for(const std::string& fieldName : timeStep.getFieldNames()) {
                result.push_back(timeStep.getVolume(fieldName));
            }
        }
    }
    return result;
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
    for(const Run& run : runs_) {
        stream << "  <tr>\n";
        stream << "    <th>" << run.getName() << "</th>\n";
        tgt::ivec3 color(run.getColor() * 255.0f);
        stream << "    <th style=\"background-color: rgb(" << color.r << ", " << color.g << ", " << color.b << ")\"></th>\n";
        stream << "    <th>" << run.getTimeSteps().size() << "</th>\n";
        stream << "    <th>" << run.getTimeSteps().front().getTime() << "</th>\n";
        stream << "    <th>" << run.getTimeSteps().back().getTime() << "</th>\n";

        for(const std::string& parameter : allParameters_) {
            const TimeStep& referenceTimeStep = run.getTimeSteps().front();
            // TODO: assumes that all fields contain the same parameters.
            const VolumeBase* referenceVolume = referenceTimeStep.getVolume(referenceTimeStep.getFieldNames().front());
            stream << "    <th>";
            if(referenceVolume->hasMetaData(parameter)) {
                stream << referenceVolume->getMetaData(parameter)->toString();
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

void EnsembleDataset::serialize(Serializer& s) const {
    s.serialize("runs", runs_);
    s.serialize("uniqueFieldNames", uniqueFieldNames_);
    s.serialize("commonFieldNames", commonFieldNames_);

    s.serialize("fieldMetaData", fieldMetaData_);
    s.serialize("allParameters", allParameters_);

    s.serialize("minNumTimeSteps", minNumTimeSteps_);
    s.serialize("maxNumTimeSteps", maxNumTimeSteps_);
    s.serialize("totalNumTimeSteps", totalNumTimeSteps_);

    s.serialize("minTimeStepDuration", minTimeStepDuration_);
    s.serialize("maxTimeStepDuration", maxTimeStepDuration_);
    s.serialize("startTime", startTime_);
    s.serialize("endTime", endTime_);
    s.serialize("commonTimeInterval",commonTimeInterval_);

    s.serialize("bounds", bounds_);
    s.serialize("commonBounds", commonBounds_);
}

void EnsembleDataset::deserialize(Deserializer& s) {
    s.deserialize("runs", runs_);
    s.deserialize("uniqueFieldNames", uniqueFieldNames_);
    s.deserialize("commonFieldNames", commonFieldNames_);

    s.deserialize("fieldMetaData", fieldMetaData_);
    s.deserialize("allParameters", allParameters_);

    s.deserialize("minNumTimeSteps", minNumTimeSteps_);
    s.deserialize("maxNumTimeSteps", maxNumTimeSteps_);
    s.deserialize("totalNumTimeSteps", totalNumTimeSteps_);

    s.deserialize("minTimeStepDuration", minTimeStepDuration_);
    s.deserialize("maxTimeStepDuration", maxTimeStepDuration_);
    s.deserialize("startTime", startTime_);
    s.deserialize("endTime", endTime_);
    s.deserialize("commonTimeInterval",commonTimeInterval_);

    s.deserialize("bounds", bounds_);
    s.deserialize("commonBounds", commonBounds_);
}

}   // namespace
