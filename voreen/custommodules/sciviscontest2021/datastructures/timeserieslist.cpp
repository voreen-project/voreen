/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "timeserieslist.h"
#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"

namespace voreen {
    TimeSeriesStep::TimeSeriesStep() : time_(-1.f) 
    {

    }
    TimeSeriesStep::TimeSeriesStep(float time, const std::vector<float>& fieldValues) 
        : time_(time), fieldValues_(fieldValues) 
    {

    }

    void TimeSeriesStep::serialize(Serializer& s) const {
        s.serialize("time", time_);
        s.serialize("fieldValues", fieldValues_);
    }
    void TimeSeriesStep::deserialize(Deserializer& s) {
        s.deserialize("time", time_);
        s.deserialize("fieldValues", fieldValues_);
    }

    TimeSeries::TimeSeries(tgt::ivec3 position, float startTime, float endTime, std::vector<TimeSeriesStep> timeSteps) 
        : position_(position), startTime_(startTime), endTime_(endTime), timeSteps_(timeSteps)
    {

    }

    tgt::ivec3 TimeSeries::getPosition() const {
        return position_;
    }
    void TimeSeries::setPosition(const tgt::ivec3& pos) {
        position_ = pos;
    }

    float TimeSeries::getStartTime() const {
        return startTime_;
    }
    void TimeSeries::setStartTime(float startTime) {
        startTime_ = startTime;
    }
    float TimeSeries::getEndTime() const {
        return endTime_;
    }
    void TimeSeries::setEndTime(float endTime) {
        endTime_ = endTime;
    }
    size_t TimeSeries::getTimeStepIndex(float time) const {
        if(timeSteps_.empty())
            return -1;

        size_t t = 0;
        while (t < getNumberOfTimeSteps()-1 && getTimeSeriesStep(t).time_ < time) t++;
        return t;
    }

    const std::vector<TimeSeriesStep>& TimeSeries::getTimeSeriesSteps() const {
        return timeSteps_;
    }
    void TimeSeries::addTimeSeriesStep(TimeSeriesStep& timeSeriesStep) {
        timeSteps_.push_back(timeSeriesStep);
    }
    const TimeSeriesStep& TimeSeries::getTimeSeriesStep(size_t idx) const {
        return timeSteps_.at(idx);
    }
    size_t TimeSeries::getNumberOfTimeSteps() const {
        return timeSteps_.size();
    }

    void TimeSeries::serialize(Serializer& s) const {
        s.serialize("position", position_);
        s.serialize("startTime", startTime_);
        s.serialize("endTime", endTime_);
        s.serialize("timeSteps", timeSteps_);
    }
    void TimeSeries::deserialize(Deserializer& s) {
        s.deserialize("position", position_);
        s.deserialize("startTime", startTime_);
        s.deserialize("endTime", endTime_);
        s.deserialize("timeSteps", timeSteps_);
    }

    TimeSeriesList::TimeSeriesList(float overallStartTime, float overallEndTime,
        std::vector<std::string> fieldNames,
        std::vector<unsigned int> componentsPerField)
        : overallStartTime_(overallStartTime), overallEndTime_(overallEndTime), fieldNames_(fieldNames),
        componentsPerField_(componentsPerField)
    {

    }

    float TimeSeriesList::getOverallStartTime() const {
        return overallStartTime_;
    }
    void TimeSeriesList::setOverallStartTime(float overallStartTime) {
        overallStartTime_ = overallStartTime;
    }
    float TimeSeriesList::getOverallEndTime() const {
        return overallEndTime_;
    }
    void TimeSeriesList::setOverallEndTime(float overallEndTime) {
        overallEndTime_ = overallEndTime;
    }

    tgt::ivec3 TimeSeriesList::getVolumeDimensions() const {
        return volumeDimensions_;
    }
    void TimeSeriesList::setVolumeDimensions(tgt::ivec3 dimensions) {
        volumeDimensions_ = dimensions;
    }

    const std::vector<std::string>& TimeSeriesList::getFieldNames() const {
        return fieldNames_;
    }
    const std::string& TimeSeriesList::getFieldName(size_t idx) const {
        return fieldNames_.at(idx);
    }
    void TimeSeriesList::setFieldNames(const std::vector<std::string>& fieldNames) {
        fieldNames_ = fieldNames;
    }
    size_t TimeSeriesList::getNumberOfFields() const {
        return fieldNames_.size();
    }
    void TimeSeriesList::setRange(std::string field, tgt::vec2 range) {
        ranges_[field] = range;
    }
    tgt::vec2 TimeSeriesList::getRange(std::string field) const {
        return ranges_.at(field);
    }

    const std::vector<unsigned int>& TimeSeriesList::getComponentsPerField() const {
        return componentsPerField_;
    }
    unsigned int TimeSeriesList::getComponentForField(size_t idx) const {
        return componentsPerField_.at(idx);
    }
    void TimeSeriesList::setComponentsPerField(const std::vector<unsigned int> componentsPerField){
        componentsPerField_ = componentsPerField;
    }
    size_t TimeSeriesList::getNumberOfComponentsPerField() const {
        return componentsPerField_.size();
    }

    const std::vector<TimeSeries>& TimeSeriesList::getAllTimeSeries() const {
        return allTimeSeries_;
    }
    const TimeSeries& TimeSeriesList::getTimeSeries(size_t idx) const {
        return allTimeSeries_.at(idx);
    }
    TimeSeries* TimeSeriesList::getTimeSeriesEditable(size_t idx) {
        return &allTimeSeries_.at(idx);
    }
    size_t TimeSeriesList::getNumberOfTimeSeries() const {
        return allTimeSeries_.size();
    }
    void TimeSeriesList::addTimeSeries(const TimeSeries& timeseries) {
        allTimeSeries_.push_back(timeseries);
    }

    void TimeSeriesList::serialize(Serializer& s) const{
        s.serialize("overallStartTime", overallStartTime_);
        s.serialize("overallEndTime", overallEndTime_);
        s.serialize("fieldNames", fieldNames_);
        s.serialize("componentsPerField", componentsPerField_);
        s.serialize("allTimeSeries", allTimeSeries_);
        s.serialize("ranges", ranges_);
        s.serialize("volumeDimensions", volumeDimensions_);
    }
    void TimeSeriesList::deserialize(Deserializer& s){
        s.deserialize("overallStartTime", overallStartTime_);
        s.deserialize("overallEndTime", overallEndTime_);
        s.deserialize("fieldNames", fieldNames_);
        s.deserialize("componentsPerField", componentsPerField_);
        s.deserialize("allTimeSeries", allTimeSeries_);
        s.deserialize("ranges", ranges_);
        s.deserialize("volumeDimensions", volumeDimensions_);
    }
}   // namespace
