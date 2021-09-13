/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_TIMESERIES_H
#define VRN_TIMESERIES_H

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/datastructures/datainvalidationobserver.h"
#include "voreen/core/io/serialization/serializable.h"

#include "tgt/vector.h"
#include <map>

namespace voreen {
    struct VRN_CORE_API TimeSeriesStep : public Serializable {
        TimeSeriesStep();
        TimeSeriesStep(float time, const std::vector<float>& fieldValues);

        virtual void serialize(Serializer& s) const;
        virtual void deserialize(Deserializer& s);

        float time_;
        std::vector<float> fieldValues_;
    };

    class VRN_CORE_API TimeSeries : public Serializable {
    public:
        TimeSeries(tgt::ivec3 position = tgt::ivec3(), float startTime = -1.f, float endTime = -1.f, 
            std::vector<TimeSeriesStep> timeSteps = std::vector<TimeSeriesStep>());

        tgt::ivec3 getPosition() const;
        void setPosition(const tgt::ivec3& pos);
        
        float getStartTime() const;
        void setStartTime(float startTime);
        float getEndTime() const;
        void setEndTime(float endTime);
        size_t getTimeStepIndex(float time) const;

        const std::vector<TimeSeriesStep>& getTimeSeriesSteps() const;
        void addTimeSeriesStep(TimeSeriesStep& timeSeriesStep);
        const TimeSeriesStep& getTimeSeriesStep(size_t idx) const;
        size_t getNumberOfTimeSteps() const;

        virtual void serialize(Serializer& s) const;
        virtual void deserialize(Deserializer& s);

    private:
        tgt::ivec3 position_;
        float startTime_;
        float endTime_;
        std::vector<TimeSeriesStep> timeSteps_;
    };

    class VRN_CORE_API TimeSeriesList : public Serializable {
    public:
        TimeSeriesList(float overallStartTime = -1.f, float overallEndTime = -1.f,
            std::vector<std::string> fieldNames = std::vector<std::string>(),
            std::vector<unsigned int> componentsPerField = std::vector<unsigned int>());

        float getOverallStartTime() const;
        void setOverallStartTime(float overallStartTime);
        float getOverallEndTime() const;
        void setOverallEndTime(float overallEndTime);

        tgt::ivec3 getVolumeDimensions() const;
        void setVolumeDimensions(tgt::ivec3 dimensions);

        const std::vector<std::string>& getFieldNames() const;
        const std::string& getFieldName(size_t idx) const;
        void setFieldNames(const std::vector<std::string>& fieldNames);
        size_t getNumberOfFields() const;
        void setRange(std::string field, tgt::vec2 range);
        tgt::vec2 getRange(std::string field) const;

        const std::vector<unsigned int>& getComponentsPerField() const;
        unsigned int getComponentForField(size_t idx) const;
        void setComponentsPerField(const std::vector<unsigned int> componentsPerField);
        size_t getNumberOfComponentsPerField() const;

        const std::vector<TimeSeries>& getAllTimeSeries() const;
        const TimeSeries& getTimeSeries(size_t idx) const;
        TimeSeries* getTimeSeriesEditable(size_t idx);
        size_t getNumberOfTimeSeries() const;
        void addTimeSeries(const TimeSeries& timeseries);

        virtual void serialize(Serializer& s) const;
        virtual void deserialize(Deserializer& s);

    private:
        float overallStartTime_;
        float overallEndTime_;
        tgt::ivec3 volumeDimensions_;
        std::vector<std::string> fieldNames_;
        std::vector<unsigned int> componentsPerField_;
        std::vector<TimeSeries> allTimeSeries_;
        // ranges for normalization
        std::map<std::string, tgt::vec2> ranges_;
    };
}   // namespace

#endif
