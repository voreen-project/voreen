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

#ifndef VRN_THRESHOLDINGFILTER_H
#define VRN_THRESHOLDINGFILTER_H

#include "parallelvolumefilter.h"

#include <functional>

namespace voreen {

enum ThresholdingStrategyType {
    LOWER_T,
    UPPER_T,
};

template<typename T>
class ThresholdingFilter : public ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>> {
public:

    ThresholdingFilter(float threshold, const T& replacement, ThresholdingStrategyType thresholdingStrategyType);
    virtual ~ThresholdingFilter();

    ParallelFilterValue<T> getValue(const typename ThresholdingFilter<T>::Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetadata, const SliceReaderMetaData& outputMetaData) const;

    ThresholdingStrategyType getThresholdingStrategyType() const;

private:

    std::function<bool(float, float)> strategy_;

    const float threshold_;
    const T replacement_;
    const ThresholdingStrategyType thresholdingStrategyType_;
};

// Implementation

template<typename T>
ThresholdingFilter<T>::ThresholdingFilter(float threshold, const T& replacement, ThresholdingStrategyType thresholdingStrategyType)
    : ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>>(0, SamplingStrategy<ParallelFilterValue<T>>::ASSERT_FALSE)
    , threshold_(threshold)
    , replacement_(replacement)
    , thresholdingStrategyType_(thresholdingStrategyType)
{
    switch (thresholdingStrategyType_) {
    case LOWER_T:
        strategy_ = std::less<T>();
        break;
    case UPPER_T:
        strategy_ = std::greater<T>();
        break;
    default:
        tgtAssert(false, "Unimplemented Thresholding Strategy");
        break;
    }
}

template<typename T>
ThresholdingFilter<T>::~ThresholdingFilter() {
}

template<typename T>
static inline float getAbsoluteMagnitudeInternal(const T& value, typename std::enable_if<std::is_same<T, float>::value>::type* = 0) {
    return value;
}

template<typename T>
static inline float getAbsoluteMagnitudeInternal(const T& value, typename std::enable_if<!std::is_same<T, float>::value>::type* = 0) {
    return tgt::length(value);
}

template<typename T>
ParallelFilterValue<T> ThresholdingFilter<T>::getValue(const typename ThresholdingFilter<T>::Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetadata, const SliceReaderMetaData& outputMetaData) const {

    T sampleValue = sample(pos);
    float value = getAbsoluteMagnitudeInternal<T>(sampleValue);
    bool replace = strategy_(value, threshold_);

    if (replace) {
        value = replacement_;
    }

    return value;
}

template<typename T>
ThresholdingStrategyType ThresholdingFilter<T>::getThresholdingStrategyType() const {
    return thresholdingStrategyType_;
}

typedef ThresholdingFilter<float>     ThresholdingFilter1D;
typedef ThresholdingFilter<tgt::vec2> ThresholdingFilter2D;
typedef ThresholdingFilter<tgt::vec3> ThresholdingFilter3D;
typedef ThresholdingFilter<tgt::vec4> ThresholdingFilter4D;

} // namespace voreen

#endif // VRN_THRESHOLDINGFILTER_H
