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

#ifndef VRN_VALUEMAPFILTER_H
#define VRN_VALUEMAPFILTER_H

#include "parallelvolumefilter.h"

#include <functional>

namespace voreen {

typedef std::vector<uint8_t> ValueMap;

template<typename T>
class ValueMapFilter : public ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>> {
public:

    ValueMapFilter(ValueMap&& valueMap, RealWorldMapping inputToLut_);
    virtual ~ValueMapFilter();

    ParallelFilterValue<T> getValue(const typename ValueMapFilter<T>::Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetadata, const SliceReaderMetaData& outputMetaData) const;
    virtual SliceReaderMetaData getMetaData(const SliceReaderMetaData& base) const;

private:
    ValueMap valueMap_;
    float minNormalized_;
    float maxNormalized_;
    RealWorldMapping inputToLut_;
};

// Implementation

static float mapValueToNormalized(uint8_t val) {
    return static_cast<float>(val) / std::numeric_limits<ValueMap::value_type>::max();
}

template<typename T>
ValueMapFilter<T>::ValueMapFilter(std::vector<uint8_t>&& valueMap, RealWorldMapping inputToLut)
    : ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>>(0, SamplingStrategy<ParallelFilterValue<T>>::ASSERT_FALSE)
    , valueMap_(std::move(valueMap))
    , minNormalized_( std::numeric_limits<float>::max())
    , maxNormalized_(-std::numeric_limits<float>::max())
    , inputToLut_(inputToLut)
{
    ValueMap::value_type min = std::numeric_limits<ValueMap::value_type>::max();
    ValueMap::value_type max = std::numeric_limits<ValueMap::value_type>::min();
    for(auto val: valueMap_) {
        min = std::min(min, val);
        max = std::max(max, val);
    }
    minNormalized_ = mapValueToNormalized(min);
    maxNormalized_ = mapValueToNormalized(max);
}

template<typename T>
ValueMapFilter<T>::~ValueMapFilter() {
}

template<typename T>
ParallelFilterValue<T> ValueMapFilter<T>::getValue(const typename ValueMapFilter<T>::Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetadata, const SliceReaderMetaData& outputMetaData) const {
    T inputNorm = sample(pos);
    return mapScalars(inputNorm, [&] (float val) {
            float lutVal = inputToLut_.normalizedToRealWorld(val);
            size_t size = valueMap_.size();
            size_t index = tgt::clamp(static_cast<size_t>(tgt::round(lutVal * (size-1))), 0UL, size);
            uint8_t mapVal = valueMap_[index];
            return mapValueToNormalized(mapVal);
            });
}

template<typename T>
SliceReaderMetaData ValueMapFilter<T>::getMetaData(const SliceReaderMetaData& base) const {
    auto md = SliceReaderMetaData::fromBase(base);
    const auto& rwm = md.getRealworldMapping();

    //TODO: check min/max via buckets in (input-)min/max range
    tgt::vec2 minmax(rwm.normalizedToRealWorld(minNormalized_), rwm.normalizedToRealWorld(maxNormalized_));
    md.setMinMaxBounds(std::vector<tgt::vec2>(base.getNumChannels(), minmax));

    return md;
}


typedef ValueMapFilter<float>     ValueMapFilter1D;
typedef ValueMapFilter<tgt::vec2> ValueMapFilter2D;
typedef ValueMapFilter<tgt::vec3> ValueMapFilter3D;
typedef ValueMapFilter<tgt::vec4> ValueMapFilter4D;

} // namespace voreen

#endif // VRN_VALUEMAPFILTER_H
