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

#ifndef VRN_MEDIANFILTER_H
#define VRN_MEDIANFILTER_H

#include "parallelvolumefilter.h"

namespace voreen {

class MedianFilter : public ParallelVolumeFilter<ParallelFilterValue1D, ParallelFilterValue1D> {
public:
    MedianFilter(const tgt::ivec3& extent, const SamplingStrategy<ParallelFilterValue1D>& samplingStrategy, const std::string& sliceBaseType);
    virtual ~MedianFilter();
    ParallelFilterValue1D getValue(const Sample& sample, const tgt::ivec3& pos) const;
private:
    tgt::ivec3 extent_;
};

template<typename  T>
class MedianFilterVector : public ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>> {
public:
    MedianFilterVector(const tgt::ivec3& extent, const SamplingStrategy<T>& samplingStrategy, const std::string& sliceBaseType);
    virtual ~MedianFilterVector();
    ParallelFilterValue<T> getValue(const typename MedianFilterVector<T>::Sample& sample, const tgt::ivec3& pos) const;
private:
    tgt::ivec3 extent_;
};

template<typename T>
MedianFilterVector<T>::MedianFilterVector(const tgt::ivec3& extent, const SamplingStrategy<T>& samplingStrategy, const std::string& sliceBaseType)
    : ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>>(extent.z, samplingStrategy, sliceBaseType)
    , extent_(extent)
{
}

template<typename T>
MedianFilterVector<T>::~MedianFilterVector() {
}

template<typename T>
ParallelFilterValue<T> MedianFilterVector<T>::getValue(const typename MedianFilterVector<T>::Sample& sample, const tgt::ivec3& pos) const {
    std::vector<T> values;
    values.reserve(tgt::hmul(extent_));

    for(int z = pos.z-extent_.z; z <= pos.z+extent_.z; ++z) {
        for(int y = pos.y-extent_.y; y <= pos.y+extent_.y; ++y) {
            for(int x = pos.x-extent_.x; x <= pos.x+extent_.x; ++x) {
                values.push_back(sample(tgt::ivec3(x,y,z)));
            }
        }
    }

    size_t argMin = 0;
    float minSum = std::numeric_limits<float>::max();
    for(size_t i = 0; i < values.size(); i++) {
        float sum = 0.0f;
        for(size_t j = 0; j < values.size(); j++) {
            sum += tgt::length(values[j] - values[i]);
        }
        if(sum < minSum) {
            minSum = sum;
            argMin = i;
        }
    }

    return values[argMin];
}

typedef MedianFilterVector<tgt::vec2> MedianFilter2D;
typedef MedianFilterVector<tgt::vec3> MedianFilter3D;
typedef MedianFilterVector<tgt::vec4> MedianFilter4D;

} // namespace voreen

#endif // VRN_MEDIANFILTER_H
