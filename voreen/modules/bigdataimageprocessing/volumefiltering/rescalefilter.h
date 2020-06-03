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

#ifndef VRN_RESCALEFILTER_H
#define VRN_RESCALEFILTER_H

#include "parallelvolumefilter.h"

#include <functional>

namespace voreen {

enum RescaleStrategyType {
    RESCALE_LOGARITHMIC_T = 0,
    RESCALE_EXPONENTIAL_T = 1,
};

template<typename T>
class RescaleFilter : public ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>> {
public:

    RescaleFilter(RescaleStrategyType strategy, const std::string& sliceBaseType);
    virtual ~RescaleFilter();

    ParallelFilterValue<T> getValue(const typename RescaleFilter<T>::Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetadata, const SliceReaderMetaData& outputMetaData) const;
    virtual SliceReaderMetaData getMetaData(const SliceReaderMetaData& base) const;

private:
    const RescaleStrategyType strategy_;
};

// Implementation

template<typename T>
RescaleFilter<T>::RescaleFilter(RescaleStrategyType strategy, const std::string& sliceBaseType)
    : ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>>(0, SamplingStrategy<ParallelFilterValue<T>>::ASSERT_FALSE, sliceBaseType)
    , strategy_(strategy)
{
}

template<typename T>
RescaleFilter<T>::~RescaleFilter() {
}

template<typename T, typename F>
static inline T mapT(const T& value, F f, typename std::enable_if<!std::is_same<T, float>::value>::type* = 0) {
    T out;
    for(int c=0; c<T::dim; ++c) {
        out[c] = f(value[c]);
    }
    return out;
}

template<typename T, typename F>
static inline T mapT(const T& value, F f, typename std::enable_if<std::is_same<T, float>::value>::type* = 0) {
    return f(value);
}

template<typename T>
static inline T rescaleInternal(const T& value, RescaleStrategyType strategy) {
    switch (strategy) {
        case RESCALE_LOGARITHMIC_T:
            return mapT(value, [] (float i) -> float {
                        float res = std::log(i);
                        if(std::isinf(res)) {
                            return -std::numeric_limits<float>::max();
                        } else {
                            return res;
                        }
                    });
        case RESCALE_EXPONENTIAL_T:
            return mapT(value, [] (float i) -> float {
                        float res = std::exp(i);
                        if(std::isinf(res)) {
                            return std::numeric_limits<float>::max();
                        } else {
                            return res;
                        }
                    });
        default:
            tgtAssert(false, "Invalid strategy");
    }
}

template<typename T>
ParallelFilterValue<T> RescaleFilter<T>::getValue(const typename RescaleFilter<T>::Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetadata, const SliceReaderMetaData& outputMetaData) const {
    T inputNorm = sample(pos);
    T inputRW = mapT(inputNorm, [&] (T i) { return inputMetadata.getRealworldMapping().normalizedToRealWorld(i);});
    T outputRW  = rescaleInternal(inputRW, strategy_);
    T outputNorm = mapT(inputNorm, [&] (T i) { return outputMetaData.getRealworldMapping().realWorldToNormalized(i);});
    return outputNorm;
}

template<typename T>
SliceReaderMetaData RescaleFilter<T>::getMetaData(const SliceReaderMetaData& base) const {
    const auto& baseRwm = base.getRealworldMapping();
    float min = rescaleInternal(baseRwm.getOffset(), strategy_);
    float max = rescaleInternal(baseRwm.getOffset() + baseRwm.getScale(), strategy_);
    SliceReaderMetaData md = [&] () {
        switch (strategy_) {
            case RESCALE_LOGARITHMIC_T:
                return SliceReaderMetaData(RealWorldMapping(
                            tgt::vec2(min, max),
                            "log(" + baseRwm.getUnit() + ")"));
            case RESCALE_EXPONENTIAL_T:
                return SliceReaderMetaData(RealWorldMapping(
                            tgt::vec2(min, max),
                            "exp(" + baseRwm.getUnit() + ")"));
            default:
                tgtAssert(false, "Invalid strategy");
        }
    } ();
    if(base.isAccurate()) {
        auto vmm = base.getVolumeMinMax();
        for(int c=0; c<base.getNumChannels(); ++c) {
            md.setMinMax(
                    rescaleInternal(vmm->getMin(c), strategy_),
                    rescaleInternal(vmm->getMax(c), strategy_),
                    c);
        }
        md.markAccurate();
    }
    return md;
}


typedef RescaleFilter<float>     RescaleFilter1D;
typedef RescaleFilter<tgt::vec2> RescaleFilter2D;
typedef RescaleFilter<tgt::vec3> RescaleFilter3D;
typedef RescaleFilter<tgt::vec4> RescaleFilter4D;

} // namespace voreen

#endif // VRN_RESCALEFILTER_H
