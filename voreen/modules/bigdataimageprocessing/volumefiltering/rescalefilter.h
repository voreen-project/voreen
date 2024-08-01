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

    RescaleFilter(RescaleStrategyType strategy);
    virtual ~RescaleFilter();

    ParallelFilterValue<T> getValue(const typename RescaleFilter<T>::Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetadata, const SliceReaderMetaData& outputMetaData) const;
    virtual SliceReaderMetaData getMetaData(const SliceReaderMetaData& base) const;

private:
    const RescaleStrategyType strategy_;
};

// Implementation

template<typename T>
RescaleFilter<T>::RescaleFilter(RescaleStrategyType strategy)
    : ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>>(0, SamplingStrategy<ParallelFilterValue<T>>::ASSERT_FALSE)
    , strategy_(strategy)
{
}

template<typename T>
RescaleFilter<T>::~RescaleFilter() {
}

template<typename T>
static inline T rescaleInternal(const T& value, RescaleStrategyType strategy) {
    switch (strategy) {
        case RESCALE_LOGARITHMIC_T:
            return mapScalars(value, [] (float i) -> float {
                        float res = std::log(i);
                        if(std::isinf(res)) {
                            return -std::numeric_limits<float>::max();
                        } else {
                            return res;
                        }
                    });
        case RESCALE_EXPONENTIAL_T:
            return mapScalars(value, [] (float i) -> float {
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
    T inputRW = mapScalars(inputNorm, [&] (float i) { return inputMetadata.getRealworldMapping().normalizedToRealWorld(i);});
    T outputRW  = rescaleInternal(inputRW, strategy_);
    T outputNorm = mapScalars(outputRW, [&] (float i) { return outputMetaData.getRealworldMapping().realWorldToNormalized(i);});
    return outputNorm;
}

template<typename T>
SliceReaderMetaData RescaleFilter<T>::getMetaData(const SliceReaderMetaData& base) const {
    const auto& baseRwm = base.getRealworldMapping();

    tgt::vec2 minmax = base.estimateMinMax();

    RealWorldMapping rwm = [&] () {
        auto baseUnit = baseRwm.getUnit().empty() ? "x" : baseRwm.getUnit();
        switch (strategy_) {
            case RESCALE_LOGARITHMIC_T:
                return RealWorldMapping(
                            rescaleInternal(minmax, strategy_),
                            "log(" + baseUnit + ")");
            case RESCALE_EXPONENTIAL_T:
                return RealWorldMapping(
                            rescaleInternal(minmax, strategy_),
                            "exp(" + baseUnit + ")");
            default:
                tgtAssert(false, "Invalid strategy");
        }
    } ();
    auto md = SliceReaderMetaData::fromBase(base);
    md.setRealWorldMapping(rwm);

    // If we have accurate min/max information from the base, we can also supply those ourselves.
    if(base.getMinMax()) {
        // Note that this mapping is ONLY possible because all supported functions are monotonically increasing!
        std::vector<tgt::vec2> minmax;
        for(auto& mm : *base.getMinMax()) {
            minmax.emplace_back(rescaleInternal(mm, strategy_));
        }
        md.setMinMax(minmax);
    }
    return md;
}


typedef RescaleFilter<float>     RescaleFilter1D;
typedef RescaleFilter<tgt::vec2> RescaleFilter2D;
typedef RescaleFilter<tgt::vec3> RescaleFilter3D;
typedef RescaleFilter<tgt::vec4> RescaleFilter4D;

} // namespace voreen

#endif // VRN_RESCALEFILTER_H
