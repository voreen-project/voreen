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

#ifndef VRN_VOLUMEFILTER_H
#define VRN_VOLUMEFILTER_H

#include "voreen/core/datastructures/volume/volumeram.h"

#include <tgt/vector.h>

#include <functional>
#include <memory>
#include <boost/optional.hpp>

namespace voreen {

// Forward declaration, defined in slicereader.h
class CachingSliceReader;
class SliceReaderMetaData;

// The VolumeFilter interface ---------------------------------------------------------------------------
class VolumeFilter {
public:
    virtual ~VolumeFilter() {}

    virtual std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const = 0;
    virtual int zExtent() const = 0;
    virtual size_t getNumInputChannels() const = 0;
    virtual size_t getNumOutputChannels() const = 0;
    virtual boost::optional<tgt::svec3> getOverwrittenDimensions() const;
    virtual SliceReaderMetaData getMetaData(const SliceReaderMetaData& base) const;
};

// ------------------------------------------------------------------------------------------------------
// Useful classes for implementation of VolumeFilters
// ------------------------------------------------------------------------------------------------------

// SamplingStrategy: Handle reading outside of volumes/slices of arbitrary types
enum SamplingStrategyType {
    MIRROR_T,
    CLAMP_T,
    PASS_THROUGH_T,
    SET_T,
    ASSERT_FALSE_T
};
template<typename T>
struct SamplingStrategy {
public:

    const static SamplingStrategy MIRROR;
    const static SamplingStrategy CLAMP;
    const static SamplingStrategy PASS_THROUGH;
    const static SamplingStrategy ASSERT_FALSE;
    static SamplingStrategy SET(T outsideVolumeValue);

    SamplingStrategy(SamplingStrategyType type, T outsideVolumeValue);

    const SamplingStrategyType type_;
    const T outsideVolumeValue_;

    typedef std::function<T(const tgt::ivec3& p)> Sampler;

    T sample(const tgt::ivec3& p, const tgt::ivec3& dim, const Sampler& getValue) const;


    // Conversion between SamplingStrategy types
    template<typename O>
    SamplingStrategy(const SamplingStrategy<O>& other);

    template<typename N>
    SamplingStrategy<N> convert(std::function<N(T)> defaultConvert) const;
};


// SimpleSlice<T>. Temporary Storage for non-float slice data
template<typename T>
class SimpleSlice {
public:
    SimpleSlice(const tgt::ivec2& dim)
        : dim_(dim)
        , data_(tgt::hmul(dim))
    {
    }

    T& at(int x, int y) {
        return data_[x+dim_.x*y];
    }

    const T& at(int x, int y) const {
        return data_[x+dim_.x*y];
    }

    typename SamplingStrategy<T>::Sampler toSampler() {
        return [this] (const tgt::ivec3& p) {
            return at(p.x, p.y);
        };
    }

private:
    const tgt::ivec2 dim_;
    std::vector<T> data_;
};



// ------------------------------------------------------------------------------------------
// Implementation ---------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------


// SamplingStrategy ---------------------------------------------------------------------------

static int mirror(int p, int dim) {
    if(p<0) {
        p = -p-1;
    }
    p %= 2*dim;
    if(p >= dim) {
        return 2*dim - p - 1;
    } else {
        return p;
    }
}

static int clamp(int p, int dim) {
    if(p < 0) {
        return 0;
    } else if(p >= dim) {
        return dim-1;
    } else {
        return p;
    }
}

template<typename T>
T SamplingStrategy<T>::sample(const tgt::ivec3& p, const tgt::ivec3& dim, const Sampler& getValue) const {
    if(p.x < 0 || p.x >= dim.x || p.y < 0 || p.y >= dim.y || p.z < 0 || p.z >= dim.z) {
        switch(type_) {
            case MIRROR_T:
                return getValue(tgt::ivec3(mirror(p.x, dim.x), mirror(p.y, dim.y), mirror(p.z, dim.z)));
            case CLAMP_T:
                return getValue(tgt::ivec3(clamp(p.x, dim.x), clamp(p.y, dim.y), clamp(p.z, dim.z)));
            case SET_T:
                return outsideVolumeValue_;
            case PASS_THROUGH_T:
                return getValue(p);
            case ASSERT_FALSE_T:
                tgtAssert(false, "Sampled outside volume");
                return T(0);
            default:
                tgtAssert(false, "Unimplemented sampling strategy type");
                return T(0);
        }
    } else {
        return getValue(p);
    }
}

template<typename T>
SamplingStrategy<T>::SamplingStrategy(SamplingStrategyType type, T outsideVolumeValue)
    : type_(type)
    , outsideVolumeValue_(outsideVolumeValue)
{
}

template<typename T>
const SamplingStrategy<T> SamplingStrategy<T>::MIRROR = {
    SamplingStrategyType::MIRROR_T, T()
};

template<typename T>
const SamplingStrategy<T> SamplingStrategy<T>::CLAMP = {
    SamplingStrategyType::CLAMP_T, T()
};

template<typename T>
const SamplingStrategy<T> SamplingStrategy<T>::PASS_THROUGH = {
    SamplingStrategyType::PASS_THROUGH_T, T()
};

template<typename T>
const SamplingStrategy<T> SamplingStrategy<T>::ASSERT_FALSE = {
    SamplingStrategyType::ASSERT_FALSE_T, T()
};

template<typename T>
SamplingStrategy<T> SamplingStrategy<T>::SET(T outsideVolumeValue) {
    return {SamplingStrategyType::SET_T, outsideVolumeValue};
}

template<typename T>
template<typename O>
SamplingStrategy<T>::SamplingStrategy(const SamplingStrategy<O>& other)
    : type_(other.type_)
    , outsideVolumeValue_(T(other.outsideVolumeValue_))
{
}

template<typename T>
template<typename N>
SamplingStrategy<N> SamplingStrategy<T>::convert(std::function<N(T)> convertOutsideVolumeValue) const {
    return SamplingStrategy<N>(
            type_,
            convertOutsideVolumeValue(outsideVolumeValue_)
            );
}

} // namespace voreen

#endif // VRN_VOLUMEFILTER_H
