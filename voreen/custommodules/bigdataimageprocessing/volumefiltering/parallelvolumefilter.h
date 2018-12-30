/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_PARALLELVOLUMEFILTER_H
#define VRN_PARALLELVOLUMEFILTER_H

#include "volumefilter.h"
#include "slicereader.h"


namespace voreen {

template<class T>
struct ParallelFilterValue {
    const static size_t dim;

    // constructor with uninitialized value
    ParallelFilterValue<T>();

    // implicit conversions to wrapped type (via constructor)
    template<class U>
    ParallelFilterValue<T>(U val);
    operator T();

    float& operator[](size_t channel);

    T val_;
};

typedef ParallelFilterValue<float> ParallelFilterValue1D;
typedef ParallelFilterValue<tgt::vec2> ParallelFilterValue2D;
typedef ParallelFilterValue<tgt::vec3> ParallelFilterValue3D;
typedef ParallelFilterValue<tgt::vec4> ParallelFilterValue4D;

// Note. InputType and OutputType need to offer the interface of ParallelFilterValue
template<typename InputType, typename OutputType>
class ParallelVolumeFilter : public VolumeFilter {
public:

    ParallelVolumeFilter(int zExtent, const SamplingStrategy<InputType>& samplingStrategy, const std::string& sliceBaseType);
    std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const;
    const std::string& getSliceBaseType() const;
    int zExtent() const;

    typedef std::function<InputType(const tgt::ivec3& pos)> Sample;
    virtual OutputType getValue(const Sample& sample, const tgt::ivec3& pos) const = 0;

    size_t getNumInputChannels() const;
    size_t getNumOutputChannels() const;

private:
    const SamplingStrategy<InputType> samplingStrategy_;
    const std::string sliceBaseType_;
    const int zExtent_;
};


// ------------------------------------------------------------------------------------------
// Implementation ---------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------

// ParallelVolumeFilter ---------------------------------------------------------------------------

template<typename InputType, typename OutputType>
ParallelVolumeFilter<InputType, OutputType>::ParallelVolumeFilter(int zExtent, const SamplingStrategy<InputType>& samplingStrategy, const std::string& sliceBaseType)
    : samplingStrategy_(samplingStrategy)
    , sliceBaseType_(sliceBaseType)
    , zExtent_(zExtent)
{
}

template<typename InputType, typename OutputType>
std::unique_ptr<VolumeRAM> ParallelVolumeFilter<InputType, OutputType>::getFilteredSlice(const CachingSliceReader* src, int z) const {
    tgtAssert(src->getZExtent() == zExtent_, "z extent mismatch");

    const tgt::ivec3& srcDim = src->getSignedDimensions();

    std::string format = (OutputType::dim == 1) ?
        sliceBaseType_ :
        "Vector" + std::to_string(OutputType::dim) + "(" + sliceBaseType_ + ")"; // only SLIGHTLY hacky

    std::unique_ptr<VolumeRAM> slice(VolumeFactory().create(format, tgt::svec3(srcDim.xy(), 1)));

    std::function<InputType(const tgt::ivec3& p)> getValueFromReader = [src] (const tgt::ivec3& p) {
        InputType out;
        for(size_t d=0; d < InputType::dim; ++d) {
            out[d] = src->getVoxelNormalized(p, d);
        }
        return out;
    };

    ParallelVolumeFilter<InputType, OutputType>::Sample sample = [this, srcDim, getValueFromReader] (const tgt::ivec3& pos) {
        return samplingStrategy_.sample(pos, srcDim, getValueFromReader);
    };

    #pragma omp parallel for
    for(int y = 0; y < srcDim.y; ++y) {
        for(int x = 0; x < srcDim.x; ++x) {
            OutputType val = getValue(sample, tgt::ivec3(x, y, z));
            for(size_t d=0; d < OutputType::dim; ++d) {
                slice->setVoxelNormalized(val[d], tgt::svec3(x,y,0), d);
            }
        }
    }
    return slice;
}

template<typename InputType, typename OutputType>
const std::string& ParallelVolumeFilter<InputType, OutputType>::getSliceBaseType() const {
    return sliceBaseType_;
}

template<typename InputType, typename OutputType>
int ParallelVolumeFilter<InputType, OutputType>::zExtent() const {
    return zExtent_;
}

template<typename InputType, typename OutputType>
size_t ParallelVolumeFilter<InputType, OutputType>::getNumInputChannels() const {
    return InputType::dim;
}

template<typename InputType, typename OutputType>
size_t ParallelVolumeFilter<InputType, OutputType>::getNumOutputChannels() const {
    return OutputType::dim;
}


// ParallelFilterValue -----------------------------------------------------------------------

template<typename T>
const size_t ParallelFilterValue<T>::dim = T::size;

// VisualStudio does not conform with the c++ standard and wants the definition of static
// member specializations as part of the declaration somehow, even though the standard says:
//
// Every program shall contain exactly one definition of every non-inline function or object that is used in that program; no diagnostic required.
//
// (See http://stackoverflow.com/questions/2342550/static-member-initialization-for-specialized-template-class)
#ifdef WIN32
template<>
const size_t ParallelFilterValue<float>::dim = 1;
#else
template<>
const size_t ParallelFilterValue<float>::dim;
#endif

template<typename T>
inline
ParallelFilterValue<T>::operator T() {
    return val_;
}

template<typename T>
inline
ParallelFilterValue<T>::ParallelFilterValue()
    : ParallelFilterValue(0)
{
}

template<typename T>
template<typename U>
inline
ParallelFilterValue<T>::ParallelFilterValue(U val)
    : val_(val)
{
}

template<typename T>
inline
float& ParallelFilterValue<T>::operator[](size_t channel) {
    //tgtAssert(channel < ParallelFilterValue<T>::dim, "invalid channel");
    return val_[channel];
}

template<>
inline
float& ParallelFilterValue<float>::operator[](size_t) {
    //tgtAssert(channel < ParallelFilterValue<float>::dim, "invalid channel");
    return val_;
}

} // namespace voreen
#endif // VRN_PARALLELVOLUMEFILTER_H

