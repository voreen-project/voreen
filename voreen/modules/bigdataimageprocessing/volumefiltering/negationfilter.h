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

#ifndef VRN_NEGATIONFILTER_H
#define VRN_NEGATIONFILTER_H

#include "parallelvolumefilter.h"

namespace voreen {

template<typename  T>
class NegationFilter : public ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>> {
public:
    NegationFilter(const tgt::bvec4& negate, const std::string& sliceBaseType);
    virtual ~NegationFilter();
    ParallelFilterValue<T> getValue(const typename NegationFilter<T>::Sample& sample, const tgt::ivec3& pos) const;
private:
    const tgt::bvec4 negate_; // TODO: We should derive the number of channels at compile time.
};

template<typename T>
NegationFilter<T>::NegationFilter(const tgt::bvec4& negate, const std::string& sliceBaseType)
    : ParallelVolumeFilter<ParallelFilterValue<T>, ParallelFilterValue<T>>(0, SamplingStrategy<ParallelFilterValue<T>>::ASSERT_FALSE, sliceBaseType)
    , negate_(negate)
{
}

template<typename T>
NegationFilter<T>::~NegationFilter() {
}

template<typename T>
ParallelFilterValue<T> NegationFilter<T>::getValue(const typename NegationFilter<T>::Sample& sample, const tgt::ivec3& pos) const {

    ParallelFilterValue<T> value = sample(pos);

    for(size_t channel=0; channel<ParallelFilterValue<T>::dim; channel++) {
        if(negate_[channel]) {
            value[channel] = -value[channel];
        }
    }

    return value;
}

typedef NegationFilter<float>     NegationFilter1D;
typedef NegationFilter<tgt::vec2> NegationFilter2D;
typedef NegationFilter<tgt::vec3> NegationFilter3D;
typedef NegationFilter<tgt::vec4> NegationFilter4D;

} // namespace voreen

#endif // VRN_MEDIANFILTER_H
