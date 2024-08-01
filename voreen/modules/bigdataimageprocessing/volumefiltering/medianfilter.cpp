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

#include "medianfilter.h"

namespace voreen {

MedianFilter::MedianFilter(const tgt::ivec3& extent, const SamplingStrategy<ParallelFilterValue1D>& samplingStrategy)
    : ParallelVolumeFilter<ParallelFilterValue1D, ParallelFilterValue1D>(extent.z, samplingStrategy)
    , extent_(extent)
{
}

MedianFilter::~MedianFilter() {
}

ParallelFilterValue1D MedianFilter::getValue(const Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetadata, const SliceReaderMetaData& outputMetaData) const {
    std::vector<float> values;
    values.reserve(tgt::hmul(extent_));

    for(int z = pos.z-extent_.z; z <= pos.z+extent_.z; ++z) {
        for(int y = pos.y-extent_.y; y <= pos.y+extent_.y; ++y) {
            for(int x = pos.x-extent_.x; x <= pos.x+extent_.x; ++x) {
                values.push_back(sample(tgt::ivec3(x,y,z)));
            }
        }
    }
    std::nth_element(values.begin(), values.begin() + values.size()/2, values.end());
    return values[values.size()/2];
}

SliceReaderMetaData MedianFilter::getMetaData(const SliceReaderMetaData& base) const {
    auto md = SliceReaderMetaData::fromBase(base);
    if(base.getMinMaxBounds()) {
        // We can only do this because the median filter never expands the
        // range of possible values compared to the input.
        md.setMinMaxBounds(*base.getMinMaxBounds());
    }
    return md;
}

} // namespace voreen
