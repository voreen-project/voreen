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

#include "binarizationfilter.h"

namespace voreen {

static const std::string BINARIZATION_SLICE_BASE_TYPE = "uint8";

BinarizationFilter::BinarizationFilter(float threshold)
    : ParallelVolumeFilter<ParallelFilterValue1D, ParallelFilterValue1D>(0, SamplingStrategy<ParallelFilterValue1D>::ASSERT_FALSE, BINARIZATION_SLICE_BASE_TYPE)
    , threshold_(threshold)
{
}

BinarizationFilter::~BinarizationFilter() {
}

ParallelFilterValue1D BinarizationFilter::getValue(const Sample& sample, const tgt::ivec3& pos, const SliceReaderMetaData& inputMetadata, const SliceReaderMetaData& outputMetaData) const {
    return sample(pos) < threshold_ ? 0.0f : 1.0f;
}

SliceReaderMetaData BinarizationFilter::getMetaData(const SliceReaderMetaData& base) const {
    auto md = SliceReaderMetaData::fromBase(base);
    md.setRealWorldMapping(RealWorldMapping(tgt::vec2(0.0f, 1.0f), ""));

    if(base.getMinMax()) {
        const auto& mm = *base.getMinMax();
        float t = base.getRealworldMapping().normalizedToRealWorld(threshold_);
        float min = 0.0f;
        float max = 1.0f;
        tgtAssert(base.getNumChannels() == 1, "Invalid number of channels");
        if(mm[0].x >= t) {
            min = 1.0f;
        } else if(mm[0].y < t) {
            max = 0.0f;
        }

        md.setMinMax({tgt::vec2(min, max)});
    }

    return md;
}

} // namespace voreen
