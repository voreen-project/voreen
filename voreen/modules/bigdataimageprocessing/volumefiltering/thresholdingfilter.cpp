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

#include "thresholdingfilter.h"

#include "slicereader.h"

namespace {
    static const auto LOWER_THRESHOLD_FUNC = [](float a, float b) { return a < b; };
    static const auto UPPER_THRESHOLD_FUNC = [](float a, float b) { return a > b; };
}

namespace voreen {

static const std::string BINARIZATION_SLICE_BASE_TYPE = "uint8";
static const float BINARIZATION_OUT_VALUE = 0.0f;
static const float BINARIZATION_IN_VALUE = 1.0f;

ThresholdingFilter::ThresholdingFilter(float threshold, float replacement, ThresholdingStrategyType thresholdingStrategyType, const std::string& sliceBaseType)
    : binarize_(false)
    , threshold_(threshold)
    , replacement_(replacement)
    , thresholdingStrategyType_(thresholdingStrategyType)
    , sliceBaseType_(sliceBaseType)
{
    switch (thresholdingStrategyType_) {
    case LOWER_T:
        strategy_ = LOWER_THRESHOLD_FUNC;
        break;
    case UPPER_T:
        strategy_ = UPPER_THRESHOLD_FUNC;
        break;
    default:
        tgtAssert(false, "Unimplemented Thresholding Strategy")
            break;
    }
}

ThresholdingFilter::ThresholdingFilter(float threshold, ThresholdingStrategyType thresholdingStrategyType)
    : binarize_(true)
    , threshold_(threshold)
    , replacement_(BINARIZATION_OUT_VALUE)
    , thresholdingStrategyType_(thresholdingStrategyType)
    , sliceBaseType_(BINARIZATION_SLICE_BASE_TYPE)
{
    switch (thresholdingStrategyType_) {
    case LOWER_T:
        strategy_ = LOWER_THRESHOLD_FUNC;
        break;
    case UPPER_T:
        strategy_ = UPPER_THRESHOLD_FUNC;
        break;
    default:
        tgtAssert(false, "Unimplemented Thresholding Strategy")
            break;
    }
}

ThresholdingFilter::~ThresholdingFilter() {
}

int ThresholdingFilter::zExtent() const {
    return 1;
}

const std::string& ThresholdingFilter::getSliceBaseType() const {
    return sliceBaseType_;
}

std::unique_ptr<VolumeRAM> ThresholdingFilter::getFilteredSlice(const CachingSliceReader* src, int z) const {
    tgtAssert(z >= 0 && z < src->getSignedDimensions().z, "Invalid z pos in slice request");

    const tgt::ivec3& dim = src->getSignedDimensions();
    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create(sliceBaseType_, tgt::svec3(dim.xy(), 1)));

    #pragma omp parallel for
    for (int y = 0; y < dim.y; ++y) {
        for (int x = 0; x < dim.x; ++x) {
            float value = src->getVoxelNormalized(tgt::ivec3(x, y, z));
            bool replace = strategy_(value, threshold_);
            if (replace) {
                value = replacement_;
            }
            else if (binarize_) {
                value = BINARIZATION_IN_VALUE;
            }
            outputSlice->setVoxelNormalized(value, x, y, 0);
        }
    }
    return outputSlice;
}

ThresholdingStrategyType ThresholdingFilter::getThresholdingStrategyType() const {
    return thresholdingStrategyType_;
}

} // namespace voreen
