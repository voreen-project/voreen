/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "binarymedianfilter.h"
#include "slicereader.h"

#include <functional>

namespace voreen {

BinaryMedianFilter::BinaryMedianFilter(const tgt::ivec3& extent, float binarizationThreshold, uint32_t objectVoxelThreshold, const SamplingStrategy<float>& samplingStrategy)
    : extent_(extent)
    , kernelDimensions_(2*extent+tgt::ivec3::one)
    , binarizationThreshold_(binarizationThreshold)
    , objectVoxelThreshold_(objectVoxelThreshold)
    , samplingStrategy_(samplingStrategy)
{
}

BinaryMedianFilter::BinaryMedianFilter(const tgt::ivec3& extent, float binarizationThreshold, const SamplingStrategy<float>& samplingStrategy)
    : BinaryMedianFilter(extent, binarizationThreshold, tgt::hmul(2*extent+tgt::ivec3::one)/2, samplingStrategy)
{
}
const std::string BinaryMedianFilter::SLICE_BASE_TYPE = "uint8";

BinaryMedianFilter::~BinaryMedianFilter() {
}

std::unique_ptr<VolumeRAM> BinaryMedianFilter::getFilteredSlice(const CachingSliceReader* src, int z) const {
    typedef SimpleSlice<uint32_t> TempSlice;
    tgtAssert(z >= 0 && z<src->getSignedDimensions().z, "Invalid z pos in slice request");

    const tgt::ivec3& dim = src->getSignedDimensions();
    TempSlice zOutput(dim.xy());

    SamplingStrategy<float>::Sampler getValueFromReader = [src] (const tgt::ivec3& p) {
        return src->getVoxelNormalized(p);
    };

    // z
    #pragma omp parallel for
    for(int y = 0; y < dim.y; ++y) {
        for(int x = 0; x < dim.x; ++x) {
            uint32_t accumulator = 0;
            for(int dz = -extent_.z; dz <= extent_.z; ++dz) {
                // Why does template parameter deduction fail here? We should be able to call .sample(...)
                if(isObjectValue(samplingStrategy_.sample(tgt::ivec3(x, y, z+dz), dim, getValueFromReader))) {
                    ++accumulator;
                }
            }
            zOutput.at(x, y) = accumulator;
        }
    }

    TempSlice yOutput(dim.xy());
    SamplingStrategy<uint32_t> ySamplingStrategy = samplingStrategy_.convert<uint32_t>([this] (float v) {
                if(isObjectValue(v)) {
                    return kernelDimensions_.z; // Outside of the volume we need to count all object voxels in the z-row
                } else {
                    return 0;
                }
            });

    // y
    #pragma omp parallel for
    for(int x = 0; x < dim.x; ++x) {
        std::vector<uint32_t> lastValueBuffer(kernelDimensions_.y);

        // Calculate number of ones at y = 0 by reading and accumulating along +-extent_.y
        // Also save those values in lastValueBuffer
        uint32_t accumulator = 0;
        for(int dy = -extent_.y; dy <= extent_.y; ++dy) {
            uint32_t val = ySamplingStrategy.sample(tgt::ivec3(x, dy, 0), dim /* wrong in z, but doesn't matter */, zOutput.toSampler());
            accumulator += val;
            lastValueBuffer.at((dy+kernelDimensions_.y)%kernelDimensions_.y) = val;
        }
        yOutput.at(x, 0) = accumulator;

        // Now run along y, subtract the value running out of the window and add the one running into the window
        // Also add the new value to the window
        for(int y = 1; y < dim.y; ++y) {
            int readPos = y + extent_.y;
            accumulator -= lastValueBuffer.at(readPos%kernelDimensions_.y);
            uint32_t val = ySamplingStrategy.sample(tgt::ivec3(x, readPos, 0), dim /* wrong in z, but doesn't matter */, zOutput.toSampler());
            accumulator += val;
            lastValueBuffer.at(readPos%kernelDimensions_.y) = val;
            yOutput.at(x, y) = accumulator;
        }
    }

    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create(SLICE_BASE_TYPE, tgt::svec3(dim.xy(), 1)));
    SamplingStrategy<uint32_t> xSamplingStrategy = samplingStrategy_.convert<uint32_t>([this] (float v) {
                if(isObjectValue(v)) {
                    return kernelDimensions_.z*kernelDimensions_.y; // Outside of the volume we need to count all object voxels in the yz-slice
                } else {
                    return 0;
                }
            });

    // x
    #pragma omp parallel for
    for(int y = 0; y < dim.y; ++y) {
        std::vector<float> lastValueBuffer(kernelDimensions_.x);

        // Calculate number of ones at x = 0 by reading and accumulating along +-extent_.x
        // Also save those values in lastValueBuffer
        uint32_t accumulator = 0;
        {
            for(int dx = -extent_.x; dx <= extent_.x; ++dx) {
                uint32_t val = xSamplingStrategy.sample(tgt::ivec3(dx, y, 0), dim /* wrong in z, but doesn't matter */, yOutput.toSampler());
                accumulator += val;
                lastValueBuffer.at((dx+kernelDimensions_.x)%kernelDimensions_.x) = val;
            }
            float binaryVal = accumulator > objectVoxelThreshold_ ? 1.0f: 0.0f;
            outputSlice->setVoxelNormalized(binaryVal, tgt::svec3(0,y,0));
        }

        // Now run along x, subtract the value running out of the window and add the one running into the window
        // Also add the new value to the window
        for(int x = 1; x < dim.x; ++x) {
            int readPos = x + extent_.x;
            accumulator -= lastValueBuffer.at(readPos%kernelDimensions_.x);
            uint32_t val = xSamplingStrategy.sample(tgt::ivec3(readPos, y, 0), dim /* wrong in z, but doesn't matter */, yOutput.toSampler());
            accumulator += val;
            lastValueBuffer.at(readPos%kernelDimensions_.x) = val;
            float binaryVal = accumulator > objectVoxelThreshold_ ? 1.0f: 0.0f;
            outputSlice->setVoxelNormalized(binaryVal, tgt::svec3(x,y,0));
        }
    }
    return outputSlice;
}

int BinaryMedianFilter::zExtent() const {
    return extent_.z;
}

bool BinaryMedianFilter::isObjectValue(float input) const {
    return input > binarizationThreshold_;
}

SliceReaderMetaData BinaryMedianFilter::getMetaData(const SliceReaderMetaData& base) const {
    auto md = SliceReaderMetaData::fromBase(base);
    md.setRealWorldMapping(RealWorldMapping(tgt::vec2(0.0, 1.0), ""));

    float t = base.getRealworldMapping().normalizedToRealWorld(binarizationThreshold_);
    float min = 0.0f;
    float max = 1.0f;
    if(base.getMinMax()) {
        const auto& mm = *base.getMinMax();
        tgtAssert(base.getNumChannels() == 1, "Invalid number of channels");
        if(mm[0].x >= t) {
            min = 1.0f;
        } else if(mm[0].y < t) {
            max = 0.0f;
        }
        md.setMinMax({tgt::vec2(min, max)});
    } else if(base.getMinMaxBounds()) {
        const auto& mm = *base.getMinMaxBounds();
        tgtAssert(base.getNumChannels() == 1, "Invalid number of channels");
        if(mm[0].x >= t) {
            min = 1.0f;
        } else if(mm[0].y < t) {
            max = 0.0f;
        }
        md.setMinMaxBounds({tgt::vec2(min, max)});
    }

    md.setBaseType(SLICE_BASE_TYPE);

    return md;
}

} // namespace voreen
