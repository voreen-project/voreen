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

#include "morphologyfilter.h"

#include "slicereader.h"

namespace {
    static const auto DILATION_FUNC = [](float a, float b) { return std::max(a, b); };
    static const auto EROSION_FUNC  = [](float a, float b) { return std::min(a, b); };
}

namespace voreen {

MorphologyFilter::MorphologyFilter(const tgt::ivec3& extent, MorphologyOperatorType type, MorphologyOperatorShape shape, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType)
    : extent_(extent)
    , type_(type)
    , shape_(shape)
    , samplingStrategy_(samplingStrategy)
    , sliceBaseType_(sliceBaseType)
{
    switch (type_) {
    case DILATION_T:
        morphFunc_ = DILATION_FUNC;
        break;
    case EROSION_T:
        morphFunc_ = EROSION_FUNC;
        break;
    default:
        tgtAssert(false, "Unimplemented morphology type");
    }
}

int MorphologyFilter::zExtent() const {
    return extent_.z;
}

const std::string& MorphologyFilter::getSliceBaseType() const {
    return sliceBaseType_;
}

std::unique_ptr<VolumeRAM> MorphologyFilter::getFilteredSlice(const CachingSliceReader* src, int z) const {
    tgtAssert(z >= 0 && z < src->getSignedDimensions().z, "Invalid z pos in slice request");

    switch (shape_) {
    case CUBE_T:
        return getFilteredSliceCubeMorphology(src, z);
    case SPHERE_T:
        return getFilteredSliceSphereMorphology(src, z);
    default:
        tgtAssert(false, "Unimplemented morphology shape");
    }
    
    return nullptr;
}

MorphologyOperatorType MorphologyFilter::getMorphologyOperatorType() const {
    return type_;
}
MorphologyOperatorShape MorphologyFilter::getMorphologyOperatorShape() const {
    return shape_;
}

std::unique_ptr<VolumeRAM> MorphologyFilter::getFilteredSliceCubeMorphology(const CachingSliceReader* src, int z) const {

    const tgt::ivec3& dim = src->getSignedDimensions();
    tgt::ivec3 halfKernelDim = extent_ / 2;
    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create(sliceBaseType_, tgt::svec3(dim.xy(), 1)));
    std::unique_ptr<VolumeRAM> srcSlice(VolumeFactory().create(sliceBaseType_, tgt::svec3(dim.xy(), 1)));

    SamplingStrategy<float>::Sampler getValueFromReader = [src](const tgt::ivec3& p) {
        return src->getVoxelNormalized(p);
    };

    // z
    #pragma omp parallel for
    for (int y = 0; y < dim.y; ++y) {
        for (int x = 0; x < dim.x; ++x) {
            float value = samplingStrategy_.sample(tgt::ivec3(x, y, z), dim, getValueFromReader);
            for (int dz = -halfKernelDim.z; dz <= halfKernelDim.z; ++dz) {
               value = morphFunc_(value, samplingStrategy_.sample(tgt::ivec3(x, y, z + dz), dim, getValueFromReader));
            }
            outputSlice->setVoxelNormalized(value, tgt::svec3(x, y, 0));
        }
    }

    std::swap(outputSlice, srcSlice);
    const VolumeRAM* srcSlicePtr = srcSlice.get();

    SamplingStrategy<float>::Sampler getValueFromSrcSlice = [&srcSlicePtr](const tgt::ivec3& p) {
        return srcSlicePtr->getVoxelNormalized(tgt::svec3(p));
    };

    // y
    #pragma omp parallel for
    for (int y = 0; y < dim.y; ++y) {
        for (int x = 0; x < dim.x; ++x) {
            float value = samplingStrategy_.sample(tgt::ivec3(x, y, 0), dim, getValueFromSrcSlice);
            for (int dy = -halfKernelDim.y; dy <= halfKernelDim.y; ++dy) {
                value = morphFunc_(value, samplingStrategy_.sample(tgt::ivec3(x, y + dy, 0), dim /* wrong in z, but doesn't matter */, getValueFromSrcSlice));
            }
            outputSlice->setVoxelNormalized(value, tgt::svec3(x, y, 0));
        }
    }

    std::swap(outputSlice, srcSlice);
    srcSlicePtr = srcSlice.get();

    // x
    #pragma omp parallel for
    for (int y = 0; y < dim.y; ++y) {
        for (int x = 0; x < dim.x; ++x) {
            float value = samplingStrategy_.sample(tgt::ivec3(x, y, 0), dim, getValueFromSrcSlice);
            for (int dx = -halfKernelDim.x; dx <= halfKernelDim.x; ++dx) {
                value = morphFunc_(value, samplingStrategy_.sample(tgt::ivec3(x+dx, y, 0), dim, getValueFromSrcSlice));
            }
            outputSlice->setVoxelNormalized(value, tgt::svec3(x, y, 0));
        }
    }
    return outputSlice;
}

std::unique_ptr<VolumeRAM> MorphologyFilter::getFilteredSliceSphereMorphology(const CachingSliceReader* src, int z) const {
    
    const tgt::ivec3& dim = src->getSignedDimensions();
    tgt::ivec3 halfKernelDim = extent_ / 2;
    tgt::vec3 kernelRadiusSq = halfKernelDim*halfKernelDim;
    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create(sliceBaseType_, tgt::svec3(dim.xy(), 1)));

    SamplingStrategy<float>::Sampler getValueFromReader = [src](const tgt::ivec3& p) {
        return src->getVoxelNormalized(p);
    };

    #pragma omp parallel for
    for (int y = 0; y < dim.y; ++y) {
        for (int x = 0; x < dim.x; ++x) {
            float value = samplingStrategy_.sample(tgt::ivec3(x, y, z), dim, getValueFromReader);
            for (int nz = z-halfKernelDim.z; nz <= z+halfKernelDim.z; nz++) {
                for (int ny = y - halfKernelDim.z; ny <= y + halfKernelDim.y; ny++) {
                    for (int nx = x - halfKernelDim.x; nx <= x + halfKernelDim.x; nx++) {

                        float dx = nx - x;
                        float dy = ny - y;
                        float dz = nz - z;

                        float d = dx*dx / kernelRadiusSq.x + dy*dy / kernelRadiusSq.y + dz*dz / kernelRadiusSq.z;
                        if(d > 1.0f) {
                            continue;
                        }

                        value = morphFunc_(value, samplingStrategy_.sample(tgt::ivec3(nx, ny, nz), dim, getValueFromReader));
                    }
                }
            }
            outputSlice->setVoxelNormalized(value, tgt::svec3(x, y, 0));
        }
    }

    return outputSlice;
}


} // namespace voreen
