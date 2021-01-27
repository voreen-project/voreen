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

#include "morphologyfilter.h"

#include "slicereader.h"

namespace {
    static const auto DILATION_FUNC = [](float a, float b) { return std::max(a, b); };
    static const auto EROSION_FUNC  = [](float a, float b) { return std::min(a, b); };
}

namespace voreen {

MorphologyFilter::MorphologyFilter(const tgt::ivec3& extent, MorphologyOperatorType type, MorphologyOperatorShape shape, const SamplingStrategy<float>& samplingStrategy)
    : extent_(extent)
    , type_(type)
    , shape_(shape)
    , samplingStrategy_(samplingStrategy)
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
    std::string type = src->getMetaData().getBaseType();
    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create(type, tgt::svec3(dim.xy(), 1)));
    std::unique_ptr<VolumeRAM> srcSlice(VolumeFactory().create(type, tgt::svec3(dim.xy(), 1)));

    SamplingStrategy<float>::Sampler getValueFromReader = [src](const tgt::ivec3& p) {
        return src->getVoxelNormalized(p);
    };

    // z
    #pragma omp parallel for
    for (int y = 0; y < dim.y; ++y) {
        for (int x = 0; x < dim.x; ++x) {
            float value = samplingStrategy_.sample(tgt::ivec3(x, y, z), dim, getValueFromReader);
            for (int dz = -extent_.z; dz <= extent_.z; ++dz) {
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
            for (int dy = -extent_.y; dy <= extent_.y; ++dy) {
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
            for (int dx = -extent_.x; dx <= extent_.x; ++dx) {
                value = morphFunc_(value, samplingStrategy_.sample(tgt::ivec3(x+dx, y, 0), dim, getValueFromSrcSlice));
            }
            outputSlice->setVoxelNormalized(value, tgt::svec3(x, y, 0));
        }
    }
    return outputSlice;
}

std::unique_ptr<VolumeRAM> MorphologyFilter::getFilteredSliceSphereMorphology(const CachingSliceReader* src, int z) const {

    const tgt::ivec3& dim = src->getSignedDimensions();
    tgt::vec3 extentSq = extent_*extent_;
    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create(src->getMetaData().getBaseType(), tgt::svec3(dim.xy(), 1)));

    SamplingStrategy<float>::Sampler getValueFromReader = [src](const tgt::ivec3& p) {
        return src->getVoxelNormalized(p);
    };

    #pragma omp parallel for
    for (int y = 0; y < dim.y; ++y) {
        for (int x = 0; x < dim.x; ++x) {
            float value = samplingStrategy_.sample(tgt::ivec3(x, y, z), dim, getValueFromReader);
            for (int dz = -extent_.z; dz <= extentSq.z; ++dz) {
                float dz2rel = extent_.z ? dz*dz/extent_.z : 0.0f;
                for (int dy = -extent_.y; dy <= extent_.y; ++dy) {
                    float dy2rel = extent_.y ? dy*dy/extentSq.y : 0.0f;
                    for (int dx = -extent_.x; dx <= extent_.x; ++dx) {
                        float dx2rel = extent_.x ? dx*dx/extentSq.x : 0.0f;

                        float d = dx2rel + dy2rel + dz2rel;
                        if(d > 1.0f) {
                            continue;
                        }

                        value = morphFunc_(value, samplingStrategy_.sample(tgt::ivec3(x+dx, y+dy, z+dz), dim, getValueFromReader));
                    }
                }
            }
            outputSlice->setVoxelNormalized(value, tgt::svec3(x, y, 0));
        }
    }

    return outputSlice;
}

SliceReaderMetaData MorphologyFilter::getMetaData(const SliceReaderMetaData& base) const {
    auto md = SliceReaderMetaData::fromBase(base);
    if(base.getMinMaxBounds()) {
        // We can only do this because the median filter never expands the
        // range of possible values compared to the input.
        md.setMinMaxBounds(*base.getMinMaxBounds());
    }
    return md;
}


} // namespace voreen
