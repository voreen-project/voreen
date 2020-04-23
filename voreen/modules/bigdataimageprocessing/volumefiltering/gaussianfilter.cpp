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

#include "gaussianfilter.h"
#include "slicereader.h"

namespace voreen {

tgt::ivec3 GaussianFilter::suitableExtent(const tgt::vec3& standardDeviation) {
    return tgt::ivec3(suitableExtent(standardDeviation.x),
                      suitableExtent(standardDeviation.y),
                      suitableExtent(standardDeviation.z));
}

float GaussianFilter::suitableStandardDeviation(int extent) {
    tgtAssert(extent >= 0, "invalid extent");
    // Chosen similar to VolumeFiltering
    return (extent+0.5f)/2.5f;
}

tgt::vec3 GaussianFilter::suitableStandardDeviation(const tgt::ivec3& extent) {
    return tgt::vec3(suitableStandardDeviation(extent.x),
                     suitableStandardDeviation(extent.y),
                     suitableStandardDeviation(extent.z));
}

int GaussianFilter::suitableExtent(float standardDeviation) {
    tgtAssert(standardDeviation >= 0, "invalid standardDeviation");
    // Chosen similar to VolumeFiltering
    return static_cast<int>(2.5*standardDeviation-0.5f);
}

GaussianFilter::GaussianFilter(float standardDeviation, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels)
    : GaussianFilter(tgt::vec3(standardDeviation), samplingStrategy, sliceBaseType, numChannels)
{
}

GaussianFilter::GaussianFilter(const tgt::vec3& standardDeviation, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels)
    : GaussianFilter(standardDeviation, GaussianFilter::suitableExtent(standardDeviation), samplingStrategy, sliceBaseType, numChannels)
{
}

GaussianFilter::GaussianFilter(int extent, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels)
    : GaussianFilter(tgt::ivec3(extent), samplingStrategy, sliceBaseType, numChannels)
{
}

GaussianFilter::GaussianFilter(const tgt::ivec3& extent, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels)
    : GaussianFilter(suitableStandardDeviation(extent), extent, samplingStrategy, sliceBaseType, numChannels)
{
}

static float* initHalfKernel(int extent, float standardDeviation) {
    tgtAssert(extent >= 0, "invalid extent");
    tgtAssert(standardDeviation > 0, "invalid standardDeviation");

    float* halfKernel = new float[extent + 1];
    std::function<float(int)> kernelFunc = [standardDeviation] (int x) { return std::exp(-0.5f * x*x / (standardDeviation*standardDeviation));};

    // find initial values
    halfKernel[0] = kernelFunc(0);
    float sum = halfKernel[0];
    for(int i=1; i < extent+1; ++i) {
        float val = kernelFunc(i);
        halfKernel[i] = val;
        sum += 2*val; //symmetric kernel => count for both sides
    }

    // normalize
    for(int i=0; i < extent+1; ++i) {
        halfKernel[i] /= sum;
    }
    return halfKernel;
}

GaussianFilter::GaussianFilter(float standardDeviation, int extent, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels)
    : GaussianFilter(tgt::vec3(standardDeviation), tgt::ivec3(extent), samplingStrategy, sliceBaseType, numChannels)
{
}

GaussianFilter::GaussianFilter(const tgt::vec3& standardDeviation, const tgt::ivec3& extent, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels)
    : neighborhoodDimensions_(extent)
    , kernelDimensions_(2*extent+tgt::ivec3::one)
    , halfKernelX_(initHalfKernel(extent.x, standardDeviation.x))
    , halfKernelY_(initHalfKernel(extent.y, standardDeviation.y))
    , halfKernelZ_(initHalfKernel(extent.z, standardDeviation.z))
    , samplingStrategy_(samplingStrategy)
    , sliceBaseType_(sliceBaseType)
    , numChannels_(numChannels)
{
    tgtAssert(tgt::hand(tgt::greaterThan(standardDeviation, tgt::vec3::zero)), "invalid standardDeviation");
    tgtAssert(tgt::hand(tgt::greaterThanEqual(extent, tgt::ivec3::zero)), "invalid extent");
}

GaussianFilter::~GaussianFilter()
{
    delete[] halfKernelX_;
    delete[] halfKernelY_;
    delete[] halfKernelZ_;
}

float GaussianFilter::getKernelValX(int centeredPos) const {
    return halfKernelX_[std::abs(centeredPos)];
}

float GaussianFilter::getKernelValY(int centeredPos) const {
    return halfKernelY_[std::abs(centeredPos)];
}

float GaussianFilter::getKernelValZ(int centeredPos) const {
    return halfKernelZ_[std::abs(centeredPos)];
}

std::unique_ptr<VolumeRAM> GaussianFilter::getFilteredSlice(const CachingSliceReader* src, int z) const {
    tgtAssert(z >= 0 && z<src->getSignedDimensions().z, "Invalid z pos in slice request");

    const tgt::ivec3& dim = src->getSignedDimensions();
    VolumeFactory volumeFactory;
    std::string format = volumeFactory.getFormat(sliceBaseType_, numChannels_);
    std::unique_ptr<VolumeRAM> outputSlice(volumeFactory.create(format, tgt::svec3(dim.xy(), 1)));
    std::unique_ptr<VolumeRAM> srcSlice(volumeFactory.create(format, tgt::svec3(dim.xy(), 1)));

    // TODO: Memory access-wise, iterating the channels first is a quite inefficient operation.
    //  Either use templates for this class or make this a nested loop inside those iterating voxels.
    //  However, for the single-channel-case being implemented before, this means practically no difference.
    for(size_t channel = 0; channel < numChannels_; channel++) {

        SamplingStrategy<float>::Sampler getValueFromReader = [src, channel](const tgt::ivec3 &p) {
            return src->getVoxelNormalized(p, channel);
        };

        // z
        #pragma omp parallel for
        for (int y = 0; y < dim.y; ++y) {
            for (int x = 0; x < dim.x; ++x) {
                float accumulator = 0;
                for (int dz = -neighborhoodDimensions_.z; dz <= neighborhoodDimensions_.z; ++dz) {
                    // Why does template parameter deduction fail here? We should be able to call .sample(...)
                    accumulator += getKernelValZ(dz) *
                                   samplingStrategy_.sample(tgt::ivec3(x, y, z + dz),
                                                            dim, getValueFromReader);
                }
                outputSlice->setVoxelNormalized(accumulator, tgt::svec3(x, y, 0), channel);
            }
        }

        std::swap(outputSlice, srcSlice);
        const VolumeRAM *srcSlicePtr = srcSlice.get();

        SamplingStrategy<float>::Sampler getValueFromSrcSlice = [&srcSlicePtr, channel](const tgt::ivec3 &p) {
            return srcSlicePtr->getVoxelNormalized(tgt::svec3(p), channel);
        };

        // y
        #pragma omp parallel for
        for (int y = 0; y < dim.y; ++y) {
            for (int x = 0; x < dim.x; ++x) {
                float accumulator = 0;
                for (int dy = -neighborhoodDimensions_.y; dy <= neighborhoodDimensions_.y; ++dy) {
                    accumulator += getKernelValY(dy) *
                                   samplingStrategy_.sample(tgt::ivec3(x, y + dy, 0),
                                                            dim /* wrong in z, but doesn't matter */,
                                                            getValueFromSrcSlice);
                }
                outputSlice->setVoxelNormalized(accumulator, tgt::svec3(x, y, 0), channel);
            }
        }

        std::swap(outputSlice, srcSlice);
        srcSlicePtr = srcSlice.get();

        // x
        #pragma omp parallel for
        for (int y = 0; y < dim.y; ++y) {
            for (int x = 0; x < dim.x; ++x) {
                float accumulator = 0;
                for (int dx = -neighborhoodDimensions_.x; dx <= neighborhoodDimensions_.x; ++dx) {
                    accumulator += getKernelValX(dx) *
                                   samplingStrategy_.sample(tgt::ivec3(x + dx, y, 0),
                                                            dim /* wrong in z, but doesn't matter */,
                                                            getValueFromSrcSlice);
                }
                outputSlice->setVoxelNormalized(accumulator, tgt::svec3(x, y, 0), channel);
            }
        }
    }

    return outputSlice;
}

int GaussianFilter::zExtent() const {
    return neighborhoodDimensions_.z;
}

const std::string& GaussianFilter::getSliceBaseType() const {
    return sliceBaseType_;
}

size_t GaussianFilter::getNumInputChannels() const {
    return numChannels_;
}

size_t GaussianFilter::getNumOutputChannels() const {
    return numChannels_;
}


} // namespace voreen
