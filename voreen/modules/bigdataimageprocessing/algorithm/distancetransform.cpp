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

#include "distancetransform.h"

#include <vector>

#include "tgt/vector.h"

namespace voreen {

template<typename T>
static inline T square(T i) {
    return i*i;
}

template<int outerDim, int innerDim, typename InputSlice, typename OutputSlice, typename InitValFunc, typename FinalValFunc>
static void dt_slice_pass(InputSlice& inputSlice, OutputSlice& outputSlice, tgt::ivec3 dim, tgt::vec3 spacingVec, InitValFunc initValFunc, FinalValFunc finalValFunc) {
    const int n = dim[innerDim];
    std::vector<int> v(n+1, 0); // Locations of parabolas in lower envelope, voxel coordinates
    std::vector<float> z(n+1, 0); // Locations of boundaries between parabolas, physical coordinates (i.e., including spacing)

    float spacing = spacingVec[innerDim];
    tgtAssert(spacing > 0, "Invalid spacing");

    for(int x = 0; x < dim[outerDim]; ++x) {
        auto f = [&] (int i) {
            tgt::svec3 slicePos(i,i,0);
            slicePos[outerDim] = x;
            float g = inputSlice->voxel(slicePos);
            return initValFunc(g);
        };

        v[0] = 0;
        z[0] = -std::numeric_limits<float>::infinity();
        z[1] = std::numeric_limits<float>::infinity();
        int k = 0;

        for(int q = 1; q < n; ++q) {
            float fq = f(q);
            if(std::isinf(fq)) {
                continue;
            }
jmp:
            int vk = v[k];
            float qs = spacing * q;
            float vks = spacing * vk;
            float s = ((fq - f(vk)) + (square(qs) - square(vks))) / (2*(qs - vks)); //note: q > vk
            tgtAssert(!std::isnan(s), "s is nan");

            if(s <= z[k]) {
                if(k > 0) {
                    k -= 1;
                    goto jmp;
                } else {
                    v[k] = q;
                    z[k] = s;
                    z[k+1] = std::numeric_limits<float>::infinity();
                }
            } else {
                k += 1;
                v[k] = q;
                z[k] = s;
                z[k+1] = std::numeric_limits<float>::infinity();
            }
        }

        k = 0;
        for(int q = 0; q < n; ++q) {
            float qs = spacing * q;
            while(z[k+1] < qs) {
                k += 1;
            }
            tgt::svec3 slicePos(q,q,0);
            slicePos[outerDim] = x;
            int vk = v[k];
            float vks = spacing * vk;
            float val = f(vk) + square(qs-vks);
            outputSlice->voxel(slicePos) = finalValFunc(val);
        }
    }
}

LZ4SliceVolume<float> compute_distance_transform(const VolumeBase& vol, float binarizationThreshold, std::string outputPath, ProgressReporter& progressReporter) {
    const tgt::svec3 dim = vol.getDimensions();
    const tgt::svec3 sliceDim(dim.x, dim.y, 1);
    const tgt::vec3 spacing = vol.getSpacing();

    SubtaskProgressReporterCollection<2> tasks(progressReporter);

    float binarizationThresholdNormalized;
    if(vol.hasMetaData(VolumeBase::META_DATA_NAME_REAL_WORLD_MAPPING)) {
        // If the input volume does not have a RealWorldMapping we need to convert the binarizationThreshold to a normalized value.
        binarizationThresholdNormalized = vol.getRealWorldMapping().realWorldToNormalized(binarizationThreshold);
    } else {
        // If the input volume does not have a RealWorldMapping we expect RW values to be normalized.
        binarizationThresholdNormalized = binarizationThreshold;
    }

    auto isBackground = [&] (float value) {
        return value < binarizationThresholdNormalized;
    };

    LZ4SliceVolumeBuilder<float> gBuilder(outputPath,
            LZ4SliceVolumeMetadata(dim)
            .withOffset(vol.getOffset())
            .withSpacing(vol.getSpacing())
            .withPhysicalToWorldTransformation(vol.getPhysicalToWorldMatrix()));

    VolumeAtomic<float> gSlice(sliceDim);
    // z-scan 1: calculate distances in forward direction
    {
        const size_t z = 0;
        std::unique_ptr<VolumeRAM> inputSlice(vol.getSlice(z));
        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                tgt::svec3 slicePos(x,y,0);
                float val = inputSlice->getVoxelNormalized(slicePos);
                float& g = gSlice.voxel(slicePos);
                if(isBackground(val)) {
                    g = 0;
                } else {
                    g = std::numeric_limits<float>::infinity();
                }
            }
        }
        gBuilder.pushSlice(gSlice);
    }
    for(size_t z = 1; z < dim.z; ++z) {
        tasks.get<0>().setProgress(static_cast<float>(z)/dim.z);
        std::unique_ptr<VolumeRAM> inputSlice(vol.getSlice(z));


        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                tgt::svec3 slicePos(x,y,0);
                float val = inputSlice->getVoxelNormalized(slicePos);
                float& g = gSlice.voxel(slicePos);
                if(isBackground(val)) {
                    g = 0;
                } else {
                    g += spacing.z;
                }
            }
        }
        gBuilder.pushSlice(gSlice);
    }
    auto gvol = std::move(gBuilder).finalize();

    auto tmpSlice = VolumeAtomic<float>(sliceDim);
    auto tmpSlicePtr = &tmpSlice;

    // z-scan 2, propagate distances in other direction
    // also, directly do y- and x- passes on the slices while they are loaded.
    {
        auto prevZSlice = gvol.loadSlice(dim.z-1);
        for(int z = dim.z-2; z >= 0; --z) {
            tasks.get<1>().setProgress(static_cast<float>(z)/dim.z);

            auto gSlice = gvol.getWriteableSlice(z);
            for(size_t y = 0; y < dim.y; ++y) {
                for(size_t x = 0; x < dim.x; ++x) {
                    tgt::svec3 slicePos(x,y,0);
                    float& gPrev = prevZSlice.voxel(slicePos);
                    float& g = gSlice->voxel(slicePos);

                    float ng = gPrev + spacing.z;
                    if(ng < g) {
                        g = ng;
                    }
                }
            }
            prevZSlice = gSlice->copy();

            // Now do x and y passes on current slice to finalize it.
            dt_slice_pass<0,1>(gSlice, tmpSlicePtr, dim, spacing, [] (float v) {return square(v);}, [] (float v) {return v;});
            dt_slice_pass<1,0>(tmpSlicePtr, gSlice, dim, spacing, [] (float v) {return v;}, [] (float v) {return std::sqrt(v);});
        }
    }

    progressReporter.setProgress(1.f);

    return gvol;
}

static bool is_peak(float l, float c, float r) {
    return l < c && c >= r || l <= c && c > r;
    //return l < c && c >= r;
}

static void count_peaks_z(LZ4WriteableSlab<uint8_t>& counts, const VolumeAtomic<float>& prev, const VolumeAtomic<float>& current, const VolumeAtomic<float>& next) {
    const tgt::svec3 dim = counts->getDimensions();
    tgtAssert(dim == prev.getDimensions(), "Dim mismatch");
    tgtAssert(dim == current.getDimensions(), "Dim mismatch");
    tgtAssert(dim == next.getDimensions(), "Dim mismatch");

    size_t n = tgt::hmul(dim);
    for(size_t i = 0; i < n; ++i) {
        float p = prev.voxel(i);
        float c = current.voxel(i);
        float n = next.voxel(i);
        if(is_peak(p, c, n)) {
            counts->voxel(i) += 1;
        }
    }
}

static void count_peaks_xy(LZ4WriteableSlab<uint8_t>& counts, const VolumeAtomic<float>& vals) {
    const tgt::ivec3 dim = counts->getDimensions();
    tgtAssert(dim == (tgt::ivec3)vals.getDimensions(), "Dim mismatch");
    tgtAssert(dim.x >= 3, "x dim too small");
    tgtAssert(dim.y >= 3, "x dim too small");


    for(int y = 1; y < dim.y-1; ++y) {
        for(int x = 1; x < dim.x-1; ++x) {
            tgt::svec3 slicePos(x,y,0);

            float c = vals.voxel(slicePos);
            bool peak_x = is_peak(vals.voxel(x-1, y, 0), c, vals.voxel(x+1, y, 0));
            bool peak_y = is_peak(vals.voxel(x, y-1, 0), c, vals.voxel(x, y+1, 0));
            counts->voxel(slicePos) += (uint8_t) peak_x + (uint8_t) peak_y;
        }
    }

    auto count_y = [&] (int x, int y) {
        tgt::svec3 slicePos(x,y,0);

        float c = vals.voxel(slicePos);
        bool peak_y = is_peak(vals.voxel(x, y-1, 0), c, vals.voxel(x, y+1, 0));
        counts->voxel(slicePos) += (uint8_t) peak_y;
    };
    for(int y = 1; y < dim.y-1; ++y) {
        count_y(0, y);
        count_y(dim.x-1, y);
    }

    auto count_x = [&] (int x, int y) {
        tgt::svec3 slicePos(x,y,0);

        float c = vals.voxel(slicePos);
        bool peak_x = is_peak(vals.voxel(x-1, y, 0), c, vals.voxel(x+1, y, 0));
        counts->voxel(slicePos) += (uint8_t) peak_x;
    };
    for(int x = 1; x < dim.x-1; ++x) {
        count_x(x, 0);
        count_x(x, dim.y-1);
    }
}

static void find_structures(LZ4WriteableSlab<uint8_t>& counts, MedialStructureType type) {
    const tgt::ivec3 dim = counts->getDimensions();
    size_t n = tgt::hmul(dim);

    for(size_t i = 0; i < n; ++i) {
        auto& c = counts->voxel(i);
        c = c >= type ? 255 : 0;
    }
}

LZ4SliceVolume<uint8_t> compute_medial_structures(const VolumeBase& vol, float binarizationThreshold, MedialStructureType structureType, std::string outputPath, ProgressReporter& progressReporter) {
    SubtaskProgressReporterCollection<2> tasks(progressReporter);

    std::string outputPathTmp = VoreenApplication::app()->getUniqueTmpFilePath();

    auto distances = compute_distance_transform(vol, binarizationThreshold, outputPathTmp, tasks.get<0>());

    const tgt::ivec3 dim = vol.getDimensions();
    const tgt::svec3 sliceDim(dim.x, dim.y, 1);

    LZ4SliceVolumeBuilder<uint8_t> builder(outputPath,
            LZ4SliceVolumeMetadata(dim)
            .withOffset(vol.getOffset())
            .withSpacing(vol.getSpacing())
            .withPhysicalToWorldTransformation(vol.getPhysicalToWorldMatrix()));


    VolumeAtomic<float> prevSlice(tgt::svec3::zero);
    VolumeAtomic<float> currentSlice(tgt::svec3::zero);
    VolumeAtomic<float> nextSlice(tgt::svec3::zero);
    for(int z = 0; z < dim.z; ++z) {
        tasks.get<1>().setProgress(static_cast<float>(z)/dim.z);
        std::unique_ptr<VolumeRAM> inputSlice(vol.getSlice(z));

        if(z < dim.z-1) {
            nextSlice = distances.loadSlice(z+1);
        }
        if(z == 0) {
            currentSlice = distances.loadSlice(z);
        }

        auto countSlice = builder.getNextWriteableSlice();
        countSlice->fill(0);

        count_peaks_xy(countSlice, currentSlice);
        if(z > 0 && z < dim.z-1) {
            count_peaks_z(countSlice, prevSlice, currentSlice, nextSlice);
        }

        find_structures(countSlice, structureType);

        prevSlice = std::move(currentSlice);
        currentSlice = std::move(nextSlice);
    }

    std::move(distances).deleteFromDisk();

    progressReporter.setProgress(1.0f);

    return std::move(builder).finalize();
}

}   // namespace
