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

#include "preprocessing.h"
#include <array>
#include <cmath>

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#define VRN_RANDOMWALKER_MEAN_NOT_MEDIAN

namespace voreen {


namespace {
struct NSmallestHeap14 {
    NSmallestHeap14()
        : data()
        , numElements(0)
    {
        data[0] = INFINITY; //sentinel for first 14 values (upper push branch)
        data[15] = -INFINITY; //sentinel for child2 of index 7

    }
    float top() {
        return data[1];
    }
    void push(float val) {
        if(numElements < 14) {
            numElements++;
            size_t i=numElements;
            size_t parent;

            parent = i>>1; if(data[parent] >= val) { goto end_not_full; } data[i] = data[parent]; i=parent; // 15->7
            parent = i>>1; if(data[parent] >= val) { goto end_not_full; } data[i] = data[parent]; i=parent; // 7->3
            parent = i>>1; if(data[parent] >= val) { goto end_not_full; } data[i] = data[parent]; i=parent; // 3->1
end_not_full:
            data[i] = val;
        } else {
            if(val < top()) {
                size_t i = 1;
                size_t child1, child2, c, childSelect;

                // 1->3
                child1 = 2;
                child2 = 3;
                childSelect = size_t(data[child1] < data[child2]);
                c = child1 + childSelect;
                if(val >= data[c]) { goto end_full; } data[i] = data[c]; i = c;

                // 3->7
                child1 = i<<1;
                child2 = child1+1;
                childSelect = size_t(data[child1] < data[child2]);
                c = child1 + childSelect;
                if(val >= data[c]) { goto end_full; } data[i] = data[c]; i = c;

                // 7->15
                child1 = i<<1;
                child2 = child1+1;
                childSelect = size_t(data[child1] < data[child2]);
                c = child1 + childSelect;
                if(val >= data[c]) { goto end_full; } data[i] = data[c]; i = c;
end_full:
                data[i] = val;
            }
        }
    }
    void clear() {
        numElements = 0;
    }
    std::array<float, 16> data;
    size_t numElements;
};
}


VolumeAtomic<float> meanFilter3x3x3(const VolumeAtomic<float>& img) {
    return meanFilter<1>(img);
}

VolumeAtomic<float> medianFilter3x3x3(const VolumeAtomic<float>& img) {
    const tgt::ivec3 start(0);
    const tgt::ivec3 end(img.getDimensions());
    const size_t numVoxels = tgt::hmul(img.getDimensions());

    const int k = 1;
    const int N=2*k+1;
    const tgt::ivec3 neighborhoodSize(k);

    // median
    tgt::ivec3 last = end - tgt::ivec3(1);
    const size_t HEAP_SIZE = N*N*N/2+1;
    tgtAssert(HEAP_SIZE == 14, "Invalid neighborhood size");
    NSmallestHeap14 heap;

    VolumeAtomic<float> res(img.getDimensions());
    VRN_FOR_EACH_VOXEL(center, start, end) {
        const tgt::ivec3 neighborhoodStart = center - neighborhoodSize;
        const tgt::ivec3 neighborhoodEnd = center + neighborhoodSize + tgt::ivec3(1);


        float pivot = img.voxel(start);
        heap.clear();
        VRN_FOR_EACH_VOXEL(pos, neighborhoodStart, neighborhoodEnd) {
            tgt::ivec3 p = tgt::clamp(pos, start, last);
            float val = img.voxel(p);
            heap.push(val);
        }
        res.voxel(center) = heap.top();
    }

    return res;
}

float estimateVariance3x3x3(const VolumeAtomic<float>& img, const VolumeAtomic<float>& mean) {
    const int k = 1;
    const tgt::ivec3 extent(k);
    const tgt::ivec3 start(extent);
    const tgt::ivec3 end(tgt::ivec3(img.getDimensions()) - extent);
    const size_t numVoxels = tgt::hmul(img.getDimensions());

    tgtAssert(img.getDimensions() == mean.getDimensions(), "Dimension mismatch");

    const int N=2*k+1;
    const int numNeighborhoodVoxels = N*N*N;

    float sumOfDifferences = 0.0f;
    VRN_FOR_EACH_VOXEL(v, start, end) {
        float estimation = mean.voxel(v);
        float val = img.voxel(v);
        float diff = estimation - val;

        sumOfDifferences += diff * diff;
    }
    float neighborhoodFactor = static_cast<float>(numNeighborhoodVoxels)/static_cast<float>(numNeighborhoodVoxels-1);

    float varianceEstimation = sumOfDifferences/numVoxels * neighborhoodFactor;

    return varianceEstimation;
}

VolumeAtomic<float> toVolumeAtomicFloat(const VolumeRAM& img) {
    tgtAssert(img.getNumChannels() == 1, "Only volumes with one channel expected");
    size_t voxels = img.getNumVoxels();
    tgt::svec3 dim = img.getDimensions();
    VolumeAtomic<float> converted(dim);
    for(size_t i=0; i<voxels; ++i) {
        converted.voxel(i) = img.getVoxelNormalized(i);
    }
    return converted;
}

VolumeAtomic<float> applyRWM(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
    tgtAssert(vol.getNumChannels() == 1, "Only volumes with one channel expected");
    size_t voxels = vol.getNumVoxels();
    VolumeAtomic<float> converted(vol.getDimensions());
    for(size_t i=0; i<voxels; ++i) {
        converted.voxel(i) = rwm.normalizedToRealWorld(vol.voxel(i));
    }
    return converted;
}

VolumeAtomic<float> applyRWM(const VolumeRAM& vol, RealWorldMapping rwm) {
    tgtAssert(vol.getNumChannels() == 1, "Only volumes with one channel expected");
    size_t voxels = vol.getNumVoxels();
    VolumeAtomic<float> converted(vol.getDimensions());
    for(size_t i=0; i<voxels; ++i) {
        converted.voxel(i) = rwm.normalizedToRealWorld(vol.getVoxelNormalized(i));
    }
    return converted;
}

}
