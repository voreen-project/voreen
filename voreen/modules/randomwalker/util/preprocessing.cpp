#include "preprocessing.h"

#include "voreen/core/utils/stringutils.h"

//#define VRN_RANDOMWALKER_MEAN_NOT_MEDIAN

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
                size_t child1, child2, c;

                // 1->3
                child1 = 2;
                child2 = 3;
                c = data[child2] < data[child1] ? child1 : child2;
                if(val >= data[c]) { goto end_full; } data[i] = data[c]; i = c;

                // 3->7
                child1 = i<<1;
                child2 = child1+1;
                c = data[child2] < data[child1] ? child1 : child2;
                if(val >= data[c]) { goto end_full; } data[i] = data[c]; i = c;

                // 7->15
                child1 = i<<1;
                child2 = child1+1;
                c = data[child2] < data[child1] ? child1 : child2;
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

VolumeAtomic<float> preprocessForAdaptiveParameterSetting(const VolumeAtomic<float>& img) {
    VolumeAtomic<float> output(img.getDimensions());
    const tgt::ivec3 start(0);
    const tgt::ivec3 end(img.getDimensions());
    const size_t numVoxels = tgt::hmul(img.getDimensions());

    const int k = 1;
    const int N=2*k+1;
    const tgt::ivec3 neighborhoodSize(k);

    clock_t tbegin = clock();
#ifdef VRN_RANDOMWALKER_MEAN_NOT_MEDIAN
    // mean
    auto conv = [&] (const VolumeAtomic<float>& input, VolumeAtomic<float>& output, int dim) {
        VRN_FOR_EACH_VOXEL(center, start, end) {
            tgt::ivec3 neigh(0);
            neigh[dim] = neighborhoodSize[dim];
            const tgt::ivec3 neighborhoodStart = tgt::max(start, center - neigh);
            const tgt::ivec3 neighborhoodEnd = tgt::min(end, center + neigh + tgt::ivec3(1));

            const int numNeighborhoodVoxels = tgt::hmul(neighborhoodEnd-neighborhoodStart);

            float sum=0.0f;
            VRN_FOR_EACH_VOXEL(pos, neighborhoodStart, neighborhoodEnd) {
                sum += input.voxel(pos);
            }
            float estimation = sum/numNeighborhoodVoxels;
            output.voxel(center) = estimation;
        }
    };
    VolumeAtomic<float> tmp(img.getDimensions());
    conv(img, output, 0);
    conv(output, tmp, 1);
    conv(tmp, output, 2);
#else
    // median
    tgt::ivec3 last = end - tgt::ivec3(1);
    const size_t HEAP_SIZE = N*N*N/2+1;
    tgtAssert(HEAP_SIZE == 14, "Invalid neighborhood size");
    NSmallestHeap14 heap;

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
        output.voxel(center) = heap.top();
    }
#endif

    float sumOfDifferences = 0.0f;
    VRN_FOR_EACH_VOXEL(center, start, end) {
        const tgt::ivec3 neighborhoodStart = tgt::max(start, center - neighborhoodSize);
        const tgt::ivec3 neighborhoodEnd = tgt::min(end, center + neighborhoodSize + tgt::ivec3(1));

        const int numNeighborhoodVoxels = tgt::hmul(neighborhoodEnd-neighborhoodStart);

        float estimation = output.voxel(center);
        float val = img.voxel(center);
        float diff = estimation - val;

        float neighborhoodFactor;
        if(numNeighborhoodVoxels > 1) {
            neighborhoodFactor = static_cast<float>(numNeighborhoodVoxels)/static_cast<float>(numNeighborhoodVoxels-1);
        } else {
            neighborhoodFactor = 1.0f;
        }

        sumOfDifferences += neighborhoodFactor * diff * diff;

        output.voxel(center) = estimation;
    }

#ifdef VRN_RANDOMWALKER_MEAN_NOT_MEDIAN
    const float varianceFactor = 2.0f/(N*N*N*N); //mean
#else
    tgtAssert(k==1, "Invalid k for variance factor");
    const float varianceFactor = 0.142; //median //TODO: this is for 2D. what about 3D?
#endif

    float rawVariance = sumOfDifferences/numVoxels;
    float varianceEstimation = rawVariance * varianceFactor;
    float stdEstimationInv = 1.0f/std::sqrt(varianceEstimation);

    VRN_FOR_EACH_VOXEL(center, start, end) {
        output.voxel(center) *= stdEstimationInv;
    }

    return output;
}

VolumeAtomic<float> preprocessForAdaptiveParameterSetting(const VolumeRAM& img) {
    if(const VolumeAtomic<float>* floatImg = dynamic_cast<const VolumeAtomic<float>*>(&img)) {
        return preprocessForAdaptiveParameterSetting(*floatImg);
    } else {
        tgt::svec3 dim = img.getDimensions();
        VolumeAtomic<float> converted(dim);
        VRN_FOR_EACH_VOXEL(center, tgt::svec3(0), dim) {
            converted.voxel(center) = img.getVoxelNormalized(center);
        }
        return preprocessForAdaptiveParameterSetting(converted);
    }
}

}
