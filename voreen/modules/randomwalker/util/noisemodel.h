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

#ifndef VRN_RANDOM_WALKER_NOISEMODEL_H
#define VRN_RANDOM_WALKER_NOISEMODEL_H

#include "tgt/assert.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/octree/octreeutils.h" // for index/position helper functions
#include "preprocessing.h"
#include <boost/math/special_functions/beta.hpp>

namespace voreen {

enum RWNoiseModel {
    RW_NOISE_GAUSSIAN,
    RW_NOISE_POISSON,
    RW_NOISE_GAUSSIAN_BIAN_MEAN,
    RW_NOISE_GAUSSIAN_BIAN_MEDIAN,
    RW_NOISE_VARIABLE_GAUSSIAN,
    RW_NOISE_TTEST,
};

// Parameter estimation according to
//
// A. Bian, X Jiang: Statistical Modeling Based Adaptive Parameter Setting for Random Walk Segmentation
// https://link.springer.com/chapter/10.1007%2F978-3-319-48680-2_61
struct RWNoiseModelGaussianBianMean {
    VolumeAtomic<float> mean;
    float diff_variance_inv;

    RWNoiseModelGaussianBianMean(RWNoiseModelGaussianBianMean&&) = default;

    static RWNoiseModelGaussianBianMean prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> mean = meanFilter3x3x3(vol);
        float variance = estimateVariance3x3x3(vol, mean);

        const int k = 1;
        const int N=2*k+1;
        const int N4=N*N*N*N;
        float correction_factor = 2.0f/N4;

        float diff_variance = variance * correction_factor;
        diff_variance = std::max(diff_variance, std::numeric_limits<float>::min());

        return RWNoiseModelGaussianBianMean {
            std::move(mean),
            1.0f/diff_variance,
        };
    }
    static RWNoiseModelGaussianBianMean prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelGaussianBianMean::prepare(toVolumeAtomicFloat(vol), rwm);
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        float voxelIntensity = mean.voxel(voxel);
        float neighborIntensity = mean.voxel(neighbor);
        float beta = 2.0f * betaBias * diff_variance_inv;
        float intDiff = (voxelIntensity - neighborIntensity);
        float intDiffSqr = intDiff*intDiff;
        float weight = exp(-beta * intDiffSqr);
        return weight;
    }
};
struct RWNoiseModelGaussianBianMedian {
    VolumeAtomic<float> mean;
    float diff_variance_inv;

    RWNoiseModelGaussianBianMedian(RWNoiseModelGaussianBianMedian&&) = default;

    static RWNoiseModelGaussianBianMedian prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> mean = medianFilter3x3x3(vol);
        float variance = estimateVariance3x3x3(vol, mean);

        // TODO: This is for 2D (value from Angs paper) and seems to work fine.
        // However, it is unclear if a better value for 3D can be found.
        float correction_factor = 0.142f;

        float diff_variance = variance * correction_factor;
        diff_variance = std::max(diff_variance, std::numeric_limits<float>::min());

        return RWNoiseModelGaussianBianMedian {
            std::move(mean),
            1.0f/diff_variance,
        };
    }
    static RWNoiseModelGaussianBianMedian prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelGaussianBianMedian::prepare(toVolumeAtomicFloat(vol), rwm);
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        float voxelIntensity = mean.voxel(voxel);
        float neighborIntensity = mean.voxel(neighbor);
        float beta = 2.0f * betaBias * diff_variance_inv;
        float intDiff = (voxelIntensity - neighborIntensity);
        float intDiffSqr = intDiff*intDiff;
        float weight = exp(-beta * intDiffSqr);
        return weight;
    }
};

VolumeAtomic<tgt::ivec3> findBestCenters(const VolumeAtomic<float>& image, const VolumeAtomic<float>& mean, const VolumeAtomic<float>& variance, int filter_extent);
float evalTTest(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent);
float evalVariableGaussian(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent);

// Parameter estimation according to
//
// A. Bian, X Jiang: T-Test Based Adaptive Random Walk Segmentation Under Multiplicative Speckle Noise Model
// https://link.springer.com/chapter/10.1007%2F978-3-319-54427-4_41
template<int filter_extent>
struct RWNoiseModelTTest {

    VolumeAtomic<float> image;
    VolumeAtomic<tgt::ivec3> best_centers;

    RWNoiseModelTTest(RWNoiseModelTTest&&) = default;

    static RWNoiseModelTTest prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {

        VolumeAtomic<float> image = vol.copy();
        VolumeAtomic<float> mean = meanFilter<filter_extent>(image);
        VolumeAtomic<float> variance = variances<filter_extent>(image, mean);

        auto bestCenters = findBestCenters(image, mean, variance, filter_extent);
        return RWNoiseModelTTest {
            std::move(image),
            std::move(bestCenters),
        };
    }
    static RWNoiseModelTTest prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelTTest::prepare(toVolumeAtomicFloat(vol), rwm);
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        return evalTTest(image, best_centers, voxel, neighbor, filter_extent);
    }
};
template<int filter_extent>
struct RWNoiseModelVariableGaussian {

    VolumeAtomic<float> image;
    VolumeAtomic<tgt::ivec3> best_centers;

    RWNoiseModelVariableGaussian(RWNoiseModelVariableGaussian&&) = default;

    static RWNoiseModelVariableGaussian prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {

        VolumeAtomic<float> image = vol.copy();
        VolumeAtomic<float> mean = meanFilter<filter_extent>(image);
        VolumeAtomic<float> variance = variances<filter_extent>(image, mean);

        auto bestCenters = findBestCenters(image, mean, variance, filter_extent);
        return RWNoiseModelVariableGaussian {
            std::move(image),
            std::move(bestCenters),
        };
    }
    static RWNoiseModelVariableGaussian prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelVariableGaussian::prepare(toVolumeAtomicFloat(vol), rwm);
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        return evalVariableGaussian(image, best_centers, voxel, neighbor, filter_extent);
    }
};

struct RWNoiseModelGaussian {
    VolumeAtomic<float> values;
    float variance_inv;

    RWNoiseModelGaussian(RWNoiseModelGaussian&&) = default;

    static RWNoiseModelGaussian prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> mean = meanFilter3x3x3(vol);
        float variance = estimateVariance3x3x3(vol, mean);

        return RWNoiseModelGaussian {
            vol.copy(),
            1.0f/variance,
        };
    }
    static RWNoiseModelGaussian prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelGaussian::prepare(toVolumeAtomicFloat(vol), rwm);
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        float voxelIntensity = values.voxel(voxel);
        float neighborIntensity = values.voxel(neighbor);

        float beta = 0.125f * betaBias * variance_inv;
        float intDiff = (voxelIntensity - neighborIntensity);
        float intDiffSqr = intDiff*intDiff;
        float weight = exp(-beta * intDiffSqr);
        return weight;
    }
};
struct RWNoiseModelPoisson {
    VolumeAtomic<float> values;

    RWNoiseModelPoisson(RWNoiseModelPoisson&&) = default;

    static RWNoiseModelPoisson prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        return RWNoiseModelPoisson {
            applyRWM(vol, rwm),
        };
    }
    static RWNoiseModelPoisson prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelPoisson {
            applyRWM(vol, rwm),
        };
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        float voxelIntensity = values.voxel(voxel);
        float neighborIntensity = values.voxel(neighbor);

        float beta = 0.5f * betaBias;
        float weight;

        float d = sqrt(voxelIntensity) - sqrt(neighborIntensity);
        weight = exp(-beta * d * d);
        return weight;
    }
};

}

#endif
