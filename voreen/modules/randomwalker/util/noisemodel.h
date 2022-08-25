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

template<RWNoiseModel N>
struct RWNoiseModelParameters {
};

template<RWNoiseModel N>
struct RWNoiseModelWeights {
};

template<typename Parameters, typename FitFun>
VolumeAtomic<tgt::ivec3> findBestCenters(const VolumeAtomic<float>& image, const Parameters* parameters, FitFun fitFun, int filter_extent);

struct GaussianParametersVariableSigma {
    float mean;
    float mul;
    float add;

    static std::vector<GaussianParametersVariableSigma> create(const VolumeAtomic<float>& image, int filter_extent);
    static float fit(GaussianParametersVariableSigma params, float f);
};

// Parameter estimation according to
//
// A. Bian, X Jiang: Statistical Modeling Based Adaptive Parameter Setting for Random Walk Segmentation
// https://link.springer.com/chapter/10.1007%2F978-3-319-48680-2_61
template<>
struct RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEAN> {
    VolumeAtomic<float> mean;
    float diff_variance_inv;

    RWNoiseModelWeights(RWNoiseModelWeights&&) = default;
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

template<>
struct RWNoiseModelParameters<RW_NOISE_GAUSSIAN_BIAN_MEAN> {
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEAN> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> mean = meanFilter3x3x3(vol);
        float variance = estimateVariance3x3x3(vol, mean);

        const int k = 1;
        const int N=2*k+1;
        const int N4=N*N*N*N;
        float correction_factor = 2.0f/N4;

        float diff_variance = variance * correction_factor;
        diff_variance = std::max(diff_variance, std::numeric_limits<float>::min());

        return RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEAN> {
            std::move(mean),
            1.0f/diff_variance,
        };
    }
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEAN> prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};

template<>
struct RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEDIAN> {
    VolumeAtomic<float> mean;
    float diff_variance_inv;

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
template<>
struct RWNoiseModelParameters<RW_NOISE_GAUSSIAN_BIAN_MEDIAN> {
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEDIAN> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> mean = medianFilter3x3x3(vol);
        float variance = estimateVariance3x3x3(vol, mean);

        // TODO: This is for 2D (value from Angs paper) and seems to work fine.
        // However, it is unclear if a better value for 3D can be found.
        float correction_factor = 0.142f;

        float diff_variance = variance * correction_factor;
        diff_variance = std::max(diff_variance, std::numeric_limits<float>::min());

        return RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEDIAN> {
            std::move(mean),
            1.0f/diff_variance,
        };
    }
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEDIAN> prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};

VolumeAtomic<tgt::ivec3> findBestCenters(const VolumeAtomic<float>& image, const VolumeAtomic<float>& mean, const VolumeAtomic<float>& variance, int filter_extent);
float evalTTest(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent);
float evalVariableGaussian(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent);
float evalPoisson(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent);

// Parameter estimation according to
//
// A. Bian, X Jiang: T-Test Based Adaptive Random Walk Segmentation Under Multiplicative Speckle Noise Model
// https://link.springer.com/chapter/10.1007%2F978-3-319-54427-4_41

template<>
struct RWNoiseModelWeights<RW_NOISE_TTEST> {
    VolumeAtomic<float> image;
    VolumeAtomic<tgt::ivec3> best_centers;
    int filter_extent;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        return evalTTest(image, best_centers, voxel, neighbor, filter_extent);
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_TTEST> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_TTEST> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> image = vol.copy();
        auto parameters = GaussianParametersVariableSigma::create(image, filter_extent);

        auto bestCenters = findBestCenters(image, parameters.data(), GaussianParametersVariableSigma::fit, filter_extent);
        return RWNoiseModelWeights<RW_NOISE_TTEST> {
            std::move(image),
            std::move(bestCenters),
            filter_extent,
        };
    }
    RWNoiseModelWeights<RW_NOISE_TTEST> prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};

template<>
struct RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN> {
    VolumeAtomic<float> image;
    VolumeAtomic<tgt::ivec3> best_centers;
    int filter_extent;
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        return evalVariableGaussian(image, best_centers, voxel, neighbor, filter_extent);
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_VARIABLE_GAUSSIAN> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> image = vol.copy();
        auto parameters = GaussianParametersVariableSigma::create(image, filter_extent);

        auto bestCenters = findBestCenters(image, parameters.data(), GaussianParametersVariableSigma::fit, filter_extent);
        return RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN> {
            std::move(image),
            std::move(bestCenters),
            filter_extent,
        };
    }
    RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN> prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};


template<>
struct RWNoiseModelWeights<RW_NOISE_GAUSSIAN> {

    VolumeAtomic<float> values;
    float variance_inv;
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
template<>
struct RWNoiseModelParameters<RW_NOISE_GAUSSIAN> {

    RWNoiseModelWeights<RW_NOISE_GAUSSIAN> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> mean = meanFilter3x3x3(vol);
        float variance = estimateVariance3x3x3(vol, mean);

        return RWNoiseModelWeights<RW_NOISE_GAUSSIAN> {
            vol.copy(),
            1.0f/variance,
        };
    }
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN> prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};


template<>
struct RWNoiseModelWeights<RW_NOISE_POISSON> {
    VolumeAtomic<float> values;
    VolumeAtomic<tgt::ivec3> best_centers;
    int filter_extent;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        return evalPoisson(values, best_centers, voxel, neighbor, filter_extent);
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_POISSON> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_POISSON> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> image = applyRWM(vol, rwm);
        VolumeAtomic<float> parameters = meanFilter(image, filter_extent);

        auto bestCenters = findBestCenters(image, parameters.voxel(), [] (float lambda, float sample) {
            return -lambda + std::log(lambda) * sample - std::lgamma(sample + 1);
        }, filter_extent);

        return RWNoiseModelWeights<RW_NOISE_POISSON> {
            std::move(image),
            std::move(bestCenters),
            filter_extent,
        };
    }
    RWNoiseModelWeights<RW_NOISE_POISSON> prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};

template<typename Parameters, typename FitFun>
VolumeAtomic<tgt::ivec3> findBestCenters(const VolumeAtomic<float>& image, const Parameters* parameters, FitFun fitFun, int filter_extent) {
    // We would take parameters as a VolumeAtomic argument, but cannot do that
    // because, for example, GaussianParametersVariableSigma is not a valid
    // VolumeElement.
    //
    // parameters MUST therefore be the dimension-order linearized storage of a
    // volume the same size as image.

    tgt::ivec3 dim = image.getDimensions();

    VolumeAtomic<tgt::ivec3> best_centers(dim);

    VRN_FOR_EACH_VOXEL(p, tgt::ivec3::zero, dim) {
        float f = image.voxel(p);

        const tgt::ivec3 extent(filter_extent);

        tgt::ivec3 begin = tgt::max(tgt::ivec3::zero, p - extent);
        tgt::ivec3 end = tgt::min(dim, p + tgt::ivec3::one + extent);

        float max = -std::numeric_limits<float>::infinity();
        size_t argmax_index = 0;

        VRN_FOR_EACH_VOXEL_INDEX(l, begin, end, dim, {
            float fit = fitFun(parameters[l], f);
            tgtAssert(std::isfinite(fit) && !std::isnan(fit), "invalid fit value");
            if(fit > max) {
                max = fit;
                argmax_index = l;
            }
        })

        tgt::ivec3 argmax = linearCoordToCubic(argmax_index, dim);
        tgtAssert(tgt::max(tgt::abs(p - argmax)) <= filter_extent, "invalid pos");
        best_centers.voxel(p) = argmax;
    }

    return best_centers;
}

}
#endif
