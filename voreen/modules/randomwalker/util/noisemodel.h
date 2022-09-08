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
    RW_NOISE_GAUSSIAN_HIERARCHICAL,
    RW_NOISE_POISSON_HIERARCHICAL,
    RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL,
    RW_NOISE_TTEST_HIERARCHICAL,
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
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        float voxelIntensity = mean.voxel(voxel);
        float neighborIntensity = mean.voxel(neighbor);
        float beta = 2.0f * diff_variance_inv;
        float intDiff = (voxelIntensity - neighborIntensity);
        float intDiffSqr = intDiff*intDiff;
        float weight = exp(-beta * intDiffSqr);
        return weight;
    }
};

template<>
struct RWNoiseModelParameters<RW_NOISE_GAUSSIAN_BIAN_MEAN> {
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEAN> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level=0, const VolumeAtomic<float>* variances=nullptr) {
        VolumeAtomic<float> mean = meanFilter3x3x3(vol);
        float variance = estimateVariance(vol, mean, 1);

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
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEAN> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level=0) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};

template<>
struct RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEDIAN> {
    VolumeAtomic<float> mean;
    float diff_variance_inv;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        float voxelIntensity = mean.voxel(voxel);
        float neighborIntensity = mean.voxel(neighbor);
        float beta = 2.0f * diff_variance_inv;
        float intDiff = (voxelIntensity - neighborIntensity);
        float intDiffSqr = intDiff*intDiff;
        float weight = exp(-beta * intDiffSqr);
        return weight;
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_GAUSSIAN_BIAN_MEDIAN> {
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEDIAN> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level=0, const VolumeAtomic<float>* variances=nullptr) {
        VolumeAtomic<float> mean = medianFilter3x3x3(vol);
        float variance = estimateVariance(vol, mean, 1);

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
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_BIAN_MEDIAN> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level=0) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};

VolumeAtomic<tgt::ivec3> findBestCenters(const VolumeAtomic<float>& image, const VolumeAtomic<float>& mean, const VolumeAtomic<float>& variance, int filter_extent);
float evalTTest(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent);
float evalVariableGaussian(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent);
float evalPoisson(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent, int level);
float evalConstGaussian(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, float variance_inv, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent);

float bhattacharyyaPoisson(float sum1, float sum2);
float bhattacharyyaConstGaussian(float mean1, float mean2, float variance_inv, float n);
float bhattacharyyaVarGaussian(float mean1, float mean2, float var1, float var2, size_t n);
float ttestFunction(float mean1, float mean2, float var1, float var2, float n1, float n2);

inline int num_voxels_in_level(int level) {
    return 1 << (3*level); //8^level
}

// Parameter estimation according to
//
// A. Bian, X Jiang: T-Test Based Adaptive Random Walk Segmentation Under Multiplicative Speckle Noise Model
// https://link.springer.com/chapter/10.1007%2F978-3-319-54427-4_41

template<>
struct RWNoiseModelWeights<RW_NOISE_TTEST> {
    VolumeAtomic<float> image;
    VolumeAtomic<tgt::ivec3> best_centers;
    int filter_extent;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        return evalTTest(image, best_centers, voxel, neighbor, filter_extent);
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_TTEST> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_TTEST> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level=0, const VolumeAtomic<float>* variances=nullptr) {
        VolumeAtomic<float> image = vol.copy();
        auto parameters = GaussianParametersVariableSigma::create(image, filter_extent);

        auto bestCenters = findBestCenters(image, parameters.data(), GaussianParametersVariableSigma::fit, filter_extent);
        return RWNoiseModelWeights<RW_NOISE_TTEST> {
            std::move(image),
            std::move(bestCenters),
            filter_extent,
        };
    }
    RWNoiseModelWeights<RW_NOISE_TTEST> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level=0) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};

template<>
struct RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN> {
    VolumeAtomic<float> image;
    VolumeAtomic<tgt::ivec3> best_centers;
    int filter_extent;
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        return evalVariableGaussian(image, best_centers, voxel, neighbor, filter_extent);
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_VARIABLE_GAUSSIAN> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level=0, const VolumeAtomic<float>* variances=nullptr) {
        VolumeAtomic<float> image = vol.copy();
        auto parameters = GaussianParametersVariableSigma::create(image, filter_extent);

        auto bestCenters = findBestCenters(image, parameters.data(), GaussianParametersVariableSigma::fit, filter_extent);
        return RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN> {
            std::move(image),
            std::move(bestCenters),
            filter_extent,
        };
    }
    RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level=0) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};


template<>
struct RWNoiseModelWeights<RW_NOISE_GAUSSIAN> {
    VolumeAtomic<float> values;
    VolumeAtomic<tgt::ivec3> best_centers;
    float variance_inv;
    int filter_extent;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        return evalConstGaussian(values, best_centers, variance_inv, voxel, neighbor, filter_extent);
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_GAUSSIAN> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_GAUSSIAN> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level=0, const VolumeAtomic<float>* variances=nullptr) {
        VolumeAtomic<float> mean = meanFilter(vol, filter_extent);
        float variance = estimateVariance(vol, mean, filter_extent);
        tgtAssert(variance > 0, "Invalid variance");

        float variance_inv = 1.0f/variance;

        auto bestCenters = findBestCenters(vol, mean.voxel(), [&] (float mu, float sample) {
            float diff = mu-sample;
            float coeff = diff*diff * variance_inv * 0.5f;
            return -coeff;
        }, filter_extent);

        return RWNoiseModelWeights<RW_NOISE_GAUSSIAN> {
            vol.copy(),
            std::move(bestCenters),
            variance_inv,
            filter_extent,
        };
    }
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level=0) {
        return prepare(toVolumeAtomicFloat(vol), rwm);
    }
};


template<>
struct RWNoiseModelWeights<RW_NOISE_POISSON> {
    VolumeAtomic<float> values;
    VolumeAtomic<tgt::ivec3> best_centers;
    int filter_extent;
    int level;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        // Theoretically, choosing level here should be better since the noise
        // distribution is not poisson otherwise, but in practice always
        // choosing 0 appears to yield better results...
        //float evalLevel = level;
        float evalLevel = 0;

        return evalPoisson(values, best_centers, voxel, neighbor, filter_extent, evalLevel);
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_POISSON> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_POISSON> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level=0, const VolumeAtomic<float>* variances=nullptr) {
        VolumeAtomic<float> image = applyRWM(vol, rwm);
        VolumeAtomic<float> parameters = meanFilter(image, filter_extent);

        auto bestCenters = findBestCenters(image, parameters.voxel(), [] (float lambda, float sample) {
            if(lambda == 0) {
                return sample == 0 ? 0.0f : -std::numeric_limits<float>::infinity();
            }
            return -lambda + std::log(lambda) * sample - std::lgamma(sample + 1);
        }, filter_extent);

        return RWNoiseModelWeights<RW_NOISE_POISSON> {
            std::move(image),
            std::move(bestCenters),
            filter_extent,
            level,
        };
    }
    RWNoiseModelWeights<RW_NOISE_POISSON> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level=0) {
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

template<>
struct RWNoiseModelWeights<RW_NOISE_GAUSSIAN_HIERARCHICAL> {
    VolumeAtomic<float> values;
    VolumeAtomic<tgt::ivec3> best_centers;
    float variance_inv;
    int filter_extent;
    int level;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        if(level==0) {
            return evalConstGaussian(values, best_centers, variance_inv, voxel, neighbor, filter_extent);
        } else {
            float mean1 = values.voxel(voxel);
            float mean2 = values.voxel(neighbor);

            int n = num_voxels_in_level(level);

            return bhattacharyyaConstGaussian(mean1, mean2, variance_inv, n);
        }
    }
};

template<>
struct RWNoiseModelParameters<RW_NOISE_GAUSSIAN_HIERARCHICAL> {
    int filter_extent;
    float variance;

    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_HIERARCHICAL> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level, const VolumeAtomic<float>* variances=nullptr) {
        float variance_inv = 1.0f/variance;

        VolumeAtomic<tgt::ivec3> bestCenters(tgt::svec3(0,0,0), false);
        if(level == 0) {
            VolumeAtomic<float> mean = meanFilter(vol, filter_extent);
            bestCenters = findBestCenters(vol, mean.voxel(), [&] (float mu, float sample) {
                float diff = mu-sample;
                float coeff = diff*diff * variance_inv * 0.5f;
                return -coeff;
            }, filter_extent);
        }

        return RWNoiseModelWeights<RW_NOISE_GAUSSIAN_HIERARCHICAL> {
            vol.copy(),
            std::move(bestCenters),
            variance_inv,
            filter_extent,
            level,
        };
    }
    RWNoiseModelWeights<RW_NOISE_GAUSSIAN_HIERARCHICAL> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level) {
        return prepare(toVolumeAtomicFloat(vol), rwm, level);
    }
};

template<>
struct RWNoiseModelWeights<RW_NOISE_POISSON_HIERARCHICAL> {
    VolumeAtomic<float> values;
    VolumeAtomic<tgt::ivec3> best_centers;
    int filter_extent;
    int level;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        if(level == 0) {
            return evalPoisson(values, best_centers, voxel, neighbor, filter_extent, level);
        } else {
            float v1 = values.voxel(voxel);
            float v2 = values.voxel(neighbor);

            int mult = num_voxels_in_level(level);
            float sum1 = v1*mult;
            float sum2 = v2*mult;

            return bhattacharyyaPoisson(sum1, sum2);
        }
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_POISSON_HIERARCHICAL> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_POISSON_HIERARCHICAL> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level, const VolumeAtomic<float>* variances=nullptr) {
        VolumeAtomic<float> image = applyRWM(vol, rwm);

        VolumeAtomic<tgt::ivec3> bestCenters(tgt::svec3(0,0,0), false);
        if(level == 0) {
            VolumeAtomic<float> parameters = meanFilter(image, filter_extent);
            bestCenters = findBestCenters(image, parameters.voxel(), [] (float lambda, float sample) {
                if(lambda == 0) {
                    return sample == 0 ? 0.0f : -std::numeric_limits<float>::infinity();
                }
                return -lambda + std::log(lambda) * sample - std::lgamma(sample + 1);
            }, filter_extent);
        }

        return RWNoiseModelWeights<RW_NOISE_POISSON_HIERARCHICAL> {
            std::move(image),
            std::move(bestCenters),
            filter_extent,
            level,
        };
    }
    RWNoiseModelWeights<RW_NOISE_POISSON_HIERARCHICAL> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level) {
        return prepare(toVolumeAtomicFloat(vol), rwm, level);
    }
};

template<>
struct RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL> {
    VolumeAtomic<float> image;
    VolumeAtomic<float> variances;
    VolumeAtomic<tgt::ivec3> best_centers;
    int filter_extent;
    int level;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        if(level == 0) {
            return evalVariableGaussian(image, best_centers, voxel, neighbor, filter_extent);
        } else {
            float mean1 = image.voxel(voxel);
            float mean2 = image.voxel(neighbor);

            float var1 = variances.voxel(voxel);
            float var2 = variances.voxel(neighbor);

            size_t n = num_voxels_in_level(level);

            return bhattacharyyaVarGaussian(mean1, mean2, var1, var2, n);
        }
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level, const VolumeAtomic<float>* variances) {
        VolumeAtomic<tgt::ivec3> bestCenters(tgt::svec3(0,0,0), false);
        if(level == 0) {
            auto parameters = GaussianParametersVariableSigma::create(vol, filter_extent);

            bestCenters = findBestCenters(vol, parameters.data(), GaussianParametersVariableSigma::fit, filter_extent);
        }

        return RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL> {
            vol.copy(),
            variances->copy(),
            std::move(bestCenters),
            filter_extent,
            level,
        };
    }
    RWNoiseModelWeights<RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level, const VolumeAtomic<float>* variances) {
        return prepare(toVolumeAtomicFloat(vol), rwm, level, variances);
    }
};

template<>
struct RWNoiseModelWeights<RW_NOISE_TTEST_HIERARCHICAL> {
    VolumeAtomic<float> image;
    VolumeAtomic<float> variances;
    VolumeAtomic<tgt::ivec3> best_centers;
    int filter_extent;
    int level;

    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor) const {
        if(level == 0) {
            return evalTTest(image, best_centers, voxel, neighbor, filter_extent);
        } else {
            float mean1 = image.voxel(voxel);
            float mean2 = image.voxel(neighbor);

            float var1 = variances.voxel(voxel);
            float var2 = variances.voxel(neighbor);

            size_t n = num_voxels_in_level(level);

            return ttestFunction(mean1, mean2, var1, var2, n, n);
        }
    }
};
template<>
struct RWNoiseModelParameters<RW_NOISE_TTEST_HIERARCHICAL> {
    int filter_extent;

    RWNoiseModelWeights<RW_NOISE_TTEST_HIERARCHICAL> prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm, int level, const VolumeAtomic<float>* variances) {
        VolumeAtomic<tgt::ivec3> bestCenters(tgt::svec3(0,0,0), false);
        if(level == 0) {
            auto parameters = GaussianParametersVariableSigma::create(vol, filter_extent);

            bestCenters = findBestCenters(vol, parameters.data(), GaussianParametersVariableSigma::fit, filter_extent);
        }

        return RWNoiseModelWeights<RW_NOISE_TTEST_HIERARCHICAL> {
            vol.copy(),
            variances->copy(),
            std::move(bestCenters),
            filter_extent,
            level,
        };
    }
    RWNoiseModelWeights<RW_NOISE_TTEST_HIERARCHICAL> prepare(const VolumeRAM& vol, RealWorldMapping rwm, int level, const VolumeAtomic<float>* variances) {
        return prepare(toVolumeAtomicFloat(vol), rwm, level, variances);
    }
};

}
#endif
