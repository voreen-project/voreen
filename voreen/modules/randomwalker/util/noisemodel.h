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
#include "preprocessing.h"
#include <boost/math/special_functions/beta.hpp>

namespace voreen {

enum RWNoiseModel {
    RW_NOISE_GAUSSIAN,
    RW_NOISE_POISSON,
    RW_NOISE_GAUSSIAN_BIAN,
    RW_NOISE_TTEST,
};

// Parameter estimation according to
//
// A. Bian, X Jiang: Statistical Modeling Based Adaptive Parameter Setting for Random Walk Segmentation
// https://link.springer.com/chapter/10.1007%2F978-3-319-48680-2_61
struct RWNoiseModelGaussianBian {
    VolumeAtomic<float> mean;
    float diff_variance_inv;

    RWNoiseModelGaussianBian(RWNoiseModelGaussianBian&&) = default;

    static RWNoiseModelGaussianBian prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {
        VolumeAtomic<float> mean = meanFilter3x3x3(vol);
        float variance = estimateVariance3x3x3(vol, mean);

        // Careful: This is _only_ valid for the mean filter as an estimator!!
        const int k = 1;
        const int N=2*k+1;
        const int N4=N*N*N*N;
        float correction_factor = 2.0f/N4;

        float diff_variance = variance * correction_factor;
        diff_variance = std::max(diff_variance, std::numeric_limits<float>::min());

        return RWNoiseModelGaussianBian {
            std::move(mean),
            1.0f/diff_variance,
        };
    }
    static RWNoiseModelGaussianBian prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelGaussianBian::prepare(toVolumeAtomicFloat(vol), rwm);
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

inline bool cmp_linear(const tgt::ivec3& p1, const tgt::ivec3& p2) {
    return p1.z < p2.z ||
        p1.z == p2.z && p1.y < p2.y ||
        p1.z == p2.z && p1.y == p2.y && p1.x < p2.x;
}

inline float square(float s) {
    return s*s;
}

inline float gaussian_pdf_exp(float val, float mean, float variance) {
    return std::exp(-0.5f/variance * square(val-mean));
}

inline float gaussian_pdf(float val, float mean, float variance) {
    return gaussian_pdf_exp(val, mean, variance)/std::sqrt(2.0f * 3.141592654f* variance);
}

// Parameter estimation according to
//
// A. Bian, X Jiang: T-Test Based Adaptive Random Walk Segmentation Under Multiplicative Speckle Noise Model
// https://link.springer.com/chapter/10.1007%2F978-3-319-54427-4_41
template<int filter_extent>
struct RWNoiseModelTTest {
    VolumeAtomic<float> image;
    VolumeAtomic<float> mean;
    VolumeAtomic<float> variance;

    RWNoiseModelTTest(RWNoiseModelTTest&&) = default;

    static RWNoiseModelTTest prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {

        VolumeAtomic<float> image = vol.copy();
        VolumeAtomic<float> mean = meanFilter<filter_extent>(image);
        VolumeAtomic<float> variance = variances<filter_extent>(image, mean);

        return RWNoiseModelTTest {
            std::move(image),
            std::move(mean),
            std::move(variance),
        };
    }
    static RWNoiseModelTTest prepare(const VolumeRAM& vol, RealWorldMapping rwm) {
        return RWNoiseModelTTest::prepare(toVolumeAtomicFloat(vol), rwm);
    }
    float getEdgeWeight(tgt::svec3 voxel, tgt::svec3 neighbor, float betaBias) const {
        tgt::ivec3 p1, p2;

        // Force weight function to be symmetric. If we don't explicitly do
        // this, there may be small differences due to numerical inaccuracies,
        // which will make the matrix asymmetric, which will 1. trigger an
        // assertion in the solver and 2. lead to undesired results from the
        // solver (values running away from the [0,1] range).
        if(cmp_linear(voxel, neighbor)) {
            p1 = voxel;
            p2 = neighbor;
        } else {
            p1 = neighbor;
            p2 = voxel;
        }

        const tgt::ivec3 dim(image.getDimensions());

        const tgt::ivec3 extent(filter_extent);

        auto best_neighborhood = [&] (tgt::ivec3 p) -> tgt::ivec3 {
            float f = image.voxel(p);

            tgt::ivec3 begin = tgt::max(tgt::ivec3::zero, p - extent);
            tgt::ivec3 end = tgt::min(dim, p + tgt::ivec3::one + extent);

            float best = 0.0f;
            auto best_p = p;
            VRN_FOR_EACH_VOXEL(n, begin, end) {
                float pdf_val = gaussian_pdf(f, mean.voxel(n), variance.voxel(n));
                if(best < pdf_val) {
                    best = pdf_val;
                    best_p = n;
                }
            }
            tgtAssert(best > 0, "invalid probability");
            return best_p;
        };

        auto neighborhood = [&] (tgt::ivec3 p) -> std::vector<tgt::ivec3> {
            tgt::ivec3 begin = tgt::max(tgt::ivec3::zero, p - extent);
            tgt::ivec3 end = tgt::min(dim, p + tgt::ivec3::one + extent);

            std::vector<tgt::ivec3> out;
            VRN_FOR_EACH_VOXEL(n, begin, end) {
                out.push_back(n);
            }
            return out;
        };

        auto mean_and_var = [&] (const std::vector<tgt::ivec3>& vals) -> std::pair<float, float> {

            float sum = 0.0f;
            for(const auto& p : vals) {
                sum += image.voxel(p);
            }
            float mean = sum/vals.size();

            float sq_sum = 0.0f;
            for(const auto& p : vals) {
                sq_sum += square(mean-image.voxel(p));
            }
            float variance = sq_sum/(vals.size()-1);

            return {mean, variance};
        };

        auto best_center1 = best_neighborhood(p1);
        auto best_center2 = best_neighborhood(p2);

        auto neighborhood1 = neighborhood(best_center1);
        auto neighborhood2 = neighborhood(best_center2);

        // Not required due to inherent order of values
        //std::sort(neighborhood1.begin(), neighborhood1.end(), cmp_linear);
        //std::sort(neighborhood2.begin(), neighborhood2.end(), cmp_linear);

        std::vector<tgt::ivec3> overlap;
        std::set_intersection (
                neighborhood1.begin(), neighborhood1.end(),
                neighborhood2.begin(), neighborhood2.end(),
                std::back_inserter(overlap), cmp_linear);

        std::sort(overlap.begin(), overlap.end(), [&](const tgt::ivec3& o1, const tgt::ivec3& o2) {
                return (tgt::distanceSq(p1, o1)-tgt::distanceSq(p2,o1)) < (tgt::distanceSq(p1, o2)-tgt::distanceSq(p2,o2));
                //return tgt::distanceSq(p1, o1) < tgt::distanceSq(p1, o2);
                //return dist_sq(p2, o1) > dist_sq(p2, o2);
                });

        auto o1begin = overlap.begin();
        auto o1end = overlap.begin()+overlap.size()/2;
        auto o2begin = o1end;
        auto o2end = overlap.end();

        std::vector<tgt::ivec3> o1(o1begin, o1end);
        std::vector<tgt::ivec3> o2(o2begin, o2end);
        std::sort(o1.begin(), o1.end(), cmp_linear);
        std::sort(o2.begin(), o2.end(), cmp_linear);

        std::vector<tgt::ivec3> n1final;
        std::vector<tgt::ivec3> n2final;

        std::set_difference(neighborhood1.begin(), neighborhood1.end(), o2.begin(), o2.end(), std::back_inserter(n1final), cmp_linear);
        std::set_difference(neighborhood2.begin(), neighborhood2.end(), o1.begin(), o1.end(), std::back_inserter(n2final), cmp_linear);

        auto mv1 = mean_and_var(n1final);
        auto mv2 = mean_and_var(n2final);
        float mean1 = mv1.first;
        float mean2 = mv2.first;
        float min_variance = 0.000001;
        float var1 = std::max(min_variance, mv1.second);
        float var2 = std::max(min_variance, mv2.second);

        float n1 = n1final.size();
        float n2 = n2final.size();

        float sn1 = var1/n1;
        float sn2 = var2/n2;

        float T_square = square(mean1 - mean2)/(sn1 + sn2);

        float m_star = square(sn1+sn2)/(square(sn1)/(n1-1.0f) + square(sn2)/(n2-1.0f));

        float m = std::round(m_star);

        float tpow = std::pow(1.0f+T_square/m, -0.5f*(m+1.0f));
        tgtAssert(!std::isnan(tpow) && std::isfinite(tpow), "Invalid tpow");

        //float gamma1 = tgamma((m+1.0f)*0.5f);
        //float gamma2 = tgamma(m*0.5f);

        //assert(!std::isnan(gamma1) && std::isfinite(gamma1));
        //assert(!std::isnan(gamma2) && std::isfinite(gamma2));

        //float w = gamma1/(std::sqrt(m)*gamma2*tpow);

        //The above sucks if m is too large, so we express the t-distribution pdf using the beta function
        float beta_term = boost::math::beta(0.5f, m*0.5f);
        tgtAssert(!std::isnan(beta_term) && std::isfinite(beta_term) && beta_term > 0, "Invalid beta");

        float w = tpow/(std::sqrt(m)*beta_term);


        tgtAssert(!std::isnan(w) && std::isfinite(w) && w >= 0, "Invalid weight");

        return w;
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
