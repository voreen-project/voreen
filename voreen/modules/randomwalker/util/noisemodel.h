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

//inline float gaussian_pdf_exp(float val, float mean, float variance) {
//    return std::exp(-0.5f/variance * square(val-mean));
//}
//
//inline float gaussian_pdf(float val, float mean, float variance) {
//    return gaussian_pdf_exp(val, mean, variance)/std::sqrt(2.0f * 3.141592654f* variance);
//}

inline float variance_of(const float* begin, const float* end, float mean) {

    float sq_sum = 0.0f;
    size_t n = std::distance(begin, end);
#ifdef VRN_MODULE_OPENMP
#pragma omp simd
#endif
    for(size_t i = 0; i<n; ++i) {
        sq_sum += square(mean-begin[i]);
    }
    float variance = sq_sum/(n-1);

    return variance;
};

#define VRN_FOR_EACH_VOXEL_INDEX(INDEX, BEGIN, END, SIZE, BODY) \
    size_t vrn_fevi_z_slice_size = (SIZE).x*(SIZE).y;\
    size_t vrn_fevi_z_begin = (BEGIN).z*vrn_fevi_z_slice_size;\
    size_t vrn_fevi_z_end = (END).z*vrn_fevi_z_slice_size;\
    for (size_t vrn_fevi_z = vrn_fevi_z_begin; vrn_fevi_z < vrn_fevi_z_end ; vrn_fevi_z += vrn_fevi_z_slice_size) {\
        size_t vrn_fevi_y_begin = vrn_fevi_z + (BEGIN).y*(SIZE).x;\
        size_t vrn_fevi_y_end = vrn_fevi_z + (END).y*(SIZE).x;\
        for (size_t vrn_fevi_y = vrn_fevi_y_begin; vrn_fevi_y < vrn_fevi_y_end ; vrn_fevi_y += (SIZE).x){\
            size_t vrn_fevi_x_begin = vrn_fevi_y + (BEGIN).x;\
            size_t vrn_fevi_x_end = vrn_fevi_y + (END).x;\
            for (size_t INDEX = vrn_fevi_x_begin; INDEX < vrn_fevi_x_end; ++INDEX){\
                BODY\
            }\
        }\
    }\

inline void process_non_overlap_single(const float* image, tgt::ivec3 imagedim, tgt::ivec3 begin, tgt::ivec3 end, tgt::ivec3 overlap_begin, tgt::ivec3 overlap_end, float*& output_cur, float& sum) {
#if 0
    for(int dim = 2; dim >= 0; --dim) {
        tgt::ivec3 part_begin = begin;
        tgt::ivec3 part_end = end;
        if(begin[dim] == overlap_begin[dim]) {
            part_begin[dim] = overlap_end[dim];
        }
        if(end[dim] == overlap_end[dim]) {
            part_end[dim] = overlap_begin[dim];
        }

        VRN_FOR_EACH_VOXEL_INDEX(l, part_begin, part_end, imagedim, {
            float val = image[l];
            *output_cur = val;
            ++output_cur;
            sum += val;
        })

        begin[dim] = overlap_begin[dim];
        end[dim] = overlap_end[dim];
    }
#else
        {
            size_t sliceSize = imagedim.y * imagedim.x;
            size_t lineSize = imagedim.x;
            for (size_t z = begin.z; z < static_cast<int>(end.z); ++z) {

                size_t zIndex = sliceSize * z;
                bool zIn = overlap_begin.z <= z && z < overlap_end.z;

                for (size_t y = begin.y; y < static_cast<int>(end.y); ++y) {

                    size_t yIndex = zIndex + lineSize * y;
                    bool yIn = zIn && overlap_begin.y <= y && y < overlap_end.y;

                    size_t iBegin = yIndex + begin.x;
                    size_t iEnd = yIndex + end.x;

                    size_t overlap_i_begin = yIndex + overlap_begin.x;
                    size_t overlap_i_end = yIndex + overlap_end.x;
                    for (size_t i = iBegin; i < iEnd; ++i) {
                        bool xIn = yIn && overlap_i_begin <= i && i < overlap_i_end;

                        if(!xIn) {
                            float val = image[i];
                            *output_cur = val;
                            ++output_cur;
                            sum += val;
                        }
                    }
                }
            }
        }
#endif
}

// Parameter estimation according to
//
// A. Bian, X Jiang: T-Test Based Adaptive Random Walk Segmentation Under Multiplicative Speckle Noise Model
// https://link.springer.com/chapter/10.1007%2F978-3-319-54427-4_41
template<int filter_extent>
struct RWNoiseModelTTest {
    struct MeanMulAdd {
        float mean;
        float mul;
        float add;
    };

    VolumeAtomic<float> image;
    std::vector<MeanMulAdd> mean_mul_add;

    RWNoiseModelTTest(RWNoiseModelTTest&&) = default;

    static RWNoiseModelTTest prepare(const VolumeAtomic<float>& vol, RealWorldMapping rwm) {

        VolumeAtomic<float> image = vol.copy();
        VolumeAtomic<float> mean = meanFilter<filter_extent>(image);
        VolumeAtomic<float> variance = variances<filter_extent>(image, mean);

        std::vector<MeanMulAdd> mean_mul_add;
        size_t n = mean.getNumVoxels();
        mean_mul_add.reserve(n);

        for(int i = 0; i<n; ++i) {
            MeanMulAdd mma;
            mma.mean = mean.voxel(i);
            float var = std::max(variance.voxel(i), 0.000001f);
            mma.add = 0.5*std::log(var);
            mma.mul = 0.5/var;
            mean_mul_add.push_back(mma);
        }

        return RWNoiseModelTTest {
            std::move(image),
            std::move(mean_mul_add),
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

        const int filter_size=2*filter_extent+1;
        const int num_voxels_max = filter_size*filter_size*filter_size;

        const tgt::ivec3 dim(image.getDimensions());

        const size_t sliceSize = dim.y*dim.x;
        const size_t lineSize = dim.x;

        const tgt::ivec3 extent(filter_extent);

        auto best_neighborhood = [&] (tgt::ivec3 p) -> tgt::ivec3 {
            float f = image.voxel(p);

            tgt::ivec3 begin = tgt::max(tgt::ivec3::zero, p - extent);
            tgt::ivec3 end = tgt::min(dim, p + tgt::ivec3::one + extent);

            float min = std::numeric_limits<float>::infinity();
            size_t argmin_index = 0;

            VRN_FOR_EACH_VOXEL_INDEX(l, begin, end, dim, {
                //float pdf_val = gaussian_pdf(f, mean.voxel(n), variance.voxel(n));
                // Instead of evaluating and maximizing the gaussian pdf (which
                // is expensive) we minimize the log instead.
                const auto& mma = mean_mul_add[l];
                float val = square(f-mma.mean) * mma.mul + mma.add;
                tgtAssert(std::isfinite(val) && !std::isnan(val), "invalid val");
                if(min > val) {
                    min = val;
                    argmin_index = l;
                }
            })

            tgt::ivec3 argmin = linearCoordToCubic(argmin_index, dim);
            tgtAssert(tgt::max(tgt::abs(p - argmin)) <= filter_extent, "invalid pos");

            return argmin;
        };

        auto best_center1 = best_neighborhood(p1);
        auto best_center2 = best_neighborhood(p2);

        tgt::ivec3 begin1 = tgt::max(tgt::ivec3::zero, best_center1 - extent);
        tgt::ivec3 end1 = tgt::min(dim, best_center1 + tgt::ivec3::one + extent);

        tgt::ivec3 begin2 = tgt::max(tgt::ivec3::zero, best_center2 - extent);
        tgt::ivec3 end2 = tgt::min(dim, best_center2 + tgt::ivec3::one + extent);

        tgt::ivec3 overlap_begin = tgt::max(begin1, begin2);
        tgt::ivec3 overlap_end = tgt::max(tgt::min(end1, end2), overlap_begin);

        float sum1 = 0.0f;
        std::vector<float> n1final;
        n1final.reserve(num_voxels_max);
        float* n1_begin = n1final.data();
        float* n1_cur = n1_begin;

        tgt::ivec3 from_b1_to_b2 = best_center2-best_center1;

        float sum2 = 0.0f;
        std::vector<float> n2final;
        n2final.reserve(num_voxels_max);
        float* n2_begin = n2final.data();
        float* n2_cur = n2_begin;

        size_t overlap_size = tgt::hmul(overlap_end - overlap_begin);
        std::vector<std::pair<float, int>> overlap;
        overlap.reserve(overlap_size);
        auto* o_begin = overlap.data();
        auto* o_cur = o_begin;

        if(end1-begin1 == end2-begin2) {
            size_t b1index = cubicCoordToLinear(begin1, dim);
            size_t e2index = cubicCoordToLinear(end2-tgt::ivec3::one, dim);
            size_t i2base = b1index + e2index;
            size_t neighborhood_counter = 0;
            {
                size_t zBegin = begin1.z*sliceSize;
                size_t zEnd = end1.z*sliceSize;
                for (size_t z = begin1.z; z < static_cast<int>(end1.z); ++z) {

                    size_t zIndex = sliceSize * z;
                    bool zIn = overlap_begin.z <= z && z < overlap_end.z;

                    for (size_t y = begin1.y; y < static_cast<int>(end1.y); ++y) {

                        size_t yIndex = zIndex + lineSize * y;
                        bool yIn = zIn && overlap_begin.y <= y && y < overlap_end.y;

                        for (size_t x = begin1.x; x < static_cast<int>(end1.x); ++x) {
                            bool xIn = yIn && overlap_begin.x <= x && x < overlap_end.x;

                            size_t i = yIndex + x;
                            float val1 = image.voxel(i);

                            if(xIn) {
                                tgt::ivec3 n(x,y,z);
                                int along_axis = tgt::dot(n, from_b1_to_b2);
                                along_axis = (along_axis << 16) + neighborhood_counter;
                                *o_cur = std::make_pair(val1,along_axis);
                                ++o_cur;
                                ++neighborhood_counter;
                            } else {
                                *n1_cur = val1;
                                ++n1_cur;
                                sum1 += val1;

                                size_t i2 = i2base - i;
                                tgtAssert(i2 < tgt::hmul(dim) ,"invalid index");
                                float val2 = image.voxel(i2);
                                *n2_cur = val2;
                                ++n2_cur;
                                sum2 += val2;
                            }
                        }
                    }
                }
            }
        } else {
            process_non_overlap_single(image.voxel(), dim, begin1, end1, overlap_begin, overlap_end, n1_cur, sum1);
            process_non_overlap_single(image.voxel(), dim, begin2, end2, overlap_begin, overlap_end, n2_cur, sum2);
            size_t neighborhood_counter = 0;
            for (size_t z = overlap_begin.z; z < static_cast<int>(overlap_end.z); ++z) {

                size_t zIndex = sliceSize * z;

                for (size_t y = overlap_begin.y; y < static_cast<int>(overlap_end.y); ++y) {

                    size_t yIndex = zIndex + lineSize * y;

                    for (size_t x = overlap_begin.x; x < static_cast<int>(overlap_end.x); ++x) {
                        size_t i = yIndex + x;
                        float val = image.voxel(i);

                        tgt::ivec3 n(x,y,z);
                        int along_axis = tgt::dot(n, from_b1_to_b2);
                        along_axis = (along_axis << 16) + neighborhood_counter;
                        *o_cur = std::make_pair(val,along_axis);
                        ++o_cur;
                        ++neighborhood_counter;
                    }
                }
            }
        }

        auto o1begin = o_begin;
        auto o1end = o_begin+std::distance(o_begin, o_cur)/2;
        auto o2end = o_cur;

        std::nth_element(o1begin, o1end, o2end, [](const std::pair<float,int>& o1, const std::pair<float,int>& o2) {
                return o1.second < o2.second;
                });


        auto o = o1begin;
        for(; o!=o1end; ++o) {
            sum1 += o->first;
            *n1_cur = o->first;
            ++n1_cur;
        }
        for(; o!=o2end; ++o) {
            sum2 += o->first;
            *n2_cur = o->first;
            ++n2_cur;
        }
        float n1 = std::distance(n1_begin, n1_cur);
        float n2 = std::distance(n2_begin, n2_cur);

        float mean1 = sum1/n1;
        float mean2 = sum2/n2;

        float min_variance = 0.000001;
        float var1 = std::max(min_variance, variance_of(n1_begin, n1_cur, mean1));
        float var2 = std::max(min_variance, variance_of(n2_begin, n2_cur, mean2));

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
