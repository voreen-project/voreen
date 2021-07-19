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

#include "noisemodel.h"

namespace voreen {

static bool cmp_linear(const tgt::ivec3& p1, const tgt::ivec3& p2) {
    return p1.z < p2.z ||
        p1.z == p2.z && p1.y < p2.y ||
        p1.z == p2.z && p1.y == p2.y && p1.x < p2.x;
}

static float square(float s) {
    return s*s;
}

static float variance_of(const float* begin, const float* end, float mean) {

    float sq_sum = 0.0f;
    size_t n = std::distance(begin, end);
#if defined(WIN32) && _MSC_VER < 1922
#pragma vector // MSVC equivalent of omp simd prior to MSVC 1922
#else
#ifdef VRN_MODULE_OPENMP
#pragma omp simd
#endif
#endif
    for(size_t i = 0; i<n; ++i) {
        sq_sum += square(mean-begin[i]);
    }
    float variance = sq_sum/(n-1);

    return variance;
};

static void process_non_overlap_single(const float* image, tgt::ivec3 imagedim, tgt::ivec3 begin, tgt::ivec3 end, tgt::ivec3 overlap_begin, tgt::ivec3 overlap_end, float*& output_cur, float& sum) {
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

VolumeAtomic<tgt::ivec3> findBestCenters(const VolumeAtomic<float>& image, const VolumeAtomic<float>& mean, const VolumeAtomic<float>& variance, int filter_extent) {
    struct MeanMulAdd {
        float mean;
        float mul;
        float add;
    };
    tgt::ivec3 dim = image.getDimensions();

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

    VolumeAtomic<tgt::ivec3> best_centers(dim);

    VRN_FOR_EACH_VOXEL(p, tgt::ivec3::zero, dim) {
        float f = image.voxel(p);

        const tgt::ivec3 extent(filter_extent);

        tgt::ivec3 begin = tgt::max(tgt::ivec3::zero, p - extent);
        tgt::ivec3 end = tgt::min(dim, p + tgt::ivec3::one + extent);

        float min = std::numeric_limits<float>::infinity();
        size_t argmin_index = 0;

        VRN_FOR_EACH_VOXEL_INDEX(l, begin, end, dim, {
            //float pdf_val = gaussian_pdf(f, mean.voxel(n), variance.voxel(n));
            // Instead of evaluating and maximizing the gaussian pdf (which
            // is expensive) we minimize the negative log instead.
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
        best_centers.voxel(p) = argmin;
    }

    return best_centers;
}

float evalTTest(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent) {
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

    auto best_center1 = best_centers.voxel(p1);
    auto best_center2 = best_centers.voxel(p2);

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

    // In most cases (except when at the border of the volume) both
    // neighborhoods are placed with rotational symmetry at opposite sides
    // of the overlap region. In this case we can sample from both regions
    // simultaneously with some index trickery: We iterate normally over
    // region 1, but do so in reverse for region 2. Then, difference
    // between the index of the first voxel of region 1 and the current
    // once is the same as the difference between the current voxel in
    // region 2 and the last voxel of the region (end-1).
    if(end1-begin1 == end2-begin2) {
        size_t b1index = cubicCoordToLinear(begin1, dim);
        size_t e2index = cubicCoordToLinear(end2-tgt::ivec3::one, dim);
        size_t i2base = b1index + e2index;
        size_t neighborhood_counter = 0;
        {
            size_t zBegin = begin1.z*sliceSize;
            size_t zEnd = end1.z*sliceSize;
            for (size_t z = begin1.z; z < end1.z; ++z) {

                size_t zIndex = sliceSize * z;
                bool zIn = overlap_begin.z <= z && z < overlap_end.z;

                for (size_t y = begin1.y; y < end1.y; ++y) {

                    size_t yIndex = zIndex + lineSize * y;
                    bool yIn = zIn && overlap_begin.y <= y && y < overlap_end.y;

                    size_t non_overlap_begin;
                    size_t non_overlap_end;
                    if(yIn) {
                        for (size_t x = overlap_begin.x; x < overlap_end.x; ++x) {
                            size_t i = yIndex + x;
                            float val = image.voxel(i);

                            tgt::ivec3 n(x,y,z);
                            int along_axis = tgt::dot(n, from_b1_to_b2);
                            along_axis = (along_axis << 16) + neighborhood_counter;
                            *o_cur = std::make_pair(val,along_axis);
                            ++o_cur;
                            ++neighborhood_counter;
                        }
                        bool overlap_at_begin = overlap_begin.x == begin1.x;
                        non_overlap_begin = overlap_at_begin ? overlap_end.x : begin1.x;
                        non_overlap_end = overlap_at_begin ? end1.x : overlap_begin.x;
                    } else {
                        non_overlap_begin = begin1.x;
                        non_overlap_end = end1.x;
                    }

                    for (size_t x = non_overlap_begin; x < non_overlap_end; ++x) {
                        size_t i = yIndex + x;
                        float val1 = image.voxel(i);
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
    } else {
        process_non_overlap_single(image.voxel(), dim, begin1, end1, overlap_begin, overlap_end, n1_cur, sum1);
        process_non_overlap_single(image.voxel(), dim, begin2, end2, overlap_begin, overlap_end, n2_cur, sum2);
        size_t neighborhood_counter = 0;
        for (size_t z = overlap_begin.z; z < overlap_end.z; ++z) {

            size_t zIndex = sliceSize * z;

            for (size_t y = overlap_begin.y; y < overlap_end.y; ++y) {

                size_t yIndex = zIndex + lineSize * y;

                for (size_t x = overlap_begin.x; x < overlap_end.x; ++x) {
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
}
