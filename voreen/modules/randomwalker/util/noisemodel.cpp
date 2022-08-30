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
#include <numeric>

namespace voreen {

static bool cmp_linear(const tgt::ivec3& p1, const tgt::ivec3& p2) {
    return p1.z < p2.z ||
        p1.z == p2.z && p1.y < p2.y ||
        p1.z == p2.z && p1.y == p2.y && p1.x < p2.x;
}

static float square(float s) {
    return s*s;
}

static float variance_of(const std::vector<float>& vec, float mean) {

    float sq_sum = 0.0f;
    size_t n = vec.size();
#if defined(WIN32) && _MSC_VER < 1922
#pragma vector // MSVC equivalent of omp simd prior to MSVC 1922
#else
#ifdef VRN_MODULE_OPENMP
#pragma omp simd
#endif
#endif
    for(size_t i = 0; i<n; ++i) {
        sq_sum += square(mean-vec[i]);
    }
    float variance = sq_sum/(n-1);

    return variance;
};

static void process_non_overlap_single(const float* image, tgt::ivec3 imagedim, tgt::ivec3 begin, tgt::ivec3 end, tgt::ivec3 overlap_begin, tgt::ivec3 overlap_end, float*& output_cur) {
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
                }
            }
        }
    }
}


static std::pair<std::vector<float>, std::vector<float>> collect_neighborhoods(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent) {
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

    std::vector<float> n1final;
    n1final.resize(num_voxels_max);
    float* n1_begin = n1final.data();
    float* n1_cur = n1_begin;

    tgt::ivec3 from_b1_to_b2 = best_center2-best_center1;

    std::vector<float> n2final;
    n2final.resize(num_voxels_max);
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

                        size_t i2 = i2base - i;
                        tgtAssert(i2 < tgt::hmul(dim) ,"invalid index");
                        float val2 = image.voxel(i2);
                        *n2_cur = val2;
                        ++n2_cur;
                    }
                }
            }
        }
    } else {
        process_non_overlap_single(image.voxel(), dim, begin1, end1, overlap_begin, overlap_end, n1_cur);
        process_non_overlap_single(image.voxel(), dim, begin2, end2, overlap_begin, overlap_end, n2_cur);
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
        *n1_cur = o->first;
        ++n1_cur;
    }
    for(; o!=o2end; ++o) {
        *n2_cur = o->first;
        ++n2_cur;
    }

    size_t n1 = std::distance(n1_begin, n1_cur);
    size_t n2 = std::distance(n2_begin, n2_cur);

    n1final.resize(n1);
    n2final.resize(n2);

    return std::make_pair(n1final, n2final);
}

float evalTTest(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent) {
    auto neighborhoods = collect_neighborhoods(image, best_centers, voxel, neighbor, filter_extent);
    auto neigh1 = neighborhoods.first;
    auto neigh2 = neighborhoods.second;


    float n1 = neigh1.size();
    float n2 = neigh2.size();

    float sum1 = std::accumulate(neigh1.begin(), neigh1.end(), 0.0f);
    float sum2 = std::accumulate(neigh2.begin(), neigh2.end(), 0.0f);

    float mean1 = sum1/n1;
    float mean2 = sum2/n2;

    float min_variance = 0.000001;
    float var1 = std::max(min_variance, variance_of(neigh1, mean1));
    float var2 = std::max(min_variance, variance_of(neigh2, mean2));

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

float evalVariableGaussian(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent) {
    auto neighborhoods = collect_neighborhoods(image, best_centers, voxel, neighbor, filter_extent);
    auto neighborhood1 = neighborhoods.first;
    auto neighborhood2 = neighborhoods.second;


    if(neighborhood1.size() > neighborhood2.size()) {
        neighborhood1.resize(neighborhood2.size());
    } else if(neighborhood1.size() < neighborhood2.size()) {
        neighborhood2.resize(neighborhood1.size());
    }
    size_t n = neighborhood1.size();
    assert(n == neighborhood2.size());

    float sum1 = std::accumulate(neighborhood1.begin(), neighborhood1.end(), 0.0f);
    float sum2 = std::accumulate(neighborhood2.begin(), neighborhood2.end(), 0.0f);

    float mean1 = sum1/n;
    float mean2 = sum2/n;

    float var1 = variance_of(neighborhood1, mean1);
    float var2 = variance_of(neighborhood2, mean2);


    float nom = std::sqrt(var1*var2);
    float denom = (var1+var2)*0.5 + square((mean1-mean2)*0.5);

    if(denom == 0.0f) {
        return 1.0f;
    }

    float quotient = nom/denom;

    assert(n>=4);
    float exponent = (n-3.0)/2;
    float w = std::pow(quotient, exponent);


    assert(!std::isnan(w) && std::isfinite(w) && w >= 0);

    return w;
}

std::vector<GaussianParametersVariableSigma> GaussianParametersVariableSigma::create(const VolumeAtomic<float>& image, int filter_extent) {
    VolumeAtomic<float> mean = meanFilter(image, filter_extent);
    VolumeAtomic<float> variance = variances(image, mean, filter_extent);

    std::vector<GaussianParametersVariableSigma> mean_mul_add;
    size_t n = mean.getNumVoxels();
    mean_mul_add.reserve(n);

    for(int i = 0; i<n; ++i) {
        GaussianParametersVariableSigma mma;
        mma.mean = mean.voxel(i);
        float var = std::max(variance.voxel(i), 0.000001f);
        mma.add = 0.5*std::log(var);
        mma.mul = 0.5/var;
        mean_mul_add.push_back(mma);
    }

    return mean_mul_add;
}
float GaussianParametersVariableSigma::fit(GaussianParametersVariableSigma params, float f) {
    // Instead of computing the exp, we maximize the log
    return - (square(f-params.mean) * params.mul + params.add);
}

float bhattacharyyaConstGaussian(float mean1, float mean2, float variance_inv, float n) {
    float diff = mean1-mean2;

    float coeff = diff*diff * variance_inv / 8.0f * n;

    float w = std::exp(-coeff);

    assert(!std::isnan(w) && std::isfinite(w) && w >= 0);

    return w;
}

float evalConstGaussian(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, float variance_inv, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent) {
    auto neighborhoods = collect_neighborhoods(image, best_centers, voxel, neighbor, filter_extent);
    auto neighborhood1 = neighborhoods.first;
    auto neighborhood2 = neighborhoods.second;


    if(neighborhood1.size() > neighborhood2.size()) {
        neighborhood1.resize(neighborhood2.size());
    } else if(neighborhood1.size() < neighborhood2.size()) {
        neighborhood2.resize(neighborhood1.size());
    }
    size_t n = neighborhood1.size();
    assert(n == neighborhood2.size());

    float sum1 = std::accumulate(neighborhood1.begin(), neighborhood1.end(), 0.0f);
    float sum2 = std::accumulate(neighborhood2.begin(), neighborhood2.end(), 0.0f);

    float mean1 = sum1/n;
    float mean2 = sum2/n;

    return bhattacharyyaConstGaussian(mean1, mean2, variance_inv, n);
}

float bhattacharyyaPoisson(float sum1, float sum2) {
    const float APPROX_THRESHOLD = 1000;

    float mean = (sum1+sum2)*0.5f;
    float exponent;

    if(sum1 < APPROX_THRESHOLD && sum2 < APPROX_THRESHOLD) {
        float t1 = std::lgamma(mean+1.0f);
        float t2 = std::lgamma(sum1+1.0f);
        float t3 = std::lgamma(sum2+1.0f);
        float t4 = (t2+t3)*0.5f;

        exponent =  t1 - t4;
    } else {
        // Use approximation if values are getting too big for lgamma
        float s1 = std::sqrt(sum1);
        float s2 = std::sqrt(sum2);
        float diff = s1-s2;
        exponent = -0.5f * (diff*diff);
    }

    float w = std::exp(exponent);
    assert(!std::isnan(w) && std::isfinite(w) && w >= 0);
    return w;
}

float evalPoisson(const VolumeAtomic<float>& image, const VolumeAtomic<tgt::ivec3>& best_centers, tgt::svec3 voxel, tgt::svec3 neighbor, int filter_extent) {
    auto neighborhoods = collect_neighborhoods(image, best_centers, voxel, neighbor, filter_extent);
    auto neighborhood1 = neighborhoods.first;
    auto neighborhood2 = neighborhoods.second;


    if(neighborhood1.size() > neighborhood2.size()) {
        neighborhood1.resize(neighborhood2.size());
    } else if(neighborhood1.size() < neighborhood2.size()) {
        neighborhood2.resize(neighborhood1.size());
    }
    size_t n = neighborhood1.size();
    assert(n == neighborhood2.size());

    float sum1 = std::accumulate(neighborhood1.begin(), neighborhood1.end(), 0.0f);
    float sum2 = std::accumulate(neighborhood2.begin(), neighborhood2.end(), 0.0f);

    return bhattacharyyaPoisson(sum1, sum2);
}

}
