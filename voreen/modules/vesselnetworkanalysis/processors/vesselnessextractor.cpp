/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "vesselnessextractor.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "modules/bigdataimageprocessing/volumefiltering/gaussianfilter.h"
#include "tgt/tgt_math.h"

#include <chrono>
#include <complex>

namespace {

using namespace voreen;

// ------------------------------------------------------------------------------------------
// Implementation of helper classes----------------------------------------------------------
// ------------------------------------------------------------------------------------------

// Space efficient symmetric 3 dimensional matrix
struct SymMat3 {

    SymMat3(float val = 0)
        : xx(val)
        , xy(val)
        , xz(val)
        , yy(val)
        , yz(val)
        , zz(val)
    {
    }

    SymMat3(float xx, float xy, float xz, float yy, float yz, float zz)
        : xx(xx)
        , xy(xy)
        , xz(xz)
        , yy(yy)
        , yz(yz)
        , zz(zz)
    {
    }

    float xx;
    float xy;
    float xz;
    float yy;
    float yz;
    float zz;
    SymMat3 operator*(const SymMat3& r) const;
    SymMat3 minusUnitTimes(float a) const;
    // Get 
    void getEigenValues(float& l1, float& l2, float& l3) const;

    // Try to get the eigenvector to the eigenvalue not supplied as a parameter.
    // Turns out, this is rather unreliable, so we use inverse iteration instead
    //tgt::vec3 getEigenvectorToTheOtherEigenValue(float notThisOne, float alsoNotThisOne) const;

    // Create a "regular" tgt matrix
    tgt::mat3 toTgtMat() const;
};

// A separable kernel that is used to convolute the image in VesselnessFeatureExtractor
// Currently, only 2nd order derivation kernels are available via static functions.
class SeparableKernel {
public:
    typedef std::function<float(float)> KernelFunc;

    // 2nd order derivation kernels (gaussian wavelets)
    static SeparableKernel deriveXX(tgt::ivec3 extent, tgt::vec3 stddev);
    static SeparableKernel deriveXY(tgt::ivec3 extent, tgt::vec3 stddev);
    static SeparableKernel deriveXZ(tgt::ivec3 extent, tgt::vec3 stddev);
    static SeparableKernel deriveYY(tgt::ivec3 extent, tgt::vec3 stddev);
    static SeparableKernel deriveYZ(tgt::ivec3 extent, tgt::vec3 stddev);
    static SeparableKernel deriveZZ(tgt::ivec3 extent, tgt::vec3 stddev);

    // Value access
    float xKernelAt(int p) const;
    float yKernelAt(int p) const;
    float zKernelAt(int p) const;


private:
    static std::vector<float> gaussKernel(float stddev, int extent);
    static std::vector<float> gaussFirstDerivativeKernel(float stddev, int extent);
    static std::vector<float> gaussSecondDerivativeKernel(float stddev, int extent);

    SeparableKernel(std::vector<float> xKernel, std::vector<float> yKernel, std::vector<float> zKernel, tgt::ivec3 extent);

    tgt::ivec3 extent_;
    std::vector<float> xKernel_;
    std::vector<float> yKernel_;
    std::vector<float> zKernel_;
};

// Extracts vesselness features (v) from the image, as well as the eigenvector to the smallest eigenvalue (e_1).
// The result is stored as a four component vector (e_1, v)
class VesselnessFeatureExtractor : public VolumeFilter {
public:

    VesselnessFeatureExtractor(float alpha, float beta, float c, const tgt::ivec3& extent, const tgt::vec3 standardDeviation, float scale);
    virtual ~VesselnessFeatureExtractor();

    float computeVesselness(const SymMat3& hessian) const;

    virtual std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const;
    virtual int zExtent() const;
    virtual size_t getNumInputChannels() const;
    virtual size_t getNumOutputChannels() const;
    virtual SliceReaderMetaData getMetaData(const SliceReaderMetaData& base) const;

    tgt::ivec3 extent_;
    //float two_alpha_squared_;
    //float two_beta_squared_;
    //float two_c_squared_;

    SeparableKernel xxKernel_;
    SeparableKernel xyKernel_;
    SeparableKernel xzKernel_;
    SeparableKernel yyKernel_;
    SeparableKernel yzKernel_;
    SeparableKernel zzKernel_;

    float scale_;
};


// ------------------------------------------------------------------------------------------
// Implementation of helper classes----------------------------------------------------------
// ------------------------------------------------------------------------------------------

SeparableKernel::SeparableKernel(std::vector<float> xKernel, std::vector<float> yKernel, std::vector<float> zKernel, tgt::ivec3 extent)
    : extent_(extent)
    , xKernel_(xKernel)
    , yKernel_(yKernel)
    , zKernel_(zKernel)
{
    tgtAssert(static_cast<int>(xKernel_.size()) == 2*extent_.x+1, "Invalid kernel size");
    tgtAssert(static_cast<int>(yKernel_.size()) == 2*extent_.y+1, "Invalid kernel size");
    tgtAssert(static_cast<int>(zKernel_.size()) == 2*extent_.z+1, "Invalid kernel size");

    //std::cout << "[";
    //for(auto val : xKernel_) {
    //    std::cout << val << ", ";
    //}
    //std::cout << "]\n";
    //std::cout.flush();
}
static SeparableKernel::KernelFunc gauss(float stdev) {
    return [stdev] (float x) {
        return std::exp(-0.5f * x*x / (stdev*stdev))/(std::sqrt(2*tgt::PIf)*stdev);
    };
}

static SeparableKernel::KernelFunc gaussFirstDerivative(float stdev) {
    return [stdev] (float x) {
        return -std::exp(-0.5f * x*x / (stdev*stdev))*x/(std::sqrt(2*tgt::PIf)*stdev*stdev*stdev);
    };
}

static SeparableKernel::KernelFunc gaussSecondDerivative(float stdev) {
    return [stdev] (float x) {
        return std::exp(-0.5f * x*x / (stdev*stdev))*(x*x-stdev*stdev)/(std::sqrt(2*tgt::PIf)*stdev*stdev*stdev*stdev*stdev);
    };
}
//Is also suitable for second derivative etc.
static tgt::ivec3 suitableExtent(tgt::vec3 stddev) {
    return tgt::ivec3(tgt::round(3.0f*stddev));
}


// Get a smooth sample by integration in the neighborhood
float smoothSample(SeparableKernel::KernelFunc f, int i) {
    int radius = 5;
    float sum = 0.0f;
    for(int d = -radius; d <= radius; ++d) {
        float df = static_cast<float>(d)/radius;
        float samplePos = static_cast<float>(i)-0.5f*df;
        sum += f(samplePos);
    }
    return sum/(2*radius + 1);
}

std::vector<float> SeparableKernel::gaussKernel(float stddev, int extent) {
    auto func = gauss(stddev);
    std::vector<float> kernel(2*extent + 1);
    float sum = 0;
    for(int i = -extent; i <= extent; ++i) {
        float val = smoothSample(func, i);
        sum += val;
        kernel[i + extent] = val;
    }
    // Normalize so that all values sum to 1
    for(int i = -extent; i <= extent; ++i) {
        kernel[i + extent] /= sum;
    }
    return kernel;
}

std::vector<float> SeparableKernel::gaussFirstDerivativeKernel(float stddev, int extent) {
    auto func = gaussFirstDerivative(stddev);
    std::vector<float> kernel(2*extent + 1);
    float positiveSum = 0;
    float negativeSum = 0;
    for(int i = -extent; i <= extent; ++i) {
        float val = smoothSample(func, i);
        if(val > 0) {
            positiveSum += val;
        } else {
            negativeSum += val;
        }
        kernel[i + extent] = val;
    }
    float positiveNormalizationFactor = std::abs(gauss(stddev)(0)/positiveSum);
    float negativeNormalizationFactor = std::abs(gauss(stddev)(0)/negativeSum);
    // Normalize so that all values sum to 0
    for(int i = -extent; i <= extent; ++i) {
        if(kernel[i + extent] > 0) {
            kernel[i + extent] *= positiveNormalizationFactor;
        } else {
            kernel[i + extent] *= negativeNormalizationFactor;
        }
    }
    return kernel;
}

std::vector<float> SeparableKernel::gaussSecondDerivativeKernel(float stddev, int extent) {
    auto func = gaussSecondDerivative(stddev);
    std::vector<float> kernel(2*extent + 1);
    float positiveSum = 0;
    float negativeSum = 0;
    for(int i = -extent; i <= extent; ++i) {
        float val = smoothSample(func, i);
        if(val > 0) {
            positiveSum += val;
        } else {
            negativeSum += val;
        }
        kernel[i + extent] = val;
    }
    float positiveNormalizationFactor = std::abs(2*gaussFirstDerivative(stddev)(stddev)/positiveSum);
    float negativeNormalizationFactor = std::abs(2*gaussFirstDerivative(stddev)(stddev)/negativeSum);
    // Normalize so that all values sum to 0
    for(int i = -extent; i <= extent; ++i) {
        if(kernel[i + extent] > 0) {
            kernel[i + extent] *= positiveNormalizationFactor;
        } else {
            kernel[i + extent] *= negativeNormalizationFactor;
        }
    }
    return kernel;
}

SeparableKernel SeparableKernel::deriveXX(tgt::ivec3 extent, tgt::vec3 stddev) {
    return SeparableKernel(
        gaussSecondDerivativeKernel(stddev.x, extent.x),
        gaussKernel(stddev.y, extent.y),
        gaussKernel(stddev.z, extent.z),
        extent
    );
}
SeparableKernel SeparableKernel::deriveXY(tgt::ivec3 extent, tgt::vec3 stddev) {
    return SeparableKernel(
        gaussFirstDerivativeKernel(stddev.x, extent.x),
        gaussFirstDerivativeKernel(stddev.y, extent.y),
        gaussKernel(stddev.z, extent.z),
        extent
    );
}
SeparableKernel SeparableKernel::deriveXZ(tgt::ivec3 extent, tgt::vec3 stddev) {
    return SeparableKernel(
        gaussFirstDerivativeKernel(stddev.x, extent.x),
        gaussKernel(stddev.y, extent.y),
        gaussFirstDerivativeKernel(stddev.z, extent.z),
        extent
    );
}
SeparableKernel SeparableKernel::deriveYY(tgt::ivec3 extent, tgt::vec3 stddev) {
    return SeparableKernel(
        gaussKernel(stddev.x, extent.x),
        gaussSecondDerivativeKernel(stddev.y, extent.y),
        gaussKernel(stddev.z, extent.z),
        extent
    );
}
SeparableKernel SeparableKernel::deriveYZ(tgt::ivec3 extent, tgt::vec3 stddev) {
    return SeparableKernel(
        gaussKernel(stddev.x, extent.x),
        gaussFirstDerivativeKernel(stddev.y, extent.y),
        gaussFirstDerivativeKernel(stddev.z, extent.z),
        extent
    );
}
SeparableKernel SeparableKernel::deriveZZ(tgt::ivec3 extent, tgt::vec3 stddev) {
    return SeparableKernel(
        gaussKernel(stddev.x, extent.x),
        gaussKernel(stddev.y, extent.y),
        gaussSecondDerivativeKernel(stddev.z, extent.z),
        extent
    );
}


float SeparableKernel::xKernelAt(int p) const {
    tgtAssert(-extent_.x <= p && p <= extent_.x, "Invalid central position");
    return xKernel_[p+extent_.x];
}

float SeparableKernel::yKernelAt(int p) const {
    tgtAssert(-extent_.y <= p && p <= extent_.y, "Invalid central position");
    return yKernel_[p+extent_.y];
}

float SeparableKernel::zKernelAt(int p) const {
    tgtAssert(-extent_.z <= p && p <= extent_.z, "Invalid central position");
    return zKernel_[p+extent_.z];
}


SymMat3 SymMat3::operator*(const SymMat3& r) const {
    SymMat3 out;
    out.xx = xx*r.xx + xy*r.xy + xz*r.xz;
    out.xy = xx*r.xy + xy*r.yy + xz*r.yz;
    out.xz = xx*r.xz + xy*r.yz + xz*r.zz;

    out.yy = xy*r.xy + yy*r.yy + yz*r.yz;
    out.yz = xy*r.xz + yy*r.yz + yz*r.zz;

    out.zz = xz*r.xz + yz*r.yz + zz*r.zz;
    return out;
}

SymMat3 SymMat3::minusUnitTimes(float a) const {
    SymMat3 out;
    out.xx = xx - a;
    out.xy = xy;
    out.xz = xz;

    out.yy = xy - a;
    out.yz = yz;

    out.zz = zz - a;
    return out;
}

void SymMat3::getEigenValues(float& l1, float& l2, float& l3) const {
    // Set up characteristic equation:   det( A - lambda I ) = 0
    //    as a cubic in lambda:  a.lambda^3 + b.lambda^2 + c.lambda + d = 0
    std::complex<double> a = -1.0;                 // -1
    std::complex<double> b = xx + yy + zz;         // trace
    std::complex<double> c = yz * yz - yy * zz     // -sum of diagonal minors
            + xz * xz - zz * xx
            + xy * xy - xx * yy;
    std::complex<double> d = xx*yy*zz + 2*xy*yz*xz - xx*yz*yz - yy*xz*xz - zz*xy*xy;

    // Solve cubic by Cardano's method (easier in complex numbers!)
    std::complex<double> p = ( b * b - 3.0 * a * c ) / ( 9.0 * a * a );
    std::complex<double> q = ( 9.0 * a * b * c - 27.0 * a * a * d - 2.0 * b * b * b ) / ( 54.0 * a * a * a );
    std::complex<double> delta = q * q - p * p * p;
    std::complex<double> deltaSqrt = sqrt(delta);
    std::complex<double> g1 = pow( q + deltaSqrt, 1.0 / 3.0 );     // warning: exponents and sqrt of complex
    std::complex<double> g2 = pow( q - deltaSqrt, 1.0 / 3.0 );
    std::complex<double> offset = -b / ( 3.0 * a );
    std::complex<double> omega = std::complex<double>( -0.5, 0.5 * sqrt( 3.0 ) );     // complex cube root of unity
    std::complex<double> omega2 = omega * omega;

    std::complex<double> cl1 = g1          + g2          + offset;
    std::complex<double> cl2 = g1 * omega  + g2 * omega2 + offset;
    std::complex<double> cl3 = g1 * omega2 + g2 * omega  + offset;

    //tgtAssert(std::abs(cl1.imag()) <= std::abs(cl1.real()), "Invalid eigenvalue for symmat");
    //tgtAssert(std::abs(cl2.imag()) <= std::abs(cl2.real()), "Invalid eigenvalue for symmat");
    //tgtAssert(std::abs(cl3.imag()) <= std::abs(cl3.real()), "Invalid eigenvalue for symmat");

    l1=cl1.real();
    l2=cl2.real();
    l3=cl3.real();
}

/*
tgt::vec3 SymMat3::getEigenvectorToTheOtherEigenValue(float notThisOne, float alsoNotThisOne) const {
    // See again:
    // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
    SymMat3 A = minusUnitTimes(notThisOne)*minusUnitTimes(alsoNotThisOne);

    // Symmetric matrix => row_n == col_n
    if(A.xx != 0.0f || A.xy != 0.0f || A.xz != 0.0f) {
        return tgt::vec3(A.xx, A.xy, A.xz);
    } else if(A.yy != 0.0f || A.yz != 0.0f) {
        return tgt::vec3(A.xy, A.yy, A.yz);
    } else if(A.yz != 0.0f) {
        return tgt::vec3(A.xz, A.yz, A.zz);
    } else {
        tgtAssert(false, "Null matrix");
        return tgt::vec3::one;
    }
}
*/
tgt::mat3 SymMat3::toTgtMat() const {
    return tgt::mat3(
            xx, xy, xz,
            xy, yy, yz,
            xz, yz, zz
            );
}

const std::string VESSELNESS_BASE_TYPE = "float";

// VesselnessFeatureExtractor ---------------------------------------------------------------------------------------------------------

VesselnessFeatureExtractor::VesselnessFeatureExtractor(float /*alpha*/, float /*beta*/, float /*c*/, const tgt::ivec3& extent, const tgt::vec3 standardDeviation, float scale)
    : extent_(extent)
    //, two_alpha_squared_(2*alpha*alpha)
    //, two_beta_squared_(2*beta*beta)
    //, two_c_squared_(2*c*c)
    , xxKernel_(SeparableKernel::deriveXX(extent, standardDeviation))
    , xyKernel_(SeparableKernel::deriveXY(extent, standardDeviation))
    , xzKernel_(SeparableKernel::deriveXZ(extent, standardDeviation))
    , yyKernel_(SeparableKernel::deriveYY(extent, standardDeviation))
    , yzKernel_(SeparableKernel::deriveYZ(extent, standardDeviation))
    , zzKernel_(SeparableKernel::deriveZZ(extent, standardDeviation))
    , scale_(scale)
{
    tgtAssert(tgt::hand(tgt::greaterThan(extent, tgt::ivec3::zero)), "Invalid extent");
}

VesselnessFeatureExtractor::~VesselnessFeatureExtractor() {
}

//#undef VRN_MODULE_OPENMP

// Filter the slice using the kernel which was separated into three 1D kernels:

template<typename T>
static void getFilteredSliceGeneric(const VesselnessFeatureExtractor& vft, const CachingSliceReader* src, int z, VolumeRAM& outputSlice) {
    tgtAssert(z >= 0 && z<src->getSignedDimensions().z, "Invalid z pos in slice request");

    const tgt::ivec3& dim = src->getSignedDimensions();

    SimpleSlice<float> z0(dim.xy(), 0);
    SimpleSlice<float> z1(dim.xy(), 0);
    SimpleSlice<float> z2(dim.xy(), 0);

    for(int dz = -vft.extent_.z; dz <= vft.extent_.z; ++dz) {
        int mz = mirror(z+dz, dim.z);
        auto sliceBase = src->getSlice(mz-z);
        auto slice = dynamic_cast<const VolumeAtomic<T>*>(sliceBase);
        const T* sliceData = slice->voxel();
        tgtAssert(slice, "somehow dispatch is broken");

        float kernelZ0 = vft.xxKernel_.zKernelAt(dz);
        float kernelZ1 = vft.xzKernel_.zKernelAt(dz);
        float kernelZ2 = vft.zzKernel_.zKernelAt(dz);
        // z
#ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for
#endif
        for(int y = 0; y < dim.y; ++y) {
            size_t posBase = y*dim.x;
#if defined(WIN32) && _MSC_VER < 1922
#pragma vector // MSVC equivalent of omp simd prior to MSVC 1922
#else
#ifdef VRN_MODULE_OPENMP
#pragma omp simd
#endif
#endif
            for(int x = 0; x < dim.x; ++x) {
                size_t pos = posBase+x;
                tgtAssert(pos < slice->getNumVoxels(), "foo");
                tgtAssert(pos < z0.size(), "foo");
                tgtAssert(pos < z1.size(), "foo");
                tgtAssert(pos < z2.size(), "foo");
                float sample = getTypeAsFloat<T>(sliceData[pos]);
                z0.at(pos) += kernelZ0*sample;
                z1.at(pos) += kernelZ1*sample;
                z2.at(pos) += kernelZ2*sample;
            }
        }
    }

    SimpleSlice<float> xx(dim.xy(), 0);
    SimpleSlice<float> xy(dim.xy(), 0);
    SimpleSlice<float> xz(dim.xy(), 0);
    SimpleSlice<float> yy(dim.xy(), 0);
    SimpleSlice<float> yz(dim.xy(), 0);
    SimpleSlice<float> zz(dim.xy(), 0);

    // y
#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
    for(int y = 0; y < dim.y; ++y) {
        size_t posBase = y*dim.x;
        for(int dy = -vft.extent_.y; dy <= vft.extent_.y; ++dy) {
            int my = mirror(y+dy, dim.y);
            size_t mPosBase = my*dim.x;

            float y0 = vft.xxKernel_.yKernelAt(dy);
            float y1 = vft.xyKernel_.yKernelAt(dy);
            float y2 = vft.yyKernel_.yKernelAt(dy);
#if defined(WIN32) && _MSC_VER < 1922
#pragma vector // MSVC equivalent of omp simd prior to MSVC 1922
#else
#ifdef VRN_MODULE_OPENMP
#pragma omp simd
#endif
#endif
            for(int x = 0; x < dim.x; ++x) {
                size_t pos = posBase+x;
                size_t mpos = mPosBase+x;
                tgtAssert(mpos < z0.size(), "foo");
                tgtAssert(pos < xx.size(), "foo");

                float sampleZ0 = z0.at(mpos);
                float sampleZ1 = z1.at(mpos);
                float sampleZ2 = z2.at(mpos);

                xx.at(pos) += y0*sampleZ0;
                xy.at(pos) += y1*sampleZ0;
                xz.at(pos) += y0*sampleZ1;
                yy.at(pos) += y2*sampleZ0;
                yz.at(pos) += y1*sampleZ1;
                zz.at(pos) += y0*sampleZ2;
            }
        }
    }


    // x
#ifdef VRN_MODULE_OPENMP
    #pragma omp parallel for
#endif
    for(int y = 0; y < dim.y; ++y) {
        size_t posBase = y*dim.x;
        std::vector<float> xxFinal(dim.x, 0);
        std::vector<float> xyFinal(dim.x, 0);
        std::vector<float> xzFinal(dim.x, 0);
        std::vector<float> yyFinal(dim.x, 0);
        std::vector<float> yzFinal(dim.x, 0);
        std::vector<float> zzFinal(dim.x, 0);
        for(int dx = -vft.extent_.x; dx <= vft.extent_.x; ++dx) {
            float x0 = vft.yyKernel_.xKernelAt(dx);
            float x1 = vft.xyKernel_.xKernelAt(dx);
            float x2 = vft.xxKernel_.xKernelAt(dx);
            for(int x = 0; x < vft.extent_.x; ++x) {
                int mx = mirror(x+dx, dim.x);
                size_t mpos = mx+posBase;
                tgtAssert(mpos < xx.size(), "foo");

                xxFinal[x] += x2*xx.at(mpos);
                xyFinal[x] += x1*xy.at(mpos);
                xzFinal[x] += x1*xz.at(mpos);
                yyFinal[x] += x0*yy.at(mpos);
                yzFinal[x] += x0*yz.at(mpos);
                zzFinal[x] += x0*zz.at(mpos);
            }
            {
#if defined(WIN32) && _MSC_VER < 1922
#pragma vector // MSVC equivalent of omp simd prior to MSVC 1922
#else
#ifdef VRN_MODULE_OPENMP
#pragma omp simd
#endif
#endif
            for(int x = vft.extent_.x; x < dim.x-vft.extent_.x; ++x) {
                int mx = x+dx;
                size_t mpos = mx+posBase;
                tgtAssert(mpos < xx.size(), "foo");

                xxFinal[x] += x2*xx.at(mpos);
                xyFinal[x] += x1*xy.at(mpos);
                xzFinal[x] += x1*xz.at(mpos);
                yyFinal[x] += x0*yy.at(mpos);
                yzFinal[x] += x0*yz.at(mpos);
                zzFinal[x] += x0*zz.at(mpos);
            }
            }
            for(int x = dim.x-vft.extent_.x; x < dim.x; ++x) {
                int mx = mirror(x+dx, dim.x);
                size_t mpos = mx+posBase;
                tgtAssert(mpos < xx.size(), "foo");

                xxFinal[x] += x2*xx.at(mpos);
                xyFinal[x] += x1*xy.at(mpos);
                xzFinal[x] += x1*xz.at(mpos);
                yyFinal[x] += x0*yy.at(mpos);
                yzFinal[x] += x0*yz.at(mpos);
                zzFinal[x] += x0*zz.at(mpos);
            }
        }

        for(int x = 0; x < dim.x; ++x) {
            SymMat3 mat {
                xxFinal.at(x),
                xyFinal.at(x),
                xzFinal.at(x),
                yyFinal.at(x),
                yzFinal.at(x),
                zzFinal.at(x),
            };
            float output = vft.computeVesselness(mat);

            outputSlice.setVoxelNormalized(output * vft.scale_, tgt::svec3(x,y,0));
        }
    }
}

std::unique_ptr<VolumeRAM> VesselnessFeatureExtractor::getFilteredSlice(const CachingSliceReader* src, int z) const {
    const tgt::ivec3& dim = src->getSignedDimensions();
    std::string inputBasetype = src->getMetaData().getBaseType();
    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create(VESSELNESS_BASE_TYPE, tgt::svec3(dim.xy(), 1)));
    DISPATCH_FOR_BASETYPE(inputBasetype, getFilteredSliceGeneric, *this, src, z, *outputSlice);
    return outputSlice;
}

int VesselnessFeatureExtractor::zExtent() const {
    return extent_.z;
}

size_t VesselnessFeatureExtractor::getNumInputChannels() const {
    return 1;
}

size_t VesselnessFeatureExtractor::getNumOutputChannels() const {
    return 1;
}


SliceReaderMetaData VesselnessFeatureExtractor::getMetaData(const SliceReaderMetaData& base) const {
    auto md = SliceReaderMetaData::fromBase(base);
    md.setRealWorldMapping(RealWorldMapping(tgt::vec2(0.0, 1.0), "Vesselness"));

    md.setBaseType(VESSELNESS_BASE_TYPE);

    return md;
}

float satoVesselness(float l1, float l2, float l3, float two_alpha1_squared, float two_alpha2_squared) {
    tgtAssert(two_alpha1_squared < two_alpha2_squared, "alpha1 >= alpha2");

    // Do a quick bubble sort so that l1 >= l2 >= l3
    if(l1 < l2) {
        std::swap(l1, l2);
    }
    if(l2 < l3) {
        std::swap(l2, l3);
    }
    if(l1 < l2) {
        std::swap(l1, l2);
    }

    //float lc = std::min(-l2, -l3);
    float lc = -l2;
    if(lc == 0) {
        return 0;
    }

    if(l1 <= 0) {
        return std::exp((-l1*l1)/(two_alpha1_squared*lc*lc))*lc;
    } else {
        return std::exp((-l1*l1)/(two_alpha2_squared*lc*lc))*lc;
    }
}

float frangiVesselness(float l1, float l2, float l3, float two_alpha_squared, float two_beta_squared, float two_c_squared) {
    // Do a quick bubble sort so that |l1| <= |l2| <= |l3|
    if(std::abs(l1) > std::abs(l2)) {
        std::swap(l1, l2);
    }
    if(std::abs(l2) > std::abs(l3)) {
        std::swap(l2, l3);
    }
    if(std::abs(l1) > std::abs(l2)) {
        std::swap(l1, l2);
    }

    double R_A = std::abs(l2) / std::abs(l3);
    double R_B = std::abs(l1) / std::sqrt(std::abs(l2*l3));
    double S = std::sqrt(l1*l1 + l2*l2 + l3*l3);

    return (l2 > 0 || l3 > 0) ? 0 : (1-std::exp(-R_A*R_A/two_alpha_squared)) * std::exp(-R_B*R_B/two_beta_squared) * (1-std::exp(-S*S/two_c_squared));
}

float erdtVesselness(float l1, float l2, float l3) {
    // Do a quick bubble sort so that l1 >= l2 >= l3
    if(l1 < l2) {
        std::swap(l1, l2);
    }
    if(l2 < l3) {
        std::swap(l2, l3);
    }
    if(l1 < l2) {
        std::swap(l1, l2);
    }

    float l2abs = std::abs(l2);
    float l3abs = std::abs(l3);
    float K = 1 - std::abs(l2abs - l3abs)/(l2abs + l3abs);
    return std::min(1.0f, std::max(0.0f, K*((2.0f/3)*l1 - l2 - l3)));
}

float VesselnessFeatureExtractor::computeVesselness(const SymMat3& H) const {
    // Declare eigenvalue variables
    float l1 = 0, l2 = 0, l3 = 0;

    H.getEigenValues(l1, l2, l3);

    // Do a quick bubble sort so that |l1| <= |l2| <= |l3|

    // Unclear in the publication of Frangi et al.: Is c supposed to be computed from global or local hessian norm. Probably global.
    //double hessian_max_norm = std::max(std::abs(H.xx), std::max(std::abs(H.xy), std::max(std::abs(H.xz), std::max(std::abs(H.yy), std::max(std::abs(H.yz), std::abs(H.zz))))));
    //double two_c_squared = 0.5*hessian_max_norm*hessian_max_norm;

    //float vesselness = frangiVesselness(l1, l2, l3, two_alpha_squared_, two_beta_squared_, two_c_squared_);
    float vesselness = satoVesselness(l1, l2, l3, 0.5 /* = 2*0.5*0.5 */, 8.0 /* = 2*2*2 */); //Turned out to yield best results
    //float vesselness = erdtVesselness(l1, l2, l3);
    return vesselness;
}


// VesselnessExtractor ---------------------------------------------------------------------------------------------------------

static tgt::vec3 getRelativeSpacing(const tgt::vec3 spacing) {
    float medianDimSpacing = 0;
    if(spacing.x > spacing.y) {
        if(spacing.x > spacing.z) {
            // x == max
            medianDimSpacing = std::max(spacing.y, spacing.z);
        } else {
            // z >= x > y
            medianDimSpacing = spacing.x;
        }
    } else {
        if(spacing.y > spacing.z) {
            // y == max
            medianDimSpacing = std::max(spacing.x, spacing.y);
        } else {
            // z >= y > x
            medianDimSpacing = spacing.y;
        }
    }
    return spacing/medianDimSpacing;
}

}

namespace voreen {

const std::string VesselnessExtractor::loggerCat_("voreen.vesselnetworkanalysis.vesselnessextractor");

VesselnessExtractor::VesselnessExtractor()
    : AsyncComputeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enabled_("enabled", "Enabled", true)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output Volume", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , vesselRadiusRangeRW_("vesselRadiusRangeRW", "Vessel Radius (mm)", 1.0, std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max())
    , scaleSpaceSteps_("scaleSpaceSteps", "Scale Space Steps", 5, 1, 10)
    , minStandardDeviationVec_("minStandardDeviationVec", "Used Min Standard Deviation (voxel)", tgt::vec3::zero, tgt::vec3::zero, tgt::vec3(std::numeric_limits<float>::max()))
    , maxStandardDeviationVec_("maxStandardDeviationVec", "Used Max Standard Deviation (voxel)", tgt::vec3::zero, tgt::vec3::zero, tgt::vec3(std::numeric_limits<float>::max()))
    , minSmoothingKernelSize_("minSmoothingKernelSize", "Min Smoothing Kernel Size", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3(std::numeric_limits<int>::max()))
    , maxSmoothingKernelSize_("maxSmoothingKernelSize", "Max Smoothing Kernel Size", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3(std::numeric_limits<int>::max()))
    , propertyDisabler_(*this)
{
    addPort(inport_);
    inport_.addCondition(new PortConditionVolumeChannelCount(1));
    addPort(outport_);

    addProperty(enabled_);
        ON_CHANGE_LAMBDA(enabled_, [this] () {
            if(enabled_.get()) {
                this->outport_.setData(nullptr);
                propertyDisabler_.restore();
            } else {
                propertyDisabler_.saveState([this] (Property* p) { return p == &enabled_; });
                propertyDisabler_.disable();
            }
        });

    addProperty(outputVolumeFilePath_);

    addProperty(vesselRadiusRangeRW_);
        ON_CHANGE(vesselRadiusRangeRW_, VesselnessExtractor, updateSmoothingProperties);

    addProperty(scaleSpaceSteps_);

    addProperty(minStandardDeviationVec_);
        minStandardDeviationVec_.setReadOnlyFlag(true);
    addProperty(maxStandardDeviationVec_);
        maxStandardDeviationVec_.setReadOnlyFlag(true);

    addProperty(minSmoothingKernelSize_);
        minSmoothingKernelSize_.setReadOnlyFlag(true);
    addProperty(maxSmoothingKernelSize_);
        maxSmoothingKernelSize_.setReadOnlyFlag(true);

    // Not relevant for sato-vesselness:
    //addProperty(blobRejectorWeight_);
    //addProperty(planeRejectorWeight_);
    //addProperty(intensityThreshold_);
        //intensityThreshold_.adaptDecimalsToRange(5);
    propertyDisabler_.saveState([this] (Property* p) { return p == &enabled_; });
}

const static float BLOB_REJECTOR_WEIGHT = 0.5;
const static float PLANE_REJECTOR_WEIGHT = 0.5;
const static float INTENSITY_THRESHOLD = 0.5;

VesselnessExtractor::~VesselnessExtractor() {}

Processor* VesselnessExtractor::create() const {
    return new VesselnessExtractor();
}

static void voxelwiseMax(SliceReader& reader, HDF5FileVolume& file, ProgressReporter* progress) {
    tgt::svec3 dim = reader.getDimensions();
    tgtAssert(dim == file.getDimensions(), "Dimension mismatch");

    reader.seek(0);
    for(size_t z = 0; z < dim.z; ++z) {
        if(progress) {
            progress->setProgress(static_cast<float>(z)/dim.z);
        }
        if(z > 0) {
            reader.advance();
        }
        const VolumeRAM* activeSlice = reader.getCurrentSlice();

        // Read the slice from the file
        std::unique_ptr<VolumeRAM> fileSlice(file.loadSlices(z, z));

        // Do a voxel wise max on fileSlice
        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                tgt::svec3 p(x, y, 0);
                float val = std::max(
                        fileSlice->getVoxelNormalized(p),
                        activeSlice->getVoxelNormalized(p));
                fileSlice->setVoxelNormalized(val, p);
            }
        }

        // Write the slice back
        file.writeSlices(fileSlice.get(), z);
    }
    if(progress) {
        progress->setProgress(1.0f);
    }
}

static std::unique_ptr<SliceReader> buildStack(const VolumeBase& input, const tgt::vec3& standardDeviationVec) {
    //tgt::ivec3 dirVesselnessExtent = tgt::ivec3::one;
    float scale = tgt::length(standardDeviationVec);

    return VolumeFilterStackBuilder(input)
        .addLayer(std::unique_ptr<VolumeFilter>(new VesselnessFeatureExtractor(PLANE_REJECTOR_WEIGHT, BLOB_REJECTOR_WEIGHT, INTENSITY_THRESHOLD, suitableExtent(standardDeviationVec), standardDeviationVec, scale)))
        //.addLayer(std::unique_ptr<VolumeFilter>(new VesselnessFinalizer(dirVesselnessExtent, SamplingStrategy<ParallelFilterValue4D>::MIRROR)))
        .build(0);
}

void VesselnessExtractor::process() {
    if(!enabled_.get()) {
        outport_.setData(inport_.getData(), false);
        return;
    }
    AsyncComputeProcessor<VesselnessExtractorInput, VesselnessExtractorOutput>::process();
}
VesselnessExtractorInput VesselnessExtractor::prepareComputeInput() {
    const VolumeBase* inputVol = inport_.getData();
    if(!inputVol) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    if(outputVolumeFilePath_.get().empty()) {
        throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
    }

    // Close the volume output file
    outport_.clear();

    const std::string baseType = VESSELNESS_BASE_TYPE;
    const std::string volumeLocation = "/vol";
    const size_t numChannels = 1;
    const int deflateLevel = 1;
    const bool truncateFile = true;
    std::unique_ptr<HDF5FileVolume> output(nullptr);
    try {
        output = HDF5FileVolume::createVolume(outputVolumeFilePath_.get(), volumeLocation, baseType, inputVol->getDimensions(), numChannels, truncateFile, deflateLevel, tgt::svec3(inputVol->getDimensions().xy(), 1), false); //May throw IOException
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }
    tgtAssert(output, "no output volume");

    if(inputVol->hasMetaData("RealWorldMapping")) {
        output->writeRealWorldMapping(inputVol->getRealWorldMapping());
    }
    output->writeOffset(inputVol->getOffset());
    output->writeSpacing(inputVol->getSpacing());
    output->writePhysicalToWorldTransformation(inputVol->getPhysicalToWorldMatrix());

    return VesselnessExtractorInput(
            *inputVol,
            std::move(output),
            scaleSpaceSteps_.get(),
            minStandardDeviationVec_.get(),
            maxStandardDeviationVec_.get()
            );
}

VesselnessExtractorOutput VesselnessExtractor::compute(VesselnessExtractorInput input, ProgressReporter& progressReporter) const {

    const float stepProgressDelta = 1.0f/scaleSpaceSteps_.get();

    progressReporter.setProgressRange(tgt::vec2(0.0, stepProgressDelta));

    // Write first slices to volume
    {
        auto stddev = input.getStandardDeviationForStep(0);
        std::unique_ptr<SliceReader> reader = buildStack(input.input, stddev);

        writeSlicesToHDF5File(*reader, *input.output, &progressReporter);
    }

    // Now build max with the following scale space steps
    for(int step=1; step < input.scaleSpaceSteps; ++step) {
        auto stddev = input.getStandardDeviationForStep(step);
        std::unique_ptr<SliceReader> reader = buildStack(input.input, stddev);

        progressReporter.setProgressRange(tgt::vec2(stepProgressDelta*step, stepProgressDelta*(step+1)));

        voxelwiseMax(*reader, *input.output, &progressReporter);
    }

    return input.output->getFileName();
}

void VesselnessExtractor::processComputeOutput(VesselnessExtractorOutput output) {
    std::unique_ptr<VolumeList> volumes(HDF5VolumeReader().read(output));
    const VolumeBase* vol = volumes->at(0);
    outport_.setData(vol);
}

void VesselnessExtractor::adjustPropertiesToInput() {
    updateSmoothingProperties();

    const VolumeBase* inputVol = inport_.getData();
    if(!inputVol) {
        return;
    }
    vesselRadiusRangeRW_.setMinValue(tgt::min(inputVol->getSpacing()));
    vesselRadiusRangeRW_.setMaxValue(tgt::min(inputVol->getBoundingBox(false).getBoundingBox().diagonal()));
    vesselRadiusRangeRW_.adaptDecimalsToRange(5);
}

void VesselnessExtractor::updateSmoothingProperties() {
    const VolumeBase* inputVol = inport_.getData();
    if(!inputVol) {
        minStandardDeviationVec_.reset();
        maxStandardDeviationVec_.reset();
        minSmoothingKernelSize_.reset();
        maxSmoothingKernelSize_.reset();
        return;
    }
    minStandardDeviationVec_.set(vesselRadiusRangeRW_.get().x/inputVol->getSpacing());
    maxStandardDeviationVec_.set(vesselRadiusRangeRW_.get().y/inputVol->getSpacing());

    tgt::ivec3 minSize = 2*suitableExtent(minStandardDeviationVec_.get()) + tgt::ivec3::one;
    tgt::ivec3 maxSize = 2*suitableExtent(maxStandardDeviationVec_.get()) + tgt::ivec3::one;
    minSmoothingKernelSize_.set(minSize);
    maxSmoothingKernelSize_.set(maxSize);
}

bool VesselnessExtractor::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    if(!inport_.isReady()) {
        setNotReadyErrorMessage("Inport not ready.");
        return false;
    }
    return true;
}
void VesselnessExtractor::initialize() {
    AsyncComputeProcessor::initialize();
}
}   // namespace
