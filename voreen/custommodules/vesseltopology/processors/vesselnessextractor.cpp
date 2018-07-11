#include "vesselnessextractor.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include "custommodules/bigdataimageprocessing/volumefiltering/gaussianfilter.h"
#include "tgt/tgt_math.h"

#include <chrono>

namespace {

using namespace voreen;

// ------------------------------------------------------------------------------------------
// Implementation of helper classes----------------------------------------------------------
// ------------------------------------------------------------------------------------------

// Space efficient symmetric 3 dimensional matrix
struct SymMat3 {
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

    VesselnessFeatureExtractor(float alpha, float beta, float c, const tgt::ivec3& extent, const tgt::vec3 standardDeviation, const SamplingStrategy<float>& samplingStrategy, const std::string sliceBaseType);
    virtual ~VesselnessFeatureExtractor();

    tgt::vec4 computeVesselnessFeatureVector(const SymMat3& hessian) const;

    virtual std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const;
    virtual int zExtent() const;
    virtual size_t getNumInputChannels() const;
    virtual size_t getNumOutputChannels() const;
    virtual const std::string& getSliceBaseType() const;

private:
    const std::string sliceBaseType_;
    SamplingStrategy<float> samplingStrategy_;

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
};

// Computes an improved vesselness from a four component vector (e_1, v) of vessel direction vector and vesselness value
// by combinining v with a value of per-voxel uniformity of vessel directions in the neighborhood of all voxels.
class VesselnessFinalizer : public ParallelVolumeFilter<ParallelFilterValue4D, ParallelFilterValue1D> {
public:
    VesselnessFinalizer(const tgt::ivec3& extent, const SamplingStrategy<ParallelFilterValue4D>& samplingStrategy, const std::string sliceBaseType);
    virtual ~VesselnessFinalizer();

    ParallelFilterValue1D getValue(const Sample& sample, const tgt::ivec3& pos) const;

private:
    tgt::ivec3 extent_;
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

    std::cout << "[";
    for(auto val : xKernel_) {
        std::cout << val << ", ";
    }
    std::cout << "]\n";
    std::cout.flush();
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
    // Find eigenvalues using the following method:
    // https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
    float p1 = xy*xy + xz*xz + yz*yz;
    if(p1 == 0) {
        l1 = xx;
        l2 = yy;
        l3 = zz;
    } else {
        float q = (xx + yy + zz) / 3; // Trace/3
        float p2 = (xx - q)*(xx - q) + (yy - q)*(yy - q) + (zz - q)*(zz - q) + 2*p1;
        float p = sqrt(p2 / 6);

        // Construct entries of symmetric matrix B = (1/p) * (A - q*E)
        float bxx = (xx - q)/p;
        float byy = (yy - q)/p;
        float bzz = (zz - q)/p;
        float bxy = xy/p;
        float bxz = xz/p;
        float byz = yz/p;

        // Now determine the determinant
        float detB = bxx*byy*bzz + 2*bxy*byz*bxz - bxx*byz*byz - byy*bxz*bxz - bzz*bxy*bxy;
        float r = detB / 2;

        float phi;
        if(r <= -1.f) {
            phi = tgt::PIf / 3;
        } else if(r >= -1.f) {
            phi = 0;
        } else {
            phi = std::acos(r) / 3;
        }

        // The eigenvalues are sorted l3 <= l2 <= l1... However we want them sorted after their absolute value :(
        l1 = q + 2*p*std::cos(phi);
        l3 = q + 2*p*std::cos(phi + (2*tgt::PIf/3));
        l2 = 3 * q - l1 - l3;
    }
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

// VesselnessFeatureExtractor ---------------------------------------------------------------------------------------------------------
//

VesselnessFeatureExtractor::VesselnessFeatureExtractor(float /*alpha*/, float /*beta*/, float /*c*/, const tgt::ivec3& extent, const tgt::vec3 standardDeviation, const SamplingStrategy<float>& samplingStrategy, const std::string sliceBaseType)
    : samplingStrategy_(samplingStrategy)
    , sliceBaseType_(sliceBaseType)
    , extent_(extent)
    //, two_alpha_squared_(2*alpha*alpha)
    //, two_beta_squared_(2*beta*beta)
    //, two_c_squared_(2*c*c)
    , xxKernel_(SeparableKernel::deriveXX(extent, standardDeviation))
    , xyKernel_(SeparableKernel::deriveXY(extent, standardDeviation))
    , xzKernel_(SeparableKernel::deriveXZ(extent, standardDeviation))
    , yyKernel_(SeparableKernel::deriveYY(extent, standardDeviation))
    , yzKernel_(SeparableKernel::deriveYZ(extent, standardDeviation))
    , zzKernel_(SeparableKernel::deriveZZ(extent, standardDeviation))
{
    tgtAssert(tgt::hand(tgt::greaterThan(extent, tgt::ivec3::zero)), "Invalid extent");
}

VesselnessFeatureExtractor::~VesselnessFeatureExtractor() {
}

// Filter the slice using the kernel which was separated into three 1D kernels:
std::unique_ptr<VolumeRAM> VesselnessFeatureExtractor::getFilteredSlice(const CachingSliceReader* src, int z) const {
    tgtAssert(z >= 0 && z<src->getSignedDimensions().z, "Invalid z pos in slice request");

    typedef SimpleSlice<SymMat3> TempSlice;

    const tgt::ivec3& dim = src->getSignedDimensions();

    SamplingStrategy<float>::Sampler getValueFromReader = [src] (const tgt::ivec3& p) {
        return src->getVoxelNormalized(p);
    };

    TempSlice zOutput(dim.xy());

    // z
    #pragma omp parallel for
    for(int y = 0; y < dim.y; ++y) {
        for(int x = 0; x < dim.x; ++x) {
            SymMat3& accumulator = zOutput.at(x, y);
            accumulator.xx = 0;
            accumulator.xy = 0;
            accumulator.xz = 0;
            accumulator.yy = 0;
            accumulator.yz = 0;
            accumulator.zz = 0;
            for(int dz = -extent_.z; dz <= extent_.z; ++dz) {
                float sample = samplingStrategy_.sample(tgt::ivec3(x, y, z+dz), dim, getValueFromReader);
                accumulator.xx += xxKernel_.zKernelAt(dz)*sample;
                accumulator.xy += xyKernel_.zKernelAt(dz)*sample;
                accumulator.xz += xzKernel_.zKernelAt(dz)*sample;
                accumulator.yy += yyKernel_.zKernelAt(dz)*sample;
                accumulator.yz += yzKernel_.zKernelAt(dz)*sample;
                accumulator.zz += zzKernel_.zKernelAt(dz)*sample;
            }
        }
    }

    TempSlice yOutput(dim.xy());
    SamplingStrategy<SymMat3> sliceSamplingStrategy = samplingStrategy_.convert<SymMat3>([] (float f) {
            return SymMat3 {
            f, f, f, f, f, f
            };
            });

    // y
    #pragma omp parallel for
    for(int y = 0; y < dim.y; ++y) {
        for(int x = 0; x < dim.x; ++x) {
            SymMat3& accumulator = yOutput.at(x, y);
            accumulator.xx = 0;
            accumulator.xy = 0;
            accumulator.xz = 0;
            accumulator.yy = 0;
            accumulator.yz = 0;
            accumulator.zz = 0;
            for(int dy = -extent_.y; dy <= extent_.y; ++dy) {
                SymMat3 sample = sliceSamplingStrategy.sample(tgt::ivec3(x, y+dy, 0), dim, zOutput.toSampler());
                accumulator.xx += xxKernel_.yKernelAt(dy)*sample.xx;
                accumulator.xy += xyKernel_.yKernelAt(dy)*sample.xy;
                accumulator.xz += xzKernel_.yKernelAt(dy)*sample.xz;
                accumulator.yy += yyKernel_.yKernelAt(dy)*sample.yy;
                accumulator.yz += yzKernel_.yKernelAt(dy)*sample.yz;
                accumulator.zz += zzKernel_.yKernelAt(dy)*sample.zz;
            }
        }
    }

    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create("Vector4(" + sliceBaseType_ + ")", tgt::svec3(dim.xy(), 1)));

    // x
    #pragma omp parallel for
    for(int y = 0; y < dim.y; ++y) {
        for(int x = 0; x < dim.x; ++x) {
            SymMat3 accumulator {
                0, 0, 0, 0, 0, 0
            };
            for(int dx = -extent_.x; dx <= extent_.x; ++dx) {
                SymMat3 sample = sliceSamplingStrategy.sample(tgt::ivec3(x+dx, y, 0), dim, yOutput.toSampler());
                accumulator.xx += xxKernel_.xKernelAt(dx)*sample.xx;
                accumulator.xy += xyKernel_.xKernelAt(dx)*sample.xy;
                accumulator.xz += xzKernel_.xKernelAt(dx)*sample.xz;
                accumulator.yy += yyKernel_.xKernelAt(dx)*sample.yy;
                accumulator.yz += yzKernel_.xKernelAt(dx)*sample.yz;
                accumulator.zz += zzKernel_.xKernelAt(dx)*sample.zz;
            }
            tgt::vec4 output = computeVesselnessFeatureVector(accumulator);
            outputSlice->setVoxelNormalized(output.x, tgt::svec3(x,y,0), 0);
            outputSlice->setVoxelNormalized(output.y, tgt::svec3(x,y,0), 1);
            outputSlice->setVoxelNormalized(output.z, tgt::svec3(x,y,0), 2);
            outputSlice->setVoxelNormalized(output.w, tgt::svec3(x,y,0), 3);
        }
    }
    return outputSlice;
}

int VesselnessFeatureExtractor::zExtent() const {
    return extent_.z;
}

size_t VesselnessFeatureExtractor::getNumInputChannels() const {
    return 1;
}

size_t VesselnessFeatureExtractor::getNumOutputChannels() const {
    return 4;
}

const std::string& VesselnessFeatureExtractor::getSliceBaseType() const {
    return sliceBaseType_;
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

tgt::vec4 VesselnessFeatureExtractor::computeVesselnessFeatureVector(const SymMat3& H) const {
    // Declare eigenvalue variables
    float l1 = 0, l2 = 0, l3 = 0;

    H.getEigenValues(l1, l2, l3);

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

    tgt::vec3 vesselDirVector(1,1,1);
    if(l2 > 0 || l3 > 0) {
        // Not a bright vessel
        vesselDirVector = tgt::vec3::zero;
    } else {
        vesselDirVector = tgt::normalize(vesselDirVector);

        // Find e_1 using inverse iteration
        tgt::mat3 A = H.toTgtMat();
        tgt::mat3 itMat;
        if((A-l1*tgt::mat3::identity).invert(itMat)) {

            for(int i=0; i<5; ++i) {
                vesselDirVector = tgt::normalize(itMat*vesselDirVector);
            }
        } else {
            // This can happen by chance, so just set e_1 to zero here
            vesselDirVector = tgt::vec3::zero;
            //LWARNINGC("voreen.vesseltopology.vesselnessextractor", "Singular Hesse Matrix, setting vessel dir to zero");
        }
    }

    // Unclear in the publication of Frangi et al.: Is c supposed to be computed from global or local hessian norm. Probably global.
    //double hessian_max_norm = std::max(std::abs(H.xx), std::max(std::abs(H.xy), std::max(std::abs(H.xz), std::max(std::abs(H.yy), std::max(std::abs(H.yz), std::abs(H.zz))))));
    //double two_c_squared = 0.5*hessian_max_norm*hessian_max_norm;

    //float vesselness = frangiVesselness(l1, l2, l3, two_alpha_squared_, two_beta_squared_, two_c_squared_);
    float vesselness = satoVesselness(l1, l2, l3, 0.5 /* = 2*0.5*0.5 */, 8.0 /* = 2*2*2 */); //Turned out to yield best results
    //float vesselness = erdtVesselness(l1, l2, l3);
    return tgt::vec4(vesselDirVector, vesselness);
}


// VesselnessFinalizer ---------------------------------------------------------------------------------------------------------



VesselnessFinalizer::VesselnessFinalizer(const tgt::ivec3& extent, const SamplingStrategy<ParallelFilterValue4D>& samplingStrategy, const std::string sliceBaseType)
    : ParallelVolumeFilter<ParallelFilterValue4D, ParallelFilterValue1D>(extent.z, samplingStrategy, sliceBaseType)
    , extent_(extent)
{
}

VesselnessFinalizer::~VesselnessFinalizer() {
}

ParallelFilterValue1D VesselnessFinalizer::getValue(const Sample& sample, const tgt::ivec3& pos) const {
    int extent = zExtent();

    tgt::vec4 thisVal = sample(pos);
    float scalarVesselness = thisVal[3];



    float dirAccumulator = 0.0f;
    tgt::vec3 thisDir = thisVal.xyz();
    for(int z = pos.z-extent; z <= pos.z+extent; ++z) {
        for(int y = pos.y-extent; y <= pos.y+extent; ++y) {
            for(int x = pos.x-extent; x <= pos.x+extent; ++x) {
                if(x != 0 || y != 0 || z != 0) {
                    tgt::vec4 val = sample(tgt::ivec3(x,y,z));
                    dirAccumulator += std::abs(tgt::dot(thisDir, val.xyz()));
                }
            }
        }
    }
    dirAccumulator /= hmul(2*extent_+ tgt::ivec3::one) - 1;


    //return scalarVesselness;
    return std::sqrt(scalarVesselness * dirAccumulator);
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

const std::string VesselnessExtractor::loggerCat_("voreen.vesseltopology.vesselnessextractor");

VesselnessExtractor::VesselnessExtractor()
    : AsyncComputeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , outputVolumeFilePath_("outputVolumeFilePath", "Output Volume", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , vesselRadiusRangeRW_("vesselRadiusRangeRW", "Vessel Radius (mm)", 1.0, 0.001, 500)
    , scaleSpaceSteps_("scaleSpaceSteps", "Scale Space Steps", 5, 1, 10)
    , minStandardDeviationVec_("minStandardDeviationVec", "Used Min Standard Deviation (voxel)", tgt::vec3::zero, tgt::vec3::zero, tgt::vec3(std::numeric_limits<float>::max()))
    , maxStandardDeviationVec_("maxStandardDeviationVec", "Used Max Standard Deviation (voxel)", tgt::vec3::zero, tgt::vec3::zero, tgt::vec3(std::numeric_limits<float>::max()))
    , minSmoothingKernelSize_("minSmoothingKernelSize", "Min Smoothing Kernel Size", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3(std::numeric_limits<int>::max()))
    , maxSmoothingKernelSize_("maxSmoothingKernelSize", "Max Smoothing Kernel Size", tgt::ivec3::zero, tgt::ivec3::zero, tgt::ivec3(std::numeric_limits<int>::max()))
{
    addPort(inport_);
    addPort(outport_);

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

static std::unique_ptr<SliceReader> buildStack(const VolumeBase& input, const tgt::vec3& standardDeviationVec, const std::string& baseType) {
    tgt::ivec3 dirVesselnessExtent = tgt::ivec3::one;

    return VolumeFilterStackBuilder(input)
        .addLayer(std::unique_ptr<VolumeFilter>(new VesselnessFeatureExtractor(PLANE_REJECTOR_WEIGHT, BLOB_REJECTOR_WEIGHT, INTENSITY_THRESHOLD, suitableExtent(standardDeviationVec), standardDeviationVec, SamplingStrategy<float>::MIRROR, baseType)))
        .addLayer(std::unique_ptr<VolumeFilter>(new VesselnessFinalizer(dirVesselnessExtent, SamplingStrategy<ParallelFilterValue4D>::MIRROR, baseType)))
        .build(0);
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

    //const std::string baseType = inputVol->getBaseType();
    const std::string baseType = "float";
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
            maxStandardDeviationVec_.get(),
            baseType
            );
}

VesselnessExtractorOutput VesselnessExtractor::compute(VesselnessExtractorInput input, ProgressReporter& progressReporter) const {

    const float stepProgressDelta = 1.0f/scaleSpaceSteps_.get();

    progressReporter.setProgressRange(tgt::vec2(0.0, stepProgressDelta));

    // Write first slices to volume
    {
        std::unique_ptr<SliceReader> reader = buildStack(input.input, input.minStandardDeviationVec, input.baseType);

        writeSlicesToHDF5File(*reader, *input.output, &progressReporter);
    }

    // Now build max with the following scale space steps
    for(int step=1; step < input.scaleSpaceSteps; ++step) {
        std::unique_ptr<SliceReader> reader = buildStack(input.input, input.getStandardDeviationForStep(step), input.baseType);

        progressReporter.setProgressRange(tgt::vec2(stepProgressDelta*step, stepProgressDelta*(step+1)));

        voxelwiseMax(*reader, *input.output, &progressReporter);
    }

    return input.output->getFileName();
}

void VesselnessExtractor::processComputeOutput(VesselnessExtractorOutput output) {
    const VolumeBase* vol = HDF5VolumeReader().read(output)->at(0);
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
}   // namespace
