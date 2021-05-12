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

#include "largetestdatagenerator.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/utils/hashing.h"
#include "../../modules/bigdataimageprocessing/algorithm/intervalwalker.h"
#include "../algorithm/boundshierarchy.h"
#include "custommodules/bigdataimageprocessingextra/ext/ziggurat.h"

#include "modules/hdf5/io/hdf5volumewriter.h"
#include "modules/hdf5/io/hdf5volumereader.h"

#include "tgt/filesystem.h"
#include "tgt/quaternion.h"

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>

#include <memory>

namespace voreen {

const std::string LargeTestDataGenerator::loggerCat_("voreen.bigdataimageprocessingextra.largetestdatagenerator");

LargeTestDataGenerator::LargeTestDataGenerator()
    : AsyncComputeProcessor<LargeTestDataGeneratorInput, LargeTestDataGeneratorOutput>()
    , outportNoisy_(Port::OUTPORT, "largetestdatagenerator.outport", "Volume Output Noisy", true, Processor::VALID)
    , outportGT_(Port::OUTPORT, "largetestdatagenerator.outportgt", "Volume Output Ground Truth", true, Processor::VALID)
    , foregroundLabelsPort_(Port::OUTPORT, "largetestdatagenerator.foregroundLabelsPort", "Foreground Seeds Outport", true, Processor::VALID)
    , backgroundLabelsPort_(Port::OUTPORT, "largetestdatagenerator.backgroundLabelsPort", "Background Seeds Outport", true, Processor::VALID)
    , volumeDimensions_("volumeDimensions", "Volume Dimensions", tgt::ivec3(2), tgt::ivec3(2), tgt::ivec3(10000))
    , structureSizeRange_("structureSizeRange", "Structure Size", 1, 1, 10000)
    , foregroundMean_("foregroundMean", "Foreground Mean Intensity", 0xffff*7/10, 0, 0xffff)
    , backgroundMean_("backgroundMean", "Background Mean Intensity", 0xffff*3/10, 0, 0xffff)
    , noiseType_("noiseType", "Noise Type")
    , gaussianNoiseSD_("gaussianNoiseSD", "Noise Standard Deviation", 0xffff/10, 0, 0xffff)
    , density_("density", "Object Density", 0.1, 0.0, 1.0)
    , seed_("seed", "RNG Seed", 0, 0, std::numeric_limits<int>::max())
    , outputVolumeNoisyFilePath_("outputVolumeFilePath", "Volume Noisy Output", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , outputVolumeGTFilePath_("outputVolumeFilePathgt", "GT Volume Output", "Path", "", "HDF5 (*.h5)", FileDialogProperty::SAVE_FILE, Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , scenario_("scenario", "Test Data Scenario")
    , retainLabel_("retainLabel", "Retain a foreground label from the output", false)
{
    addPort(outportNoisy_);
    addPort(outportGT_);
    addPort(foregroundLabelsPort_);
    addPort(backgroundLabelsPort_);

    addProperty(volumeDimensions_);
        ON_CHANGE_LAMBDA(volumeDimensions_, [this] () {
            structureSizeRange_.setMaxValue(tgt::min(volumeDimensions_.get()));
        });
        volumeDimensions_.setTracking(false);
    addProperty(structureSizeRange_);
        structureSizeRange_.setMaxValue(tgt::min(volumeDimensions_.get()));
    addProperty(foregroundMean_);
        foregroundMean_.setTracking(false);
    addProperty(backgroundMean_);
        backgroundMean_.setTracking(false);
    addProperty(noiseType_);
        noiseType_.addOption("gaussian", "Gaussian", LargeTestDataGeneratorInput::GAUSSIAN);
        noiseType_.addOption("poisson", "Poisson", LargeTestDataGeneratorInput::POISSON);
        ON_CHANGE_LAMBDA(noiseType_, [this] () {
            gaussianNoiseSD_.setVisibleFlag(noiseType_.getValue() == LargeTestDataGeneratorInput::GAUSSIAN);
        });
        noiseType_.selectByValue(LargeTestDataGeneratorInput::GAUSSIAN);
    addProperty(gaussianNoiseSD_);
        gaussianNoiseSD_.setTracking(false);
    addProperty(density_);
        density_.setTracking(false);
    addProperty(seed_);
        seed_.setTracking(false);
    addProperty(scenario_);
        scenario_.addOption("cells", "Cells", LargeTestDataGeneratorInput::CELLS);
        scenario_.addOption("vessels", "Vessels", LargeTestDataGeneratorInput::VESSELS);
        scenario_.selectByValue(LargeTestDataGeneratorInput::CELLS);
    addProperty(outputVolumeNoisyFilePath_);
    addProperty(outputVolumeGTFilePath_);
    addProperty(retainLabel_);
}

LargeTestDataGeneratorInput LargeTestDataGenerator::prepareComputeInput() {
    // Reset output volume to make sure it (and the hdf5filevolume) are not used any more
    outportGT_.setData(nullptr);
    outportNoisy_.setData(nullptr);

    const std::string volumeNoisyFilePath = outputVolumeNoisyFilePath_.get();
    const std::string volumeGTFilePath = outputVolumeGTFilePath_.get();
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;

    tgt::svec3 dim = volumeDimensions_.get();

    float gaussianNoiseSD = gaussianNoiseSD_.get();

    LargeTestDataGeneratorInput::random_engine_type randomEngine {};
    randomEngine.seed(seed_.get());

    const std::string baseTypeNoisy = "uint16";
    const std::string baseTypeGT = "uint8";
    const tgt::vec3 spacing = tgt::vec3::one;
    const tgt::vec3 offset = tgt::vec3::zero;
    const RealWorldMapping rwmNoisy = RealWorldMapping::createDenormalizingMapping(baseTypeNoisy);
    const RealWorldMapping rwmGT = RealWorldMapping(1,0,"");
    const int deflateLevelNoisy = 0;
    const int deflateLevelGT = 1;

    if(volumeNoisyFilePath.empty() || volumeGTFilePath.empty()) {
        throw InvalidInputException("No volume file path specified!", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<HDF5FileVolume> outputVolumeNoisy = nullptr;
    std::unique_ptr<HDF5FileVolume> outputVolumeGT = nullptr;
    try {
        outputVolumeNoisy = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeNoisyFilePath, volumeLocation, baseTypeNoisy, dim, 1, true, deflateLevelNoisy, tgt::svec3(dim.x, dim.y, 1), false));
        outputVolumeGT = std::unique_ptr<HDF5FileVolume>(HDF5FileVolume::createVolume(volumeGTFilePath, volumeLocation, baseTypeGT, dim, 1, true, deflateLevelGT, tgt::svec3(dim.x, dim.y, 1), false));
    } catch(tgt::IOException& e) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    outputVolumeNoisy->writeSpacing(spacing);
    outputVolumeNoisy->writeOffset(offset);
    outputVolumeNoisy->writeRealWorldMapping(rwmNoisy);

    outputVolumeGT->writeSpacing(spacing);
    outputVolumeGT->writeOffset(offset);
    outputVolumeGT->writeRealWorldMapping(rwmGT);

    LINFO("Using structure size range: " << structureSizeRange_.get());
    LINFO("Using voldim: " << dim);
    LINFO("Using noise SD: " << gaussianNoiseSD);
    LINFO("Using seed: " << seed_.get());

    return LargeTestDataGeneratorInput(
        scenario_.getValue(),
        std::move(outputVolumeNoisy),
        std::move(outputVolumeGT),
        randomEngine,
        foregroundMean_.get(),
        backgroundMean_.get(),
        noiseType_.getValue(),
        gaussianNoiseSD,
        density_.get(),
        structureSizeRange_.get(),
        retainLabel_.get()
    );
}

struct Balls {
    Balls()
        : center_()
        , radius_()
    {
    }

    void add(tgt::ivec3 center, int radius) {
        center_.push_back(center);
        radius_.push_back(radius);
    }
    size_t size() const {
        tgtAssert(center_.size() == radius_.size(), "ball component sizes mismatch");
        return center_.size();
    }
    bool insideSingle(tgt::ivec3 p, size_t i) const {
        int radius = radius_[i];
        tgt::ivec3 center = center_[i];

        return tgt::distanceSq(p, center) < radius * radius;
    }

    bool inside(tgt::ivec3 p) const {
        int size = this->size();
        int ret = 0;

        auto rad = radius_.begin();
        auto cen = center_.begin();
#ifdef VRN_MODULE_OPENMP
#pragma omp simd
#endif
        for(int i=0; i < size; ++i) {
            int radius = rad[i];
            tgt::ivec3 center = cen[i];
            if(tgt::distanceSq(p, center) < radius * radius) {
                ret += 1;
            }
        }
        return ret > 0;
    }
    inline tgt::IntBounds bounds(size_t i) {
        tgt::ivec3 r(radius(i));
        tgt::ivec3 p(center(i));
        tgt::IntBounds full_bounds(p-r, p+r);
        return full_bounds;
    }
    tgt::ivec3 center(size_t i) {
        tgtAssert(i < center_.size(), "Invalid index");
        return center_[i];
    }
    int radius(size_t i) {
        tgtAssert(i < radius_.size(), "Invalid index");
        return radius_[i];
    }
    void clear() {
        center_.clear();
        radius_.clear();
    }

private:
    std::vector<tgt::ivec3> center_;
    std::vector<int> radius_;
};


static inline tgt::vec3 findOrthogonal(tgt::vec3 v) {
    tgt::vec3 orthBase1(1,0,0);
    tgt::vec3 orthBase2(0,1,0);
    tgt::vec3 orthBase = std::abs(tgt::dot(orthBase1, v)) < std::abs(tgt::dot(orthBase2, v)) ? orthBase1 : orthBase2;
    return tgt::cross(orthBase, v);
}
static inline float simdFriendlyInverseSqrt(float f) {
    // Taken from https://en.wikipedia.org/wiki/Fast_inverse_square_root, but
    // adapted to rule out undefined behaviour: Use memcpy instead of pointer/
    // union/reinterpret casts.
    const float x2 = f * 0.5F;
	const float threehalfs = 1.5F;

    uint32_t i;
    std::memcpy(&i, &f, sizeof i);
	i  = 0x5f3759df - ( i >> 1 );
    std::memcpy(&f, &i, sizeof i);
	f  *= ( threehalfs - ( x2 * f * f ) );
	return f;
}

struct Cylinders {
    Cylinders()
        : start_()
        , end_()
        , radius_()
    {
    }

    void add(tgt::ivec3 start, tgt::ivec3 end, float radius) {
        start_.push_back(start);
        end_.push_back(end);
        radius_.push_back(radius);
    }
    size_t size() const {
        tgtAssert(start_.size() == radius_.size(), "cylinder sizes mismatch");
        tgtAssert(start_.size() == end_.size(), "cylinder sizes mismatch");
        return start_.size();
    }

    bool inside(tgt::ivec3 pint) const {
        return distant(pint, 1.0);
    }
    bool insideSingle(tgt::ivec3 p, size_t i) const {
        return distantSingle(p, 1.0, i);
    }

    bool distantSingle(tgt::ivec3 pint, float radiusMult, size_t i) const {
        float radius = radiusMult * radius_[i];
        tgt::vec3 s = start_[i];
        tgt::vec3 e = end_[i];
        tgt::vec3 p = pint;

        tgt::vec3 pnorm = p - s;
        tgt::vec3 ax = e-s;
        float lenInv = simdFriendlyInverseSqrt(tgt::lengthSq(ax));
        ax = ax*lenInv;
        float len = 1.0/lenInv;

        float proj = tgt::dot(ax, pnorm);

        tgt::vec3 onDisk = pnorm - ax*proj;

        return proj >= 0 && proj <= len && tgt::lengthSq(onDisk) < radius*radius;
    }
    bool distant(tgt::ivec3 pint, float radiusMult) const {
        int size = this->size();
        int ret = 0;
        tgt::vec3 p = pint;

        auto rad = radius_.begin();
        auto start = start_.begin();
        auto end = end_.begin();
#ifdef VRN_MODULE_OPENMP
#pragma omp simd
#endif
        for(int i=0; i < size; ++i) {
            float radius = radiusMult * rad[i];
            tgt::vec3 s = start[i];
            tgt::vec3 e = end[i];

            tgt::vec3 pnorm = p - s;
            tgt::vec3 ax = e-s;
            float lenInv = simdFriendlyInverseSqrt(tgt::lengthSq(ax));
            ax = ax*lenInv;
            float len = 1.0/lenInv;

            float proj = tgt::dot(ax, pnorm);

            tgt::vec3 onDisk = pnorm - ax*proj;

            if(proj >= 0 && proj <= len && tgt::lengthSq(onDisk) < radius*radius) {
                ret += 1;
            }
        }
        return ret > 0;
    }
    inline tgt::IntBounds bounds(size_t i) {
        tgt::ivec3 r(radius(i));
        tgt::ivec3 s(start(i));
        tgt::ivec3 e(end(i));
        tgt::IntBounds full_bounds(s-r, s+r);
        full_bounds.addPoint(e-r);
        full_bounds.addPoint(e+r);
        return full_bounds;
    }
    void clear() {
        start_.clear();
        end_.clear();
        radius_.clear();
    }
    tgt::vec3 start(size_t i) {
        return start_[i];
    }
    tgt::vec3 end(size_t i) {
        return end_[i];
    }
    float radius(size_t i) {
        return radius_[i];
    }

private:
    std::vector<tgt::vec3> start_;
    std::vector<tgt::vec3> end_;
    std::vector<float> radius_;
};

static bool clipLineToBB(const tgt::Bounds& bb, tgt::vec3& p0, tgt::vec3& p1) {
    // Using the Liang–Barsky algorithm:
    tgt::vec3 delta = p1-p0;
    tgt::vec3 min = bb.getLLF();
    tgt::vec3 max = bb.getURB();

    tgt::vec3 p[] = { -delta, delta };
    tgt::vec3 q[] = { p0 - min, max - p0 };

    float umin = 0.0;
    float umax = 1.0;
    for(int j=0; j<2; ++j) {
        for(int i=0; i<3; ++i) {
            float& pi = p[j][i];
            float& qi = q[j][i];
            if(pi == 0) {
                if(qi < 0) {
                    return false;
                }
            } else if(pi < 0) {
                umin = std::max(umin, qi/pi);
            } else {
                umax = std::min(umax, qi/pi);
            }
        }
    }

    if(umin > umax) {
        return false;
    }
    p1 = p0 + umax * delta;
    p0 = p0 + umin * delta;
    return true;
}

static void initCylinders(LargeTestDataGeneratorInput& input, Balls& balls, Cylinders& cylinders, std::vector<std::vector<tgt::vec3>>& foregroundLabels, std::vector<std::vector<tgt::vec3>>& backgroundLabels, ProgressReporter& progress) {
    SubtaskProgressReporterCollection<2> progressSteps(progress, {0.5, 0.5});

    const tgt::svec3 dim = input.outputVolumeNoisy->getDimensions();
    float avgDim = tgt::hadd(dim)/3.0f;
    int minRadius = std::max(1, input.structureSizeRange.x/2);
    int maxRadius = std::max(1, input.structureSizeRange.y/2);

    const float maxRelativeLen = 2.0f;
    const float deviationAngleStd = tgt::PIf * 0.1;
    const size_t maxTries = 10;
    const float radiusToBaseLenFactor = 5.0f;
    std::uniform_real_distribution<> lDistr(0.5f, 2.0f);

    std::uniform_int_distribution<> xDistr(0, dim.x-1);
    std::uniform_int_distribution<> yDistr(0, dim.y-1);
    std::uniform_int_distribution<> zDistr(0, dim.z-1);

    size_t totalVolume = tgt::hmul(dim);
    float volumeToFill = totalVolume * input.density;
    float remainingVolume = volumeToFill;


    tgt::Bounds volumeBounds(tgt::vec3(-0.5), tgt::vec3(dim) - tgt::vec3(0.5));

    struct BranchSeed {
        tgt::vec3 begin;
        tgt::vec3 dir;
        float radius;
        float baseLen;
    };
    std::deque<BranchSeed> queue;


    while(remainingVolume > 0) {
        progressSteps.get<0>().setProgress(remainingVolume/volumeToFill);

        tgt::vec3 start;
        tgt::vec3 dir;
        size_t tr = 0;
        do {
            ++tr;
            if(tr > maxTries) {
                break;
            }

            start = tgt::vec3(xDistr(input.randomEngine), yDistr(input.randomEngine), zDistr(input.randomEngine));
            dir = tgt::vec3::zero;
            size_t dimension = std::uniform_int_distribution<>(0, 2)(input.randomEngine);
            size_t side = std::uniform_int_distribution<>(0, 1)(input.randomEngine);
            if(side == 0) {
                start[dimension] = 0.0f;
                dir[dimension] = 1.0;
            } else {
                start[dimension] = dim[dimension]-1;
                dir[dimension] = -1.0;
            }
        } while(cylinders.distant(start, 2.0f));
        if(tr > maxTries) {
            LWARNINGC(LargeTestDataGenerator::loggerCat_, "Failed to fulfill density request");
            break;
        }
        float rootRadiusLog = std::uniform_real_distribution<>(std::log(minRadius), std::log(maxRadius))(input.randomEngine);
        float rootRadius = tgt::clamp(std::exp(rootRadiusLog), float(minRadius), float(maxRadius));
        float rootBaseLen = std::uniform_real_distribution<>(0, rootRadius * radiusToBaseLenFactor)(input.randomEngine);
        balls.add(start, rootRadius);
        queue.push_back(BranchSeed {
            start,
            dir,
            rootRadius,
            rootBaseLen,
        });

        while(!queue.empty()) {
            auto seed = queue.front();
            queue.pop_front();

            tgt::vec3 end;
            tgt::vec3 axis;
            size_t tr = 0;
            do {
                ++tr;
                if(tr > maxTries) {
                    break;
                }
                float deviationAngle = std::normal_distribution<>(tgt::PIf * 0.1, deviationAngleStd)(input.randomEngine);
                float rotationAngle = std::uniform_real_distribution<>(0.0, 2 * tgt::PIf)(input.randomEngine);
                float len = std::max(lDistr(input.randomEngine) * seed.baseLen, 2.0*seed.radius);

                axis = seed.dir;
                tgt::vec3 orth = tgt::normalize(findOrthogonal(axis));

                orth = tgt::Quaternion<float>::rotate(orth, rotationAngle, axis);
                axis = tgt::Quaternion<float>::rotate(axis, deviationAngle, orth);
                end = seed.begin + axis * len;
                tgtAssert(!std::isnan(end.x) && !std::isnan(end.x) && !std::isnan(end.x), "Invalid end calculated");
            } while(cylinders.distant(end, 2.0f));
            if(tr > maxTries) {
                continue;
            }

            tgt::vec3 p1 = seed.begin;
            tgt::vec3 p2 = end;
            float vol;
            if(clipLineToBB(volumeBounds, p1, p2)) {
                if(tgt::distance(p1, p2) < 2.0f) {
                    continue;
                }
                vol = tgt::distance(p1, p2) * seed.radius * seed.radius * tgt::PIf;
            } else {
                LWARNINGC(LargeTestDataGenerator::loggerCat_, "centerline completely outside of volume");
                continue;
            }
            balls.add(end, seed.radius);
            cylinders.add(seed.begin, end, seed.radius);
            remainingVolume -= vol;

            if(remainingVolume < 0) {
                break;
            }

            if(!volumeBounds.containsPoint(end)) {
                continue;
            }

            float crosssection = seed.radius*seed.radius;
            std::uniform_real_distribution<> cDistr(crosssection/10.0, crosssection/2.0);

            float c1 = cDistr(input.randomEngine);
            float c2 = crosssection - c1;
            float childRadii[] = {std::sqrt(c1), std::sqrt(c2)};
            for(float childRadius : childRadii) {
                if(childRadius > minRadius) {
                    queue.push_back(BranchSeed {
                        end,
                        axis,
                        childRadius,
                        float(childRadius * radiusToBaseLenFactor),
                    });
                }
            }
        }
    }
    LINFOC(LargeTestDataGenerator::loggerCat_, "Placed " << cylinders.size() << " Cylinders");

    std::vector<std::pair<size_t, tgt::IntBounds>> ballBounds;
    for(size_t i=0; i < balls.size(); ++i) {
        ballBounds.emplace_back(i, balls.bounds(i));
    }
    auto maybeBallHierarchy = BoundsHierarchy<int, size_t>::buildTopDown(std::move(ballBounds));
    if(!maybeBallHierarchy) {
        return;
    }
    auto& ballHierarchy = *maybeBallHierarchy;

    std::vector<std::pair<size_t, tgt::IntBounds>> cylinderBounds;
    for(size_t i=0; i < cylinders.size(); ++i) {
        cylinderBounds.emplace_back(i, cylinders.bounds(i));
    }
    auto maybecylinderHierarchy = BoundsHierarchy<int, size_t>::buildTopDown(std::move(cylinderBounds));
    if(!maybecylinderHierarchy) {
        return;
    }
    auto& cylinderHierarchy = *maybecylinderHierarchy;

    auto invalid = [&] (const tgt::vec3& p) {
        if(!volumeBounds.containsPoint(p)) {
            return true;
        }
        tgt::ivec3 pi = tgt::round(p);
        {
            auto indices = ballHierarchy.findBounds(p);
            for(auto i : indices) {
                if(balls.insideSingle(p, i)) {
                    return true;
                }
            }
        }
        {
            auto indices = cylinderHierarchy.findBounds(p);
            for(auto i : indices) {
                if(cylinders.insideSingle(p, i)) {
                    return true;
                }
            }
        }
        return false;
        //return cylinders.inside(p) || balls.inside(p) || !volumeBounds.containsPoint(p);
    };

    size_t numElements = cylinders.size();
    for(int i=0; i<numElements; ++i) {
        progressSteps.get<1>().setProgress(static_cast<float>(i)/numElements);
        auto s = cylinders.start(i);
        auto e = cylinders.end(i);

        bool clipSuccess = clipLineToBB(volumeBounds, s, e);
        tgtAssert(clipSuccess, "clip must have already been successful earlier");

        float radius = cylinders.radius(i);

        int seedsAlongLine = 2;
        int seedsPerPointFg = 1;
        int seedsPerPointBg = 1;

        int max_tries = 10;
        int tries = max_tries;

        tgt::vec3 axis = e - s;
        tgt::vec3 orth = tgt::normalize(findOrthogonal(axis));
        for(int j=0; j<seedsAlongLine; ++j) {
            float alpha = static_cast<float>(j+1)/(seedsAlongLine+1);

            const tgt::vec3 p = alpha * s + (1.0f - alpha) * e;

            for(int j=0; j<seedsPerPointFg; ++j) {
                std::vector<tgt::vec3> segmentFg;
                segmentFg.push_back(p);
                segmentFg.push_back(tgt::vec3(p)+tgt::vec3(0.001));
                foregroundLabels.push_back(segmentFg);
            }

            for(int j=0; j<seedsPerPointBg; ++j) {
                for(int tr=0; tr<max_tries; ++tr) {

                    float rotationAngle = std::uniform_real_distribution<>(0.0, 2 * tgt::PIf)(input.randomEngine);
                    tgt::vec3 offset = 2.0f * radius * tgt::Quaternion<float>::rotate(orth, rotationAngle, axis);

                    tgt::vec3 po = p + offset;

                    if(invalid(po)) {
                        continue;
                    }

                    std::vector<tgt::vec3> segmentBg;
                    segmentBg.push_back(po);
                    segmentBg.push_back(tgt::vec3(po)+tgt::vec3(0.001));
                    backgroundLabels.push_back(segmentBg);
                    break;
                }
            }
        }
    }
}

static void initCells(LargeTestDataGeneratorInput& input, Balls& balls, Cylinders& cylinders, std::vector<std::vector<tgt::vec3>>& foregroundLabels, std::vector<std::vector<tgt::vec3>>& backgroundLabels, ProgressReporter& progress) {
    SubtaskProgressReporterCollection<2> progressSteps(progress, {0.5, 0.5});

    const tgt::svec3 dim = input.outputVolumeNoisy->getDimensions();
    int minRadius = std::max(1, input.structureSizeRange.x/2);
    int maxRadius = std::max(1, input.structureSizeRange.y/2);

    size_t elementVolumeEstimate;

    float a = minRadius;
    float b = maxRadius;
    if(a == b) {
        elementVolumeEstimate = tgt::round(4.0/3.0 * tgt::PIf * a*a*a);
    } else {
        elementVolumeEstimate = tgt::round(tgt::PIf*(b*b*b*b-a*a*a*a)/(3 * (b - a)));
    }

    size_t totalVolume = tgt::hmul(dim);
    size_t numElements = std::round(totalVolume/elementVolumeEstimate * input.density);
    LINFOC(LargeTestDataGenerator::loggerCat_, "Placing " << numElements << " Spheres");

    std::uniform_int_distribution<> xDistr(0, dim.x-1);
    std::uniform_int_distribution<> yDistr(0, dim.y-1);
    std::uniform_int_distribution<> zDistr(0, dim.z-1);
    std::uniform_int_distribution<> rDistr(minRadius, maxRadius);

    for(int i=0; i<numElements; ++i) {
        progressSteps.get<0>().setProgress(static_cast<float>(i)/numElements);
        tgt::ivec3 p(xDistr(input.randomEngine), yDistr(input.randomEngine), zDistr(input.randomEngine));
        balls.add(p, rDistr(input.randomEngine));

        std::vector<tgt::vec3> segment;
        segment.push_back(p);
        segment.push_back(tgt::vec3(p)+tgt::vec3(0.001));
        foregroundLabels.push_back(segment);
    }

    std::vector<std::pair<size_t, tgt::IntBounds>> bounds;
    for(size_t i=0; i < balls.size(); ++i) {
        bounds.emplace_back(i, balls.bounds(i));
    }
    auto maybeHierarchy = BoundsHierarchy<int, size_t>::buildTopDown(std::move(bounds));
    if(!maybeHierarchy) {
        return;
    }
    auto& hierarchy = *maybeHierarchy;

    tgt::Bounds volumeBounds(tgt::vec3(-0.5), tgt::vec3(dim) - tgt::vec3(0.5));

    auto invalid = [&] (tgt::ivec3 p) {
        if(!volumeBounds.containsPoint(p)) {
            return true;
        }
        auto indices = hierarchy.findBounds(p);
        for(auto i : indices) {
            if(balls.insideSingle(p, i)) {
                return true;
            }
        }
        return false;
    };

    int max_tries = 10;
    std::uniform_real_distribution<float> dirDistr(-1.0f, 1.0f);

    size_t failures = 0;
    for(int i=0; i<numElements; ++i) {
        progressSteps.get<1>().setProgress(static_cast<float>(i)/numElements);
        tgt::vec3 p = balls.center(i);
        float radius = balls.radius(i);

        int seedsPerCell = 1;

        for(int j=0; j<seedsPerCell; ++j) {
            int tr;
            for(tr=0; tr<max_tries; ++tr) {
                tgt::vec3 dir(
                        dirDistr(input.randomEngine),
                        dirDistr(input.randomEngine),
                        dirDistr(input.randomEngine)
                        );
                if(tgt::lengthSq(dir) == 0) {
                    continue;
                }
                dir = tgt::normalize(dir);
                tgt::vec3 offset = 2.0f * radius * dir;

                p += offset;

                if(invalid(p)) {
                    continue;
                }

                std::vector<tgt::vec3> segmentBg;
                segmentBg.push_back(p);
                segmentBg.push_back(tgt::vec3(p)+tgt::vec3(0.001));
                backgroundLabels.push_back(segmentBg);
                break;
            }
            if(tr == max_tries) {
                ++failures;
            }
        }
    }

    if(failures > 0) {
        LWARNINGC(LargeTestDataGenerator::loggerCat_, "Failed to position " << failures << "background seeds");
    }
}

template<typename NoiseTypeImpl>
static void rasterize(
        LargeTestDataGeneratorInput& input,
        NoiseTypeImpl noise,
        Balls& balls,
        Cylinders& cylinders,
        ProgressReporter& progress
        ) {
    const tgt::svec3 dim = input.outputVolumeNoisy->getDimensions();
    tgtAssert(input.outputVolumeGT->getDimensions() == dim, "Dimension mismatch");

    std::mutex intervalWalkerMutex;
    std::vector<Interval<int, size_t>> ballIntervals;
    for(size_t i=0; i < balls.size(); ++i) {
        auto full_bounds = balls.bounds(i);
        ballIntervals.emplace_back(full_bounds.getLLF().z, full_bounds.getURB().z, i);
    }
    IntervalWalker<int, size_t> ballWalker(0, std::move(ballIntervals));

    std::vector<Interval<int, size_t>> cylinderIntervals;
    for(size_t i=0; i < cylinders.size(); ++i) {
        auto full_bounds = cylinders.bounds(i);
        cylinderIntervals.emplace_back(full_bounds.getLLF().z, full_bounds.getURB().z, i);
    }
    IntervalWalker<int, size_t> cylinderWalker(0, std::move(cylinderIntervals));

    std::uniform_int_distribution<uint64_t> sliceBaseSeedDistr(0, std::numeric_limits<uint64_t>::max());
    uint64_t baseSeed = sliceBaseSeedDistr(input.randomEngine);

    ThreadedTaskProgressReporter parallelProgress(progress, dim.z);

//#undef VRN_MODULE_OPENMP

    bool aborted = false;
#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for(int zDimWork=0; zDimWork<dim.z; ++zDimWork) {
#ifdef VRN_MODULE_OPENMP
        if (aborted) {
            continue;
        }
#endif
        std::vector<Interval<int, size_t>> ballIntervals;
        std::vector<Interval<int, size_t>> cylinderIntervals;

        int z;
        {
            std::lock_guard<std::mutex> guard(intervalWalkerMutex);

            auto ballIt = ballWalker.next();
            for(auto it = ballIt.next(); it != ballIt.end(); it = ballIt.next()) {
                auto& i = it->value;
                auto full_bounds = balls.bounds(i);
                ballIntervals.emplace_back(full_bounds.getLLF().y, full_bounds.getURB().y, i);
            }

            auto cylinderIt = cylinderWalker.next();
            for(auto it = cylinderIt.next(); it != cylinderIt.end(); it = cylinderIt.next()) {
                auto& i = it->value;
                auto full_bounds = cylinders.bounds(i);
                cylinderIntervals.emplace_back(full_bounds.getLLF().y, full_bounds.getURB().y, i);
            }
            z = ballIt.currentPos();
            tgtAssert(z == cylinderIt.currentPos(), "IntervalWalker pos mismatch");
        }

        IntervalWalker<int, size_t> ballWalker(0, std::move(ballIntervals));
        IntervalWalker<int, size_t> cylinderWalker(0, std::move(cylinderIntervals));

        pcg32_state randomEngine = pcg32_init(baseSeed + z);

        VolumeAtomic<uint16_t> sliceNoisy(tgt::vec3(dim.x, dim.y, 1));
        VolumeAtomic<uint8_t> sliceGT(tgt::vec3(dim.x, dim.y, 1));


        for(int y=0; y<dim.y; ++y) {

            auto ballIt = ballWalker.next();
            std::vector<Interval<int, size_t>> ballIntervals;
            for(auto it = ballIt.next(); it != ballIt.end(); it = ballIt.next()) {
                auto& i = it->value;
                auto full_bounds = balls.bounds(i);
                ballIntervals.emplace_back(full_bounds.getLLF().x, full_bounds.getURB().x, i);
            }
            IntervalWalker<int, size_t> ballWalker(0, std::move(ballIntervals));

            auto cylinderIt = cylinderWalker.next();
            std::vector<Interval<int, size_t>> cylinderIntervals;
            for(auto it = cylinderIt.next(); it != cylinderIt.end(); it = cylinderIt.next()) {
                auto& i = it->value;
                auto full_bounds = cylinders.bounds(i);
                cylinderIntervals.emplace_back(full_bounds.getLLF().x, full_bounds.getURB().x, i);
            }
            IntervalWalker<int, size_t> cylinderWalker(0, std::move(cylinderIntervals));

            for(int x=0; x<dim.x; ++x) {
                tgt::ivec3 p(x,y,z);
                size_t slicePos = y*dim.x+x;

                bool inside = false;
                auto ballIt = ballWalker.next();
                for(auto it = ballIt.next(); it != ballIt.end(); it = ballIt.next()) {
                    auto& i = it->value;
                    inside |= balls.insideSingle(p, i);
                }
                auto cylinderIt = cylinderWalker.next();
                for(auto it = cylinderIt.next(); it != cylinderIt.end(); it = cylinderIt.next()) {
                    auto& i = it->value;
                    inside |= cylinders.insideSingle(p, i);
                }
                //bool inside = (balls.inside(p)) || cylinders.inside(p);
                uint16_t val = inside ? input.foregroundMean : input.backgroundMean;

                sliceNoisy.voxel(slicePos) = noise.apply(val, randomEngine);

                sliceGT.voxel(slicePos) = inside ? 255 : 0;
            }
        }
        input.outputVolumeNoisy->writeSlices(&sliceNoisy, z);
        input.outputVolumeGT->writeSlices(&sliceGT, z);

        if(parallelProgress.reportStepDone()) {
#ifdef VRN_MODULE_OPENMP
            #pragma omp critical
            aborted = true;
#else
            aborted = true;
            break;
#endif
        }
    }
    if(aborted) {
        throw boost::thread_interrupted();
    }
}

struct GaussianNoiseGenerator {
    float standardDeviation_;

    inline uint16_t apply(uint16_t input, pcg32_state& randomEngine) const {
        int val = input;
        val += tgt::fastround(gsl_ran_gaussian_ziggurat(randomEngine, standardDeviation_));
        return tgt::clamp(val, 0, 0xffff);
    };
};
struct PoissonNoiseGenerator {
    inline uint16_t apply(uint16_t input, pcg32_state& randomEngine) const {
        if(input < 50) {
            if(input == 0) {
                return 0;
            }
            int x = 0;
            float lambda = static_cast<float>(input);
            float p = std::exp(-lambda);
            float s = p;

            uint32_t rint = pcg32(randomEngine);
            float u = static_cast<float>(rint);
            while(u > s * static_cast<float>(0xffffffff)) {
                ++x;
                p *= lambda/x;
                s += p;
            }
            return x;
        } else {
            int val = input;
            float std = std::sqrt(static_cast<float>(input));
            val += tgt::fastround(gsl_ran_gaussian_ziggurat(randomEngine, std));
            return tgt::clamp(val, 0, 0xffff);
        }
    };
};


LargeTestDataGeneratorOutput LargeTestDataGenerator::compute(LargeTestDataGeneratorInput input, ProgressReporter& progressReporter) const {
    SubtaskProgressReporterCollection<2> globalProgressSteps(progressReporter, {0.02, 0.98});

    tgtAssert(input.outputVolumeNoisy, "No outputVolume");
    tgtAssert(input.outputVolumeGT, "No outputVolume");
    const tgt::svec3 dim = input.outputVolumeNoisy->getDimensions();
    tgtAssert(input.outputVolumeGT->getDimensions() == dim, "Dimension mismatch");
    Balls balls{};
    Cylinders cylinders{};


    std::vector<std::vector<tgt::vec3>> foregroundLabels{};
    std::vector<std::vector<tgt::vec3>> backgroundLabels{};

    auto& initProgress = globalProgressSteps.template get<0>();
    switch(input.scenario) {
        case LargeTestDataGeneratorInput::CELLS: {
            initCells(input, balls, cylinders, foregroundLabels, backgroundLabels, initProgress);
            break;
        }
        case LargeTestDataGeneratorInput::VESSELS: {
            initCylinders(input, balls, cylinders, foregroundLabels, backgroundLabels, initProgress);
            break;
        }
        default: {
            tgtAssert(false, "Invalid scenario");
        }
    }

    auto& p = globalProgressSteps.template get<1>();
    switch(input.noiseType) {
        case LargeTestDataGeneratorInput::GAUSSIAN:
            rasterize(input, GaussianNoiseGenerator { input.gaussianNoiseSD }, balls, cylinders, p);
            break;
        case LargeTestDataGeneratorInput::POISSON:
            rasterize(input, PoissonNoiseGenerator {}, balls, cylinders, p);
            break;
    }

    if(input.retainLabel) {
        int tries = 100;
        bool retained = false;
        auto& labels = foregroundLabels;
        std::uniform_int_distribution<uint64_t> indexDistr(0, labels.size()-1);
        while(tries > 0) {
            uint64_t index = indexDistr(input.randomEngine);
            auto& segment = labels.at(index);
            if(!segment.empty()) {
                segment.pop_back();
                break;
            } else {
                --tries;
            }
        }
        if(tries == 0) {
            LWARNINGC("LargeTestDataGenerator", "Failed to retain a background label!");
        }
    }


    std::unique_ptr<PointSegmentListGeometryVec3> fg(new PointSegmentListGeometryVec3());
    fg->setData(foregroundLabels);
    std::unique_ptr<PointSegmentListGeometryVec3> bg(new PointSegmentListGeometryVec3());
    bg->setData(backgroundLabels);
    return {
        std::move(fg),
        std::move(bg),
        input.outputVolumeNoisy->getFileName(),
        input.outputVolumeGT->getFileName(),
    };
    //outputVolume will be destroyed and thus closed now.
}

void LargeTestDataGenerator::processComputeOutput(LargeTestDataGeneratorOutput output) {
    // outputVolume has been destroyed and thus closed by now.
    // So we can open it again (and use HDF5VolumeReader's implementation to read all the metadata with the file)
    const VolumeBase* volNoisy = HDF5VolumeReader().read(output.outputVolumeNoisyFilePath)->at(0);
    const VolumeBase* volGT = HDF5VolumeReader().read(output.outputVolumeGTFilePath)->at(0);

    outportNoisy_.setData(volNoisy);
    outportGT_.setData(volGT);

    foregroundLabelsPort_.setData(output.foregroundLabels.release());
    backgroundLabelsPort_.setData(output.backgroundLabels.release());
}

bool LargeTestDataGenerator::isReady() const {
    if(!isInitialized()) {
        setNotReadyErrorMessage("Not initialized.");
        return false;
    }
    return true;
}

LargeTestDataGenerator::~LargeTestDataGenerator() {
}
VoreenSerializableObject* LargeTestDataGenerator::create() const {
    return new LargeTestDataGenerator();
}

} // namespace voreen
