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

#include "octreewalker.h"

#include "../solver/randomwalkerweights.h"
#include "../util/preprocessing.h"
#include "../util/seeds.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormorphology.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresample.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatornumsignificant.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanagermmap.h"
#include "voreen/core/datastructures/octree/volumeoctreenodegeneric.h"
#include "voreen/core/utils/hashing.h"
#include "tgt/vector.h"
#include "tgt/memory.h"

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include <climits>
#include <random>

#ifdef VRN_MODULE_OPENMP
#include <omp.h>
#endif

namespace voreen {

#if defined(VRN_MODULE_OPENMP) && 1
#define VRN_OCTREEWALKER_USE_OMP
#endif

namespace {

const std::string PREV_RESULT_FILE_NAME = "prev_result.vvod";
const std::string PREV_RESULT_VALID_FILE_NAME = "prev_result_valid.token";
const std::string PREV_RESULT_OCTREE_KEY = "Octree";
const std::string PREV_RESULT_FOREGROUND_KEY = "foregroundSeeds";
const std::string PREV_RESULT_BACKGROUND_KEY = "backgroundSeeds";
const std::string BRICK_BUFFER_SUBDIR =      "brickBuffer";
const std::string BRICK_BUFFER_FILE_PREFIX = "buffer_";

static inline size_t volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

static inline size_t volumeCoordsToIndex(const tgt::ivec3& coords, const tgt::ivec3& dim) {
    return coords.z*dim.y*dim.x + coords.y*dim.x + coords.x;
}

static void freeTreeComponents(VolumeOctreeNode* root, std::unordered_set<const VolumeOctreeNode*>& nodesToSave, OctreeBrickPoolManagerBase& brickPoolManager) {
    if(root && nodesToSave.find(root) == nodesToSave.end()) {
        if(root->hasBrick()) {
            brickPoolManager.deleteBrick(root->getBrickAddress());
        }
        for(int i=0; i<8; ++i) {
            freeTreeComponents(root->children_[i], nodesToSave, brickPoolManager);
        }
        delete root;
    }
}

static void freeNodes(VolumeOctreeNode* root) {
    if(root) {
        for(int i=0; i<8; ++i) {
            freeNodes(root->children_[i]);
        }
    }
    delete root;
}

static uint16_t absdiff(uint16_t v1, uint16_t v2) {
    return std::max(v1, v2) - std::min(v1, v2);
}
static float absdiffRel(uint16_t v1, uint16_t v2) {
    return static_cast<float>(absdiff(v1, v2))/0xffff;
}

static uint16_t normToBrick(float val) {
    return tgt::clamp(val, 0.0f, 1.0f) * 0xffff;
}
static float brickToNorm(uint16_t val) {
    const float NORM_TO_BRICK_FACTOR = 1.0f/0xffff;
    return static_cast<float>(val)*NORM_TO_BRICK_FACTOR;
}

static constexpr size_t mask(size_t width) {
    return (1 << width) - 1;
}

static const size_t PACKING_EXPONENT_BITS = 7;
static const size_t PACKING_DIGIT_BITS = 16 - PACKING_EXPONENT_BITS;

static uint16_t packNormalizedFloat(float value) {
    tgtAssert(0.0f <= value && value < 1.0f, "Value not normalized");
    if(!std::isnormal(value)) {
        return 0;
    }

    uint32_t v;
    std::memcpy(&v, &value, sizeof(v));

    uint16_t exponent_bits = (v >> 23) & mask(PACKING_EXPONENT_BITS);
    uint16_t digit_bits = (v >> (23-PACKING_DIGIT_BITS)) & mask(PACKING_DIGIT_BITS);

    uint16_t exp = exponent_bits << PACKING_DIGIT_BITS;
    return exp | digit_bits;
}

static float unpackNormalizedFloat(uint16_t bits) {
    uint16_t exponent_bits = (bits >> PACKING_DIGIT_BITS) & mask(PACKING_EXPONENT_BITS);
    uint16_t digit_bits = bits & mask(PACKING_DIGIT_BITS);

    uint32_t ieee_bits = digit_bits << (23-PACKING_DIGIT_BITS) | exponent_bits << 23;
    float v;
    std::memcpy(&v, &ieee_bits, sizeof(v));

    tgtAssert(0 <= v && v < 1, "Unpacked invalid value");
    return v;
}

static bool hasExtentParameter(RWNoiseModel n) {
    switch(n) {
        case RW_NOISE_GAUSSIAN:
        case RW_NOISE_POISSON:
        case RW_NOISE_VARIABLE_GAUSSIAN:
        case RW_NOISE_TTEST:
        case RW_NOISE_GAUSSIAN_HIERARCHICAL:
        case RW_NOISE_POISSON_HIERARCHICAL:
        case RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL:
        case RW_NOISE_TTEST_HIERARCHICAL:
            return true;
        case RW_NOISE_GAUSSIAN_BIAN_MEAN:
        case RW_NOISE_GAUSSIAN_BIAN_MEDIAN:
            return false;
    }
}

static bool requiresVarianceTree(RWNoiseModel n) {
    switch(n) {
        case RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL:
        case RW_NOISE_TTEST_HIERARCHICAL:
            return true;
        case RW_NOISE_GAUSSIAN_BIAN_MEAN:
        case RW_NOISE_GAUSSIAN_BIAN_MEDIAN:
        case RW_NOISE_GAUSSIAN:
        case RW_NOISE_POISSON:
        case RW_NOISE_VARIABLE_GAUSSIAN:
        case RW_NOISE_TTEST:
        case RW_NOISE_GAUSSIAN_HIERARCHICAL:
        case RW_NOISE_POISSON_HIERARCHICAL:
            return false;
    }
}

static bool requiresGlobalVarianceEstimate(RWNoiseModel n) {
    switch(n) {
        case RW_NOISE_GAUSSIAN_HIERARCHICAL:
            return true;
        case RW_NOISE_GAUSSIAN_BIAN_MEAN:
        case RW_NOISE_GAUSSIAN_BIAN_MEDIAN:
        case RW_NOISE_GAUSSIAN:
        case RW_NOISE_POISSON:
        case RW_NOISE_VARIABLE_GAUSSIAN:
        case RW_NOISE_TTEST:
        case RW_NOISE_POISSON_HIERARCHICAL:
        case RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL:
        case RW_NOISE_TTEST_HIERARCHICAL:
            return false;
    }
}
}

OctreeWalkerOutput::OctreeWalkerOutput(
    OctreeWalkerPreviousResult&& result,
    std::unordered_set<const VolumeOctreeNode*>&& sharedNodes,
    VarianceTree&& varianceTree,
    std::chrono::duration<float> duration
)   : result_(std::move(result))
    , sharedNodes_(std::move(sharedNodes))
    , varianceTree_(std::move(varianceTree))
    , duration_(duration)
{ }

OctreeWalkerOutput::OctreeWalkerOutput(OctreeWalkerOutput&& other)
    : result_(std::move(other.result_))
    , sharedNodes_(std::move(other.sharedNodes_))
    , duration_(other.duration_)
    , varianceTree_(std::move(other.varianceTree_))
    , movedOut_(other.movedOut_)
{
}

OctreeWalkerOutput::~OctreeWalkerOutput() {
    // If OctreeWalkerOutput is dropped prematurely (e.g., if the computation
    // is interrupted after the final result was created, we need to make sure
    // NOT to drop the Brickpoolmanager, BUT free all (new) nodes of the (now
    // discarded) tree.
    if(result_.isPresent()) {
        std::move(result_).destroyButRetainNodes(sharedNodes_);
    }
}


const std::string OctreeWalker::loggerCat_("voreen.RandomWalker.OctreeWalker");

OctreeWalker::OctreeWalker()
    : AsyncComputeProcessor<OctreeWalkerInput, OctreeWalkerOutput>()
    , inportVolume_(Port::INPORT, "volume.input")
    , inportForegroundSeeds_(Port::INPORT, "geometry.seedsForeground", "geometry.seedsForeground", true)
    , inportBackgroundSeeds_(Port::INPORT, "geometry.seedsBackground", "geometry.seedsBackground", true)
    , outportProbabilities_(Port::OUTPORT, "volume.probabilities", "volume.probabilities", false)
    , noiseModel_("noiseModel", "Noise Model")
    , minEdgeWeight_("minEdgeWeight", "Min Edge Weight: 10^(-t)", 5, 0, 10)
    , parameterEstimationNeighborhoodExtent_("parameterEstimationNeighborhoodExtent", "Extent", 1, 1, 3)
    , preconditioner_("preconditioner", "Preconditioner")
    , errorThreshold_("errorThreshold", "Error Threshold: 10^(-t)", 2, 0, 10)
    , maxIterations_("conjGradIterations", "Max Iterations", 1000, 1, 5000)
    , conjGradImplementation_("conjGradImplementation", "Implementation")
    , homogeneityThreshold_("homogeneityThreshold", "Homogeneity Threshold", 0.01, 0.0, 1.0)
    , binaryPruningDelta_("binaryPruningDelta", "Binary Pruning Delta", 0.01, 0.0, 0.5)
    , incrementalSimilarityThreshold_("incrementalSimilarityThreshold", "Incremental Similarity Treshold", 0.01, 0.0, 1.0)
    , clearResult_("clearResult", "Clear Result", Processor::INVALID_RESULT, Property::LOD_DEBUG)
    , resultPath_("resultPath", "Result Cache Path", "Result Cache Path", "", "", FileDialogProperty::DIRECTORY, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , prevResultPath_("")
    , previousResult_(boost::none)
    , varianceTree_(VarianceTree::none())
    , brickPoolManager_(nullptr)
{
    // ports
    addPort(inportVolume_);
    addPort(inportForegroundSeeds_);
    addPort(inportBackgroundSeeds_);
    addPort(outportProbabilities_);

    // random walker properties
    addProperty(noiseModel_);
        noiseModel_.addOption("gaussian", "Gaussian with constant variance (Drees 2022)", RW_NOISE_GAUSSIAN);
        noiseModel_.addOption("shot", "Poisson (Drees 2022)", RW_NOISE_POISSON);
        noiseModel_.addOption("variable_gaussian", "Gaussian with signal-dependent variance (Drees 2022)", RW_NOISE_VARIABLE_GAUSSIAN);
        noiseModel_.addOption("gaussian_bian", "Gaussian with mean filter (Bian 2016)", RW_NOISE_GAUSSIAN_BIAN_MEAN);
        noiseModel_.addOption("gaussian_bian_median", "Gaussian with median filter (Bian 2016)", RW_NOISE_GAUSSIAN_BIAN_MEDIAN);
        noiseModel_.addOption("ttest", "T-Test (Bian 2016)", RW_NOISE_TTEST);
        //noiseModel_.addOption("gaussian_hierarchical", "Gaussian with constant variance (Drees 2022), hierarchical", RW_NOISE_GAUSSIAN_HIERARCHICAL); // This really does not work too well :(
        noiseModel_.addOption("shot_hierarchical", "Poisson (Drees 2022), hierarchical", RW_NOISE_POISSON_HIERARCHICAL);
        noiseModel_.addOption("variable_gaussian_hierarchical", "Gaussian with signal-dependent variance (Drees 2022), hierarchical", RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL);
        noiseModel_.addOption("ttest_hierarchical", "T-Test (Bian 2016), hierarchical", RW_NOISE_TTEST_HIERARCHICAL);
        noiseModel_.selectByValue(RW_NOISE_GAUSSIAN);
        noiseModel_.setGroupID("rwparam");
        ON_CHANGE_LAMBDA(noiseModel_, [this] () {
            RWNoiseModel m = noiseModel_.getValue();
            parameterEstimationNeighborhoodExtent_.setVisibleFlag(hasExtentParameter(m));
        });
    addProperty(minEdgeWeight_);
        minEdgeWeight_.setGroupID("rwparam");
        minEdgeWeight_.setTracking(false);
    addProperty(parameterEstimationNeighborhoodExtent_);
        parameterEstimationNeighborhoodExtent_.setGroupID("rwparam");
        parameterEstimationNeighborhoodExtent_.setTracking(false);
    addProperty(homogeneityThreshold_);
        homogeneityThreshold_.setGroupID("rwparam");
        homogeneityThreshold_.adaptDecimalsToRange(5);
        homogeneityThreshold_.setTracking(false);
    addProperty(binaryPruningDelta_);
        binaryPruningDelta_.setGroupID("rwparam");
        binaryPruningDelta_.adaptDecimalsToRange(5);
        binaryPruningDelta_.setTracking(false);
    addProperty(incrementalSimilarityThreshold_);
        incrementalSimilarityThreshold_.setGroupID("rwparam");
        incrementalSimilarityThreshold_.adaptDecimalsToRange(5);
        incrementalSimilarityThreshold_.setTracking(false);
    setPropertyGroupGuiName("rwparam", "Random Walker Parametrization");

    // conjugate gradient solver
    addProperty(preconditioner_);
        preconditioner_.addOption("none", "None");
        preconditioner_.addOption("jacobi", "Jacobi");
        preconditioner_.select("jacobi");
        preconditioner_.setGroupID("conjGrad");
    addProperty(errorThreshold_);
        errorThreshold_.setGroupID("conjGrad");
        errorThreshold_.setTracking(false);
    addProperty(maxIterations_);
        maxIterations_.setGroupID("conjGrad");
        maxIterations_.setTracking(false);
    addProperty(conjGradImplementation_);
        conjGradImplementation_.addOption("blasCPU", "CPU");
#ifdef VRN_MODULE_OPENMP
        conjGradImplementation_.addOption("blasMP", "OpenMP");
        conjGradImplementation_.select("blasMP");
#endif
#ifdef VRN_MODULE_OPENCL
        conjGradImplementation_.addOption("blasCL", "OpenCL");
        conjGradImplementation_.select("blasCL");
#endif
        conjGradImplementation_.setGroupID("conjGrad");
    setPropertyGroupGuiName("conjGrad", "Conjugate Gradient Solver");

    // conjugate gradient solver
    addProperty(clearResult_);
        ON_CHANGE_LAMBDA(clearResult_, [this] () {
                interruptComputation();
                clearPreviousResults();
                });
        clearResult_.setGroupID("resultcache");
    addProperty(resultPath_);
        ON_CHANGE_LAMBDA(resultPath_, [this] () {
                if(resultPath_.get() != prevResultPath_) {
                    interruptComputation();
                    clearPreviousResults();
                }
                });
        resultPath_.setGroupID("resultcache");
    setPropertyGroupGuiName("resultcache", "Result Cache");
}

OctreeWalker::~OctreeWalker() {
}

Processor* OctreeWalker::create() const {
    return new OctreeWalker();
}

void OctreeWalker::initialize() {
    AsyncComputeProcessor::initialize();

#ifdef VRN_MODULE_OPENCL
    voreenBlasCL_.initialize();
#endif
}

void OctreeWalker::deinitialize() {
    interruptComputation();
    clearPreviousResults();

    AsyncComputeProcessor::deinitialize();
}

bool OctreeWalker::isReady() const {
    bool ready = false;
    ready |= outportProbabilities_.isConnected();
    ready &= inportVolume_.isReady();
    ready &= inportForegroundSeeds_.isReady();
    ready &= inportBackgroundSeeds_.isReady();
    return ready;
}

void OctreeWalker::serialize(Serializer& ser) const {
    AsyncComputeProcessor<OctreeWalkerInput, OctreeWalkerOutput>::serialize(ser);
    ser.serialize("prevResultPath", prevResultPath_);
}

void OctreeWalker::deserialize(Deserializer& s) {
    AsyncComputeProcessor<OctreeWalkerInput, OctreeWalkerOutput>::deserialize(s);
    try {
        s.deserialize("prevResultPath", prevResultPath_);
    } catch (SerializationException& e) {
        s.removeLastError();
    }

    std::unique_ptr<VolumeOctree> previousOctree(nullptr);


    std::string previousResultFile = tgt::FileSystem::cleanupPath(prevResultPath_ + "/" + PREV_RESULT_FILE_NAME);

    PointSegmentListGeometryVec3 foregroundSeeds;
    PointSegmentListGeometryVec3 backgroundSeeds;
    if(!prevResultPath_.empty()) {
        if(tgt::FileSystem::fileExists(prevResultPath_ + "/" + PREV_RESULT_VALID_FILE_NAME)) {
            if(tgt::FileSystem::fileExists(previousResultFile)) {
                XmlDeserializer d(prevResultPath_);

                previousOctree.reset(new VolumeOctree());
                try {
                    std::ifstream fs(previousResultFile);
                    d.read(fs);
                    d.deserialize(PREV_RESULT_OCTREE_KEY, *previousOctree);

                    resultPath_.set(prevResultPath_);

                    try {
                        d.deserialize(PREV_RESULT_FOREGROUND_KEY, foregroundSeeds);
                        d.deserialize(PREV_RESULT_BACKGROUND_KEY, backgroundSeeds);
                    } catch (SerializationException& e) {
                        d.removeLastError();
                        LWARNING("No foreground/background seed present from previous result." << e.what());
                    }
                } catch (std::exception& e) {
                    LERROR("Failed to deserialize previous solution: " << e.what());
                    previousOctree.reset(nullptr);
                } catch (SerializationException& e) {
                    d.removeLastError();
                    LERROR("Failed to deserialize previous solution: " << e.what());
                    previousOctree.reset(nullptr);
                }
            }
        } else {
            // This may happen if Voreen crashed or was killed during execution.
            LERROR("Skipping deserialization of invalid previous solution");
        }
    }


    previousResult_ = boost::none;
    if(previousOctree) {
        tgt::vec3 spacing(1.0f);
        tgt::vec3 offset(1.0f);

        OctreeBrickPoolManagerMmap* obpmm = dynamic_cast<OctreeBrickPoolManagerMmap*>(previousOctree->getBrickPoolManager());
        if(obpmm) {
            brickPoolManager_.reset(obpmm);
            previousResult_ = OctreeWalkerPreviousResult(
                *previousOctree,
                std::unique_ptr<VolumeBase>(new Volume(previousOctree.release(), spacing, offset)),
                std::move(foregroundSeeds),
                std::move(backgroundSeeds)
            );
        }
    }
}

template<typename MeanFunc, typename VarFunc>
static void compute_variance(tgt::svec3 basePos, tgt::svec3 halfBrickSize, VolumeAtomic<uint16_t>& outBrick, MeanFunc meanfunc, VarFunc varfunc) {
    VRN_FOR_EACH_VOXEL(i, tgt::svec3(0), halfBrickSize) {
        double varsum = 0;
        double meansum = 0;
        double meansqsum = 0;
        VRN_FOR_EACH_VOXEL(c, tgt::svec3(0), tgt::svec3(2)) {
            tgt::svec3 p = i*static_cast<size_t>(2)+c;

            double mean = brickToNorm(meanfunc(p));
            double var = unpackNormalizedFloat(varfunc(p));

            varsum += var;
            meansum += mean;
            meansqsum += mean*mean;
        }
        double mean = meansum / 8;
        double meansq = mean*mean;

        double squaredterm = (varsum + meansqsum) / 8;

        double var = squaredterm - meansq;

        uint16_t val = packNormalizedFloat(var);
        outBrick.voxel(basePos+i) = val;
    }
}


static VolumeOctreeNodeGeneric<1>* buildVarianceTreeRecursively(const VolumeOctreeNodeGeneric<1>& inputNode, const OctreeBrickPoolManagerBase& inputBrickPoolManager, OctreeBrickPoolManagerMmap& varianceBrickPoolManager, int level, const tgt::svec3& brickDataSize) {
    if(level == 0 || inputNode.isLeaf() || inputNode.isHomogeneous()) {
        return new VolumeOctreeNodeGeneric<1>(OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS, true, 0);
    } else {
        auto outputNode = new VolumeOctreeNodeGeneric<1>(varianceBrickPoolManager.allocateBrick(), true, 0);
        auto outputBrickAddr = outputNode->getBrickAddress();
        tgtAssert(outputBrickAddr != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS, "No brick");
        BrickPoolBrick outputBrick(outputBrickAddr, brickDataSize, varianceBrickPoolManager);

        const tgt::svec3 halfBrickSize = brickDataSize/static_cast<size_t>(2);

        for(int i=0; i<8; ++i) {
            auto inputChildBase = inputNode.children_[i];
            if(!inputChildBase) {
                outputNode->children_[i] = nullptr;
                continue;
            }
            if(!inputChildBase->inVolume()) {
                outputNode->children_[i] = new VolumeOctreeNodeGeneric<1>(OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS, false);
                continue;
            }
            auto inputChild = dynamic_cast<const VolumeOctreeNodeGeneric<1>*>(inputChildBase);
            auto inputChildBrickAddr = inputChild->getBrickAddress();
            tgtAssert(inputChild, "Too many channels");

            VolumeOctreeNodeGeneric<1>* varianceChild = buildVarianceTreeRecursively(*inputChild, inputBrickPoolManager, varianceBrickPoolManager, level-1, brickDataSize);
            auto varianceChildBrickAddr = varianceChild->getBrickAddress();
            outputNode->children_[i] = varianceChild;

            tgt::svec3 basePos = tgt::svec3::zero;
            basePos.x += (i & 1) == 0 ? 0 : halfBrickSize.x;
            basePos.y += (i & 2) == 0 ? 0 : halfBrickSize.y;
            basePos.z += (i & 4) == 0 ? 0 : halfBrickSize.z;
            if(varianceChildBrickAddr != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS) {
                BrickPoolBrickConst varianceChildBrick(varianceChildBrickAddr, brickDataSize, varianceBrickPoolManager);

                if(inputChildBrickAddr != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS) {
                    BrickPoolBrickConst inputChildBrick(inputChildBrickAddr, brickDataSize, inputBrickPoolManager);

                    compute_variance(basePos, halfBrickSize, outputBrick.data(), [&] (tgt::svec3 p) { return inputChildBrick.data().voxel(p); }, [&] (tgt::svec3 p) { return varianceChildBrick.data().voxel(p); });
                } else {
                    compute_variance(basePos, halfBrickSize, outputBrick.data(), [&] (tgt::svec3 p) { return inputChild->getAvgValue(); }, [&] (tgt::svec3 p) { return varianceChildBrick.data().voxel(p); });
                }
            } else {
                if(inputChildBrickAddr != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS) {
                    BrickPoolBrickConst inputChildBrick(inputChildBrickAddr, brickDataSize, inputBrickPoolManager);
                    compute_variance(basePos, halfBrickSize, outputBrick.data(), [&] (tgt::svec3 p) { return inputChildBrick.data().voxel(p); }, [&] (tgt::svec3 p) { return varianceChild->getAvgValue(); });
                } else {
                    compute_variance(basePos, halfBrickSize, outputBrick.data(), [&] (tgt::svec3 p) { return inputChild->getAvgValue(); }, [&] (tgt::svec3 p) { return varianceChild->getAvgValue(); });
                }
            }
        }

        return outputNode;
    }
}

VarianceTree VarianceTree::none() {
    return VarianceTree {
        nullptr,
        LocatedVolumeOctreeNode(nullptr, 0, tgt::svec3::zero, tgt::svec3::zero),
    };
}

bool VarianceTree::isNone() const {
    bool ret = brickPoolManager_.get() == nullptr;
    tgtAssert((root_.node_ == nullptr) == ret, "invalid state");
    return ret;
}

VarianceTree::VarianceTree(std::unique_ptr<OctreeBrickPoolManagerMmap>&& brickPoolManager, LocatedVolumeOctreeNode root)
    : brickPoolManager_(std::move(brickPoolManager))
    , root_(root)
{
}

VarianceTree::VarianceTree(VarianceTree&& other)
    : brickPoolManager_(std::move(other.brickPoolManager_))
    , root_(other.root_)
{
    other.root_.node_ = nullptr;
}
VarianceTree& VarianceTree::operator=(VarianceTree&& other) {
    if(this != &other) {
        // Destruct the current object, but keep the memory.
        this->~VarianceTree();
        // Call the move constructor on the memory region of the current object.
        new(this) VarianceTree(std::move(other));
    }

    return *this;
}
VarianceTree::~VarianceTree() {
    freeNodes(root_.node_);
    if(brickPoolManager_) {
        brickPoolManager_->deinitialize();
        tgt::FileSystem::deleteDirectoryRecursive(brickPoolManager_->getBrickPoolPath());
    }
}

static VarianceTree computeVarianceTree(const OctreeWalkerInput& input) {
    const OctreeBrickPoolManagerBase& inputPoolManager = *input.octree_.getBrickPoolManager();
    const tgt::svec3 brickDataSize = input.octree_.getBrickDim();
    const tgt::svec3 dim = input.octree_.getDimensions();

    auto inputRoot = dynamic_cast<const VolumeOctreeNodeGeneric<1>*>(input.octree_.getRootNode());
    tgtAssert(inputRoot, "Too many channels");

    std::string brickpooltmppath = VoreenApplication::app()->getUniqueTmpFilePath();
    tgt::FileSystem::createDirectoryRecursive(brickpooltmppath);
    auto varianceBrickPoolManager = tgt::make_unique<OctreeBrickPoolManagerMmap>(brickpooltmppath, BRICK_BUFFER_FILE_PREFIX);
    varianceBrickPoolManager->initialize(inputPoolManager.getBrickMemorySizeInByte());

    int rootLevel = input.octree_.getNumLevels()-1;
    VolumeOctreeNodeGeneric<1>* varianceRoot = buildVarianceTreeRecursively(*inputRoot, inputPoolManager, *varianceBrickPoolManager, rootLevel, brickDataSize);

    return VarianceTree {
        std::move(varianceBrickPoolManager),
        LocatedVolumeOctreeNode(varianceRoot, rootLevel, tgt::svec3::zero, dim),
    };
}

OctreeWalker::ComputeInput OctreeWalker::prepareComputeInput() {
    tgtAssert(inportVolume_.hasData(), "no input volume");

    // clear previous results and update property ranges, if input volume has changed
    if (inportVolume_.hasChanged()) {
        outportProbabilities_.setData(0);
        varianceTree_ = VarianceTree::none();
    }
    auto vol = inportVolume_.getThreadSafeData();

    // clear previous results and update property ranges, if input volume has changed
    if (!vol->hasRepresentation<VolumeOctree>()) {
        throw InvalidInputException("No octree Representation", InvalidInputException::S_ERROR);
    }
    const VolumeOctree* octreePtr = vol->getRepresentation<VolumeOctree>();
    tgtAssert(octreePtr, "No octree");
    const VolumeOctree& octree = *octreePtr;

    if (octree.getNumChannels() != 1) {
        throw InvalidInputException("Only single channel volumes are supported.", InvalidInputException::S_ERROR);
    }


    // select BLAS implementation and preconditioner
    const VoreenBlas* voreenBlas = getVoreenBlasFromProperties();
    VoreenBlas::ConjGradPreconditioner precond = VoreenBlas::NoPreconditioner;
    if (preconditioner_.isSelected("jacobi"))
        precond = VoreenBlas::Jacobi;

    float errorThresh = 1.f / pow(10.f, static_cast<float>(errorThreshold_.get()));
    int maxIterations = maxIterations_.get();

    prevResultPath_ = resultPath_.get();
    if (prevResultPath_.empty() ||  !tgt::FileSystem::dirExists(tgt::FileSystem::parentDir(prevResultPath_))) {
        throw InvalidInputException("Parent directory of previous result path " + prevResultPath_ + " does not exist.", InvalidInputException::S_ERROR);
    }
    std::string brickPoolPath = tgt::FileSystem::cleanupPath(prevResultPath_ + "/" + BRICK_BUFFER_SUBDIR);
    if (!tgt::FileSystem::dirExists(brickPoolPath)) {
        tgt::FileSystem::createDirectoryRecursive(brickPoolPath);
    }

    const tgt::svec3 volumeDim = octree.getDimensions();
    const tgt::svec3 brickDim = octree.getBrickDim();
    const size_t brickSize = tgt::hmul(brickDim);
    const size_t numChannels = 1;
    const size_t maxLevel = octree.getNumLevels()-1;
    size_t brickSizeInBytes = brickSize * sizeof(uint16_t);

    if(previousResult_) {
        auto& oldTree = previousResult_->octree();
        if (incrementalSimilarityThreshold_.get() == 0.0f) {
            LINFO("Not reusing previous solution due to incremental similarity threshold of 0");
            clearPreviousResults();
        } else if(oldTree.getDimensions() != volumeDim
               || oldTree.getBrickDim() != brickDim
               || oldTree.getNumLevels() != octree.getNumLevels()
               || brickPoolManager_->getBrickMemorySizeInByte() != brickSizeInBytes) {

            LINFO("Not reusing previous solution for incompatible input volume");
            clearPreviousResults();
        }
    }

    if(previousResult_) {
        LINFO("Reusing previous solution for compatible input volume");
    }

    // Create a new brickpool if we need a new one
    if(!brickPoolManager_ || brickPoolManager_->getBrickMemorySizeInByte() != brickSizeInBytes) {
        brickPoolManager_.reset(new OctreeBrickPoolManagerMmap(brickPoolPath, BRICK_BUFFER_FILE_PREFIX));
        brickPoolManager_->initialize(brickSizeInBytes);
    }

    return ComputeInput {
        previousResult_ ? &previousResult_->octree() : nullptr,
        previousResult_ ? &previousResult_->foregroundSeeds_ : nullptr,
        previousResult_ ? &previousResult_->backgroundSeeds_ : nullptr,
        brickPoolManager_,
        *vol,
        *octreePtr,
        inportForegroundSeeds_.getThreadSafeAllData(),
        inportBackgroundSeeds_.getThreadSafeAllData(),
        minEdgeWeight_.get(),
        parameterEstimationNeighborhoodExtent_.get(),
        voreenBlas,
        precond,
        errorThresh,
        maxIterations,
        homogeneityThreshold_.get(),
        binaryPruningDelta_.get(),
        incrementalSimilarityThreshold_.get(),
        noiseModel_.getValue(),
        varianceTree_,
    };
}
struct BrickNeighborhood {
    BrickNeighborhood() = delete;
    BrickNeighborhood(const BrickNeighborhood& other) = delete;
    BrickNeighborhood& operator=(const BrickNeighborhood& other) = delete;
    BrickNeighborhood(BrickNeighborhood&& other) = default;
    BrickNeighborhood& operator=(BrickNeighborhood&& other) = default;

    tgt::mat4 centerBrickToNeighborhood() const {
        return tgt::mat4::createTranslation(tgt::vec3(centerBrickLlf_));
    }
    tgt::mat4 neighborhoodToCenterBrick() const {
        return tgt::mat4::createTranslation(-tgt::vec3(centerBrickLlf_));
    }
    tgt::mat4 voxelToNeighborhood() const {
        return centerBrickToNeighborhood() * voxelToCenterBrick_;
    }
    bool isEmpty() const {
        return data_.getDimensions() == tgt::svec3(0);
    }

    VolumeAtomic<float> data_;
    tgt::svec3 centerBrickLlf_; // In coordinate system ...
    tgt::svec3 centerBrickUrb_; // ... of neighborhood buffer
    tgt::svec3 dimensions_;
    tgt::mat4 voxelToCenterBrick_;
    float min_;
    float max_;
    float avg_;

    template<typename ToFloat>
    static BrickNeighborhood fromNode(const VolumeOctreeNodeLocation& current, size_t sampleLevel, const LocatedVolumeOctreeNodeConst& root, const tgt::svec3& brickBaseSize, const OctreeBrickPoolManagerBase& brickPoolManager, ToFloat toFloat = brickToNorm) {
        const tgt::svec3 volumeDim = root.location().voxelDimensions();

        const tgt::mat4 brickToVoxel = current.brickToVoxel();
        const tgt::mat4 voxelToBrick = current.voxelToBrick();

        //const tgt::ivec3 neighborhoodSize = brickBaseSize;
        //const tgt::ivec3 neighborhoodSize = tgt::ivec3(2);
        const tgt::vec3 neighborhoodSize = brickBaseSize/size_t(8);
        const tgt::vec3 neighborhoodSizeGlobal = brickToVoxel.getRotationalPart().transform(neighborhoodSize);

        const tgt::vec3 voxelLlf = tgt::max(tgt::vec3(0.0),         tgt::vec3(current.llf_) - neighborhoodSizeGlobal);
        const tgt::vec3 voxelUrb = tgt::min(tgt::vec3(volumeDim), tgt::vec3(current.urb_) + neighborhoodSizeGlobal);

        const tgt::vec3 regionLlf = voxelToBrick.transform(voxelLlf);
        const tgt::vec3 regionUrb = voxelToBrick.transform(voxelUrb);

        // Note: Here in and in other cases: We use `ceil` for the URB-Corner
        // of brick bounding boxes because we must include "overhanging" bricks
        // from lower levels. A minimal example with a volume of dimension 3
        // and a brick size of 2:
        //
        // |. .|. _|
        // | .   X |
        //
        // Without `ceil`, the root level voxel marked with X would not be
        // processed.
        // Afaik (this may be incorrect in general, but is at least true for
        // the `VolumeOctreeLevelExtractor`), Voreen would not include X when
        // rendering the octree, but we most still process it here, in order to
        // avoid ignoring labels at the (URB) border of the volume.
        const tgt::svec3 regionDim = tgt::ceil(regionUrb - regionLlf);
        const tgt::svec3 brickDim = tgt::ceil(voxelToBrick.getRotationalPart().transform(current.urb_ - current.llf_));
        const tgt::svec3 neighborhoodOffset = tgt::round(voxelToBrick.getRotationalPart().transform(tgt::vec3(current.llf_) - voxelLlf));
        tgtAssert(tgt::hand(tgt::lessThanEqual(brickDim+neighborhoodOffset, regionDim)), "Center brick region with offset should never be larger than the entire region.");

        VolumeAtomic<float> output(regionDim);

        float min = std::numeric_limits<float>::infinity();
        float max = -std::numeric_limits<float>::infinity();
        float sum = 0.0f;

        // In output buffer space
        const tgt::svec3 regionBegin(0);
        const tgt::svec3 regionEnd = regionDim;
        const tgt::svec3 centerBegin = neighborhoodOffset;
        const tgt::svec3 centerEnd = neighborhoodOffset+brickDim;

        VRN_FOR_EACH_VOXEL(blockIndex, tgt::svec3(0), tgt::svec3(3)) {
            tgt::svec3 blockBegin;
            tgt::svec3 blockEnd;
            for(int dim=0; dim<3; ++dim) {
                switch(blockIndex[dim]) {
                    case 0: {
                        blockBegin[dim] = regionBegin[dim];
                        blockEnd[dim] = centerBegin[dim];
                        break;
                    }
                    case 1: {
                        blockBegin[dim] = centerBegin[dim];
                        blockEnd[dim] = centerEnd[dim];
                        break;
                    }
                    case 2: {
                        blockBegin[dim] = centerEnd[dim];
                        blockEnd[dim] = regionEnd[dim];
                        break;
                    }
                }
            }
            tgt::svec3 blockDimensions = blockEnd - blockBegin;
            if(tgt::hor(tgt::equal(blockDimensions, tgt::svec3(0)))) {
                continue;
            }

            tgt::mat4 bufferToBrick = tgt::mat4::createTranslation(-tgt::vec3(centerBegin));

            tgt::mat4 bufferToVoxel = brickToVoxel * bufferToBrick;
            tgt::svec3 samplePoint = tgt::fastround(bufferToVoxel.transform(blockBegin));
            auto node = root.findChildNode(samplePoint, brickBaseSize, sampleLevel, true);
            if(node.node().hasBrick()) {
                tgt::mat4 bufferToSampleBrick = node.location().voxelToBrick() * bufferToVoxel;
                tgt::ivec3 maxpos = node.location().brickDimensions() - tgt::svec3(1);


                // The following is equivalent to this, but faster.
                //
                //VRN_FOR_EACH_VOXEL(point, blockBegin, blockEnd) {
                //    tgt::svec3 samplePos = tgt::round(bufferToSampleBrick.transform(point));
                //    samplePos = tgt::clamp(samplePos, tgt::svec3(0), maxpos);
                //    float val = brick.getVoxelNormalized(samplePos);

                //    output.setVoxelNormalized(val, point);
                //    min = std::min(val, min);
                //    max = std::max(val, max);
                //    sum += val;
                //}

                BrickPoolBrickConst brick(node.node().getBrickAddress(), brickBaseSize, brickPoolManager);
                const uint16_t* brickData = brick.data().voxel();
                tgt::ivec3 brickDim = brick.data().getDimensions();

                float* outputData = output.voxel();
                tgt::ivec3 outputDim = output.getDimensions();

                // Note: We know (and thus assume and use) that bufferToSampleBrick is only a scale with translation
                for(int z=blockBegin.z; z!=blockEnd.z; ++z) {
                    int sampleZ = tgt::fastround(bufferToSampleBrick.elemRowCol[2][2] * z + bufferToSampleBrick.elemRowCol[2][3]);
                    sampleZ = tgt::clamp(sampleZ, 0, maxpos.z);
                    int brickIndexBaseY = sampleZ * brickDim.y;

                    int outputBaseY = z * outputDim.y;

                    for(int y=blockBegin.y; y!=blockEnd.y; ++y) {
                        int sampleY = tgt::fastround(bufferToSampleBrick.elemRowCol[1][1] * y + bufferToSampleBrick.elemRowCol[1][3]);
                        sampleY = tgt::clamp(sampleY, 0, maxpos.y);
                        int brickIndexBaseX = brickDim.x * (sampleY + brickIndexBaseY);

                        int outputBaseX = outputDim.x * (y + outputBaseY);
                        int outputIndex = outputBaseX + blockBegin.x;

                        for(int x=blockBegin.x; x!=blockEnd.x; ++x) {
                            int sampleX = tgt::fastround(bufferToSampleBrick.elemRowCol[0][0] * x + bufferToSampleBrick.elemRowCol[0][3]);
                            sampleX = tgt::clamp(sampleX, 0, maxpos.x);
                            int brickIndex = brickIndexBaseX + sampleX;


                            uint16_t valuint = brickData[brickIndex];
                            float val = toFloat(valuint);

                            outputData[outputIndex] = val;

                            min = std::min(val, min);
                            max = std::max(val, max);
                            sum += val;

                            ++outputIndex;
                        }
                    }
                }
            } else {
                float val = static_cast<float>(node.node().getAvgValue())/0xffff;
                min = std::min(val, min);
                max = std::max(val, max);
                sum += val * tgt::hmul(blockEnd - blockBegin);
                VRN_FOR_EACH_VOXEL(point, blockBegin, blockEnd) {
                    output.setVoxelNormalized(val, point);
                }
            }
        }
        float avg = sum / output.getNumVoxels();
        return BrickNeighborhood {
            std::move(output),
            centerBegin,
            centerEnd,
            regionDim,
            voxelToBrick,
            min,
            max,
            avg
        };
    }
};

class RandomWalkerSeedsBrick : public RandomWalkerSeeds {
    static const float PREVIOUS_CONFLICT;  // Conflicting labels in this voxel, already present in previous iteration
    static const float CURRENT_CONFLICT;   // Conflicting labels in this voxel with new labels taking part
    static const float UNLABELED;          // No label in this voxel
    static const float FOREGROUND;
    static const float BACKGROUND;
    static const float FOREGROUND_NEW;
    static const float BACKGROUND_NEW;
public:
    RandomWalkerSeedsBrick(tgt::svec3 bufferDimensions, tgt::mat4 voxelToSeeds, const PointSegmentListGeometryVec3& foregroundSeedList, const PointSegmentListGeometryVec3& backgroundSeedList, const PointSegmentListGeometryVec3* fgslOld, const PointSegmentListGeometryVec3* bgslOld)
        : seedBuffer_(bufferDimensions)
        , conflicts_(false)
        , newConflicts_(false)
        , minSeed_(1.0)
        , maxSeed_(0.0)
    {
        numSeeds_ = 0;
        seedBuffer_.fill(UNLABELED);

        tgt::IntBounds brickBounds(tgt::ivec3::zero, tgt::ivec3(bufferDimensions)-tgt::ivec3::one);

        auto collectLabelsFromGeometry = [&] (const PointSegmentListGeometryVec3& seedList, const PointSegmentListGeometryVec3* oldList, float label) {
            for (int m=0; m<seedList.getNumSegments(); m++) {
                const std::vector<tgt::vec3>& points = seedList.getData()[m];
                const std::vector<tgt::vec3>* old = nullptr;
                if (oldList && m < oldList->getData().size()) {
                    old = &oldList->getData()[m];
                }
                if (points.empty())
                    continue;
                for (size_t i=0; i<points.size()-1; i++) {
                    bool new_seeds = true;
                    if(old && i+1 < old->size()) {
                        new_seeds = points[i] != (*old)[i] || points[i+1] != (*old)[i+1];
                    }

                    tgt::vec3 left = voxelToSeeds*points[i];
                    tgt::vec3 right = voxelToSeeds*points[i+1];
                    tgt::vec3 dir = tgt::normalize(right - left);

                    for (float t=0.f; t<tgt::length(right-left); t += 1.f) {
                        tgt::ivec3 point = tgt::iround(left + t*dir);
                        if(!brickBounds.containsPoint(point)) {
                            continue;
                        }
                        float& seedVal = seedBuffer_.voxel(point);
                        minSeed_ = std::min(label, minSeed_);
                        maxSeed_ = std::max(label, maxSeed_);

                        float lbl = label;
                        if(new_seeds) {
                            if(lbl == FOREGROUND) {
                                lbl = FOREGROUND_NEW;
                            } else {
                                lbl = BACKGROUND_NEW;
                            }
                        }
                        if (seedVal == UNLABELED) {
                            seedVal = lbl;
                            ++numSeeds_;
                        } else if(seedVal == CURRENT_CONFLICT) {
                            // ignore, already processed
                        } else if(seedVal == PREVIOUS_CONFLICT) {
                            if(new_seeds) {
                                seedVal = CURRENT_CONFLICT;
                                newConflicts_ = true;
                            } else {
                                // Ignore, still a previous conflict
                            }
                        } else { // Seedval is either foreground or background
                            float seedClass = tgt::clamp(seedVal, 0.0f, 1.0f);
                            bool alreadyNew = seedClass != seedVal;
                            if(label == seedClass) {
                                // Same label as before, but possibly new
                                if(new_seeds) {
                                    seedVal = lbl;
                                } else {
                                    // Ignore, seed value is either equal or NEW_* already
                                }
                            } else {
                                // Conflict, either previous or caused by new labels
                                conflicts_ = true;
                                if(new_seeds || alreadyNew) {
                                    seedVal = CURRENT_CONFLICT;
                                    newConflicts_ = true;
                                } else {
                                    seedVal = PREVIOUS_CONFLICT;
                                }
                                // Either way the current voxel will not act as
                                // a seed for the current brick:
                                --numSeeds_;
                            }
                        }
                    }
                }
                if(old) {
                    for(size_t i = points.size(); i < old->size(); ++i) {
                        // If the previous label array is longer than this one
                        // and one of these old points is in this brick, we
                        // definitely want to recompute this brick as its label
                        // situation has changed.
                        auto& p = (*old)[i];
                        if(brickBounds.containsPoint(voxelToSeeds*p)) {
                            conflicts_ = true;
                            newConflicts_ = true;
                            break;
                        }
                    }
                }
            }
            // Similarly, if there were previous label point lists in this
            // brick that are not present anymore, we definitely have to
            // recompute this brick:
            if(oldList) {
                size_t oldListSize = oldList->getData().size();
                for(size_t i=seedList.getNumSegments(); i<oldListSize; ++i) {
                    for (auto& p : oldList->getData()[i]) {
                        if(brickBounds.containsPoint(voxelToSeeds*p)) {
                            conflicts_ = true;
                            newConflicts_ = true;
                            break;
                        }
                    }
                }
            }
        };
        collectLabelsFromGeometry(foregroundSeedList, fgslOld, FOREGROUND);
        collectLabelsFromGeometry(backgroundSeedList, bgslOld, BACKGROUND);
    }
    RandomWalkerSeedsBrick(RandomWalkerSeedsBrick&& other) = default;

    void addNeighborhoodBorderSeeds(const BrickNeighborhood& neighborhood) {
        tgtAssert(neighborhood.data_.getDimensions() == neighborhood.dimensions_, "Invalid buffer dimensions");

        auto collectLabelsFromNeighbor = [&] (size_t dim, size_t sliceIndex) {
            tgt::svec3 begin(0);
            tgt::svec3 end(neighborhood.dimensions_);

            begin[dim] = sliceIndex;
            end[dim] = sliceIndex+1;

            // Do not collect parent level border labels at the border of the
            // volume. There is no additional information in this case.
            //
            // We are at the border iff their is no neighborhood around the centerbrick,
            // i.e., if the center brick llf/urb is at 0/neighborhood dimensions
            if(begin[dim] == neighborhood.centerBrickLlf_[dim] || end[dim] == neighborhood.centerBrickUrb_[dim]) {
                return;
            }

            VRN_FOR_EACH_VOXEL(seed, begin, end) {
                float& seedVal = seedBuffer_.voxel(seed);
                if (seedVal == UNLABELED || seedVal == CURRENT_CONFLICT || seedVal == PREVIOUS_CONFLICT) {
                    // Only exising, user-defined labels are not overwritten by
                    // a border seed.
                    // Unlabeled values are fine either way, conflicting voxels
                    // have been counted as conflicts before, but can now
                    // provide a sensible value from the parent node/brick.
                    float val = neighborhood.data_.voxel(seed);
                    seedVal = val;
                    ++numSeeds_;
                    minSeed_ = std::min(val, minSeed_);
                    maxSeed_ = std::max(val, maxSeed_);
                }
            }
        };

        collectLabelsFromNeighbor(0, 0);
        collectLabelsFromNeighbor(0, neighborhood.dimensions_.x-1);
        collectLabelsFromNeighbor(1, 0);
        collectLabelsFromNeighbor(1, neighborhood.dimensions_.y-1);
        collectLabelsFromNeighbor(2, 0);
        collectLabelsFromNeighbor(2, neighborhood.dimensions_.z-1);
    }
    virtual ~RandomWalkerSeedsBrick() {}
    virtual void initialize() {};

    virtual bool isSeedPoint(size_t index) const {
        float val = seedBuffer_.voxel(index);
        return val != UNLABELED && val != CURRENT_CONFLICT && val != PREVIOUS_CONFLICT;
    }
    virtual bool isSeedPoint(const tgt::ivec3& voxel) const {
        float val = seedBuffer_.voxel(voxel);
        return val != UNLABELED && val != CURRENT_CONFLICT && val != PREVIOUS_CONFLICT;
    }
    virtual float getSeedValue(size_t index) const {
        tgtAssert(isSeedPoint(index), "Getting seed value from non-seed");
        return tgt::clamp(seedBuffer_.voxel(index), 0.0f, 1.0f);
    }
    virtual float getSeedValue(const tgt::ivec3& voxel) const {
        tgtAssert(isSeedPoint(voxel), "Getting seed value from non-seed");
        return tgt::clamp(seedBuffer_.voxel(voxel), 0.0f, 1.0f);
    }

    std::vector<size_t> generateVolumeToRowsTable() {
        size_t numVoxels = tgt::hmul(seedBuffer_.getDimensions());
        std::vector<size_t> volIndexToRow(numVoxels, -1);

        // compute volIndexToRow values
        size_t curRow = 0;
        for (size_t i=0; i<numVoxels; i++) {
            if (isSeedPoint(i)) {
                volIndexToRow[i] = -1;
            } else {
                volIndexToRow[i] = curRow;
                curRow++;
            }
        }
        return volIndexToRow;
    }

    tgt::svec3 bufferDimensions() const {
        return seedBuffer_.getDimensions();
    }
    bool hasConflicts() const {
        return conflicts_;
    }
    bool hasNewConflicts() const {
        return newConflicts_;
    }
    float minSeed() {
        return minSeed_;
    }
    float maxSeed() {
        return maxSeed_;
    }

    VolumeAtomic<float> seedBuffer_;
    bool conflicts_;
    bool newConflicts_;
    float minSeed_;
    float maxSeed_;
};

const float RandomWalkerSeedsBrick::FOREGROUND        =  1.0f;
const float RandomWalkerSeedsBrick::BACKGROUND        =  0.0f;
const float RandomWalkerSeedsBrick::FOREGROUND_NEW    =  2.0f;
const float RandomWalkerSeedsBrick::BACKGROUND_NEW    = -1.0f;
const float RandomWalkerSeedsBrick::PREVIOUS_CONFLICT = -4.0f;
const float RandomWalkerSeedsBrick::CURRENT_CONFLICT  = -3.0f;
const float RandomWalkerSeedsBrick::UNLABELED         = -2.0f;

struct RandomWalkerVoxelAccessorBrick final : public RandomWalkerVoxelAccessor {
    RandomWalkerVoxelAccessorBrick(const VolumeAtomic<float>& brick)
        : brick_(brick)
    {
    }
    inline float voxel(const tgt::svec3& pos) {
        return brick_.voxel(pos);
    }
private:
    const VolumeAtomic<float>& brick_;
};

template<RWNoiseModel NoiseModel>
static void processVoxelWeights(const RandomWalkerSeedsBrick& seeds, EllpackMatrix<float>& mat, float* vec, size_t* volumeIndexToRowTable, RWNoiseModelWeights<NoiseModel>& model, const tgt::svec3& volDim, float minWeight, tgt::vec3 spacing) {

    float minSpacing = tgt::min(spacing);
    auto edgeWeight = [&] (tgt::vec3 voxel, tgt::vec3 neighbor) {
        float weight = model.getEdgeWeight(voxel, neighbor);
        weight = std::max(weight, minWeight);

        return weight;
    };
    tgtAssert(volumeIndexToRowTable, "no volumeIndexToRowTable passed");
    tgtAssert(mat.isInitialized(), "matrix not initialized");

    VRN_FOR_EACH_VOXEL(voxel, tgt::ivec3(0), tgt::ivec3(volDim)) {

        size_t index = volumeCoordsToIndex(voxel, volDim);

        //float curIntensity = voxelFun.voxel(voxel);

        float weightSum = 0;

        bool currentIsSeedpoint = seeds.isSeedPoint(index);

        for(int dim=0; dim<3; ++dim) {
            if(voxel[dim] > 0) {
                tgt::ivec3 neighbor = voxel;
                neighbor[dim] -= 1;

                size_t neighborIndex = volumeCoordsToIndex(neighbor, volDim);
                //float neighborIntensity = voxelFun.voxel(neighbor);

                float weight = edgeWeight(voxel, neighbor);

                if(seeds.isSeedPoint(neighbor)) {
                    if(!currentIsSeedpoint) {
                        size_t curRow = volumeIndexToRowTable[index];
                        vec[curRow] += weight * seeds.getSeedValue(neighbor);
                    }
                } else {
                    size_t nRow = volumeIndexToRowTable[neighborIndex];
                    if(!currentIsSeedpoint) {
                        size_t curRow = volumeIndexToRowTable[index];
                        //tgtAssert(nRow >= 0 && nRow < numUnseeded_, "Invalid row");
                        tgtAssert(mat.getIndex(curRow, nRow) == -1, "foo");
                        tgtAssert(mat.getIndex(nRow, curRow) == -1, "foo");
                        mat.getWritableValue(curRow, nRow) = -weight;
                        mat.getWritableValue(nRow, curRow) = -weight;
                    } else {
                        vec[nRow] += weight * seeds.getSeedValue(voxel);
                    }

                    tgtAssert(mat.getIndex(nRow, nRow) != -1, "foo");
                    // Update weight sum of neighbor with smaller index.
                    mat.getWritableValue(nRow, nRow) += weight;
                }

                weightSum += weight;
            }
        }

        if(!currentIsSeedpoint) {
            // This is the first time writing to mat at this location, so overwriting is fine.
            size_t curRow = volumeIndexToRowTable[index];
            tgtAssert(mat.getIndex(curRow, curRow) == -1, "foo");
            mat.getWritableValue(curRow, curRow) = weightSum;
        }
    }
}

template<RWNoiseModel NoiseModel>
static uint64_t processOctreeBrick(RWNoiseModelParameters<NoiseModel> parameters, OctreeWalkerInput& input, VolumeOctreeNodeLocation& outputNodeGeometry, Histogram1D& histogram, uint16_t& min, uint16_t& max, uint16_t& avg, bool& hasNewSeedConflicts, OctreeBrickPoolManagerBase& outputPoolManager, LocatedVolumeOctreeNode* outputRoot, const LocatedVolumeOctreeNodeConst& inputRoot, boost::optional<LocatedVolumeOctreeNode> prevRoot, const VarianceTree& varianceTree, PointSegmentListGeometryVec3& foregroundSeeds, PointSegmentListGeometryVec3& backgroundSeeds, std::mutex& clMutex, const ProfileDataCollector& ramProfiler, const ProfileDataCollector& vramProfiler) {
    auto canSkipChildren = [&] (float min, float max) {
        float parentValueRange = max-min;
        bool minMaxSkip = max < 0.5-input.binaryPruningDelta_ || min > 0.5+input.binaryPruningDelta_;
        return minMaxSkip || parentValueRange < input.homogeneityThreshold_;
    };

    const OctreeBrickPoolManagerBase& inputPoolManager = *input.octree_.getBrickPoolManager();
    const tgt::svec3 brickDataSize = input.octree_.getBrickDim();

    boost::optional<BrickNeighborhood> seedsNeighborhood = boost::none;

    bool stop = false;
    RandomWalkerSeedsBrick seeds = [&] () {
        if(outputRoot) {
            seedsNeighborhood = BrickNeighborhood::fromNode(outputNodeGeometry, outputNodeGeometry.level_+1, *outputRoot, brickDataSize, outputPoolManager, brickToNorm);
            tgt::svec3 seedBufferDimensions = seedsNeighborhood->data_.getDimensions();
            tgt::mat4 voxelToSeedTransform = seedsNeighborhood->voxelToNeighborhood();

            RandomWalkerSeedsBrick seeds(seedBufferDimensions, voxelToSeedTransform, foregroundSeeds, backgroundSeeds, input.previousForegroundSeeds_, input.previousBackgroundSeeds_);
            seeds.addNeighborhoodBorderSeeds(*seedsNeighborhood);

            if(!seeds.hasConflicts() && canSkipChildren(seeds.minSeed(), seeds.maxSeed())) {
                //LINFOC(OctreeWalker::loggerCat_, "skip block early");
                stop = true;
                avg = normToBrick(seedsNeighborhood->avg_);
                min = normToBrick(seeds.minSeed());
                max = normToBrick(seeds.maxSeed());
            }

            return seeds;
        } else {
            tgt::mat4 voxelToSeedTransform = outputNodeGeometry.voxelToBrick();
            // Note: For the reason to use `ceil` here see BrickNeighborhood::fromNode above.
            tgt::svec3 seedBufferDimensions = tgt::ceil(voxelToSeedTransform.getRotationalPart().transform(outputNodeGeometry.voxelDimensions()));
            RandomWalkerSeedsBrick seeds(seedBufferDimensions, voxelToSeedTransform, foregroundSeeds, backgroundSeeds, input.previousForegroundSeeds_, input.previousBackgroundSeeds_);
            return seeds;
        }
    }();
    ProfileAllocation seedsNeighborhoodAllocation(ramProfiler, seedsNeighborhood ? seedsNeighborhood->data_.getNumBytes() : 0);
    ProfileAllocation seedsAllocation(ramProfiler, seeds.seedBuffer_.getNumBytes());

    hasNewSeedConflicts = seeds.hasNewConflicts();
    if(stop) {
        return OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS;
    }

    tgt::svec3 walkerBlockDim = seeds.bufferDimensions();

    size_t numVoxels = tgt::hmul(walkerBlockDim);
    size_t numSeeds = seeds.getNumSeeds();
    size_t systemSize = numVoxels - numSeeds;

    BrickNeighborhood inputNeighborhood = BrickNeighborhood::fromNode(outputNodeGeometry, outputNodeGeometry.level_, inputRoot, brickDataSize, inputPoolManager, brickToNorm);
    ProfileAllocation neighborhoodAllocation(ramProfiler, inputNeighborhood.data_.getNumBytes());

    VolumeAtomic<float> varianceNeighborhood(tgt::svec3::zero, false);
    if(requiresVarianceTree(input.noiseModel_) && outputNodeGeometry.level_ > 0) {
        auto n = BrickNeighborhood::fromNode(outputNodeGeometry, outputNodeGeometry.level_, varianceTree.root_, brickDataSize, *varianceTree.brickPoolManager_, [&] (uint16_t v) { return unpackNormalizedFloat(v); });
        varianceNeighborhood = std::move(n.data_);
    }

    if(numSeeds == 0) {
        // No way to decide between foreground and background
        if(!seeds.hasConflicts()) {
            avg = 0xffff/2;
            min = avg;
            max = avg;
            for(int i=0; i<numVoxels; ++i) {
                histogram.addSample(0.5f);
            }
            return OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS;
        } else {
            // There are seeds, but they all overlap (=> conflicts). We need to try again in lower levels.
            uint64_t outputBrickAddr = outputPoolManager.allocateBrick();
            BrickPoolBrick outputBrick(outputBrickAddr, brickDataSize, outputPoolManager);
            ProfileAllocation outputBrickAllocation(ramProfiler, outputBrick.data().getNumBytes());

            const tgt::svec3 brickStart = inputNeighborhood.centerBrickLlf_;
            const tgt::svec3 brickEnd = inputNeighborhood.centerBrickUrb_;
            tgt::svec3 centerBrickSize = brickEnd - brickStart;
            if(seedsNeighborhood) {
                auto& parentProbs = *seedsNeighborhood;
                float sum = 0.5f;
                VRN_FOR_EACH_VOXEL(pos, brickStart, brickEnd) {
                    uint16_t val = parentProbs.data_.voxel(pos);

                    outputBrick.data().voxel(pos - brickStart) = val;

                    min = std::min(val, min);
                    max = std::max(val, max);
                    sum += val;

                    histogram.addSample(brickToNorm(val));
                }
                avg = sum/tgt::hmul(brickEnd);
            } else {
                avg = 0xffff/2;
                min = avg;
                max = avg;
                float avgf = brickToNorm(avg);
                VRN_FOR_EACH_VOXEL(pos, brickStart, brickEnd) {
                    outputBrick.data().voxel(pos) = avg;
                    histogram.addSample(avgf);
                }
            }
            return outputBrickAddr;
        }
    }

    auto volIndexToRow = seeds.generateVolumeToRowsTable();
    ProfileAllocation volIndexToRowAllocation(ramProfiler, volIndexToRow.size() * sizeof(size_t));

    std::vector<float> initialization(systemSize, 0.5f);
    ProfileAllocation initializationAllocation(ramProfiler, systemSize * sizeof(float));
    if(seedsNeighborhood) {
        auto& neighborhood = *seedsNeighborhood;
        VRN_FOR_EACH_VOXEL(pos, tgt::svec3(0), neighborhood.data_.getDimensions()) {
            size_t logicalIndex = volumeCoordsToIndex(pos, walkerBlockDim);
            if (!seeds.isSeedPoint(logicalIndex)) {
                initialization[volIndexToRow[logicalIndex]] = neighborhood.data_.voxel(pos);
            }
        }
    }

    ProfileAllocation solutionAllocation(ramProfiler, systemSize * sizeof(float));
    auto solution = tgt::make_unique<float[]>(systemSize);
    std::fill_n(solution.get(), systemSize, 0.5f);

    EllpackMatrix<float> mat(systemSize, systemSize, 7);
    mat.initializeBuffers();
    ProfileAllocation profmatvalues(ramProfiler, 7 * systemSize * sizeof(float));
    ProfileAllocation profmatrows(ramProfiler, 7 * systemSize * sizeof(size_t));

    float minWeight = 1.f / pow(10.f, static_cast<float>(input.minWeight_));

    auto rwm = input.volume_.getRealWorldMapping();

    auto model = parameters.prepare(inputNeighborhood.data_, rwm, outputNodeGeometry.level_, &varianceNeighborhood);
    ProfileAllocation noiseModelSize(ramProfiler, volIndexToRow.size() * sizeof(float) * 4); //Assuming worst case: ttest

    auto vec = std::vector<float>(systemSize, 0.0f);
    ProfileAllocation vecAllocation(ramProfiler, vec.size() * sizeof(float));

    tgt::vec3 spacing = input.volume_.getSpacing();

    processVoxelWeights<NoiseModel>(seeds, mat, vec.data(), volIndexToRow.data(), model, walkerBlockDim, minWeight, spacing);

    for(int i=0; i<10; ++i) {
        int iterations;
        {
            //std::lock_guard<std::mutex> guard(clMutex);
            ProfileAllocation systemAllocation = cgSystemFloatEllpack(vramProfiler, systemSize);
            iterations = input.blas_->sSpConjGradEll(mat, vec.data(), solution.get(), initialization.data(),
                input.precond_, input.errorThreshold_, input.maxIterations_);
        }
        if(iterations < input.maxIterations_) {
            break;
        }
        LWARNINGC(OctreeWalker::loggerCat_, "MAX ITER NOT SUFFICIENT: " << i);
    }

    const tgt::svec3 brickStart = inputNeighborhood.centerBrickLlf_;
    const tgt::svec3 brickEnd = inputNeighborhood.centerBrickUrb_;
    tgt::svec3 centerBrickSize = brickEnd - brickStart;

    uint64_t sum = 0;

    uint64_t outputBrickAddr = outputPoolManager.allocateBrick();
    {
        BrickPoolBrick outputBrick(outputBrickAddr, brickDataSize, outputPoolManager);
        ProfileAllocation outputBrickAllocation(ramProfiler, outputBrick.data().getNumBytes());

        VRN_FOR_EACH_VOXEL(pos, brickStart, brickEnd) {
            size_t logicalIndex = volumeCoordsToIndex(pos, walkerBlockDim);
            float valf;
            if (seeds.isSeedPoint(logicalIndex)) {
                valf = seeds.getSeedValue(logicalIndex);
            } else {
                valf = solution[volIndexToRow[logicalIndex]];
            }
            valf = tgt::clamp(valf, 0.0f, 1.0f);
            uint16_t val = valf*0xffff;

            outputBrick.data().voxel(pos - brickStart) = val;

            min = std::min(val, min);
            max = std::max(val, max);
            sum += val;

            histogram.addSample(valf);
        }
        avg = sum/tgt::hmul(centerBrickSize);
    }

    return outputBrickAddr;
}

const tgt::svec3 OCTREEWALKER_CHILD_POSITIONS[] = {
    tgt::svec3(0,0,0),
    tgt::svec3(1,0,0),
    tgt::svec3(1,1,0),
    tgt::svec3(0,1,0),
    tgt::svec3(0,1,1),
    tgt::svec3(1,1,1),
    tgt::svec3(1,0,1),
    tgt::svec3(0,0,1),
};

OctreeWalker::ComputeOutput OctreeWalker::compute(ComputeInput input, ProgressReporter& progressReporter) const {
    progressReporter.setProgress(0.0);

    auto start = clock::now();

    const tgt::svec3 volumeDim = input.octree_.getDimensions();
    const tgt::svec3 brickDim = input.octree_.getBrickDim();
    const size_t brickSize = tgt::hmul(brickDim);
    const size_t numChannels = 1;
    const size_t maxLevel = input.octree_.getNumLevels()-1;

    auto& outputBrickPoolManager = *input.brickPoolManager_;

    boost::optional<LocatedVolumeOctreeNode> prevRoot = [&] () -> boost::optional<LocatedVolumeOctreeNode> {
        if(input.previousResult_) {
            auto& prev = *input.previousResult_;
            tgtAssert(volumeDim == prev.getDimensions(), "prev result: Dimension mismatch");
            tgtAssert(prev.getRootNode(), "prev result: No root");
            tgtAssert(input.octree_.getNumLevels() == prev.getNumLevels(), "prev result: num levels mismatch");
            return LocatedVolumeOctreeNode(prev.getRootNode(), maxLevel, tgt::svec3(0), prev.getDimensions());
        } else {
            return boost::none;
        }
    }();

    ProfileAllocation inputTreeNodes(ramProfiler_, input.octree_.getNumNodes() * sizeof(VolumeOctreeNodeGeneric<1>));

    struct NodeToProcess {
        VolumeOctreeNode** outputNodeSlot; // Never null
        tgt::svec3 llf;
        tgt::svec3 urb;

        VolumeOctreeNode*& outputNode() {
            return *outputNodeSlot;
        }
    };

    // These contain nodes whose CHILDREN AND BRICKS are still referenced and thus must not be deleted.
    std::unordered_set<const VolumeOctreeNode*> nodesToSave;

    LocatedVolumeOctreeNode outputRootNode(nullptr, maxLevel, tgt::svec3(0), volumeDim);

    // Mark previous solution as invalid during computation in case Voreen is killed or crashes
    std::string valid_token_path = prevResultPath_ + "/" + PREV_RESULT_VALID_FILE_NAME;
    tgt::FileSystem::deleteFile(valid_token_path);

    // If computation is canceled: Delete new nodes, but spare nodes that are part of another tree.
    tgt::ScopeGuard nodeCleanup([&] () {
        freeTreeComponents(&outputRootNode.node(), nodesToSave, *brickPoolManager_);

        // The previous result is still valid after interuption _and_ cleanup
        std::ofstream f(valid_token_path);
        f.flush();
    });

    uint16_t globalMin = 0xffff;
    uint16_t globalMax = 0;

    const size_t HISTOGRAM_BUCKETS = 1 << 12;
#ifdef VRN_OCTREEWALKER_USE_OMP
    size_t num_threads = omp_get_max_threads();
    std::vector<Histogram1D> histograms(num_threads, Histogram1D(0.0, 1.0, HISTOGRAM_BUCKETS));
    LINFO("Using parallel octree walker variant.");
#else
    Histogram1D histogram(0.0, 1.0, HISTOGRAM_BUCKETS);
    LINFO("Using sequential octree walker variant.");
#endif

    std::vector<NodeToProcess> nodesToProcess;
    if(input.octree_.getRootNode()->isHomogeneous()) {
        LWARNING("Input octree consists of single, homogeneous node");

        uint16_t intensity = 0xffff/2;
        globalMin = intensity;
        globalMax = intensity;

#ifdef VRN_OCTREEWALKER_USE_OMP
        Histogram1D& histogram = histograms.at(0);
#endif
        histogram.addSample(0.5);

        tgtAssert(!outputRootNode.node_, "Rootnode should be null");
        outputRootNode.node_ = new VolumeOctreeNodeGeneric<1>(OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS, true, intensity);
    } else {
        nodesToProcess.push_back(
                NodeToProcess {
                &outputRootNode.node_,
                tgt::svec3::zero,
                volumeDim,
                }
            );
    }
    LocatedVolumeOctreeNodeConst inputRoot = input.octree_.getLocatedRootNode();

    float variance;
    if(requiresGlobalVarianceEstimate(input.noiseModel_)) {
        const size_t num_sample_bricks = 10;
        const size_t max_tries = 100;

        size_t num_sampled = 0;
        size_t try_ = 0;
        variance = 0.0f;

        auto dimInBricks = volumeDim/brickDim;
        size_t numBricks = tgt::hmul(dimInBricks);

        auto rng_seed = 1;
        std::mt19937 prng(rng_seed);
        std::uniform_int_distribution<size_t> distr(0, numBricks+1);

        const OctreeBrickPoolManagerBase& inputPoolManager = *input.octree_.getBrickPoolManager();

        while(try_ < max_tries && num_sampled < num_sample_bricks) {
            try_ += 1;

            size_t nodeIndex = distr(prng);
            tgt::svec3 posInBricks = linearCoordToCubic(nodeIndex, dimInBricks);
            tgt::svec3 samplePos = posInBricks * brickDim;

            auto node = inputRoot.findChildNode(samplePos, brickDim, 0, false);
            if(node.node().isHomogeneous()) {
                continue;
            }
            BrickPoolBrickConst brick(node.node().getBrickAddress(), brickDim, inputPoolManager);

            // Crop brick data to valid values
            auto nodeBrickDim = node.location().voxelURB() - node.location().voxelLLF();
            std::unique_ptr<VolumeAtomic<uint16_t>> brickData(brick.data().getSubVolume(nodeBrickDim));

            auto brickFloat = toVolumeAtomicFloat(*brickData);

            const int extent = 1;
            auto mean = meanFilter(brickFloat, extent);
            float var = estimateVariance(brickFloat, mean, extent);

            variance += var;
            num_sampled += 1;
        }

        variance /= num_sampled;

        LINFO("Estimated variance: " << variance);
    }

    VarianceTree newVarianceTree = VarianceTree::none();
    const VarianceTree* varianceTreePtr;
    if(input.varianceTree_.isNone() && requiresVarianceTree(input.noiseModel_)) {
        LINFO("Computing variance octree...");
        newVarianceTree = computeVarianceTree(input);
        varianceTreePtr = &newVarianceTree;
    } else {
        varianceTreePtr = &input.varianceTree_;
    }
    const VarianceTree& varianceTree = *varianceTreePtr;

    PointSegmentListGeometryVec3 foregroundSeeds;
    PointSegmentListGeometryVec3 backgroundSeeds;
    getSeedListsFromPorts(input.foregroundGeomSeeds_, foregroundSeeds);
    getSeedListsFromPorts(input.backgroundGeomSeeds_, backgroundSeeds);

    std::mutex clMutex;

    // Level order iteration => Previos level is always available
    for(int level = maxLevel; level >=0; --level) {
        // Note in the following that 1/4 seems to better represent the actual progress (rather than 1/8).
        // This may be due to the fact that the actual work we have to do happens on the (2D!) _surface_
        // of objects in the volume.
        float progressBegin = 1.0/(1 << (2 * (level + 1))); // 1/4 ^ (level+1 => next level)
        float progressEnd = 1.0/(1 << (2 * (level))); // 1/4 ^ (level)
        SubtaskProgressReporter levelProgress(progressReporter, tgt::vec2(progressBegin, progressEnd));

        LINFO("Level " << level << ": " << nodesToProcess.size() << " Nodes to process.");

        std::vector<NodeToProcess> nextNodesToProcess;
        ProfileAllocation nextNodesToProcessAllocation(ramProfiler_, 0);


        const int numNodes = nodesToProcess.size();
        ThreadedTaskProgressReporter parallelProgress(levelProgress, numNodes);
        bool aborted = false;

#ifdef VRN_OCTREEWALKER_USE_OMP
        std::atomic<size_t> histogramIDCounter(0);
#pragma omp parallel
        {
        size_t histogramID = histogramIDCounter.fetch_add(1);
        Histogram1D& histogram = histograms.at(histogramID);

#pragma omp for schedule(dynamic, 1) // Schedule: Process nodes/bricks locally to utilize brick cache
#endif
        for (int nodeId = 0; nodeId < numNodes; ++nodeId) {

#ifdef VRN_OCTREEWALKER_USE_OMP
            if(aborted) {
                continue;
            }
            if(parallelProgress.reportStepDone()) {
                #pragma omp critical
                aborted = true;
                continue;
            }
#else
            if(parallelProgress.reportStepDone()) {
                aborted = true;
                break;
            }
#endif

            // Make sure to hit LRU cache: Go from back to front
            auto& node = nodesToProcess[numNodes-nodeId-1];
            ProfileAllocation inputTreeNodes(ramProfiler_, nodesToProcess.size() * sizeof(NodeToProcess));

            bool inVolume = tgt::hand(tgt::lessThan(node.llf, volumeDim));
            if(!inVolume) {
                ramProfiler_.allocate(sizeof(VolumeOctreeNodeGeneric<1>));
                node.outputNode() = new VolumeOctreeNodeGeneric<1>(OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS, false);
                continue;
            }

            uint16_t min = 0xffff;
            uint16_t max = 0;
            uint16_t avg = 0xffff/2;
            uint64_t newBrickAddr;
            bool hasNewSeedsConflicts;
            {
                VolumeOctreeNodeLocation outputNodeGeometry(level, node.llf, node.urb);
                switch (input.noiseModel_) {
                    case RW_NOISE_GAUSSIAN_BIAN_MEAN:
                        newBrickAddr = processOctreeBrick<RW_NOISE_GAUSSIAN_BIAN_MEAN>({}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                        break;
                    case RW_NOISE_GAUSSIAN_BIAN_MEDIAN:
                        newBrickAddr = processOctreeBrick<RW_NOISE_GAUSSIAN_BIAN_MEDIAN>({}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                        break;
                    case RW_NOISE_TTEST:
                        newBrickAddr = processOctreeBrick<RW_NOISE_TTEST>({input.parameterEstimationNeighborhoodExtent_}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                        break;
                    case RW_NOISE_GAUSSIAN:
                        newBrickAddr = processOctreeBrick<RW_NOISE_GAUSSIAN>({input.parameterEstimationNeighborhoodExtent_}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                        break;
                    case RW_NOISE_VARIABLE_GAUSSIAN:
                        newBrickAddr = processOctreeBrick<RW_NOISE_VARIABLE_GAUSSIAN>({input.parameterEstimationNeighborhoodExtent_}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                        break;
                    case RW_NOISE_POISSON:
                        newBrickAddr = processOctreeBrick<RW_NOISE_POISSON>({input.parameterEstimationNeighborhoodExtent_}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                        break;
                    case RW_NOISE_GAUSSIAN_HIERARCHICAL:
                        newBrickAddr = processOctreeBrick<RW_NOISE_GAUSSIAN_HIERARCHICAL>({input.parameterEstimationNeighborhoodExtent_, variance}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                    case RW_NOISE_POISSON_HIERARCHICAL:
                        newBrickAddr = processOctreeBrick<RW_NOISE_POISSON_HIERARCHICAL>({input.parameterEstimationNeighborhoodExtent_}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                        break;
                    case RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL:
                        newBrickAddr = processOctreeBrick<RW_NOISE_VARIABLE_GAUSSIAN_HIERARCHICAL>({input.parameterEstimationNeighborhoodExtent_}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                        break;
                    case RW_NOISE_TTEST_HIERARCHICAL:
                        newBrickAddr = processOctreeBrick<RW_NOISE_TTEST_HIERARCHICAL>({input.parameterEstimationNeighborhoodExtent_}, input, outputNodeGeometry, histogram, min, max, avg, hasNewSeedsConflicts, outputBrickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, varianceTree, foregroundSeeds, backgroundSeeds, clMutex, ramProfiler_, vramProfiler_);
                        break;
                    default:
                        tgtAssert(false, "Invalid noise model selected");
                }
            }

            #pragma omp critical
            {
                globalMin = std::min(globalMin, min);
            }
            #pragma omp critical
            {
                globalMax = std::max(globalMax, max);
            }

            bool childrenToProcess = newBrickAddr != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS && level > 0;
            VolumeOctreeNode* newNode = nullptr;

            // Check if previous result is sufficiently close to current one
            // and if there are no (new) labels which have not been processed
            // in this branch due to conflicts.
            // If this is the case: Take branch from old octree.
            if(childrenToProcess && prevRoot && !hasNewSeedsConflicts) {
                auto prevNode = prevRoot->findChildNode(node.llf, brickDim, level, false);

                if(prevNode.node().hasBrick()
                    && prevNode.location().level() == level
                    && absdiffRel(prevNode.node().getMinValue(), min) < input.incrementalSimilarityThreshold_
                    && absdiffRel(prevNode.node().getMaxValue(), max) < input.incrementalSimilarityThreshold_
                    && absdiffRel(prevNode.node().getAvgValue(), avg) < input.incrementalSimilarityThreshold_
                ) {
                    uint16_t maxDiff = 0;
                    const tgt::svec3 begin(0);
                    const tgt::svec3 end(prevNode.geometry_.brickDimensions());
                    {
                        BrickPoolBrickConst prevBrick(prevNode.node().getBrickAddress(), brickDim, outputBrickPoolManager);
                        BrickPoolBrickConst currentBrick(newBrickAddr, brickDim, outputBrickPoolManager);

                        VRN_FOR_EACH_VOXEL(pos, begin, end) {
                            uint16_t current = currentBrick.data().voxel(pos);
                            uint16_t prev = prevBrick.data().voxel(pos);
                            uint16_t diff = absdiff(current, prev);
                            maxDiff = std::max(maxDiff, diff);
                        }
                    }
                    float reldiff = static_cast<float>(maxDiff) / 0xffff;
                    if(reldiff < input.incrementalSimilarityThreshold_) {
                        outputBrickPoolManager.deleteBrick(newBrickAddr);
                        tgtAssert((level == 0) == prevNode.node_->isLeaf(), "Previous node is leaf, but not at bottom");

                        newNode = prevNode.node_;
                        childrenToProcess = false;

                        #pragma omp critical
                        {
                            nodesToSave.insert(newNode);
                            ramProfiler_.allocate(sizeof(size_t));
                        }
                        //LINFO("Using similar branch from previous iteration");
                    }
                }
            }

            if(!newNode) {
                auto newNodeGeneric = new VolumeOctreeNodeGeneric<1>(newBrickAddr, true);
                ramProfiler_.allocate(sizeof(VolumeOctreeNodeGeneric<1>));
                newNodeGeneric->avgValues_[0] = avg;
                newNodeGeneric->minValues_[0] = min;
                newNodeGeneric->maxValues_[0] = max;

                newNode = newNodeGeneric;
            }

            if(childrenToProcess) {

                tgt::svec3 childBrickSize = brickDim * (size_t(1) << (level-1));
                for(auto child : OCTREEWALKER_CHILD_POSITIONS) {
                    const size_t childId = volumeCoordsToIndex(child, tgt::svec3::two);

                    tgt::svec3 start = node.llf + childBrickSize * child;
                    tgt::svec3 end = tgt::min(start + childBrickSize, volumeDim);
                    #pragma omp critical
                    {
                        nextNodesToProcess.push_back(
                                NodeToProcess {
                                &newNode->children_[childId],
                                start,
                                end,
                            }
                        );
                    }
                    nextNodesToProcessAllocation.increaseBy(sizeof(NodeToProcess));
                }
            }

            node.outputNode() = newNode;
        }
#ifdef VRN_OCTREEWALKER_USE_OMP
        }
#endif
        if(aborted) {
            throw boost::thread_interrupted();
        }

        nodesToProcess = nextNodesToProcess;
    }

    nodeCleanup.dismiss();

#ifdef VRN_OCTREEWALKER_USE_OMP
    Histogram1D histogram(0.0, 1.0, HISTOGRAM_BUCKETS);
    for(auto& h : histograms) {
        for(size_t i = 0; i < HISTOGRAM_BUCKETS; ++i) {
            histogram.increaseBucket(i, h.getBucket(i));
        }
    }
#endif

    std::vector<Histogram1D*> octreeHistograms;
    octreeHistograms.push_back(new Histogram1D(histogram));
    auto octree = new VolumeOctree(outputRootNode.node_, &outputBrickPoolManager, std::move(octreeHistograms), brickDim, input.octree_.getDimensions(), numChannels);
    auto output = tgt::make_unique<Volume>(octree, &input.volume_);

    float min = static_cast<float>(globalMin)/0xffff;
    float max = static_cast<float>(globalMax)/0xffff;
    output->addDerivedData(new VolumeMinMax(min, max, min, max));
    output->addDerivedData(new VolumeHistogramIntensity(histogram));
    output->setRealWorldMapping(RealWorldMapping(1.0, 0.0f, ""));

    // It's not really feasible, but also not important to compute an actual hash
    boost::uuids::basic_random_generator<boost::mt19937> uuidGenerator{};
    boost::uuids::uuid uuid = uuidGenerator();
    std::string uuidstr = boost::lexical_cast<std::string>(uuid);
    std::string hash = VoreenHash::getHash(uuidstr.data(), uuidstr.size());
    output->setHash(hash);

    auto finish = clock::now();
    return ComputeOutput (
        OctreeWalkerPreviousResult(
            *octree,
            std::move(output),
            std::move(foregroundSeeds),
            std::move(backgroundSeeds)
            ),
        std::move(nodesToSave),
        std::move(newVarianceTree),
        finish - start
    );
}

void OctreeWalker::processComputeOutput(ComputeOutput output) {
    LINFO("Total runtime: " << output.duration_.count() << " sec");
#ifdef VRN_MODULE_RANDOMWALKER_PROFILING_ENABLED
    LINFO("Total ram: " << ramProfiler_.peak());
    LINFO("Total vram: " << vramProfiler_.peak());
    ramProfiler_.reset();
    vramProfiler_.reset();
#endif

    // Set new output
    outportProbabilities_.setData(&output.result_.volume(), false);

    // previousOctree_ is now not referenced anymore, so we are free to clean up.
    if(previousResult_) {
        std::move(*previousResult_).destroyButRetainNodes(output.sharedNodes_);
    }
    previousResult_ = std::move(output.result_);

    // Update origin property of newly created volume.
    std::string previousResultDir = resultPath_.get();
    std::string previousResultFile = tgt::FileSystem::cleanupPath(previousResultDir + "/" + PREV_RESULT_FILE_NAME);
    previousResult_->volume().setOrigin(previousResultFile);

    // Write out updated octree (metadata) cache file itself.
    XmlSerializer s(previousResultDir);
    std::ofstream fs(previousResultFile);
    s.serialize(PREV_RESULT_OCTREE_KEY, previousResult_->octree());
    s.serialize(PREV_RESULT_FOREGROUND_KEY, previousResult_->foregroundSeeds_);
    s.serialize(PREV_RESULT_BACKGROUND_KEY, previousResult_->backgroundSeeds_);
    s.write(fs);

    if(!output.varianceTree_.isNone()) {
        varianceTree_ = std::move(output.varianceTree_);
    }

    // Finally, the cache result can be marked as valid
    std::ofstream f(prevResultPath_ + "/" + PREV_RESULT_VALID_FILE_NAME);
    f.flush();
}
void OctreeWalker::clearPreviousResults() {
    // First: Reset output
    outportProbabilities_.setData(nullptr, false);

    // previousOctree_ is now not referenced anymore, so we are free to clean up.
    previousResult_ = boost::none;

    if(brickPoolManager_ && brickPoolManager_->isInitialized()) {
        brickPoolManager_->deinitialize();
    }
    brickPoolManager_.reset(nullptr);
    prevResultPath_ = "";
}

const VoreenBlas* OctreeWalker::getVoreenBlasFromProperties() const {

#ifdef VRN_MODULE_OPENMP
    if (conjGradImplementation_.isSelected("blasMP")) {
        return &voreenBlasMP_;
    }
#endif
#ifdef VRN_MODULE_OPENCL
    if (conjGradImplementation_.isSelected("blasCL")) {
        return &voreenBlasCL_;
    }
#endif

    return &voreenBlasCPU_;
}

VolumeOctree& OctreeWalkerPreviousResult::octree() {
    return *octree_;
}
VolumeBase& OctreeWalkerPreviousResult::volume() {
    return *volume_;
}
OctreeWalkerPreviousResult::OctreeWalkerPreviousResult(OctreeWalkerPreviousResult&& other)
{
    *this = std::move(other);
}
OctreeWalkerPreviousResult::OctreeWalkerPreviousResult(VolumeOctree& octree, std::unique_ptr<VolumeBase>&& volume, PointSegmentListGeometryVec3 foregroundSeeds, PointSegmentListGeometryVec3 backgroundSeeds)
    : octree_(&octree)
    , volume_(std::move(volume))
    , foregroundSeeds_(std::move(foregroundSeeds))
    , backgroundSeeds_(std::move(backgroundSeeds))
{
}

OctreeWalkerPreviousResult& OctreeWalkerPreviousResult::operator=(OctreeWalkerPreviousResult&& other) {
    volume_ = std::move(other.volume_);
    octree_ = other.octree_;
    foregroundSeeds_ = std::move(other.foregroundSeeds_);
    backgroundSeeds_ = std::move(other.backgroundSeeds_);

    other.octree_ = nullptr;
    return *this;
}

OctreeWalkerPreviousResult::~OctreeWalkerPreviousResult() {
    if(volume_) {
        tgtAssert(octree_, "Octree should be present if volume is");
        auto res = std::move(*octree_).decompose();
        // Brickpoolmanager reference is not required here. The important thing
        // is that the previous result does not deconstruct the
        // brickPoolManager (which is owned by the propessor).

        // Clean up (entire) old tree
        // If this is not desired, destroyButRetainNodes has to be called instead.
        freeNodes(res.second);
    }
}
void OctreeWalkerPreviousResult::destroyButRetainNodes(std::unordered_set<const VolumeOctreeNode*>& nodesToSave) && {
    auto tmp = std::move(*this);
    tgtAssert(tmp.octree_, "No octree");
    tgtAssert(tmp.volume_, "No volume");

    auto res = std::move(*tmp.octree_).decompose();

    // Clean up old tree
    freeTreeComponents(res.second, nodesToSave, *res.first);

    tmp.volume_.reset(nullptr); // Free (now empty) octree itself, a representation of the volume
                                // -- sans nodes and brick pool manager: Nodes are already freed
                                // and the brick pool manager is owned by the Processor.
    tmp.octree_ = nullptr;          // Unset pointer to invalid now memory
}
bool OctreeWalkerPreviousResult::isPresent() const {
    tgtAssert((octree_ == nullptr) == (volume_ == nullptr), "Octree/Volume present mismatch");
    return volume_ != nullptr;
}

}   // namespace
