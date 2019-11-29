/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "../solver/randomwalkerseeds.h"
#include "../solver/randomwalkerweights.h"

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
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "tgt/vector.h"
#include "tgt/memory.h"

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include <climits>

#ifdef VRN_MODULE_OPENMP
#include <omp.h>
#endif

namespace voreen {

#if defined(VRN_MODULE_OPENMP) && 1
#define VRN_OCTREEWALKER_USE_OMP
#endif

//#define VRN_OCTREEWALKER_MEAN_NOT_MEDIAN

namespace {

inline size_t volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

inline size_t volumeCoordsToIndex(const tgt::ivec3& coords, const tgt::ivec3& dim) {
    return coords.z*dim.y*dim.x + coords.y*dim.x + coords.x;
}

static void freeTreeComponents(VolumeOctreeNode* root, std::unordered_set<const VolumeOctreeNode*>& nodesToSave, OctreeBrickPoolManagerMmap& brickPoolManager) {
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
    return static_cast<float>(val)/0xffff;
}
}


const std::string OctreeWalker::loggerCat_("voreen.RandomWalker.OctreeWalker");

OctreeWalker::OctreeWalker()
    : AsyncComputeProcessor<OctreeWalkerInput, OctreeWalkerOutput>()
    , inportVolume_(Port::INPORT, "volume.input")
    , inportForegroundSeeds_(Port::INPORT, "geometry.seedsForeground", "geometry.seedsForeground", true)
    , inportBackgroundSeeds_(Port::INPORT, "geometry.seedsBackground", "geometry.seedsBackground", true)
    , outportProbabilities_(Port::OUTPORT, "volume.probabilities", "volume.probabilities", false)
    , usePrevProbAsInitialization_("usePrevProbAsInitialization", "Use Previous Probabilities as Initialization", false, Processor::VALID, Property::LOD_ADVANCED)
    , minEdgeWeight_("minEdgeWeight", "Min Edge Weight: 10^(-t)", 5, 0, 10)
    , preconditioner_("preconditioner", "Preconditioner")
    , errorThreshold_("errorThreshold", "Error Threshold: 10^(-t)", 2, 0, 10)
    , maxIterations_("conjGradIterations", "Max Iterations", 1000, 1, 5000)
    , conjGradImplementation_("conjGradImplementation", "Implementation")
    , homogeneityThreshold_("homogeneityThreshold", "Homogeneity Threshold", 0.01, 0.0, 1.0)
    , incrementalSimilarityThreshold_("incrementalSimilarityThreshold", "Incremental Similarity Treshold", 0.01, 0.0, 1.0)
    , previousOctree_(nullptr)
    , previousVolume_(nullptr)
    , brickPoolManager_(nullptr)
{
    // ports
    addPort(inportVolume_);
        ON_CHANGE(inportVolume_, OctreeWalker, clearPreviousResults);
    addPort(inportForegroundSeeds_);
    addPort(inportBackgroundSeeds_);
    addPort(outportProbabilities_);

    addProperty(usePrevProbAsInitialization_);

    // random walker properties
    setPropertyGroupGuiName("rwparam", "Random Walker Parametrization");
    addProperty(minEdgeWeight_);
        minEdgeWeight_.setGroupID("rwparam");
    addProperty(homogeneityThreshold_);
        homogeneityThreshold_.setGroupID("rwparam");
        homogeneityThreshold_.adaptDecimalsToRange(5);
    addProperty(incrementalSimilarityThreshold_);
        incrementalSimilarityThreshold_.setGroupID("rwparam");
        incrementalSimilarityThreshold_.adaptDecimalsToRange(5);

    // conjugate gradient solver
    setPropertyGroupGuiName("conjGrad", "Conjugate Gradient Solver");
    addProperty(preconditioner_);
        preconditioner_.addOption("none", "None");
        preconditioner_.addOption("jacobi", "Jacobi");
        preconditioner_.select("jacobi");
        preconditioner_.setGroupID("conjGrad");
    addProperty(errorThreshold_);
        errorThreshold_.setGroupID("conjGrad");
    addProperty(maxIterations_);
        maxIterations_.setGroupID("conjGrad");
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

    updateGuiState();
}

void OctreeWalker::deinitialize() {
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

OctreeWalker::ComputeInput OctreeWalker::prepareComputeInput() {
    //edgeWeightTransFunc_.setVolume(inportVolume_.getData());

    tgtAssert(inportVolume_.hasData(), "no input volume");

    // clear previous results and update property ranges, if input volume has changed
    if (inportVolume_.hasChanged()) {
        outportProbabilities_.setData(0);
    }
    auto vol = inportVolume_.getThreadSafeData();

    // clear previous results and update property ranges, if input volume has changed
    if (!vol->hasRepresentation<VolumeOctree>()) {
        throw InvalidInputException("No octree Representation", InvalidInputException::S_ERROR);
    }
    const VolumeOctree* octreePtr = vol->getRepresentation<VolumeOctree>();
    tgtAssert(octreePtr, "No octree");


    // select BLAS implementation and preconditioner
    const VoreenBlas* voreenBlas = getVoreenBlasFromProperties();
    VoreenBlas::ConjGradPreconditioner precond = VoreenBlas::NoPreconditioner;
    if (preconditioner_.isSelected("jacobi"))
        precond = VoreenBlas::Jacobi;

    float errorThresh = 1.f / pow(10.f, static_cast<float>(errorThreshold_.get()));
    int maxIterations = maxIterations_.get();

    return ComputeInput {
        previousOctree_,
        brickPoolManager_,
        *vol,
        *octreePtr,
        inportForegroundSeeds_.getThreadSafeAllData(),
        inportBackgroundSeeds_.getThreadSafeAllData(),
        minEdgeWeight_.get(),
        voreenBlas,
        precond,
        errorThresh,
        maxIterations,
        homogeneityThreshold_.get(),
        incrementalSimilarityThreshold_.get(),
    };
}

static void getSeedListsFromPorts(std::vector<PortDataPointer<Geometry>>& geom, PointSegmentListGeometry<tgt::vec3>& seeds) {

    for (size_t i=0; i<geom.size(); i++) {
        const PointSegmentListGeometry<tgt::vec3>* seedList = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(geom.at(i).get());
        if (!seedList)
            LWARNINGC(OctreeWalker::loggerCat_, "Invalid geometry. PointSegmentListGeometry<vec3> expected.");
        else {
            auto transformMat = seedList->getTransformationMatrix();
            for (int j=0; j<seedList->getNumSegments(); j++) {
                std::vector<tgt::vec3> points;
                for(auto& vox : seedList->getSegment(j)) {
                    points.push_back(transformMat.transform(vox));
                }
                seeds.addSegment(points);
            }
        }
    }
}

static VolumeOctreeNode* findLeafNodeFor(VolumeOctreeNode* root, tgt::svec3& llf, tgt::svec3& urb, size_t& level, const tgt::svec3& point, const tgt::svec3& brickDataSize, size_t targetLevel) {
    tgtAssert(tgt::hand(tgt::lessThanEqual(llf, point)) && tgt::hand(tgt::lessThan(point, urb)), "Invalid point pos");
    tgtAssert(root, "No root");

    if(root->isLeaf() || level == targetLevel) {
        return root;
    }

    tgtAssert(level >= 0, "Invalid level");
    tgt::svec3 newLlf = llf;
    tgt::svec3 newUrb = urb;
    size_t newLevel = level - 1;
    tgt::svec3 brickSize = brickDataSize * (1UL << newLevel);
    size_t index = 0;
    if(point.x >= llf.x+brickSize.x) {
        index += 1;
        newLlf.x = llf.x+brickSize.x;
    } else {
        newUrb.x = llf.x+brickSize.x;
    }
    if(point.y >= llf.y+brickSize.y) {
        index += 2;
        newLlf.y = llf.y+brickSize.y;
    } else {
        newUrb.y = llf.y+brickSize.y;
    }
    if(point.z >= llf.z+brickSize.z) {
        index += 4;
        newLlf.z = llf.z+brickSize.z;
    } else {
        newUrb.z = llf.z+brickSize.z;
    }

    VolumeOctreeNode* child = root->children_[index];
    tgtAssert(child, "No child in non leaf node");

    if(child->isHomogeneous()) {
        // Parent has better resolution
        return root;
    }
    level = newLevel;
    urb = newUrb;
    llf = newLlf;
    return findLeafNodeFor(child, llf, urb, level, point, brickDataSize, targetLevel);
}
struct OctreeWalkerNodeGeometry {
    OctreeWalkerNodeGeometry(size_t level, tgt::svec3 llf, tgt::svec3 urb)
        : level_(level)
        , llf_(llf)
        , urb_(urb)
    {
    }

    tgt::svec3 voxelDimensions() const {
        return urb_ - llf_;
    }

    tgt::svec3 brickDimensions() const {
        return voxelDimensions() / scale();
    }

    size_t scale() const {
        return 1 << level_;
    }

    tgt::mat4 voxelToBrick() const {
        return tgt::mat4::createScale(tgt::vec3(1.0f/scale())) * tgt::mat4::createTranslation(-tgt::vec3(llf_));
    }

    tgt::mat4 brickToVoxel() const {
        return tgt::mat4::createTranslation(llf_) * tgt::mat4::createScale(tgt::vec3(scale()));
    }

    size_t level_;
    tgt::svec3 llf_;
    tgt::svec3 urb_;
};
struct OctreeWalkerNode {
    OctreeWalkerNode(VolumeOctreeNode* node, size_t level, tgt::svec3 llf, tgt::svec3 urb)
        : node_(node)
        , geometry_(level, llf, urb)
    {
    }
    OctreeWalkerNode findChildNode(const tgt::svec3& point, const tgt::svec3& brickDataSize, size_t targetLevel) const {
        size_t level = geometry_.level_;
        tgt::svec3 llf = geometry_.llf_;
        tgt::svec3 urb = geometry_.urb_;

        tgtAssert(level >= targetLevel, "Invalid target level");

        VolumeOctreeNode* node = findLeafNodeFor(node_, llf, urb, level, point, brickDataSize, targetLevel);
        return OctreeWalkerNode(node, level, llf, urb);
    }

    VolumeOctreeNode& node() {
        return *node_;
    }

    VolumeOctreeNode* node_; // Never null
    OctreeWalkerNodeGeometry geometry_;
};

// TODO Move OctreeBrickPoolManager and generate from member function => Const and mutable variant
struct OctreeWalkerNodeBrick {
    OctreeWalkerNodeBrick() = delete;
    OctreeWalkerNodeBrick(const OctreeWalkerNodeBrick&) = delete;
    OctreeWalkerNodeBrick& operator=(const OctreeWalkerNodeBrick&) = delete;

    OctreeWalkerNodeBrick(uint64_t addr, const tgt::svec3& brickDataSize, const OctreeBrickPoolManagerBase& pool)
        : addr_(addr)
        , data_(pool.getWritableBrick(addr_), brickDataSize, false) // data is not owned!
        , pool_(pool)
    {
    }
    ~OctreeWalkerNodeBrick() {
        pool_.releaseBrick(addr_, OctreeBrickPoolManagerBase::WRITE);
    }
    float getVoxelNormalized(const tgt::svec3& pos) const {
        return brickToNorm(data_.voxel(pos));
    }
    uint64_t addr_;
    VolumeAtomic<uint16_t> data_;
    const OctreeBrickPoolManagerBase& pool_;
};
struct OctreeWalkerNodeBrickConst {
    OctreeWalkerNodeBrickConst() = delete;
    OctreeWalkerNodeBrickConst(const OctreeWalkerNodeBrickConst&) = delete;
    OctreeWalkerNodeBrickConst& operator=(const OctreeWalkerNodeBrickConst&) = delete;

    OctreeWalkerNodeBrickConst(uint64_t addr, const tgt::svec3& brickDataSize, const OctreeBrickPoolManagerBase& pool)
        : addr_(addr)
        , data_(const_cast<uint16_t*>(pool.getBrick(addr_)), brickDataSize, false) // data is not owned!
        , pool_(pool)
    {
    }
    ~OctreeWalkerNodeBrickConst() {
        pool_.releaseBrick(addr_, OctreeBrickPoolManagerBase::READ);
    }

    float getVoxelNormalized(const tgt::svec3& pos) const {
        return brickToNorm(data_.voxel(pos));
    }
    uint64_t addr_;
    const VolumeAtomic<uint16_t> data_;
    const OctreeBrickPoolManagerBase& pool_;
};
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
    tgt::svec3 centerBrickUrb_; // ... of seed buffer
    tgt::svec3 dimensions_;
    tgt::mat4 voxelToCenterBrick_;
    float min_;
    float max_;
    float avg_;

    static BrickNeighborhood empty(tgt::svec3 dimensions, int scale) {
        return BrickNeighborhood {
            VolumeAtomic<float>(tgt::svec3(0)),
            tgt::svec3(0),
            dimensions,
            dimensions,
            tgt::mat4::identity,
            0.0f,
            0.0f,
            0.0f,
        };
    }
    static BrickNeighborhood fromNode(const OctreeWalkerNodeGeometry& current, size_t sampleLevel, const OctreeWalkerNode& root, const tgt::svec3& brickBaseSize, const OctreeBrickPoolManagerBase& brickPoolManager) {
        const tgt::svec3 volumeDim = root.geometry_.voxelDimensions();

        const tgt::mat4 brickToVoxel = current.brickToVoxel();
        const tgt::mat4 voxelToBrick = current.voxelToBrick();

        //const tgt::ivec3 neighborhoodSize = brickBaseSize;
        //const tgt::ivec3 neighborhoodSize = tgt::ivec3(2);
        const tgt::ivec3 neighborhoodSize = brickBaseSize/8UL;

        const tgt::ivec3 brickLlf(0);
        const tgt::ivec3 brickUrb = voxelToBrick.transform(current.urb_);

        const tgt::svec3 voxelLlf = tgt::max(tgt::vec3(0),         brickToVoxel.transform(brickLlf - neighborhoodSize));
        const tgt::svec3 voxelUrb = tgt::min(tgt::vec3(volumeDim), brickToVoxel.transform(brickUrb + neighborhoodSize));

        const tgt::ivec3 regionLlf = voxelToBrick.transform(voxelLlf);
        const tgt::ivec3 regionUrb = voxelToBrick.transform(voxelUrb);

        const tgt::svec3 regionDim = regionUrb - regionLlf;

        VolumeAtomic<float> output(regionDim);

        float min = std::numeric_limits<float>::infinity();
        float max = -std::numeric_limits<float>::infinity();
        float sum = 0.0f;

        VRN_FOR_EACH_VOXEL(blockIndex, tgt::svec3(0), tgt::svec3(3)) {
            tgt::ivec3 blockLlf;
            tgt::ivec3 blockUrb;
            for(int dim=0; dim<3; ++dim) {
                switch(blockIndex[dim]) {
                    case 0: {
                        blockLlf[dim] = regionLlf[dim];
                        blockUrb[dim] = brickLlf[dim];
                        break;
                    }
                    case 1: {
                        blockLlf[dim] = brickLlf[dim];
                        blockUrb[dim] = brickUrb[dim];
                        break;
                    }
                    case 2: {
                        blockLlf[dim] = brickUrb[dim];
                        blockUrb[dim] = regionUrb[dim];
                        break;
                    }
                }
            }
            tgt::svec3 blockDimensions = blockUrb - blockLlf;
            if(tgt::hor(tgt::equal(blockDimensions, tgt::svec3(0)))) {
                continue;
            }
            tgt::svec3 samplePoint = brickToVoxel.transform(blockLlf);
            auto node = root.findChildNode(samplePoint, brickBaseSize, sampleLevel);
            if(node.node().hasBrick()) {
                OctreeWalkerNodeBrickConst brick(node.node().getBrickAddress(), brickBaseSize, brickPoolManager);

                tgt::mat4 centerToSampleBrick = node.geometry_.voxelToBrick() * brickToVoxel;
                VRN_FOR_EACH_VOXEL(point, blockLlf, blockUrb) {
                    tgt::vec3 samplePos = centerToSampleBrick.transform(point);
                    samplePos = tgt::clamp(samplePos, tgt::vec3(0), tgt::vec3(node.geometry_.brickDimensions() - tgt::svec3(1)));
                    float val = brick.getVoxelNormalized(samplePos);
                    tgt::vec3 neighborhoodBufferPos = point - regionLlf;
                    output.setVoxelNormalized(val, neighborhoodBufferPos);
                    min = std::min(val, min);
                    max = std::max(val, max);
                    sum += val;
                }
            } else {
                float val = static_cast<float>(node.node().getAvgValue())/0xffff;
                min = std::min(val, min);
                max = std::max(val, max);
                sum += val * tgt::hmul(blockUrb - blockLlf);
                VRN_FOR_EACH_VOXEL(point, blockLlf, blockUrb) {
                    tgt::vec3 neighborhoodBufferPos = point - regionLlf;
                    output.setVoxelNormalized(val, neighborhoodBufferPos);
                }
            }
        }
        float avg = sum / output.getNumVoxels();
        return BrickNeighborhood {
            std::move(output),
            -regionLlf,
            -regionLlf+brickUrb,
            regionDim,
            voxelToBrick,
            min,
            max,
            avg
        };
    }
};

struct NSmallestHeap14 {
    NSmallestHeap14()
        : data()
        , numElements(0)
    {
        data[0] = FLT_MAX; //sentinel for first 14 values (upper push branch)
        data[15] = FLT_MIN; //sentinel for child2 of index 7
    }
    float nthlargest() {
        return data[1];
    }
    void push(float val) {
        if(numElements < 14) {
            numElements++;
            uint8_t i=numElements;
            data[i] = val;
            uint8_t parent;

            parent = i>>1; if(data[parent] >= val) { return; } std::swap(data[parent], data[i]); i=parent; // 15->7
            parent = i>>1; if(data[parent] >= val) { return; } std::swap(data[parent], data[i]); i=parent; // 7->3
            parent = i>>1; if(data[parent] >= val) { return; } std::swap(data[parent], data[i]);           // 3->1
        } else {
            if(val < nthlargest()) {
                uint8_t i = 1;
                data[i] = val;
                uint8_t child1, child2, c;

                // 1->3
                child1 = i<<1;
                child2 = child1+1;
                c = data[child1] > data[child2] ? child1 : child2;
                if(data[c] > val) { std::swap(data[c], data[i]); i = c; } else { return; }

                // 3->7
                child1 = i<<1;
                child2 = child1+1;
                c = data[child1] > data[child2] ? child1 : child2;
                if(data[c] > val) { std::swap(data[c], data[i]); i = c; } else { return; }

                // 7->15
                child1 = i<<1;
                child2 = child1+1;
                c = data[child1] > data[child2] ? child1 : child2;
                if(data[c] > val) { std::swap(data[c], data[i]); }
            }
        }
    }
    void clear() {
        numElements = 0;
    }
    std::array<float, 16> data;
    uint8_t numElements;
};

class RandomWalkerSeedsBrick : public RandomWalkerSeeds {
    static const float CONFLICT;
    static const float UNLABELED;
    static const float FOREGROUND;
    static const float BACKGROUND;
public:
    RandomWalkerSeedsBrick(tgt::svec3 bufferDimensions, tgt::mat4 voxelToSeeds, const PointSegmentListGeometryVec3& foregroundSeedList, const PointSegmentListGeometryVec3& backgroundSeedList)
        : seedBuffer_(bufferDimensions)
        , numConflicts_(0)
    {
        numSeeds_ = 0;
        seedBuffer_.fill(UNLABELED);

        // foreground geometry seeds
        auto collectLabelsFromGeometry = [&] (const PointSegmentListGeometryVec3& seedList, float label) {
            for (int m=0; m<seedList.getNumSegments(); m++) {
                const std::vector<tgt::vec3>& foregroundPoints = seedList.getData()[m];
                if (foregroundPoints.empty())
                    continue;
                for (size_t i=0; i<foregroundPoints.size()-1; i++) {
                    tgt::vec3 left = voxelToSeeds*foregroundPoints[i];
                    tgt::vec3 right = voxelToSeeds*foregroundPoints[i+1];
                    tgt::vec3 dir = tgt::normalize(right - left);
                    for (float t=0.f; t<tgt::length(right-left); t += 1.f) {
                        tgt::ivec3 point = tgt::iround(left + t*dir);
                        if(tgt::hor(tgt::lessThan(point, tgt::ivec3::zero)) || tgt::hor(tgt::greaterThanEqual(point, tgt::ivec3(bufferDimensions)))) {
                            continue;
                        }
                        float& seedVal = seedBuffer_.voxel(point);
                        if (seedVal == UNLABELED) {
                            seedVal = label;
                            ++numSeeds_;
                        } else if(seedVal != label && seedVal != CONFLICT) {
                            seedVal = CONFLICT;
                            --numSeeds_;
                            ++numConflicts_;
                        }
                    }
                }
            }
        };
        collectLabelsFromGeometry(foregroundSeedList, FOREGROUND);
        collectLabelsFromGeometry(backgroundSeedList, BACKGROUND);
    }
    RandomWalkerSeedsBrick(RandomWalkerSeedsBrick&& other) = default;

    void addNeighborhoodBorderSeeds(const BrickNeighborhood& neighborhood, tgt::svec3 volumeDimensions) {
        tgtAssert(neighborhood.data_.getDimensions() == neighborhood.dimensions_, "Invalid buffer dimensions");

        tgt::ivec3 volumeLlfSeeds = neighborhood.voxelToNeighborhood().transform(tgt::vec3(0.0));
        tgt::ivec3 volumeUrbSeeds = neighborhood.voxelToNeighborhood().transform(volumeDimensions);
        auto collectLabelsFromNeighbor = [&] (size_t dim, size_t sliceIndex) {
            tgt::svec3 begin(0);
            tgt::svec3 end(neighborhood.dimensions_);

            begin[dim] = sliceIndex;
            end[dim] = sliceIndex+1;

            // Do not collect parent level border labels at the border of the volume. There is no additional information in this case.
            if(begin[dim] == volumeLlfSeeds[dim] || end[dim] == volumeUrbSeeds[dim]) {
                return;
            }

            VRN_FOR_EACH_VOXEL(seed, begin, end) {
                float& seedVal = seedBuffer_.voxel(seed);
                if (seedVal == UNLABELED) {
                    float val = neighborhood.data_.voxel(seed);
                    seedVal = val;
                    ++numSeeds_;
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
        return val != UNLABELED && val != CONFLICT;
    }
    virtual bool isSeedPoint(const tgt::ivec3& voxel) const {
        float val = seedBuffer_.voxel(voxel);
        return val != UNLABELED && val != CONFLICT;
    }
    virtual float getSeedValue(size_t index) const {
        tgtAssert(isSeedPoint(index), "Getting seed value from non-seed");
        return seedBuffer_.voxel(index);
    }
    virtual float getSeedValue(const tgt::ivec3& voxel) const {
        tgtAssert(isSeedPoint(voxel), "Getting seed value from non-seed");
        return seedBuffer_.voxel(voxel);
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
    uint64_t numConflicts() const {
        return numConflicts_;
    }

private:
    VolumeAtomic<float> seedBuffer_;
    uint64_t numConflicts_;
};

const float RandomWalkerSeedsBrick::CONFLICT = -2.0f;
const float RandomWalkerSeedsBrick::UNLABELED = -1.0f;
const float RandomWalkerSeedsBrick::FOREGROUND = 1.0f;
const float RandomWalkerSeedsBrick::BACKGROUND = 0.0f;

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

static VolumeAtomic<float> preprocessImageForRandomWalker(const VolumeAtomic<float>& img) {
    VolumeAtomic<float> output(img.getDimensions());
    const tgt::ivec3 start(0);
    const tgt::ivec3 end(img.getDimensions());
    const size_t numVoxels = tgt::hmul(img.getDimensions());

    const int k = 1;
    const int N=2*k+1;
    const tgt::ivec3 neighborhoodSize(k);

#ifdef VRN_OCTREEWALKER_MEAN_NOT_MEDIAN
    // mean
    auto conv = [&] (const VolumeAtomic<float>& input, VolumeAtomic<float>& output, int dim) {
        VRN_FOR_EACH_VOXEL(center, start, end) {
            tgt::ivec3 neigh(0);
            neigh[dim] = neighborhoodSize[dim];
            const tgt::ivec3 neighborhoodStart = tgt::max(start, center - neigh);
            const tgt::ivec3 neighborhoodEnd = tgt::min(end, center + neigh + tgt::ivec3(1));

            const int numNeighborhoodVoxels = tgt::hmul(neighborhoodEnd-neighborhoodStart);

            float sum=0.0f;
            VRN_FOR_EACH_VOXEL(pos, neighborhoodStart, neighborhoodEnd) {
                sum += input.voxel(pos);
            }
            float estimation = sum/numNeighborhoodVoxels;
            output.voxel(center) = estimation;
        }
    };
    VolumeAtomic<float> tmp(img.getDimensions());
    conv(img, output, 0);
    conv(output, tmp, 1);
    conv(tmp, output, 2);
#else
    // median
#if 0
    tgt::ivec3 last = end - tgt::ivec3(1);
    VRN_FOR_EACH_VOXEL(center, start, end) {
        const tgt::ivec3 neighborhoodStart = center - neighborhoodSize;
        const tgt::ivec3 neighborhoodEnd = center + neighborhoodSize + tgt::ivec3(1);

        std::array<float, N*N*N> vals;
        int i=0;
        VRN_FOR_EACH_VOXEL(pos, neighborhoodStart, neighborhoodEnd) {
            tgt::ivec3 p = tgt::clamp(pos, start, last);
            vals[i++] = img.voxel(p);
        }
        tgtAssert(i==N*N*N, "OI");
        int centerIndex = i/2;
        std::nth_element(vals.begin(), vals.begin()+centerIndex, vals.end());
        //std::sort(vals.begin(), vals.begin()+i);
        output.voxel(center) = vals[centerIndex];
    }
#else
    tgt::ivec3 last = end - tgt::ivec3(1);
    const size_t HEAP_SIZE = N*N*N/2+1;
    tgtAssert(HEAP_SIZE == 14, "Invalid neighborhood size");
    NSmallestHeap14 heap;
    VRN_FOR_EACH_VOXEL(center, start, end) {
        const tgt::ivec3 neighborhoodStart = center - neighborhoodSize;
        const tgt::ivec3 neighborhoodEnd = center + neighborhoodSize + tgt::ivec3(1);

        heap.clear();
        VRN_FOR_EACH_VOXEL(pos, neighborhoodStart, neighborhoodEnd) {
            tgt::ivec3 p = tgt::clamp(pos, start, last);
            float val = img.voxel(p);
            heap.push(val);
        }
        output.voxel(center) = heap.nthlargest();
    }
#endif
#endif

    float sumOfDifferences = 0.0f;
    VRN_FOR_EACH_VOXEL(center, start, end) {
        const tgt::ivec3 neighborhoodStart = tgt::max(start, center - neighborhoodSize);
        const tgt::ivec3 neighborhoodEnd = tgt::min(end, center + neighborhoodSize + tgt::ivec3(1));

        const int numNeighborhoodVoxels = tgt::hmul(neighborhoodEnd-neighborhoodStart);

        float estimation = output.voxel(center);
        float val = img.voxel(center);
        float diff = estimation - val;

        float neighborhoodFactor;
        if(numNeighborhoodVoxels > 1) {
            neighborhoodFactor = static_cast<float>(numNeighborhoodVoxels)/static_cast<float>(numNeighborhoodVoxels-1);
        } else {
            neighborhoodFactor = 1.0f;
        }

        sumOfDifferences += neighborhoodFactor * diff * diff;

        output.voxel(center) = estimation;
    }

#ifdef VRN_OCTREEWALKER_MEAN_NOT_MEDIAN
    const float varianceFactor = 2.0f/(N*N*N*N); //mean
#else
    tgtAssert(k==1, "Invalid k for variance factor");
    const float varianceFactor = 0.142; //median //TODO: this is for 2D. what about 3d?
#endif

    float rawVariance = sumOfDifferences/numVoxels;
    float varianceEstimation = rawVariance * varianceFactor;
    float stdEstimationInv = 1.0f/std::sqrt(varianceEstimation);

    VRN_FOR_EACH_VOXEL(center, start, end) {
        output.voxel(center) *= stdEstimationInv;
    }

    return output;
}
template<typename Accessor>
static void processVoxelWeights(const tgt::ivec3& voxel, const RandomWalkerSeedsBrick& seeds, EllpackMatrix<float>& mat, float* vec, size_t* volumeIndexToRowTable, Accessor& voxelFun, const tgt::svec3& volDim, float minWeight) {
    auto edgeWeight = [minWeight] (float voxelIntensity, float neighborIntensity) {
        float beta = 0.5f;
        float intDiff = (voxelIntensity - neighborIntensity);
        float intDiffSqr = intDiff*intDiff;
        float weight = exp(-beta * intDiffSqr);
        weight = std::max(weight, minWeight);

        return weight;
    };
    tgtAssert(volumeIndexToRowTable, "no volumeIndexToRowTable passed");
    tgtAssert(mat.isInitialized(), "matrix not initialized");

    size_t index = volumeCoordsToIndex(voxel, volDim);

    float curIntensity = voxelFun.voxel(voxel);

    float weightSum = 0;

    bool currentIsSeedpoint = seeds.isSeedPoint(index);

    for(int dim=0; dim<3; ++dim) {
        if(voxel[dim] > 0) {
            tgt::ivec3 neighbor = voxel;
            neighbor[dim] -= 1;

            size_t neighborIndex = volumeCoordsToIndex(neighbor, volDim);
            float neighborIntensity = voxelFun.voxel(neighbor);

            float weight = edgeWeight(curIntensity, neighborIntensity);

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
                    mat.setValue(curRow, nRow, -weight);
                    mat.setValue(nRow, curRow, -weight);
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
        mat.setValue(curRow, curRow, weightSum);
    }
}

static uint64_t processOctreeBrick(OctreeWalkerInput& input, OctreeWalkerNodeGeometry& outputNodeGeometry, Histogram1D& histogram, uint16_t& min, uint16_t& max, uint16_t& avg, bool& hasSeedConflicts, bool parentHadSeedsConflicts, OctreeBrickPoolManagerBase& outputPoolManager, OctreeWalkerNode* outputRoot, const OctreeWalkerNode inputRoot, boost::optional<OctreeWalkerNode> prevRoot, PointSegmentListGeometryVec3& foregroundSeeds, PointSegmentListGeometryVec3& backgroundSeeds, std::mutex& clMutex) {
    auto canSkipChildren = [&] (float min, float max) {
        float parentValueRange = max-min;
        const float delta = 0.01;
        bool minMaxSkip = max < 0.5-delta || min > 0.5+delta;
        return parentValueRange < input.homogeneityThreshold_ || minMaxSkip;
    };

    const OctreeBrickPoolManagerBase& inputPoolManager = *input.octree_.getBrickPoolManager();
    const tgt::svec3 brickDataSize = input.octree_.getBrickDim();
    const tgt::svec3 volumeDim = input.octree_.getDimensions();

    boost::optional<BrickNeighborhood> seedsNeighborhood = boost::none;

    bool stop = false;
    RandomWalkerSeedsBrick seeds = [&] () {
        if(outputRoot) {
            seedsNeighborhood = BrickNeighborhood::fromNode(outputNodeGeometry, outputNodeGeometry.level_+1, *outputRoot, brickDataSize, outputPoolManager);
            tgt::svec3 seedBufferDimensions = seedsNeighborhood->data_.getDimensions();
            tgt::mat4 voxelToSeedTransform = seedsNeighborhood->voxelToNeighborhood();

            if(canSkipChildren(seedsNeighborhood->min_, seedsNeighborhood->max_) && !parentHadSeedsConflicts) {
                //LINFOC(OctreeWalker::loggerCat_, "skip block early");
                stop = true;
                avg = normToBrick(seedsNeighborhood->avg_);
                min = normToBrick(seedsNeighborhood->min_);
                max = normToBrick(seedsNeighborhood->max_);
            }

            RandomWalkerSeedsBrick seeds(seedBufferDimensions, voxelToSeedTransform, foregroundSeeds, backgroundSeeds);
            seeds.addNeighborhoodBorderSeeds(*seedsNeighborhood, volumeDim);
            return seeds;
        } else {
            tgt::svec3 seedBufferDimensions = outputNodeGeometry.voxelDimensions() / outputNodeGeometry.scale();
            tgt::mat4 voxelToSeedTransform = tgt::mat4::createScale(tgt::vec3(1.0f / outputNodeGeometry.scale()));
            RandomWalkerSeedsBrick seeds(seedBufferDimensions, voxelToSeedTransform, foregroundSeeds, backgroundSeeds);
            return seeds;
        }
    }();
    hasSeedConflicts = seeds.numConflicts() != 0;
    if(stop) {
        return OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS;
    }

    tgt::svec3 walkerBlockDim = seeds.bufferDimensions();

    size_t numVoxels = tgt::hmul(walkerBlockDim);
    size_t numSeeds = seeds.getNumSeeds();
    size_t systemSize = numVoxels - numSeeds;

    BrickNeighborhood inputNeighborhood = BrickNeighborhood::fromNode(outputNodeGeometry, outputNodeGeometry.level_, inputRoot, brickDataSize, inputPoolManager);

    if(numSeeds == 0) {
        // No way to decide between foreground and background
        if(!hasSeedConflicts) {
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
            OctreeWalkerNodeBrick outputBrick(outputBrickAddr, brickDataSize, outputPoolManager);

            const tgt::svec3 brickStart = inputNeighborhood.centerBrickLlf_;
            const tgt::svec3 brickEnd = inputNeighborhood.centerBrickUrb_;
            tgt::svec3 centerBrickSize = brickEnd - brickStart;
            if(seedsNeighborhood) {
                auto& parentProbs = *seedsNeighborhood;
                float sum = 0.5f;
                VRN_FOR_EACH_VOXEL(pos, brickStart, brickEnd) {
                    uint16_t val = parentProbs.data_.voxel(pos);

                    outputBrick.data_.voxel(pos - brickStart) = val;

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
                    outputBrick.data_.voxel(pos) = avg;
                    histogram.addSample(avgf);
                }
            }
            return outputBrickAddr;
        }
    }

    auto volIndexToRow = seeds.generateVolumeToRowsTable();

    std::vector<float> initialization(systemSize, 0.0f);
    if(seedsNeighborhood) {
        auto& neighborhood = *seedsNeighborhood;
        VRN_FOR_EACH_VOXEL(pos, tgt::svec3(0), neighborhood.data_.getDimensions()) {
            size_t logicalIndex = volumeCoordsToIndex(pos, walkerBlockDim);
            if (!seeds.isSeedPoint(logicalIndex)) {
                initialization[volIndexToRow[logicalIndex]] = neighborhood.data_.voxel(pos);
            }
        }
    }

    auto solution = tgt::make_unique<float[]>(systemSize);
    std::fill_n(solution.get(), systemSize, 0.5f);

    EllpackMatrix<float> mat(systemSize, systemSize, 7);
    mat.initializeBuffers();

    float beta = 0.5f;
    float minWeight = 1.f / pow(10.f, static_cast<float>(input.minWeight_));

    RandomWalkerEdgeWeightIntensity edgeWeightFun(tgt::vec2(0.0f, 1.0f), beta, minWeight);

    auto rwInput = preprocessImageForRandomWalker(inputNeighborhood.data_);
    RandomWalkerVoxelAccessorBrick voxelAccessor(rwInput);

    auto vec = std::vector<float>(systemSize, 0.0f);

    VRN_FOR_EACH_VOXEL(pos, tgt::ivec3(0), tgt::ivec3(walkerBlockDim)) {
        processVoxelWeights(pos, seeds, mat, vec.data(), volIndexToRow.data(), voxelAccessor, walkerBlockDim, minWeight);
    }

    for(int i=0; i<10; ++i) {
        int iterations;
        {
            //std::lock_guard<std::mutex> guard(clMutex);
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
        OctreeWalkerNodeBrick outputBrick(outputBrickAddr, brickDataSize, outputPoolManager);

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

            outputBrick.data_.voxel(pos - brickStart) = val;

            min = std::min(val, min);
            max = std::max(val, max);
            sum += val;

            histogram.addSample(valf);
        }
        avg = sum/tgt::hmul(centerBrickSize);
    }

    if(canSkipChildren(brickToNorm(min), brickToNorm(max)) && !hasSeedConflicts) {
        outputPoolManager.deleteBrick(outputBrickAddr);
        return OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS;
    }

    return outputBrickAddr;
}

const std::string BRICK_BUFFER_SUBDIR =      "brickBuffer";
const std::string BRICK_BUFFER_FILE_PREFIX = "buffer_";

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

    std::string octreePath = "/home/dominik/nosnapshot/tmp/octreewalkertest/";

    std::string brickPoolPath = tgt::FileSystem::cleanupPath(octreePath + "/" + BRICK_BUFFER_SUBDIR);
    if (!tgt::FileSystem::dirExists(brickPoolPath)) {
        tgt::FileSystem::createDirectoryRecursive(brickPoolPath);
    }

    size_t brickSizeInBytes = brickSize * sizeof(uint16_t);
    if(!input.brickPoolManager_ || input.brickPoolManager_->getBrickMemorySizeInByte() != brickSizeInBytes) {
        if(input.brickPoolManager_) {
            input.brickPoolManager_->deinitialize();
        }
        input.brickPoolManager_.reset(new OctreeBrickPoolManagerMmap(brickPoolPath, BRICK_BUFFER_FILE_PREFIX));
        input.brickPoolManager_->initialize(brickSizeInBytes);
    }

    auto& brickPoolManager = *input.brickPoolManager_;

    boost::optional<OctreeWalkerNode> prevRoot = [&] () -> boost::optional<OctreeWalkerNode> {
        if(input.previousResult_) {
            auto& prev = *input.previousResult_;
            tgtAssert(volumeDim == prev.getDimensions(), "prev result: Dimension mismatch");
            tgtAssert(prev.getRootNode(), "prev result: No root");
            tgtAssert(input.octree_.getNumLevels() == prev.getNumLevels(), "prev result: num levels mismatch");
            return OctreeWalkerNode(prev.getRootNode(), maxLevel, tgt::svec3(0), prev.getDimensions());
        } else {
            return boost::none;
        }
    }();

    struct NodeToProcess {
        const VolumeOctreeNode* inputNode;
        VolumeOctreeNode** outputNodeSlot; // Never null
        tgt::svec3 llf;
        tgt::svec3 urb;
        bool parentHadSeedsConflicts;

        VolumeOctreeNode*& outputNode() {
            return *outputNodeSlot;
        }
    };

    // These contain nodes whose CHILDREN AND BRICKS are still referenced and thus must not be deleted.
    std::unordered_set<const VolumeOctreeNode*> nodesToSave;

    OctreeWalkerNode outputRootNode(nullptr, maxLevel, tgt::svec3(0), volumeDim);
    // If computation is canceled: Delete new nodes, but spare nodes that are part of another tree.
    tgt::ScopeGuard nodeCleanup([&] () {
        freeTreeComponents(outputRootNode.node_, nodesToSave, *brickPoolManager_);
    });

    std::vector<NodeToProcess> nodesToProcess;
    nodesToProcess.push_back(
        NodeToProcess {
            input.octree_.getRootNode(),
            &outputRootNode.node_,
            tgt::svec3::zero,
            volumeDim,
            false,
        }
    );

    const OctreeWalkerNode inputRoot(const_cast<VolumeOctreeNode*>(input.octree_.getRootNode()), input.octree_.getActualTreeDepth()-1, tgt::svec3(0), input.octree_.getDimensions());

    uint16_t globalMin = 0xffff;
    uint16_t globalMax = 0;

    Histogram1D histogram(0.0, 1.0, 256);

    auto rwm = input.volume_.getRealWorldMapping();

    PointSegmentListGeometryVec3 foregroundSeeds;
    PointSegmentListGeometryVec3 backgroundSeeds;
    getSeedListsFromPorts(input.foregroundGeomSeeds_, foregroundSeeds);
    getSeedListsFromPorts(input.backgroundGeomSeeds_, backgroundSeeds);

    std::mutex clMutex;
#ifdef VRN_OCTREEWALKER_USE_OMP
    LINFO("Using parallel octree walker variant.");
#else
    LINFO("Using sequential octree walker variant.");
#endif

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

        const int numNodes = nodesToProcess.size();
        ThreadedTaskProgressReporter parallelProgress(levelProgress, numNodes);
        bool aborted = false;

#ifdef VRN_OCTREEWALKER_USE_OMP
#pragma omp parallel for schedule(dynamic, 1) // Schedule: Process nodes/bricks locally to utilize brick cache
#endif
        for (int nodeId = 0; nodeId < numNodes; ++nodeId) {

#ifdef VRN_OCTREEWALKER_USE_OMP
            if(aborted) {
                continue;
            }
            if(parallelProgress.reportStepDone()) {
                #pragma omp critical
                aborted = true;
            }
#else
            if(parallelProgress.reportStepDone()) {
                aborted = true;
                break;
            }
#endif

            // Make sure to hit LRU cache: Go from back to front
            auto& node = nodesToProcess[numNodes-nodeId-1];

            tgtAssert(node.inputNode, "No input node");
            if(!node.inputNode->inVolume()) {
                node.outputNode() = new VolumeOctreeNodeGeneric<1>(OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS, false);
                continue;
            }

            uint16_t min = 0xffff;
            uint16_t max = 0;
            uint16_t avg = 0xffff/2;
            uint64_t newBrickAddr;
            bool hasSeedsConflicts;
            {
                tgtAssert(node.inputNode->hasBrick(), "No Brick");

                OctreeWalkerNodeGeometry outputNodeGeometry(level, node.llf, node.urb);
                newBrickAddr = processOctreeBrick(input, outputNodeGeometry, histogram, min, max, avg, hasSeedsConflicts, node.parentHadSeedsConflicts, brickPoolManager, level == maxLevel ? nullptr : &outputRootNode, inputRoot, prevRoot, foregroundSeeds, backgroundSeeds, clMutex);
            }

            globalMin = std::min(globalMin, min);
            globalMax = std::max(globalMax, max);

            bool childrenToProcess = newBrickAddr != OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS && !node.inputNode->isLeaf();
            VolumeOctreeNode* newNode = nullptr;

            // Check if previous result is sufficiently close to current one.
            // If this is the case: copy branch from old octree.
            if(childrenToProcess && prevRoot && !hasSeedsConflicts) {
                auto prevNode = prevRoot->findChildNode(node.llf, brickDim, level);

                if(prevNode.node().hasBrick()
                    && absdiffRel(prevNode.node().getMinValue(), min) < input.incrementalSimilarityThreshold_
                    && absdiffRel(prevNode.node().getMaxValue(), max) < input.incrementalSimilarityThreshold_
                    && absdiffRel(prevNode.node().getAvgValue(), avg) < input.incrementalSimilarityThreshold_
                ) {
                    uint16_t maxDiff = 0;
                    const tgt::svec3 begin(0);
                    const tgt::svec3 end(prevNode.geometry_.brickDimensions());
                    {
                        OctreeWalkerNodeBrickConst prevBrick(prevNode.node().getBrickAddress(), brickDim, brickPoolManager);
                        OctreeWalkerNodeBrickConst currentBrick(newBrickAddr, brickDim, brickPoolManager);

                        VRN_FOR_EACH_VOXEL(pos, begin, end) {
                            uint16_t current = currentBrick.data_.voxel(pos);
                            uint16_t prev = prevBrick.data_.voxel(pos);
                            uint16_t diff = absdiff(current, prev);
                            maxDiff = std::max(maxDiff, diff);
                        }
                    }
                    float reldiff = static_cast<float>(maxDiff) / 0xffff;
                    if(reldiff < input.incrementalSimilarityThreshold_) {
                        brickPoolManager.deleteBrick(newBrickAddr);

                        newNode = prevNode.node_;
                        childrenToProcess = false;

                        #pragma omp critical
                        {
                            nodesToSave.insert(newNode);
                        }
                        //LINFO("Using similar branch from previous iteration");
                    }
                }
            }

            if(!newNode) {
                auto newNodeGeneric = new VolumeOctreeNodeGeneric<1>(newBrickAddr, true);
                newNodeGeneric->avgValues_[0] = avg;
                newNodeGeneric->minValues_[0] = min;
                newNodeGeneric->maxValues_[0] = max;

                newNode = newNodeGeneric;
            }

            if(childrenToProcess) {

                tgt::svec3 childBrickSize = brickDim * (1UL << (level-1));
                for(auto child : OCTREEWALKER_CHILD_POSITIONS) {
                    const size_t childId = volumeCoordsToIndex(child, tgt::svec3::two);
                    VolumeOctreeNode* inputChildNode = node.inputNode->children_[childId];
                    tgtAssert(inputChildNode, "No child node");

                    tgt::svec3 start = node.llf + childBrickSize * child;
                    tgt::svec3 end = tgt::min(start + childBrickSize, volumeDim);
                    #pragma omp critical
                    {
                        nextNodesToProcess.push_back(
                                NodeToProcess {
                                inputChildNode,
                                &newNode->children_[childId],
                                start,
                                end,
                                hasSeedsConflicts,
                            }
                        );
                    }
                }
            }

            node.outputNode() = newNode;
        }
        if(aborted) {
            throw boost::thread_interrupted();
        }

        nodesToProcess = nextNodesToProcess;
    }

    nodeCleanup.dismiss();
    auto octree = new VolumeOctree(outputRootNode.node_, &brickPoolManager, brickDim, input.octree_.getDimensions(), numChannels);
    auto output = tgt::make_unique<Volume>(octree, &input.volume_);

    float min = static_cast<float>(globalMin)/0xffff;
    float max = static_cast<float>(globalMax)/0xffff;
    output->addDerivedData(new VolumeMinMax(min, max, min, max));
    output->addDerivedData(new VolumeHistogramIntensity(histogram));
    output->setRealWorldMapping(RealWorldMapping(1.0, 0.0f, "Probability"));
    auto finish = clock::now();
    return ComputeOutput {
        octree,
        std::move(output),
        nodesToSave,
        finish - start,
    };
}

void OctreeWalker::processComputeOutput(ComputeOutput output) {
    if (!output.volume_) {
        LERROR("Failed to compute Random Walker solution");
        return;
    }
    LINFO("Total runtime: " << output.duration_.count() << " sec");

    // Set new output
    outportProbabilities_.setData(output.volume_.get(), false);

    // previousOctree_ is now not referenced anymore, so we are free to clean up.
    if(previousVolume_) {
        tgtAssert(previousOctree_, "Previous result volume without octree");

        auto res = std::move(*previousOctree_).decompose();
        // Brickpoolmanager reference is not required here. The important thing is that the previous result does not deconstruct the brickPoolManager

        // Clean up old tree
        freeTreeComponents(res.second, output.previousNodesToSave, *brickPoolManager_);
    }

    previousOctree_ = output.octree_;
    previousVolume_ = std::move(output.volume_);
}
void OctreeWalker::clearPreviousResults() {
    // First: Reset output
    outportProbabilities_.setData(nullptr, false);

    // previousOctree_ is now not referenced anymore, so we are free to clean up.
    if(previousOctree_) {
        tgtAssert(previousVolume_, "Previous result octree without volume");

        auto res = std::move(*previousOctree_).decompose();
        // Brickpoolmanager reference is not required here. The important thing is that the previous result does not deconstruct the brickPoolManager

        // Clean up old tree
        freeNodes(res.second);
    }
    previousOctree_ = nullptr;
    previousVolume_.reset(nullptr);

    if(brickPoolManager_) {
        brickPoolManager_->deinitialize();
    }
    brickPoolManager_.reset(nullptr);
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


void OctreeWalker::updateGuiState() {
    //bool useTransFunc = enableTransFunc_.get();
    //edgeWeightTransFunc_.setVisibleFlag(useTransFunc);
    //edgeWeightBalance_.setVisibleFlag(useTransFunc);
}

}   // namespace
