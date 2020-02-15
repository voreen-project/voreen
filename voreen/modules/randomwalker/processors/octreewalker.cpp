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

#include "../solver/randomwalkerweights.h"
#include "../util/preprocessing.h"

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
#include "voreen/core/utils/hashing.h"
#include "tgt/vector.h"
#include "tgt/memory.h"

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include <climits>

#ifdef VRN_MODULE_OPENMP
#include <omp.h>
#endif

namespace voreen {

#if defined(VRN_MODULE_OPENMP) && 1
#define VRN_OCTREEWALKER_USE_OMP
#endif

namespace {

const std::string RESULT_OCTREE_FILE_NAME = "prev_result_octree.vvod";
const std::string BRICK_BUFFER_SUBDIR =      "brickBuffer";
const std::string BRICK_BUFFER_FILE_PREFIX = "buffer_";

inline size_t volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

inline size_t volumeCoordsToIndex(const tgt::ivec3& coords, const tgt::ivec3& dim) {
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
    return static_cast<float>(val)/0xffff;
}

}

OctreeWalkerOutput::OctreeWalkerOutput(
    VolumeOctree* octree,
    std::unique_ptr<VolumeBase>&& volume,
    std::unordered_set<const VolumeOctreeNode*>&& sharedNodes,
    std::chrono::duration<float> duration
)   : octree_(octree)
    , volume_(std::move(volume))
    , sharedNodes_(std::move(sharedNodes))
    , duration_(duration)
{ }

OctreeWalkerOutput::OctreeWalkerOutput(OctreeWalkerOutput&& other)
    : octree_(other.octree_)
    , volume_(std::move(other.volume_))
    , sharedNodes_(std::move(other.sharedNodes_))
    , duration_(other.duration_)
{
    other.octree_ = nullptr;
    other.volume_ = nullptr;
}

OctreeWalkerOutput::~OctreeWalkerOutput() {
    // If OctreeWalkerOutput is dropped prematurely (e.g., if the computation
    // is interrupted after the final result was created, we need to make sure
    // NOT to drop the Brickpoolmanager, BUT free all nodes of the (now
    // discarded) tree.
    if(volume_) {
        tgtAssert(octree_, "Previous result octree without volume");

        auto res = std::move(*octree_).decompose();
        // Brickpoolmanager reference is not required here. The important thing is that the previous result does not deconstruct the brickPoolManager

        // Clean up old tree
        freeTreeComponents(res.second, sharedNodes_, *res.first);
    }
}


const std::string OctreeWalker::loggerCat_("voreen.RandomWalker.OctreeWalker");

OctreeWalker::OctreeWalker()
    : AsyncComputeProcessor<OctreeWalkerInput, OctreeWalkerOutput>()
    , inportVolume_(Port::INPORT, "volume.input")
    , inportForegroundSeeds_(Port::INPORT, "geometry.seedsForeground", "geometry.seedsForeground", true)
    , inportBackgroundSeeds_(Port::INPORT, "geometry.seedsBackground", "geometry.seedsBackground", true)
    , outportProbabilities_(Port::OUTPORT, "volume.probabilities", "volume.probabilities", false)
    , minEdgeWeight_("minEdgeWeight", "Min Edge Weight: 10^(-t)", 5, 0, 10)
    , preconditioner_("preconditioner", "Preconditioner")
    , errorThreshold_("errorThreshold", "Error Threshold: 10^(-t)", 2, 0, 10)
    , maxIterations_("conjGradIterations", "Max Iterations", 1000, 1, 5000)
    , conjGradImplementation_("conjGradImplementation", "Implementation")
    , homogeneityThreshold_("homogeneityThreshold", "Homogeneity Threshold", 0.01, 0.0, 1.0)
    , incrementalSimilarityThreshold_("incrementalSimilarityThreshold", "Incremental Similarity Treshold", 0.01, 0.0, 1.0)
    , resultPath_("resultPath", "Result Cache Path", "Result Cache Path", "", "", FileDialogProperty::DIRECTORY)
    , prevResultPath_("")
    , previousOctree_(nullptr)
    , previousVolume_(nullptr)
    , brickPoolManager_(nullptr)
{
    // ports
    addPort(inportVolume_);
    addPort(inportForegroundSeeds_);
    addPort(inportBackgroundSeeds_);
    addPort(outportProbabilities_);

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

    // conjugate gradient solver
    addProperty(resultPath_);
        ON_CHANGE_LAMBDA(resultPath_, [this] () {
                if(resultPath_.get() != prevResultPath_) {
                    interruptComputation();
                    clearPreviousResults();
                }
                });
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

    std::string previousResultFile = tgt::FileSystem::cleanupPath(prevResultPath_ + "/" + RESULT_OCTREE_FILE_NAME);

    if(!prevResultPath_.empty() && tgt::FileSystem::fileExists(previousResultFile)) {
        XmlDeserializer d(prevResultPath_);

        previousOctree_ = new VolumeOctree();
        try {
            std::ifstream fs(previousResultFile);
            d.read(fs);
            d.deserialize("Octree", *previousOctree_);

            resultPath_.set(prevResultPath_);
        } catch (std::exception& e) {
            LERROR("Failed to deserialize previous solution: " << e.what());
            delete previousOctree_;
            previousOctree_ = nullptr;
        } catch (SerializationException& e) {
            d.removeLastError();
            LERROR("Failed to deserialize previous solution: " << e.what());
            delete previousOctree_;
            previousOctree_ = nullptr;
        }
    }

    if(previousOctree_) {
        tgt::vec3 spacing(1.0f);
        tgt::vec3 offset(1.0f);

        OctreeBrickPoolManagerMmap* obpmm = dynamic_cast<OctreeBrickPoolManagerMmap*>(previousOctree_->getBrickPoolManager());
        if(obpmm) {
            previousVolume_.reset(new Volume(previousOctree_, spacing, offset));
            brickPoolManager_.reset(obpmm);
        } else {
            delete previousOctree_;
            previousOctree_ = nullptr;
        }
    }
}


OctreeWalker::ComputeInput OctreeWalker::prepareComputeInput() {
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

    // Check if the previous result is compatible with the current input
    if(previousOctree_ && (
              previousOctree_->getDimensions() != volumeDim
           || previousOctree_->getBrickDim() != brickDim
           || previousOctree_->getNumLevels() != octree.getNumLevels()
           || brickPoolManager_->getBrickMemorySizeInByte() != brickSizeInBytes)) {

        // If not, clear previous results
        clearPreviousResults();
    }

    if(previousOctree_) {
        LINFO("Reusing previous solution for compatible input volume");
    } else {
        LINFO("Not reusing previous solution for incompatible input volume");
    }

    // Create a new brickpool if we need a new one
    if(!brickPoolManager_) {
        brickPoolManager_.reset(new OctreeBrickPoolManagerMmap(brickPoolPath, BRICK_BUFFER_FILE_PREFIX));
        brickPoolManager_->initialize(brickSizeInBytes);
    }

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

    static BrickNeighborhood fromNode(const VolumeOctreeNodeLocation& current, size_t sampleLevel, const LocatedVolumeOctreeNodeConst& root, const tgt::svec3& brickBaseSize, const OctreeBrickPoolManagerBase& brickPoolManager) {
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
            tgt::svec3 samplePoint = tgt::round(bufferToVoxel.transform(blockBegin));
            auto node = root.findChildNode(samplePoint, brickBaseSize, sampleLevel);
            if(node.node().hasBrick()) {
                BrickPoolBrickConst brick(node.node().getBrickAddress(), brickBaseSize, brickPoolManager);

                tgt::mat4 bufferToSampleBrick = node.location().voxelToBrick() * bufferToVoxel;
                VRN_FOR_EACH_VOXEL(point, blockBegin, blockEnd) {
                    tgt::svec3 samplePos = tgt::round(bufferToSampleBrick.transform(point));
                    samplePos = tgt::clamp(samplePos, tgt::svec3(0), node.location().brickDimensions() - tgt::svec3(1));
                    float val = brick.getVoxelNormalized(samplePos);

                    output.setVoxelNormalized(val, point);
                    min = std::min(val, min);
                    max = std::max(val, max);
                    sum += val;
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
                        minSeed_ = std::min(label, minSeed_);
                        maxSeed_ = std::max(label, maxSeed_);
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
                if (seedVal == UNLABELED) {
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
    float minSeed() {
        return minSeed_;
    }
    float maxSeed() {
        return maxSeed_;
    }

private:
    VolumeAtomic<float> seedBuffer_;
    uint64_t numConflicts_;
    float minSeed_;
    float maxSeed_;
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

static uint64_t processOctreeBrick(OctreeWalkerInput& input, VolumeOctreeNodeLocation& outputNodeGeometry, Histogram1D& histogram, uint16_t& min, uint16_t& max, uint16_t& avg, bool& hasSeedConflicts, bool parentHadSeedsConflicts, OctreeBrickPoolManagerBase& outputPoolManager, LocatedVolumeOctreeNode* outputRoot, const LocatedVolumeOctreeNodeConst& inputRoot, boost::optional<LocatedVolumeOctreeNode> prevRoot, PointSegmentListGeometryVec3& foregroundSeeds, PointSegmentListGeometryVec3& backgroundSeeds, std::mutex& clMutex) {
    auto canSkipChildren = [&] (float min, float max) {
        float parentValueRange = max-min;
        const float delta = 0.01;
        bool minMaxSkip = max < 0.5-delta || min > 0.5+delta;
        return parentValueRange < input.homogeneityThreshold_ || minMaxSkip;
    };

    const OctreeBrickPoolManagerBase& inputPoolManager = *input.octree_.getBrickPoolManager();
    const tgt::svec3 brickDataSize = input.octree_.getBrickDim();

    boost::optional<BrickNeighborhood> seedsNeighborhood = boost::none;

    bool stop = false;
    RandomWalkerSeedsBrick seeds = [&] () {
        if(outputRoot) {
            seedsNeighborhood = BrickNeighborhood::fromNode(outputNodeGeometry, outputNodeGeometry.level_+1, *outputRoot, brickDataSize, outputPoolManager);
            tgt::svec3 seedBufferDimensions = seedsNeighborhood->data_.getDimensions();
            tgt::mat4 voxelToSeedTransform = seedsNeighborhood->voxelToNeighborhood();

            RandomWalkerSeedsBrick seeds(seedBufferDimensions, voxelToSeedTransform, foregroundSeeds, backgroundSeeds);
            seeds.addNeighborhoodBorderSeeds(*seedsNeighborhood);

            if(!parentHadSeedsConflicts && canSkipChildren(std::min(seedsNeighborhood->min_, seeds.minSeed()), std::max(seeds.maxSeed(), seedsNeighborhood->max_))) {
                //LINFOC(OctreeWalker::loggerCat_, "skip block early");
                stop = true;
                avg = normToBrick(seedsNeighborhood->avg_);
                min = normToBrick(seedsNeighborhood->min_);
                max = normToBrick(seedsNeighborhood->max_);
            }

            return seeds;
        } else {
            tgt::mat4 voxelToSeedTransform = outputNodeGeometry.voxelToBrick();
            // Note: For the reason to use `ceil` here see BrickNeighborhood::fromNode above.
            tgt::svec3 seedBufferDimensions = tgt::ceil(voxelToSeedTransform.getRotationalPart().transform(outputNodeGeometry.voxelDimensions()));
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
            BrickPoolBrick outputBrick(outputBrickAddr, brickDataSize, outputPoolManager);

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

    std::vector<float> initialization(systemSize, 0.5f);
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

    auto rwInput = preprocessForAdaptiveParameterSetting(inputNeighborhood.data_);
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
        BrickPoolBrick outputBrick(outputBrickAddr, brickDataSize, outputPoolManager);

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

    if(canSkipChildren(brickToNorm(min), brickToNorm(max)) && !hasSeedConflicts) {
        outputPoolManager.deleteBrick(outputBrickAddr);
        return OctreeBrickPoolManagerBase::NO_BRICK_ADDRESS;
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

    auto& brickPoolManager = *input.brickPoolManager_;

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

    LocatedVolumeOctreeNode outputRootNode(nullptr, maxLevel, tgt::svec3(0), volumeDim);
    // If computation is canceled: Delete new nodes, but spare nodes that are part of another tree.
    tgt::ScopeGuard nodeCleanup([&] () {
        freeTreeComponents(&outputRootNode.node(), nodesToSave, *brickPoolManager_);
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

    LocatedVolumeOctreeNodeConst inputRoot = input.octree_.getLocatedRootNode();

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

                VolumeOctreeNodeLocation outputNodeGeometry(level, node.llf, node.urb);
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
                        BrickPoolBrickConst prevBrick(prevNode.node().getBrickAddress(), brickDim, brickPoolManager);
                        BrickPoolBrickConst currentBrick(newBrickAddr, brickDim, brickPoolManager);

                        VRN_FOR_EACH_VOXEL(pos, begin, end) {
                            uint16_t current = currentBrick.data().voxel(pos);
                            uint16_t prev = prevBrick.data().voxel(pos);
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

                tgt::svec3 childBrickSize = brickDim * (size_t(1) << (level-1));
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
    std::vector<Histogram1D*> octreeHistograms;
    octreeHistograms.push_back(new Histogram1D(histogram));
    auto octree = new VolumeOctree(outputRootNode.node_, &brickPoolManager, std::move(octreeHistograms), brickDim, input.octree_.getDimensions(), numChannels);
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
        octree,
        std::move(output),
        std::move(nodesToSave),
        finish - start
    );
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
        freeTreeComponents(res.second, output.sharedNodes_, *brickPoolManager_);
    }

    previousOctree_ = output.octree_;
    previousVolume_ = std::move(output.volume_);

    output.octree_ = nullptr;

    if(previousOctree_) {
        std::string previousResultDir = resultPath_.get();
        std::string previousResultFile = tgt::FileSystem::cleanupPath(previousResultDir + "/" + RESULT_OCTREE_FILE_NAME);
        previousVolume_->setOrigin(previousResultFile);

        XmlSerializer s(previousResultDir);
        std::ofstream fs(previousResultFile);
        s.serialize("Octree", *previousOctree_);
        s.write(fs);

        tgtAssert(previousVolume_, "Previous result octree without volume");
        previousVolume_->setOrigin(previousResultFile);
    }
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

}   // namespace
