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

#include "../solver/randomwalkersolver.h"
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
#include "voreen/core/datastructures/octree/octreebrickpoolmanagerdisk.h"
#include "voreen/core/datastructures/octree/volumeoctreenodegeneric.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "tgt/vector.h"
#include "tgt/memory.h"

#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"

#include <climits>

namespace voreen {

const std::string OctreeWalker::loggerCat_("voreen.RandomWalker.OctreeWalker");
using tgt::vec3;

OctreeWalker::OctreeWalker()
    : AsyncComputeProcessor<OctreeWalkerInput, OctreeWalkerOutput>(),
    inportVolume_(Port::INPORT, "volume.input"),
    inportForegroundSeeds_(Port::INPORT, "geometry.seedsForeground", "geometry.seedsForeground", true),
    inportBackgroundSeeds_(Port::INPORT, "geometry.seedsBackground", "geometry.seedsBackground", true),
    outportProbabilities_(Port::OUTPORT, "volume.probabilities", "volume.probabilities", false),
    usePrevProbAsInitialization_("usePrevProbAsInitialization", "Use Previous Probabilities as Initialization", false, Processor::VALID, Property::LOD_ADVANCED),
    beta_("beta", "Edge Weight Scale: 2^beta", 12, 0, 20),
    minEdgeWeight_("minEdgeWeight", "Min Edge Weight: 10^(-t)", 5, 0, 10),
    preconditioner_("preconditioner", "Preconditioner"),
    errorThreshold_("errorThreshold", "Error Threshold: 10^(-t)", 2, 0, 10),
    maxIterations_("conjGradIterations", "Max Iterations", 1000, 1, 5000),
    conjGradImplementation_("conjGradImplementation", "Implementation"),
    currentInputVolume_(0)
{
    // ports
    addPort(inportVolume_);
    addPort(inportForegroundSeeds_);
    addPort(inportBackgroundSeeds_);
    addPort(outportProbabilities_);

    addProperty(usePrevProbAsInitialization_);

    // random walker properties
    addProperty(beta_);
    addProperty(minEdgeWeight_);
    beta_.setGroupID("rwparam");
    minEdgeWeight_.setGroupID("rwparam");
    setPropertyGroupGuiName("rwparam", "Random Walker Parametrization");

    // conjugate gradient solver
    preconditioner_.addOption("none", "None");
    preconditioner_.addOption("jacobi", "Jacobi");
    preconditioner_.select("jacobi");
    addProperty(preconditioner_);
    addProperty(errorThreshold_);
    addProperty(maxIterations_);
    conjGradImplementation_.addOption("blasCPU", "CPU");
#ifdef VRN_MODULE_OPENMP
    conjGradImplementation_.addOption("blasMP", "OpenMP");
    conjGradImplementation_.select("blasMP");
#endif
#ifdef VRN_MODULE_OPENCL
    conjGradImplementation_.addOption("blasCL", "OpenCL");
    conjGradImplementation_.select("blasCL");
#endif
    addProperty(conjGradImplementation_);
    preconditioner_.setGroupID("conjGrad");
    errorThreshold_.setGroupID("conjGrad");
    maxIterations_.setGroupID("conjGrad");
    conjGradImplementation_.setGroupID("conjGrad");
    setPropertyGroupGuiName("conjGrad", "Conjugate Gradient Solver");
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


    std::vector<float> prevProbs;

    float errorThresh = 1.f / pow(10.f, static_cast<float>(errorThreshold_.get()));
    int maxIterations = maxIterations_.get();

    return ComputeInput {
        *vol,
        *octreePtr,
        inportForegroundSeeds_.getThreadSafeAllData(),
        inportBackgroundSeeds_.getThreadSafeAllData(),
        beta_.get(),
        minEdgeWeight_.get(),
        voreenBlas,
        precond,
        errorThresh,
        maxIterations,
    };
}

namespace {
    struct LoopRecord {
        int iteration;
        int level;
        float scaleFactor;
        tgt::ivec3 workDim;
        size_t numSeeds;
        size_t numForegroundSeeds;
        size_t numBackgroundSeeds;
        tgt::vec2 probabilityRange;
        int numIterations;
        std::chrono::duration<float> timeIteration;
        std::chrono::duration<float> timeSetup;
        std::chrono::duration<float> timeSolving;
        std::chrono::duration<float> timeSeedAnalysis;

        void print() {
            size_t numVoxels = tgt::hmul(workDim);
            std::string cat = "voreen.RandomWalker.OctreeWalker";
            LINFOC(cat, iteration << ". Iteration: level=" << level << ", scaleFactor=" << scaleFactor
                << ", dim=" << workDim);
            LINFOC(cat, "* num voxels: " << numVoxels << ", num seeds:  " << numSeeds << " (ratio: " << (float)numSeeds/tgt::hmul(workDim) << ")");
            LINFOC(cat, "* num unseeded: " << numVoxels-numSeeds);
            LINFOC(cat, "* probability range: " << probabilityRange);
            LINFOC(cat, "* runtime: " << timeIteration.count() << " sec");
            LINFOC(cat, "  - system setup: " << timeSetup.count() << " sec");
            LINFOC(cat, "  - solving: " << timeSolving.count() << " sec"
                << " (iterations: " << numIterations << ")");
            LINFOC(cat, "  - seed analysis: " << timeSeedAnalysis.count() << " sec");
        }
    };
} // namespace anonymous

static void getSeedListsFromPorts(std::vector<PortDataPointer<Geometry>>& geom, PointSegmentListGeometry<tgt::vec3>& seeds) {

    for (size_t i=0; i<geom.size(); i++) {
        const PointSegmentListGeometry<tgt::vec3>* seedList = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(geom.at(i).get());
        if (!seedList)
            LWARNINGC("voreen.RandomWalker.OctreeWalker", "Invalid geometry. PointSegmentListGeometry<vec3> expected.");
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

/*
struct OctreeBrickNeighbor {
    const OctreeBrick* brick_;
    tgt::mat4 toNeighborMat_;
    bool isPresent() {
        return brick_;
    }
    float readValue(tgt::ivec3 pos) {
        tgtAssert(isPresent(), "Brick not present");
        tgt::vec3 p = toNeighborMat_.transform(pos);
        return brick_->data_.getVoxelNormalized(p);
    }
};
*/

struct RandomWalkerVoxelAccessorBrick : public RandomWalkerVoxelAccessor {
    RandomWalkerVoxelAccessorBrick(const OctreeBrick& brick)
        : brick_(brick)
    {
    }
    virtual float voxel(const tgt::svec3& pos) {
        return brick_.data_.getVoxelNormalized(pos);
    }
private:
    const OctreeBrick& brick_;
};

namespace {

inline size_t volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

inline size_t volumeCoordsToIndex(const tgt::ivec3& coords, const tgt::ivec3& dim) {
    return coords.z*dim.y*dim.x + coords.y*dim.x + coords.x;
}
}

class RandomWalkerSeedsBrick : public RandomWalkerSeeds {
    const uint8_t UNLABELED = 0xff;
    const uint8_t FOREGROUND = 1;
    const uint8_t BACKGROUND = 0;
public:
    RandomWalkerSeedsBrick(const OctreeBrick& brick, const PointSegmentListGeometryVec3& foregroundSeedList, const PointSegmentListGeometryVec3& backgroundSeedList)
        : brick_(brick)
        , seedBuffer_(tgt::hmul(brick_.dim_), UNLABELED)
        , numForegroundSeeds_(0)
        , numBackgroundSeeds_(0)
    {
        tgt::ivec3 volDim = brick_.dim_;
        size_t numVoxels = tgt::hmul(brick_.dim_);

        tgt::mat4 voxelToBrick = brick_.voxelToBrick();
        // foreground geometry seeds
        auto collectLabels = [&] (const PointSegmentListGeometryVec3& seedList, uint8_t label, size_t& counter) {
            for (int m=0; m<seedList.getNumSegments(); m++) {
                const std::vector<tgt::vec3>& foregroundPoints = seedList.getData()[m];
                if (foregroundPoints.empty())
                    continue;
                for (size_t i=0; i<foregroundPoints.size()-1; i++) {
                    tgt::vec3 left = voxelToBrick*foregroundPoints[i];
                    tgt::vec3 right = voxelToBrick*foregroundPoints[i+1];
                    tgt::vec3 dir = tgt::normalize(right - left);
                    for (float t=0.f; t<tgt::length(right-left); t += 1.f) {
                        tgt::ivec3 point = tgt::iround(left + t*dir);
                        if(tgt::hor(tgt::lessThan(point, tgt::ivec3::zero)) || tgt::hor(tgt::greaterThanEqual(point, volDim))) {
                            continue;
                        }
                        size_t index = volumeCoordsToIndex(point, volDim);
                        tgtAssert(index < numVoxels, "Invalid index");
                        if (seedBuffer_[index] == UNLABELED) {
                            seedBuffer_[index] = label;
                            counter++;
                        }
                    }
                }
            }
        };
        collectLabels(foregroundSeedList, FOREGROUND, numForegroundSeeds_);
        collectLabels(backgroundSeedList, BACKGROUND, numBackgroundSeeds_);
        numSeeds_ = numForegroundSeeds_ + numBackgroundSeeds_;
    }
    virtual ~RandomWalkerSeedsBrick() {}
    virtual void initialize() {};

    virtual bool isSeedPoint(size_t index) const {
        return seedBuffer_[index];
    }
    virtual bool isSeedPoint(const tgt::ivec3& voxel) const {
        return isSeedPoint(volumeCoordsToIndex(voxel, brick_.dim_));
    }
    virtual float getSeedValue(size_t index) const {
        if (seedBuffer_[index] == FOREGROUND)
            return 1.f;
        else if (seedBuffer_[index] == BACKGROUND)
            return 0.f;
        else
            return -1.f;
    }
    virtual float getSeedValue(const tgt::ivec3& voxel) const {
        return getSeedValue(volumeCoordsToIndex(voxel, brick_.dim_));
    }

    std::vector<size_t> generateVolumeToRowsTable() {
        size_t numVoxels = tgt::hmul(brick_.dim_);
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

private:
    size_t numForegroundSeeds_;
    size_t numBackgroundSeeds_;
    const OctreeBrick& brick_;
    std::vector<uint8_t> seedBuffer_;
};


void OctreeWalker::processOctreeBrick(ComputeInput& input, OctreeBrick& brick, uint16_t* outputBrick, ProgressReporter& progressReporter, Histogram1D& histogram, uint16_t& min, uint16_t& max, uint16_t& avg) const {
    //TODO: catch out of memory

    tgt::svec3 brickDim = brick.dim_;

    //TODO: Only copy once
    PointSegmentListGeometryVec3 foregroundSeeds;
    PointSegmentListGeometryVec3 backgroundSeeds;
    getSeedListsFromPorts(input.foregroundGeomSeeds_, foregroundSeeds);
    getSeedListsFromPorts(input.backgroundGeomSeeds_, backgroundSeeds);

    RandomWalkerSeedsBrick seeds(brick, foregroundSeeds, backgroundSeeds);

    float* oldSystemSolution = nullptr;
    size_t numVoxels = tgt::hmul(brickDim);
    size_t numSeeds = seeds.getNumSeeds(); //TODO revisit when adding additional border seeds
    size_t systemSize = numVoxels - numSeeds;

    // No way to decide between foreground and background
    if(numSeeds == 0) {
        min = 0xffff/2;
        max = 0xffff/2;
        avg = 0xffff/2;
        for(int i=0; i<numVoxels; ++i) {
            histogram.addSample(0.5f);
        }
        std::fill_n(outputBrick, numVoxels, 0xffff/2);
        return;
    }

    auto solution = tgt::make_unique<float[]>(systemSize);

    EllpackMatrix<float> mat;
    mat.setDimensions(systemSize, systemSize, 7);
    mat.initializeBuffers();

    auto vmm = input.volume_.getDerivedData<VolumeMinMax>();
    tgt::vec2 intensityRange(vmm->getMin(), vmm->getMax());
    std::unique_ptr<RandomWalkerEdgeWeight> edgeWeightFun;

    float beta = static_cast<float>(1<<input.beta_);
    float minWeight = 1.f / pow(10.f, static_cast<float>(input.minWeight_));

    edgeWeightFun.reset(new RandomWalkerEdgeWeightIntensity(intensityRange, beta, minWeight));

    std::unique_ptr<RandomWalkerVoxelAccessor> voxelAccessor(new RandomWalkerVoxelAccessorBrick(brick));

    RandomWalkerWeights edgeWeights(std::move(voxelAccessor), std::move(edgeWeightFun), brickDim);

    auto vec = std::vector<float>(systemSize, 0.0f);

    auto volIndexToRow = seeds.generateVolumeToRowsTable();

    for (int z=0; z<brickDim.z; z++) {
        for (int y=0; y<brickDim.y; y++) {
            for (int x=0; x<brickDim.x; x++) {
                edgeWeights.processVoxel(tgt::ivec3(x, y, z), &seeds, mat, vec.data(), volIndexToRow.data());
            }
        }
    }

    int iterations = input.blas_->sSpConjGradEll(mat, vec.data(), solution.get(), oldSystemSolution,
        input.precond_, input.errorThreshold_, input.maxIterations_, progressReporter);

    tgt::svec3 physicalBrickDim = brick.data_.getDimensions();
    size_t physicalBrickSize = tgt::hmul(physicalBrickDim);
    std::fill_n(outputBrick, physicalBrickSize, 0.5f);

    uint64_t sum = 0;

    for (int z=0; z<brickDim.z; z++) {
        for (int y=0; y<brickDim.y; y++) {
            for (int x=0; x<brickDim.x; x++) {
                size_t physicalIndex = volumeCoordsToIndex(x, y, z, physicalBrickDim);
                size_t logicalIndex = volumeCoordsToIndex(x, y, z, brickDim);
                float valf;
                if (seeds.isSeedPoint(logicalIndex)) {
                    valf = seeds.getSeedValue(logicalIndex);
                } else {
                    valf = solution[volIndexToRow[logicalIndex]];
                }
                uint16_t val = valf*0xffff;

                outputBrick[physicalIndex] = val;

                min = std::min(val, min);
                max = std::max(val, max);
                sum += val;

                histogram.addSample(valf);
            }
        }
    }
    avg = sum/numVoxels;
}

const std::string BRICK_BUFFER_SUBDIR =      "brickBuffer";
const std::string BRICK_BUFFER_FILE_PREFIX = "buffer_";

OctreeWalker::ComputeOutput OctreeWalker::compute(ComputeInput input, ProgressReporter& progressReporter) const {
    OctreeWalkerOutput invalidResult = OctreeWalkerOutput {
        std::unique_ptr<Volume>(nullptr),
        std::chrono::duration<float>(0),
    };

    progressReporter.setProgress(0.0);

    auto start = clock::now();

    const tgt::svec3 volumeDim = input.octree_.getDimensions();
    const tgt::svec3 brickDim = input.octree_.getBrickDim();
    const size_t brickSize = tgt::hmul(brickDim);
    const size_t numChannels = 1;
    const size_t maxLevel = input.octree_.getNumLevels()-1;

    const float homogeneityThreshold = 0.1; //TODO property

    std::string octreePath = "/home/dominik/nosnapshot/tmp/octreewalkertest/";

    std::string brickPoolPath = tgt::FileSystem::cleanupPath(octreePath + "/" + BRICK_BUFFER_SUBDIR);
    if (!tgt::FileSystem::dirExists(brickPoolPath))
        tgt::FileSystem::createDirectoryRecursive(brickPoolPath);

    size_t brickSizeInBytes = brickSize * sizeof(uint16_t);
    OctreeBrickPoolManagerBase* brickPoolManager = new OctreeBrickPoolManagerDisk(brickSizeInBytes,
            VoreenApplication::app()->getCpuRamLimit(), brickPoolPath, BRICK_BUFFER_FILE_PREFIX);
    brickPoolManager->initialize(brickSizeInBytes);

    struct NodeToProcess {
        const VolumeOctreeNode* inputNode;
        VolumeOctreeNode* outputNode;
        tgt::svec3 llf;
        tgt::svec3 urb;
    };

    std::vector<NodeToProcess> nodesToProcess;
    VolumeOctreeNode* newRootNode = new VolumeOctreeNodeGeneric<1>(brickPoolManager->allocateBrick(), true);
    nodesToProcess.push_back(
        NodeToProcess {
            input.octree_.getRootNode(),
            newRootNode,
            tgt::svec3::zero,
            volumeDim
        }
    );

    uint16_t globalMin = 0xffff;
    uint16_t globalMax = 0;

    Histogram1D histogram(0.0, 1.0, 256);

    // Level order iteration => Previos level is always available
    for(int level = maxLevel; level >=0; --level) {
        float progressBegin = 1.0/(1 << (3 * (level + 1))); // 1/8 ^ (level+1 => next level)
        float progressEnd = 1.0/(1 << (3 * (level))); // 1/8 ^ (level)
        SubtaskProgressReporter levelProgress(progressReporter, tgt::vec2(progressBegin, progressEnd));

        std::vector<NodeToProcess> nextNodesToProcess;
        int i=0;
        float perNodeProgress = 1.0f/nodesToProcess.size();
        for(auto& node : nodesToProcess) {
            SubtaskProgressReporter progress(levelProgress, tgt::vec2(i*perNodeProgress, (i+1)*perNodeProgress));
            tgtAssert(node.inputNode, "No input node");
            const uint16_t* brickData = input.octree_.getNodeBrick(node.inputNode);

            tgtAssert(node.inputNode->hasBrick(), "No Brick");

            OctreeBrick brick(brickData, brickDim, node.llf, node.urb, level);

            auto newBrick = brickPoolManager->getWritableBrick(node.outputNode->getBrickAddress());

            uint16_t min = 0xffff;
            uint16_t max = 0;
            uint16_t avg = 0xffff/2;
            processOctreeBrick(input, brick, newBrick, progress, histogram, min, max, avg);

            globalMin = std::min(globalMin, min);
            globalMax = std::max(globalMax, max);

            auto genericNode = dynamic_cast<VolumeOctreeNodeGeneric<1>*>(node.outputNode);
            tgtAssert(genericNode, "Failed dynamic_cast");
            genericNode->avgValues_[0] = avg;
            genericNode->minValues_[0] = min;
            genericNode->maxValues_[0] = max;

            input.octree_.releaseNodeBrick(node.inputNode);
            brickPoolManager->releaseBrick(node.outputNode->getBrickAddress(), OctreeBrickPoolManagerBase::WRITE);

            float range = static_cast<float>(max-min)/0xffff;

            if(max-min < homogeneityThreshold) {
                brickPoolManager->deleteBrick(node.outputNode->getBrickAddress());
                node.outputNode->setBrickAddress(NO_BRICK_ADDRESS);
            } else if(!node.inputNode->isLeaf() /* TODO handle early leaf */) {
                tgt::svec3 childBrickSize = brickDim * (1UL << (level-1));
                VRN_FOR_EACH_VOXEL(child, tgt::svec3::zero, tgt::svec3::two) {
                    const size_t childId = volumeCoordsToIndex(child, tgt::svec3::two);
                    VolumeOctreeNode* inputChildNode = node.inputNode->children_[childId];
                    tgtAssert(inputChildNode, "No child node");

                    VolumeOctreeNode* outputChildNode;
                    if(inputChildNode->inVolume()) {
                        outputChildNode = new VolumeOctreeNodeGeneric<1>(brickPoolManager->allocateBrick(), true);

                        tgt::svec3 start = node.llf + childBrickSize * child;
                        tgt::svec3 end = tgt::min(start + childBrickSize, volumeDim);
                        nextNodesToProcess.push_back(
                                NodeToProcess {
                                inputChildNode,
                                outputChildNode,
                                start,
                                end
                            }
                        );
                    } else {
                        outputChildNode = new VolumeOctreeNodeGeneric<1>(NO_BRICK_ADDRESS, false);
                    }
                    node.outputNode->children_[childId] = outputChildNode;
                }
            }
            progress.setProgress(1.0f);
            ++i;
        }
        nodesToProcess = nextNodesToProcess;
    }


    auto output = tgt::make_unique<Volume>(new VolumeOctree(newRootNode, brickPoolManager, brickDim, input.octree_.getDimensions(), numChannels), &input.volume_);

    float min = static_cast<float>(globalMin)/0xffff;
    float max = static_cast<float>(globalMax)/0xffff;
    output->addDerivedData(new VolumeMinMax(min, max, min, max));
    output->addDerivedData(new VolumeHistogramIntensity(histogram));
    output->setRealWorldMapping(RealWorldMapping(1.0, 0.0f, "Probability"));
    auto finish = clock::now();
    return ComputeOutput {
        std::move(output), finish - start,
    };
}

void OctreeWalker::processComputeOutput(ComputeOutput output) {
    if (output.volume_) {
        LINFO("Total runtime: " << output.duration_.count() << " sec");
    } else {
        LERROR("Failed to compute Random Walker solution");
    }
    outportProbabilities_.setData(output.volume_.release());
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
