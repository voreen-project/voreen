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

#include "voreen/core/datastructures/octree/volumeoctree.h"

#include "voreen/core/datastructures/octree/octreeutils.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/utils/memoryinfo.h"

#include "voreen/core/utils/voreenfilepathhelper.h"

#include "tgt/assert.h"
#include "tgt/logmanager.h"
#include "tgt/tgt_math.h"
#include "tgt/stopwatch.h"
#include "tgt/filesystem.h"

#include <sstream>
#include <queue>

#ifdef VRN_MODULE_OPENMP
#include "omp.h"
#endif

using tgt::svec3;
using tgt::vec3;

namespace {
    const size_t MAX_CHANNELS = 4; //< maximum number of channels that are supported
    const std::string NODE_BUFFER_FILE_NAME = "nodebuffer.raw";
    const size_t NUM_HISTOGRAM_BUCKETS_BITS = 12;
    const size_t NUM_HISTOGRAM_BUCKETS = 1 << NUM_HISTOGRAM_BUCKETS_BITS; //< This may buckets will be used during octree construction

    static void deleteSubTree(voreen::VolumeOctreeNode* root) {
        if (!root)
            return;

        deleteSubTree(root->children_[0]);
        root->children_[0] = 0;

        deleteSubTree(root->children_[1]);
        root->children_[1] = 0;

        deleteSubTree(root->children_[2]);
        root->children_[2] = 0;

        deleteSubTree(root->children_[3]);
        root->children_[3] = 0;

        deleteSubTree(root->children_[4]);
        root->children_[4] = 0;

        deleteSubTree(root->children_[5]);
        root->children_[5] = 0;

        deleteSubTree(root->children_[6]);
        root->children_[6] = 0;

        deleteSubTree(root->children_[7]);
        root->children_[7] = 0;

        delete root;
    }


    /**
     * Helper class representing a 3D grid of nodes. Only used during iterative construction.
     */
    class NodeGrid3D {
    public:
        NodeGrid3D(tgt::svec3 gridDimensions)
            : dim_(gridDimensions)
            , numNodes_(tgt::hmul(dim_))
            , grid_(numNodes_)
        {
            tgtAssert(tgt::hand(tgt::greaterThan(dim_, tgt::svec3::zero)), "invalid grid dim");
        }
        ~NodeGrid3D() {
            for(auto& node : grid_) {
                deleteSubTree(node);
            }
        }

        voreen::VolumeOctreeNode* getNode(const tgt::svec3& pos) const {
            return grid_[cubicCoordToLinear(pos, dim_)];
        }
        voreen::VolumeOctreeNode* takeNode(const tgt::svec3& pos) {
            auto& current = grid_[cubicCoordToLinear(pos, dim_)];
            auto ret = current;
            current = nullptr;
            return ret;
        }
        void setNode(voreen::VolumeOctreeNode* node, const tgt::svec3 pos) {
            auto& current = grid_[cubicCoordToLinear(pos, dim_)];
            tgtAssert(!current, "Overwriting node in NodeGrid3D");
            current = node;
        }

        tgt::svec3 getDim() const {
            return dim_;
        }
        bool isComplete() const {
            for (size_t i=0; i<numNodes_; i++) {
                if (!grid_[i])
                    return false;
            }
            return true;
        }
        bool isEmpty() const {
            for (size_t i=0; i<numNodes_; i++) {
                if (grid_[i])
                    return false;
            }
            return true;
        }

    private:
        tgt::svec3 dim_;
        size_t numNodes_;

        // It does not really make sense to use unique_ptr here, since
        // VolumeOctreeNode (in most cases?) in fact owns its children (in a
        // way), but does not delete them on destruction...
        //
        // However, we own the nodes (including their children) here, so we
        // have to clean them up in our destructor.
        std::vector<voreen::VolumeOctreeNode*> grid_;
    };

    // voxel value conversion templates
    template<class T>
    inline uint16_t convertVoxelValueToUInt16(T value) {
        float normValue;
        if (voreen::VolumeElement<T>::isInteger()) {
            normValue = (static_cast<float>(value) - static_cast<float>(voreen::VolumeElement<T>::rangeMin())) /
                static_cast<float>(voreen::VolumeElement<T>::rangeMax()-voreen::VolumeElement<T>::rangeMin());
        }
        else {
            normValue = static_cast<float>(value);
        }

        return static_cast<uint16_t>(normValue*65535.f);
    }
    template<>
    inline uint16_t convertVoxelValueToUInt16(uint8_t value) {
        return static_cast<uint16_t>(value << 8);
    }
    template<>
    inline uint16_t convertVoxelValueToUInt16(uint16_t value) {
        return value;
    }
    template<>
    inline uint16_t convertVoxelValueToUInt16(uint32_t value) {
        return static_cast<uint16_t>(value >> 16);
    }
    template<>
    inline uint16_t convertVoxelValueToUInt16(float value) {
        return static_cast<uint16_t>(value*65535.f);
    }


    struct HalfSampleMean {
        static inline uint16_t halfsample(const uint16_t& lll, const uint16_t& llh, const uint16_t& lhl, const uint16_t& lhh, const uint16_t& hll, const uint16_t& hlh, const uint16_t& hhl, const uint16_t& hhh) {
            uint64_t halfValue =
                (uint64_t) lll +
                (uint64_t) llh +
                (uint64_t) lhl +
                (uint64_t) lhh +
                (uint64_t) hll +
                (uint64_t) hlh +
                (uint64_t) hhl +
                (uint64_t) hhh ;
            halfValue /= 8;
            return halfValue;
        }
    };

    struct HalfSampleMin {
        static inline uint16_t halfsample(const uint16_t& lll, const uint16_t& llh, const uint16_t& lhl, const uint16_t& lhh, const uint16_t& hll, const uint16_t& hlh, const uint16_t& hhl, const uint16_t& hhh) {
            uint16_t llx = std::min(lll, llh);
            uint16_t lhx = std::min(lhl, lhh);
            uint16_t hlx = std::min(hll, hlh);
            uint16_t hhx = std::min(hhl, hhh);
            uint16_t lxx = std::min(lhx, llx);
            uint16_t hxx = std::min(hhx, hlx);
            return std::min(hxx, lxx);
        }
    };

    struct HalfSampleMax {
        static inline uint16_t halfsample(const uint16_t& lll, const uint16_t& llh, const uint16_t& lhl, const uint16_t& lhh, const uint16_t& hll, const uint16_t& hlh, const uint16_t& hhl, const uint16_t& hhh) {
            uint16_t llx = std::max(lll, llh);
            uint16_t lhx = std::max(lhl, lhh);
            uint16_t hlx = std::max(hll, hlh);
            uint16_t hhx = std::max(hhl, hhh);
            uint16_t lxx = std::max(lhx, llx);
            uint16_t hxx = std::max(hhx, hlx);
            return std::max(hxx, lxx);
        }
    };
    struct HalfSampleMedian {
        static inline uint16_t halfsample(const uint16_t& lll, const uint16_t& llh, const uint16_t& lhl, const uint16_t& lhh, const uint16_t& hll, const uint16_t& hlh, const uint16_t& hhl, const uint16_t& hhh) {
            std::array<uint16_t, 8> vals {{lll, llh, lhl, lhh, hll, hlh, hhl, hhh}};
            const int MEDPOS = 4; // This introduces a bias towards higher values, but what can you do...
            std::nth_element(vals.begin(), vals.begin()+MEDPOS, vals.end());
            return vals[MEDPOS];
        }
    };
} // namespace anonymous

namespace voreen {

const std::string VolumeOctree::loggerCat_("voreen.VolumeOctree");

VolumeOctree::VolumeOctree(const std::vector<const VolumeBase*>& channelVolumes, size_t brickDim, float homogeneityThreshold,
            HalfSampleAggregateFunction halfSampleFn, OctreeBrickPoolManagerBase* brickPoolManager,
            size_t numThreads, ProgressReporter* progressReporter)
    : VolumeOctreeBase(tgt::svec3(brickDim), !channelVolumes.empty() ? channelVolumes.front()->getDimensions() : svec3(brickDim), channelVolumes.size())
    , rootNode_(nullptr)
    //, tempBrickBufferUsed_(false)
{
    if (channelVolumes.empty())
        throw VoreenException("No input volumes passed");
    if (channelVolumes.size() > MAX_CHANNELS)
        throw VoreenException("More than " + itos(MAX_CHANNELS) + " not supported");

    if (!brickPoolManager)
        throw VoreenException("No brick pool manager passed");
    brickPoolManager_ = brickPoolManager;

    // check whether input channel volumes have equal properties (dimension, voxel format)
    svec3 volumeDimensions = channelVolumes.front()->getDimensions();
    std::string format = channelVolumes.front()->getFormat();
    for (size_t i=1; i<channelVolumes.size(); i++) {
        if (channelVolumes.at(i)->getDimensions() != volumeDimensions)
            throw VoreenException("Dimensions of input volumes do not match: " +
                genericToString(channelVolumes.at(i)->getDimensions()) + " != " + genericToString(volumeDimensions));
        if (channelVolumes.at(i)->getFormat() != format)
            throw VoreenException("Formats of input volumes do not match: [" + channelVolumes.at(i)->getFormat() + " != " + format + "]");
    }

    tgtAssert(getNumLevels() > 0 && getNumLevels() < tgt::max(channelVolumes.front()->getDimensions()), "invalid level count");

    // convert normalized homogeneity threshold to uint16_t
    bool octreeOptimization = (homogeneityThreshold >= 0.f);
    uint16_t homogeneityThresholdUInt16 = tgt::clamp(tgt::iround(homogeneityThreshold * 65535), 0, 65535);

    // construct tree
    try {
        buildOctreeIteratively(channelVolumes, octreeOptimization, homogeneityThresholdUInt16, halfSampleFn, numThreads, progressReporter);

        LDEBUG("Flushing brick pool to disk...");
        if (progressReporter)
            progressReporter->setProgressRange(tgt::vec2(progressReporter->getProgress(), 1.f));
        brickPoolManager_->flushPoolToDisk(progressReporter);
    }
    catch (...) {
        if (brickPoolManager_ && brickPoolManager_->isInitialized()) {
            brickPoolManager_->deinitialize();
        }
        delete brickPoolManager_;
        brickPoolManager_ = nullptr;

        deleteSubTree(rootNode_);
        rootNode_ = nullptr;

        throw;
    }
    tgtAssert(rootNode_, "no root node after octree construction");
    updateTreeMetaDataCache();

    if (progressReporter) {
        progressReporter->setProgressRange(tgt::vec2(0.f, 1.f));
        progressReporter->setProgress(1.f);
    }
}

VolumeOctree::VolumeOctree(const VolumeBase* volume, size_t brickDim, float homogeneityThreshold,
            HalfSampleAggregateFunction halfSampleFn, OctreeBrickPoolManagerBase* brickPoolManager,
            size_t numThreads, ProgressReporter* progressReporter)
    : VolumeOctree(std::vector<const VolumeBase*> {volume}, brickDim, homogeneityThreshold, halfSampleFn, brickPoolManager, numThreads, progressReporter)
{ }

/**
 * Construct a VolumeOctree from preprocessed parts, i.e., an existing hierarchy of nodes whose bricks are stored
 * in the passed brickPoolManager.
 *
 * Note that the VolumeOctree takes ownership of the passed histograms.
 */
VolumeOctree::VolumeOctree(VolumeOctreeNode* root, OctreeBrickPoolManagerBase* brickPoolManager, std::vector<Histogram1D*>&& histograms, const tgt::svec3& brickDim, const tgt::svec3& volumeDim, size_t numChannels)
    : VolumeOctreeBase(brickDim, volumeDim, numChannels)
    , rootNode_(root)
    , brickPoolManager_(brickPoolManager)
    , histograms_(std::move(histograms))
{
    tgtAssert(tgt::hmul(brickDim) * sizeof(uint16_t) == brickPoolManager_->getBrickMemorySizeInByte(), "Brick size mismatch");
    updateTreeMetaDataCache();
}

std::pair<OctreeBrickPoolManagerBase*, VolumeOctreeNode*> VolumeOctree::decompose() && {
    auto tmp = std::move(*this);
    auto ret = std::make_pair(tmp.brickPoolManager_, tmp.rootNode_);
    tmp.rootNode_ = nullptr;
    tmp.brickPoolManager_ = nullptr;
    return ret;
}


VolumeOctree::VolumeOctree(VolumeOctree&& other)
    : rootNode_(other.rootNode_)
    , brickPoolManager_(other.brickPoolManager_)
    , numNodes_(0)
    , numBricks_(0)
    , actualDepth_(0)
{
    other.rootNode_ = nullptr;
    other.brickPoolManager_ = nullptr;
}
VolumeOctree& VolumeOctree::operator=(VolumeOctree&& other) {
    if(&other != this) {
        this->~VolumeOctree();
        new(this) VolumeOctree(std::move(other));
    }
    return *this;
}

// default constructor for serialization (private)
VolumeOctree::VolumeOctree()
    : rootNode_(nullptr)
    , brickPoolManager_(0)
{}

VolumeOctree::~VolumeOctree() {
    if (rootNode_)
        deleteSubTree(rootNode_);
    rootNode_ = 0;

    // deinitialize brick pool manager
    //tgtAssert(brickPoolManager_, "no brick pool manager");
    if (brickPoolManager_)
        brickPoolManager_->deinitialize();
    delete brickPoolManager_;
    brickPoolManager_ = 0;

    for (size_t i=0; i<histograms_.size(); i++)
        delete histograms_.at(i);
    histograms_.clear();
}

void VolumeOctree::updateTreeMetaDataCache() {
    tgtAssert(rootNode_, "no root node");
    numNodes_ = rootNode_->getNodeCount();
    numBricks_ = rootNode_->getNumBricks();
    actualDepth_ = rootNode_->getDepth();
}

VolumeOctree* VolumeOctree::create() const {
    return new VolumeOctree();
}

size_t VolumeOctree::getNumNodes() const {
    tgtAssert(rootNode_, "no root node");
    return numNodes_;
}

size_t VolumeOctree::getNumBricks() const {
    tgtAssert(rootNode_, "no root node");
    return numBricks_;
}

size_t VolumeOctree::getActualTreeDepth() const {
    tgtAssert(rootNode_, "no root node");
    return actualDepth_;
}

uint64_t VolumeOctree::getBrickPoolMemoryAllocated() const {
    return brickPoolManager_->getBrickPoolMemoryAllocated();
}

uint64_t VolumeOctree::getBrickPoolMemoryUsed() const {
    tgtAssert(rootNode_, "no root node");
    return static_cast<uint64_t>(rootNode_->getNumBricks())*static_cast<uint64_t>(getBrickMemorySize());
}

std::string VolumeOctree::getDescription() const {
    std::ostringstream desc;

    desc << "Volume dim: " << getVolumeDim() << std::endl;
    desc << "Channels:   " << getNumChannels() << std::endl;
    desc << "Octree dim: " << getOctreeDim() << std::endl;
    desc << "Brick dim:  " << getBrickDim() << std::endl;
    desc << "Num levels: " << getNumLevels() << std::endl;
    desc << "Tree depth: " << getActualTreeDepth() << std::endl;
    desc << "Num nodes/bricks: \t" << getNumNodes() << "/" << rootNode_->getNumBricks() << std::endl;

    size_t brickMemSize = getBrickMemorySize();
    desc << "Brick memory size: \t" << formatMemorySize(brickMemSize) << " \t(" << brickMemSize << " Bytes)" << std::endl;

    uint64_t brickBufferMemUsed = getBrickPoolMemoryUsed();
    uint64_t brickBufferMemAllocated = getBrickPoolMemoryAllocated();
    desc << "Brick pool mem used:\t" << formatMemorySize(brickBufferMemUsed) << " \t(" << brickBufferMemUsed << " Bytes)" << std::endl;
    desc << "Brick pool mem alloc:\t" << formatMemorySize(brickBufferMemAllocated) << " \t(" << brickBufferMemAllocated << " Bytes)" << std::endl;

    uint64_t volumeMemSize = static_cast<uint64_t>(tgt::hmul(getVolumeDim()))*getBytesPerVoxel();
    desc << "Volume memory size: \t" << formatMemorySize(volumeMemSize) << " \t(" << volumeMemSize << " Bytes)" /*<< std::endl()*/;

    return desc.str();
}

const Histogram1D* VolumeOctree::getHistogram(size_t channel /*= 0*/) const {
    tgtAssert(channel < getNumChannels(), "invalid channel");
    tgtAssert(histograms_.size() == getNumChannels(), "number of histograms does not match channel count");
    return histograms_.at(channel);
}

uint16_t VolumeOctree::getVoxel(const tgt::svec3& pos, size_t channel /*= 0*/, size_t nodeLevel /* = 0*/) const {
    tgtAssert(channel < getNumChannels(), "invalid channel");

    if (tgt::hor(tgt::greaterThanEqual(pos, getVolumeDim())))
        throw std::invalid_argument("Voxel outside volume dimensions: " + genericToString(pos));

    svec3 nodeLLF, nodeURB;
    vec3 dummy;
    const VolumeOctreeNode* node = getNode(vec3(pos) / vec3(getVolumeDim()), nodeLevel, nodeLLF, nodeURB /*just outside the node*/, dummy, dummy);
    tgtAssert(node, "null pointer returned as node");

    if (!node->hasBrick()) {
        return node->getAvgValue(channel);
    }
    else { // sample brick
        const uint16_t* brick = getNodeBrick(node);
        tgtAssert(brick, "no brick returned");

        // if the node has a brick it must be at level 0 (our target level), we can thus just sample using pos - nodeLLF within the brick
        //svec3 brickPos = pos - nodeLLF;
        svec3 brickPos = svec3(tgt::ifloor((vec3(pos-nodeLLF)/vec3(nodeURB-nodeLLF)) * vec3(getBrickDim())));
        brickPos = tgt::clamp(brickPos, svec3::zero, getBrickDim()-svec3(1));
        tgtAssert(tgt::hand(tgt::greaterThanEqual(brickPos, svec3::zero)) && tgt::hand(lessThan(brickPos, getBrickDim())), "Position within brick outside of brick dimensions!");

        size_t voxelOffset = cubicCoordToLinear(brickPos, getBrickDim()) * getNumChannels() + channel;
        tgtAssert(voxelOffset < (getBrickMemorySize() / getBytesPerVoxel() * getNumChannels()), "invalid voxel offset");
        uint16_t voxelValue = brick[voxelOffset];

        releaseNodeBrick(node);

        return voxelValue;
    }
}

const VolumeOctreeNode* VolumeOctree::getRootNode() const {
    return rootNode_;
}
VolumeOctreeNode* VolumeOctree::getRootNode() {
    return rootNode_;
}

const VolumeOctreeNode* VolumeOctree::getNode(const tgt::vec3& point, size_t& level,
    tgt::svec3& voxelLLF, tgt::svec3& voxelURB /* just OUTSIDE */, tgt::vec3& normLLF, tgt::vec3& normURB) const
{
    tgtAssert(rootNode_, "no root node");
    if (!inRange(point, vec3(0.f), vec3(1.f)))
        throw std::invalid_argument("Point coordinates outside unit cube: " + genericToString(point));
    if (level >= getNumLevels())
        throw std::invalid_argument("Passed level larger than octree depth: " + itos(level) + " (octree depth: " + itos(getNumLevels()) + ")");

    svec3 voxel = samplePosToVoxel(point, getVolumeDim());
    tgtAssert(tgt::hand(tgt::lessThan(voxel, getVolumeDim())), "invalid voxel");

    size_t resultLevel;
    const VolumeOctreeNode* node = getNodeAtVoxel(voxel, getNumLevels()-1, level, rootNode_, svec3::zero, getOctreeDim(), resultLevel, voxelLLF, voxelURB);
    tgtAssert(resultLevel >= level, "invalid return level");
    tgtAssert(tgt::hand(tgt::lessThan(voxelLLF, voxelURB)) && tgt::hand(tgt::lessThanEqual(voxelURB, getOctreeDim())), "invalid llf/urb");

    level = resultLevel;
    normLLF = vec3(voxelLLF) / vec3(getVolumeDim());
    normURB = vec3(voxelURB) / vec3(getVolumeDim());

    return node;
}

const uint16_t* VolumeOctree::getNodeBrick(const VolumeOctreeNode* node) const {
    tgtAssert(brickPoolManager_, "no brick pool manager");
    tgtAssert(node, "null pointer");

    return brickPoolManager_->getBrick(node->getBrickAddress());
}

void VolumeOctree::releaseNodeBrick(const VolumeOctreeNode* node) const {
    tgtAssert(brickPoolManager_, "no brick pool manager");
    tgtAssert(node, "null pointer");

    return brickPoolManager_->releaseBrick(node->getBrickAddress());

}

VolumeRAM* VolumeOctree::createVolume(size_t level /*= 0*/, clock_t timeLimit /*= 0*/, bool* complete /*= 0*/) const {
    if (level >= getNumLevels())
        throw std::invalid_argument("Passed level larger than octree depth: " + itos(level) + " (octree depth: " + itos(getNumLevels()) + ")");

    // start clock
    tgt::Stopwatch runtimeWatch;
    runtimeWatch.start();

    svec3 levelVolumeDim = tgt::ceil(tgt::vec3(getVolumeDim()) / vec3(1 << (level)));
    VolumeRAM* levelVolumeRam = 0;
    try {
        switch (getNumChannels()) {
        case 1:
            levelVolumeRam = new VolumeRAM_UInt16(levelVolumeDim, true);
            break;
        case 2:
            levelVolumeRam = new VolumeRAM_2xUInt16(levelVolumeDim, true);
            break;
        case 3:
            levelVolumeRam = new VolumeRAM_3xUInt16(levelVolumeDim, true);
            break;
        case 4:
            levelVolumeRam = new VolumeRAM_4xUInt16(levelVolumeDim, true);
            break;
        default:
            tgtAssert(false, "more than 4 channels");
        }
    }
    catch (std::bad_alloc&) {
        LERROR("Failed to allocate texture buffer");
        return 0;
    }
    tgtAssert(levelVolumeRam, "output volume not created");
    uint16_t* levelVolumeBuffer = reinterpret_cast<uint16_t*>(levelVolumeRam->getData());

    bool volumeComplete;
    composeNodeTexture(rootNode_, svec3(0,0,0), getNumLevels()-1, level, levelVolumeBuffer, levelVolumeDim, timeLimit, runtimeWatch, volumeComplete);
    if (complete)
        *complete = volumeComplete;

    return levelVolumeRam;
}

VolumeRAM* VolumeOctree::createSlice(SliceAlignment sliceAlignment, size_t sliceIndex, size_t level /*= 0*/,
    clock_t timeLimit /*= 0*/, bool* complete /*= 0*/, tgt::svec3 begin /*= tgt::svec3::zero*/, tgt::svec3 end /*= tgt::svec3(-1)*/) const {
    if (level >= getNumLevels())
        throw std::invalid_argument("Passed level larger than octree depth: " + itos(level) + " (octree depth: " + itos(getNumLevels()) + ")");
    if (sliceIndex >= getVolumeDim()[sliceAlignment])
        throw std::invalid_argument("Invalid slice index:" + itos(sliceIndex));
    if(end == tgt::svec3(-1))
        end = getVolumeDim() - tgt::svec3::one;
    if(!tgt::hand(tgt::lessThanEqual(begin,end)))
        throw std::invalid_argument("End is smaller than begin: "+ genericToString(end) + " < " + genericToString(begin));
    if(!tgt::hand(tgt::lessThan(begin,getVolumeDim())))
        throw std::invalid_argument("Begin must be smaller than getVolumeDim(): "+ genericToString(getVolumeDim()) + " < " + genericToString(begin));
    if(tgt::clamp(tgt::svec3(sliceIndex),begin,end)[sliceAlignment] != sliceIndex)
        throw std::invalid_argument("SliceIndex must between begin and end");

    // start clock
    tgt::Stopwatch runtimeWatch;
    if(timeLimit > 0)
        runtimeWatch.start();

    // determine dimension of output slice with regard to selected level
    const svec3 levelVolumeDim = tgt::ceil(tgt::vec3(end + tgt::svec3::one - begin) / vec3(1 << level));

    svec3 sliceVolumeDim = levelVolumeDim;
    sliceVolumeDim[sliceAlignment] = 1;

    // transform requested slice index/begin/end to selected level
    const tgt::svec3 beginAtLevel = begin / static_cast<size_t>(1 << level);
    const tgt::svec3 endAtLevel = beginAtLevel + levelVolumeDim;
    const size_t sliceIndexAtLevel = std::min(sliceIndex / static_cast<size_t>(1 << level), (beginAtLevel+levelVolumeDim)[sliceAlignment]-1);

    // allocate output slice volume
    VolumeRAM* sliceVolumeRam = 0;
    try {
        switch (getNumChannels()) {
        case 1:
            sliceVolumeRam = new VolumeRAM_UInt16(sliceVolumeDim, true);
            break;
        case 2:
            sliceVolumeRam = new VolumeRAM_2xUInt16(sliceVolumeDim, true);
            break;
        case 3:
            sliceVolumeRam = new VolumeRAM_3xUInt16(sliceVolumeDim, true);
            break;
        case 4:
            sliceVolumeRam = new VolumeRAM_4xUInt16(sliceVolumeDim, true);
            break;
        default:
            tgtAssert(false, "more than 4 channels");
        }
    }
    catch (std::bad_alloc&) {
        LERROR("Failed to allocate texture buffer");
        return 0;
    }
    tgtAssert(sliceVolumeRam, "output volume not created");

    // recursively copy voxel data from bricks to output slice
    uint16_t* sliceVolumeBuffer = reinterpret_cast<uint16_t*>(sliceVolumeRam->getData());
    bool sliceComplete;
    composeNodeSliceTexture(sliceAlignment, rootNode_, svec3(0,0,0), sliceIndexAtLevel, getNumLevels()-1, level,
        sliceVolumeBuffer, sliceVolumeDim, timeLimit, runtimeWatch, sliceComplete, beginAtLevel, endAtLevel);
    if (complete)
        *complete = sliceComplete;

    return sliceVolumeRam;
}


//------------------
// private functions

void VolumeOctree::buildOctreeIteratively(const std::vector<const VolumeBase*>& volumes, bool octreeOptimization,
            uint16_t homogeneityThreshold, HalfSampleAggregateFunction halfSampleFn, size_t numThreads,
            ProgressReporter* progressReporter) {
    tgtAssert(volumes.size() > 0, "no channel volumes passed");
    tgtAssert(brickPoolManager_, "no brick pool manager");
    tgtAssert(numThreads > 0, "num threads must no be zero");

    // initialize brick pool manager
    tgtAssert(brickPoolManager_, "no brick pool manager");
    brickPoolManager_->initialize(getBrickMemorySize());

    // setup multi-threading
#ifndef VRN_MODULE_OPENMP
    if (numThreads > 1) {
        LWARNING("Multi-threaded octree construction requires the OpenMP module");
        numThreads = 1;
    }
#endif
    numThreads = std::max<size_t>(numThreads, 1);

    // log construction / memory usage info
    LDEBUG("Creating octree iteratively (" <<
        "Volume dim: "      << getVolumeDim()   << ", " <<
        "Volume mem size: " << formatMemorySize(tgt::hmul(getVolumeDim())*getBytesPerVoxel()) << ", " <<
        "Channels:   "      << getNumChannels() << ", " <<
        "Octree dim: "      << getOctreeDim()   << ", " <<
        "Octree depth: "    << getNumLevels()   << ", " <<
        "Brick dim:  "      << getBrickDim()    << ", " <<
        "Brick mem size: "  << formatMemorySize(getBrickMemorySize()) << ", " <<
        "Brick pool manager: " << brickPoolManager_->getClassName() << " [" << brickPoolManager_->getDescription() << "], " <<
        "Num threads: "     << numThreads << ")";
    );
    LDEBUG(MemoryInfo::getProcessMemoryUsageAsString());
    LDEBUG(MemoryInfo::getAvailableMemoryAsString());

    // determine whether tree is built from VolumeRAM or VolumeDisk representations, and determine input volume type
    bool ramMode = true;
    const std::string inputDataFormat = volumes.front()->getFormat();
    for (size_t i=0; i<volumes.size(); i++) {
        ramMode &= volumes.at(i)->hasRepresentation<VolumeRAM>();
        tgtAssert(volumes.at(i)->getFormat() == inputDataFormat, "format of input volumes does not match"); //< checked by constructor
    }
    if (!ramMode) { // volumes not present in RAM => check that disk representations are available
        for (size_t i=0; i<volumes.size(); i++) {
            if (!volumes.at(i)->getRepresentation<VolumeDisk>())
                throw VoreenException("Neither VolumeRAM nor VolumeDisk representation(s) available");
        }
    }

    if (progressReporter)
        progressReporter->setProgress(0.f);

    // initialize histograms
    std::vector< std::vector< std::vector<uint64_t> > > histogramBuffers; //< one buffer per thread per channel
    for (size_t threadID=0; threadID<numThreads; threadID++) {
        std::vector< std::vector<uint64_t> > threadBuffers;
        for (size_t ch=0; ch<getNumChannels(); ch++)
            threadBuffers.push_back(std::vector<uint64_t>(NUM_HISTOGRAM_BUCKETS, 0));
        histogramBuffers.push_back(threadBuffers);
    }

    //
    // 1. Create nodes at level 0 (full resolution bricks) from input volumes
    //
    LDEBUG("- Creating level 0 nodes");
    tgt::svec3 numNodesPerDim = tgt::ceil(tgt::vec3(getDimensions()) / tgt::vec3(getBrickDim()));
    std::unique_ptr<NodeGrid3D> level0Grid(new NodeGrid3D(numNodesPerDim));
    for (size_t nodeIndexZ = 0; nodeIndexZ < numNodesPerDim.z; nodeIndexZ++) {
        size_t startSlice = nodeIndexZ*getBrickDim().z;

        // if nodes completely outside volume (first z slice index >= volumeDim.z) => create empty dummy nodes
        tgtAssert(startSlice < getVolumeDim().z, "Invalid start slice");
        size_t endSlice = std::min(startSlice + getBrickDim().z-1, getVolumeDim().z-1);

        // load slices for current brick plate [startSlice;endSlice]
        size_t sliceRangeMemSize = (endSlice-startSlice+1)*tgt::hmul(getVolumeDim().xy())*getBytesPerVoxel() / getNumChannels();   // mem size per channel !!
        LDEBUG("-- Loading slice range [" << startSlice << "," << endSlice << "] (memory size: " << itos(getNumChannels()) << "x" << formatMemorySize(sliceRangeMemSize) << ")");
        LDEBUG("--- Before: " << MemoryInfo::getProcessMemoryUsageAsString());
        LDEBUG("--- Before: " << MemoryInfo::getAvailableMemoryAsString());

        // RAII helper struct holding the temporary slice range volumes
        struct SliceVolumes {
            std::vector<const VolumeRAM*> channelVolumes;
            std::vector<const void*> channelDataBuffers;
            ~SliceVolumes() {
                LDEBUG("-- Deleting loaded/extracted slices");
                for (size_t i=0; i<channelVolumes.size(); i++)
                    delete channelVolumes.at(i);

                LDEBUG("--- After: " << MemoryInfo::getProcessMemoryUsageAsString());
                LDEBUG("--- After: " << MemoryInfo::getAvailableMemoryAsString());
            }
        };
        SliceVolumes sliceVolumes;

        // create slice volumes
        for (size_t ch=0; ch<getNumChannels(); ch++) {
            VolumeRAM* sliceVolume = 0;
            if (ramMode) { // extract current slice range from RAM volume
                const VolumeRAM* channelVolume = volumes.at(ch)->getRepresentation<VolumeRAM>();
                tgtAssert(channelVolume, "no ram representation"); //< checked above
                tgt::svec3 slicePackageLLF(0, 0, startSlice);
                tgt::svec3 slicePackageDim(channelVolume->getDimensions().x, channelVolume->getDimensions().y, endSlice-startSlice+1);
                sliceVolume = channelVolume->getSubVolume(slicePackageDim, slicePackageLLF);
                if (!sliceVolume)
                    throw VoreenException("Failed to extract slices [" + itos(startSlice) + "," + itos(endSlice) + "] from input RAM volume for channel " + itos(ch));
            }
            else { // load current slice range from disk
                const VolumeDisk* channelDiskVolume = volumes.at(ch)->getRepresentation<VolumeDisk>();
                tgtAssert(channelDiskVolume, "no disk representation"); //< checked above
                sliceVolume = channelDiskVolume->loadSlices(startSlice, endSlice);
                if (!sliceVolume)
                    throw VoreenException("Failed to extract slices [" + itos(startSlice) + "," + itos(endSlice) + "] from input disk volume for channel " + itos(ch));
            }
            tgtAssert(sliceVolume, "slice volume not created"); //< should have thrown exception
            sliceVolumes.channelVolumes.push_back(sliceVolume);
            sliceVolumes.channelDataBuffers.push_back(sliceVolume->getData());
        }
        tgtAssert(sliceVolumes.channelVolumes.size() == getNumChannels() && sliceVolumes.channelDataBuffers.size() == getNumChannels(),
            "invalid number of slice volumes");
        LDEBUG("--- After: " << MemoryInfo::getProcessMemoryUsageAsString());
        LDEBUG("--- After: " << MemoryInfo::getAvailableMemoryAsString());

        // create nodes for current slice plate [startSlice;endSlice] (multi-threaded)
        LDEBUG("-- Creating nodes for slice range [" << startSlice << "," << endSlice << "]");
        LDEBUG("--- Before: " << MemoryInfo::getProcessMemoryUsageAsString());
        LDEBUG("--- Before: " << MemoryInfo::getAvailableMemoryAsString());

        #ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for
        #endif
        for (int threadID = 0; threadID < (int)numThreads; threadID++) {
            tgtAssert((int)histogramBuffers.size() > threadID, "missing histogram buffer");
            uint16_t avgValues[MAX_CHANNELS], minValues[MAX_CHANNELS], maxValues[MAX_CHANNELS];

            tgt::svec3 nodeIndex(0, 0, nodeIndexZ);
            for (nodeIndex.y = threadID; nodeIndex.y < numNodesPerDim.y; nodeIndex.y += numThreads) {
                for (nodeIndex.x = 0; nodeIndex.x < numNodesPerDim.x; nodeIndex.x++) {
                    VolumeOctreeNode* node = 0;
                    tgt::svec3 nodeLLF = nodeIndex*getBrickDim();
                    tgt::svec3 nodeURB = nodeLLF + getBrickDim();
                    if (tgt::hor(tgt::greaterThanEqual(nodeLLF, getVolumeDim()))) { // outside volume in x or y direction => create empty node
                        node = VolumeOctreeBase::createNode(getNumChannels());
                    }
                    else { // create node from slice volume data buffer(s)
                        nodeLLF.z = 0;
                        nodeURB.z = getBrickDim().z;
                        tgt::svec3 textureDim(getVolumeDim().x, getVolumeDim().y, endSlice-startSlice+1);
                        tgtAssert(textureDim == sliceVolumes.channelVolumes.front()->getDimensions(), "invalid texture dim");
                        if (inputDataFormat == "uint8")
                            node = createTreeNodeFromTexture<uint8_t>(nodeLLF, nodeURB,
                                sliceVolumes.channelDataBuffers, textureDim,
                                octreeOptimization, homogeneityThreshold, avgValues, minValues, maxValues, histogramBuffers.at(threadID));
                        else if (inputDataFormat == "int8")
                            node = createTreeNodeFromTexture<int8_t>(nodeLLF, nodeURB,
                                sliceVolumes.channelDataBuffers, textureDim,
                                octreeOptimization, homogeneityThreshold, avgValues, minValues, maxValues, histogramBuffers.at(threadID));
                        else if (inputDataFormat == "uint16")
                            node = createTreeNodeFromTexture<uint16_t>(nodeLLF, nodeURB,
                                sliceVolumes.channelDataBuffers, textureDim,
                                octreeOptimization, homogeneityThreshold, avgValues, minValues, maxValues, histogramBuffers.at(threadID));
                        else if (inputDataFormat == "int16")
                            node = createTreeNodeFromTexture<int16_t>(nodeLLF, nodeURB,
                                sliceVolumes.channelDataBuffers, textureDim,
                                octreeOptimization, homogeneityThreshold, avgValues, minValues, maxValues, histogramBuffers.at(threadID));
                        else if (inputDataFormat == "uint32")
                            node = createTreeNodeFromTexture<uint32_t>(nodeLLF, nodeURB,
                                sliceVolumes.channelDataBuffers, textureDim,
                                octreeOptimization, homogeneityThreshold, avgValues, minValues, maxValues, histogramBuffers.at(threadID));
                        else if (inputDataFormat == "int32")
                            node = createTreeNodeFromTexture<int32_t>(nodeLLF, nodeURB,
                                sliceVolumes.channelDataBuffers, textureDim,
                                octreeOptimization, homogeneityThreshold, avgValues, minValues, maxValues, histogramBuffers.at(threadID));
                        else if (inputDataFormat == "float")
                            node = createTreeNodeFromTexture<float>(nodeLLF, nodeURB,
                                sliceVolumes.channelDataBuffers, textureDim,
                                octreeOptimization, homogeneityThreshold, avgValues, minValues, maxValues, histogramBuffers.at(threadID));
                        else if (inputDataFormat == "double")
                            node = createTreeNodeFromTexture<double>(nodeLLF, nodeURB,
                                sliceVolumes.channelDataBuffers, textureDim,
                                octreeOptimization, homogeneityThreshold, avgValues, minValues, maxValues, histogramBuffers.at(threadID));
                        else
                            throw VoreenException("Unknown/unsupported input data format: " + inputDataFormat);

                        tgtAssert(node, "no node created");
                        tgtAssert(minValues[0] <= avgValues[0] && avgValues[0] <= maxValues[0], "invalid avg/min/max values");
                        tgtAssert(node->getAvgValue() == avgValues[0] && node->getMinValue() == minValues[0] && node->getMaxValue() == maxValues[0],
                            "avg/min/max values of returned node differ from returned avg/min/max values");

                    } // node creation
                    tgtAssert(node, "no node created");

                    #ifdef VRN_MODULE_OPENMP
                    #pragma omp critical
                    #endif
                    {
                    tgtAssert(level0Grid->getNode(nodeIndex) == 0, "node already created");
                    level0Grid->setNode(node, nodeIndex);
                    }

                } // nodeIndex.x
            } // nodeIndex.y
        } // multi-threading loop

        // update progress bar
        if (progressReporter) {
            float level0Progress = (float)(endSlice-1) / getVolumeDim().z;
            progressReporter->setProgress(level0Progress*0.7f);
        }

    } // nodeIndexZ (slice plate)

    tgtAssert(level0Grid->isComplete(), "level 0 node grid not completely constructed");


    //
    // 2. iteratively create next level nodes from current level
    //

    size_t currentLevel = 1;
    // brickDim == octreeDim => tree has only one level => already finished
    if (level0Grid->getDim() == tgt::svec3::one) {
        rootNode_ = level0Grid->takeNode(tgt::svec3(0, 0, 0));
    }
    else { // iterative construction of upper levels
        std::unique_ptr<NodeGrid3D> currentLevelGrid = std::move(level0Grid);
        while (tgt::hor(tgt::greaterThan(currentLevelGrid->getDim(), tgt::svec3::two))) {
            // check current grid
            tgtAssert(currentLevelGrid->isComplete(), "current level grid is not complete");

            // log memory usage
            LDEBUG("- Creating level " << currentLevel << " nodes");
            LDEBUG("-- Before: " << MemoryInfo::getProcessMemoryUsageAsString());
            LDEBUG("-- Before: " << MemoryInfo::getAvailableMemoryAsString());

            const tgt::svec3 parentLevelGridDim = tgt::ceil(tgt::vec3(currentLevelGrid->getDim()) / tgt::vec3::two);
            const tgt::svec3 childLevelVolumeDim = tgt::ceil(tgt::vec3(getVolumeDim()) / tgt::vec3(1 << (currentLevel-1)));
            const svec3 brickDim = getBrickDim();

            // create parent level grid
            std::unique_ptr<NodeGrid3D> parentLevelGrid(new NodeGrid3D(parentLevelGridDim));

            // create parent nodes
            for (size_t parentNodeZ=0; parentNodeZ<parentLevelGridDim.z; parentNodeZ++) {

                #ifdef VRN_MODULE_OPENMP
                #pragma omp parallel for
                #endif
                for (int threadID = 0; threadID < (int)numThreads; threadID++) {
                    tgt::svec3 parentNodeID(0, 0, parentNodeZ);
                    for (parentNodeID.y=threadID; parentNodeID.y<parentLevelGridDim.y; parentNodeID.y += numThreads) {
                        for (parentNodeID.x=0; parentNodeID.x<parentLevelGridDim.x; parentNodeID.x++) {
                            VolumeOctreeNode* childNodes[8];
                            for (size_t i=0; i<8; i++) {
                                tgt::svec3 childNodeID = linearCoordToCubic(i, tgt::svec3::two);
                                tgt::svec3 childPos = parentNodeID*tgt::svec3::two + childNodeID;
                                VolumeOctreeNode* childNode;
                                if(tgt::hand(tgt::lessThan(childPos, currentLevelGrid->getDim()))) {
                                    childNode = currentLevelGrid->takeNode(childPos);
                                } else {
                                    // Create dummy node.
                                    childNode = VolumeOctreeBase::createNode(getNumChannels());
                                }
                                childNodes[i] = childNode;
                            }
                            tgt::svec3 childLevelVolumeLlf = tgt::min(parentNodeID * brickDim * tgt::svec3::two, childLevelVolumeDim);
                            tgt::svec3 childLevelVolumeUrb = tgt::min(childLevelVolumeLlf + (brickDim * tgt::svec3::two), childLevelVolumeDim);
                            tgt::svec3 inBrickUrb = childLevelVolumeUrb - childLevelVolumeLlf;

                            uint16_t avgValues[MAX_CHANNELS], minValues[MAX_CHANNELS], maxValues[MAX_CHANNELS];
                            VolumeOctreeNode* parentNode = createParentNode(childNodes, octreeOptimization, homogeneityThreshold,
                                inBrickUrb, avgValues, minValues, maxValues, halfSampleFn);
                            tgtAssert(minValues[0] <= avgValues[0] && avgValues[0] <= maxValues[0], "invalid avg/min/max values");
                            tgtAssert(parentNode->getAvgValue() == avgValues[0] && parentNode->getMinValue() == minValues[0] && parentNode->getMaxValue() == maxValues[0],
                                "avg/min/max values of returned node differ from returned avg/min/max values");

                            #ifdef VRN_MODULE_OPENMP
                            #pragma omp critical
                            #endif
                            {
                            tgtAssert(parentLevelGrid->getNode(parentNodeID) == 0, "parent node already exists");
                            parentLevelGrid->setNode(parentNode, parentNodeID);
                            }

                        } // parentNodeIndex.x
                    } // parentNodeIndex.y
                } // multi-threading loop

                // update progress bar
                if (progressReporter && currentLevel == 1) {
                    float inLevelProgress = (float)(parentNodeZ+1) / (float)parentLevelGridDim.z;
                    progressReporter->setProgress(0.7f + inLevelProgress*0.2f);
                }

            } // parentNodeIndexZ

            // log memory usage
            LDEBUG("-- After: " << MemoryInfo::getProcessMemoryUsageAsString());
            LDEBUG("-- After: " << MemoryInfo::getAvailableMemoryAsString());

            tgtAssert(currentLevelGrid->isEmpty(), "Grid not empty");
            // advance to parent level
            currentLevelGrid = std::move(parentLevelGrid);
            currentLevel++;
        }
        tgtAssert(tgt::hand(tgt::lessThanEqual(currentLevelGrid->getDim(), tgt::svec3::two)), "level grid dimensions of [2 2 2] or smaller expected");
        tgtAssert(tgt::hand(tgt::notEqual(currentLevelGrid->getDim(), tgt::svec3::zero)), "level grid dimensions contain zero");

        LDEBUG("- Creating root node");
        uint16_t avgValues[MAX_CHANNELS], minValues[MAX_CHANNELS], maxValues[MAX_CHANNELS];

        tgt::svec3 inBrickUrb = tgt::ceil(tgt::vec3(getVolumeDim()) / tgt::vec3(1 << (currentLevel - 1)));

        VolumeOctreeNode* childNodes[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        VRN_FOR_EACH_VOXEL(childPos, tgt::svec3::zero, currentLevelGrid->getDim()) {
            size_t childIndex = cubicCoordToLinear(childPos, tgt::svec3::two);
            childNodes[childIndex] = currentLevelGrid->takeNode(childPos);
        }
        for(auto& node : childNodes) {
            if(!node) {
                // Create dummy node
                node = VolumeOctreeBase::createNode(getNumChannels());
            }
        }
        rootNode_ = createParentNode(childNodes, octreeOptimization, homogeneityThreshold,
            inBrickUrb, avgValues, minValues, maxValues, halfSampleFn);
        tgtAssert(minValues[0] <= avgValues[0] && avgValues[0] <= maxValues[0], "invalid avg/min/max values");
    }

    tgtAssert(rootNode_, "no root node");
    tgtAssert(rootNode_->getNumBricks() <= rootNode_->getNodeCount(), "number of bricks larger than number of nodes");

    // create histograms from buffers
    tgtAssert(histogramBuffers.size() == numThreads, "invalid histogram buffer");
    tgtAssert(histogramBuffers.front().size() == getNumChannels(), "invalid histogram buffer");
    const size_t numBuckets = 1<<std::min<size_t>(NUM_HISTOGRAM_BUCKETS_BITS, volumes.front()->getBytesPerVoxel() * 8); //< 4096 for 16 bit and more, 256 for 8 bit
    const size_t bucketSize = NUM_HISTOGRAM_BUCKETS / numBuckets;
    const float realWorldMin = volumes.front()->getRealWorldMapping().normalizedToRealWorld(0.f);
    const float realWorldMax = volumes.front()->getRealWorldMapping().normalizedToRealWorld(1.f);
    for (size_t ch=0; ch<getNumChannels(); ch++) {
        Histogram1D* channelHistogram = new Histogram1D(realWorldMin, realWorldMax, (int)numBuckets);
        for (size_t threadID=0; threadID<numThreads; threadID++) {
            tgtAssert(histogramBuffers.at(threadID).size() > ch, "invalid histogram buffer");
            std::vector<uint64_t>& threadHistBuffer = histogramBuffers.at(threadID).at(ch);
            tgtAssert(threadHistBuffer.size() == NUM_HISTOGRAM_BUCKETS, "invalid histogram buffer size");
            for (size_t bucket=0; bucket<numBuckets; bucket++) {
                uint64_t bucketValue = 0;
                for (size_t i=bucket*bucketSize; i<(bucket+1)*bucketSize; i++)
                    bucketValue += threadHistBuffer[i];
                channelHistogram->increaseBucket(bucket, bucketValue);
            }
        }
        histograms_.push_back(channelHistogram);
    }

    // log memory usage info
    LDEBUG("Finished octree creation (" << brickPoolManager_->getClassName() << " [" << brickPoolManager_->getDescription() << "]" << ")");
    LDEBUG("- " << MemoryInfo::getProcessMemoryUsageAsString());
    LDEBUG("- " << MemoryInfo::getAvailableMemoryAsString());

}

const VolumeOctreeNode* VolumeOctree::getNodeAtVoxel(const svec3& voxel, const size_t curLevel, const size_t targetLevel,
    const VolumeOctreeNode* node, const svec3& nodeLlf, const svec3& nodeUrb /* just OUTSIDE */,
    size_t& resultLevel, svec3& resultLlf, svec3& resultUrb /* just OUTSIDE */) const
{
    tgtAssert(node, "null pointer passed");
    tgtAssert(curLevel >= targetLevel && curLevel < getNumLevels(), "invalid current level");
    tgtAssert(inRange(voxel, nodeLlf, nodeUrb-svec3(1)), "point coords outside node dimensions");

    if (curLevel == targetLevel || !node->children_[0]) { // current node level requested, or current node is leaf => stop descent
        resultLevel = curLevel;
        resultLlf = nodeLlf;
        resultUrb = nodeUrb;
        return node;
    }
    else {
        svec3 nodeDim = nodeUrb - nodeLlf;
        svec3 nodeHalfDim = nodeDim/svec3(2);
        svec3 voxelOffset = voxel - nodeLlf;
        tgtAssert(inRange(voxelOffset, svec3::zero, nodeDim-svec3(1)), "invalid voxel offset");
        svec3 childNodeID = voxelOffset / nodeHalfDim;
        tgtAssert(inRange(childNodeID, svec3::zero, svec3::one), "invalid child node id");

        const VolumeOctreeNode* childNode = node->children_[cubicCoordToLinear(childNodeID, svec3::two)];
        tgtAssert(childNode, "child node is null");
        svec3 childLlf = nodeLlf + childNodeID*nodeHalfDim;
        svec3 childUrb = childLlf + nodeHalfDim;
        const VolumeOctreeNode* resultNode = getNodeAtVoxel(voxel, curLevel-1, targetLevel,
            childNode, childLlf, childUrb, resultLevel, resultLlf, resultUrb);
        tgtAssert(resultNode, "null pointer returned as node");
        return resultNode;
    }
}

void VolumeOctree::composeNodeTexture(const VolumeOctreeNode* node, const svec3& nodeOffset, size_t curLevel, size_t targetLevel,
    uint16_t* textureBuffer, const svec3& textureDim, clock_t timeLimit, tgt::Stopwatch& runtimeWatch, bool& complete) const
{
    tgtAssert(brickPoolManager_, "no brick pool manager");
    tgtAssert(node, "null pointer passed");
    tgtAssert(curLevel >= targetLevel && curLevel < getNumLevels(), "invalid current level");
    tgtAssert(textureBuffer, "null pointer passed");
    tgtAssert(tgt::hand(tgt::lessThanEqual(nodeOffset, getOctreeDim())), "node offset larger than texture dimensions");

    const size_t numChannels = getNumChannels();
    svec3 nodeDim = getBrickDim() * svec3(1<<(curLevel - targetLevel));

    // skip leaf brick, if time limit for the volume creation has been reached
    bool timeLimitReached = (timeLimit > 0 && runtimeWatch.getRuntime() >= timeLimit);
    bool skipBrickDueToTimeLimit =
        (timeLimitReached && curLevel == targetLevel && node->hasBrick() && !brickPoolManager_->isBrickInRAM(node->getBrickAddress()) );

    if (node->isHomogeneous() || skipBrickDueToTimeLimit) { // homogeneous => no brick and no children => use avg value
        uint16_t avgValues[MAX_CHANNELS];
        for (size_t channel=0; channel<numChannels; channel++)
            avgValues[channel] = node->getAvgValue(channel);
        VRN_FOR_EACH_VOXEL(voxel, nodeOffset, nodeOffset+nodeDim) {
            if (tgt::hand(tgt::lessThan(voxel, textureDim))) {
                const size_t textureLinearCoord = cubicCoordToLinear(voxel, textureDim)*numChannels;
                tgtAssert(textureLinearCoord+numChannels-1 < tgt::hmul(textureDim)*numChannels, "invalid texture coord");
                std::copy(&avgValues[0], &avgValues[numChannels], textureBuffer+textureLinearCoord);
            }
        }
        complete = !skipBrickDueToTimeLimit;
    }
    else if (curLevel == targetLevel) { // final level => copy brick texture to target buffer (or use avg value)
        tgtAssert(node->hasBrick(), "node expected to have a brick");
        copyBrickToTexture(brickPoolManager_->getBrick(node->getBrickAddress()), getBrickDim(), textureBuffer, textureDim, nodeOffset);
        brickPoolManager_->releaseBrick(node->getBrickAddress());
        complete = true;
    }
    else { // higher level => let child nodes copy their sub-node textures to target texture
        tgtAssert(node->hasBrick(), "node has no brick");
        tgtAssert(!node->isLeaf(), "node not expected to be leaf"); //< higher level leaves have no brick (see above)
        svec3 subNodeDim = nodeDim / svec3(2);
        complete = true;
        VRN_FOR_EACH_VOXEL(childCoord, svec3::zero, svec3::two) {
            const VolumeOctreeNode* child = node->children_[cubicCoordToLinear(childCoord, svec3::two)];
            bool childComplete;
            tgtAssert(child, "no child node");
            composeNodeTexture(child, nodeOffset + childCoord*subNodeDim,
                curLevel-1, targetLevel, textureBuffer, textureDim, timeLimit, runtimeWatch, childComplete);
            complete &= childComplete;
        }
    }
}

void VolumeOctree::composeNodeSliceTexture(SliceAlignment sliceAlignment, const VolumeOctreeNode* node,
    const tgt::svec3& nodeOffsetInTexture, size_t sliceIndexInNode, size_t curLevel, size_t targetLevel,
    uint16_t* textureBuffer, const tgt::svec3& textureDim, clock_t timeLimit, tgt::Stopwatch& runtimeWatch,
    bool& complete, const tgt::svec3 begin, const tgt::svec3 end) const
{
    tgtAssert(brickPoolManager_, "no brick pool manager");
    tgtAssert(node, "null pointer passed");
    tgtAssert(curLevel >= targetLevel && curLevel < getNumLevels(), "invalid current level");
    tgtAssert(textureBuffer, "null pointer passed");
    tgtAssert(tgt::hand(tgt::lessThan(nodeOffsetInTexture, getOctreeDim())), "node offset larger than octree dimensions");
    tgtAssert(tgt::hor(tgt::lessThan(nodeOffsetInTexture, getVolumeDim())), "node offset outside volume dimensions");

    const size_t numChannels = getNumChannels();
    const svec3 brickDim = getBrickDim();
    const svec3 nodeDimInTexture = brickDim * svec3(1<<(curLevel-targetLevel));
    tgtAssert(sliceIndexInNode < nodeDimInTexture[sliceAlignment], "invalid slice index");

    // skip leaf brick, if time limit for the slice creation has been reached
    bool timeLimitReached = (timeLimit > 0 && runtimeWatch.getRuntime() >= timeLimit);
    bool skipBrickDueToTimeLimit =
        (timeLimitReached && curLevel == targetLevel && node->hasBrick() && !brickPoolManager_->isBrickInRAM(node->getBrickAddress()) );

    if (node->isHomogeneous() || skipBrickDueToTimeLimit) { // homogeneous => no brick and no children => use avg value
        svec3 textureVoxelStart = tgt::max(nodeOffsetInTexture, begin);
        svec3 textureVoxelEnd = tgt::min(nodeOffsetInTexture + nodeDimInTexture, end);
        textureVoxelStart[sliceAlignment] = 0;
        textureVoxelEnd[sliceAlignment] = 1;
        tgt::svec3 textureBegin = begin;
        textureBegin[sliceAlignment] = 0;
        //textureVoxelEnd = tgt::min(textureVoxelEnd, textureDim); //< node might lie partially outside target texture (NPOT) [no longer because of end]
        VRN_FOR_EACH_VOXEL(textureVoxel, textureVoxelStart, textureVoxelEnd) {
            for (size_t channel=0; channel<numChannels; channel++)
                textureBuffer[cubicCoordToLinear(textureVoxel-textureBegin, textureDim)*numChannels + channel] = node->getAvgValue(channel);
        }
        complete = !skipBrickDueToTimeLimit;
    }
    else if (curLevel == targetLevel) { // final level => copy brick slice to target texture
        tgtAssert(sliceIndexInNode < getBrickDim().z, "slice index outside brick");
        tgtAssert(node->hasBrick(), "node has no brick (should not get here)");
        const uint16_t* brick = brickPoolManager_->getBrick(node->getBrickAddress());

        svec3 textureVoxelStart = tgt::max(nodeOffsetInTexture, begin);
        svec3 textureVoxelEnd = tgt::min(nodeOffsetInTexture + brickDim,end);
        textureVoxelStart[sliceAlignment] = 0;
        textureVoxelEnd[sliceAlignment] = 1;
        tgt::svec3 textureBegin = begin;
        textureBegin[sliceAlignment] = 0;
        tgt::svec3 brickNodeOffsetInTexture = nodeOffsetInTexture;
        brickNodeOffsetInTexture[sliceAlignment] = 0;
        //textureVoxelEnd = tgt::min(textureVoxelEnd, textureDim); //< node might lie partially outside target texture (NPOT) [no longer because of end]
        VRN_FOR_EACH_VOXEL(textureVoxel, textureVoxelStart, textureVoxelEnd) {
            for (size_t channel=0; channel<numChannels; channel++) {
                size_t textureLinearCoord = cubicCoordToLinear(textureVoxel-textureBegin, textureDim)*numChannels + channel;

                svec3 brickVoxel = textureVoxel - brickNodeOffsetInTexture;
                brickVoxel[sliceAlignment] += sliceIndexInNode;
                tgtAssert(tgt::hand(tgt::lessThan(brickVoxel, brickDim)), "invalid brick voxel");
                size_t brickLinearCoord = cubicCoordToLinear(brickVoxel, brickDim)*numChannels + channel;

                textureBuffer[textureLinearCoord] = brick[brickLinearCoord];
            }
        }

        brickPoolManager_->releaseBrick(node->getBrickAddress());

        complete = true;
    }
    else { // inner node => let child nodes recursively copy their sub-node slice textures to target texture
        // note: child nodes are zyx ordered
        tgtAssert(!node->isLeaf() && node->children_[0], "node has no children (should not get here)");

        svec3 childNodeDim = nodeDimInTexture / svec3(2);
        size_t childLayer;
        if (sliceIndexInNode < childNodeDim[sliceAlignment]) { //< slice lies in lower four child nodes
            childLayer = 0;
        }
        else { //< slice lies in upper four child nodes
            childLayer = 1;
            sliceIndexInNode -= childNodeDim[sliceAlignment];
        }
        svec3 childStart = svec3(0, 0, 0);
        svec3 childEnd = svec3(2, 2, 2);
        childStart[sliceAlignment] = childLayer;
        childEnd[sliceAlignment] = childLayer+1;
        complete = true;
        VRN_FOR_EACH_VOXEL(childCoord, childStart, childEnd) {
            const VolumeOctreeNode* child = node->children_[cubicCoordToLinear(childCoord, svec3::two)];
            tgtAssert(child, "no child node");
            svec3 childNodeOffset = nodeOffsetInTexture + childCoord*childNodeDim;
            bool childComplete;
            composeNodeSliceTexture(sliceAlignment, child, childNodeOffset, sliceIndexInNode, curLevel-1, targetLevel,
                textureBuffer, textureDim, timeLimit, runtimeWatch, childComplete, begin, end);
            complete &= childComplete;
        }
    }
}


//------------------
// low-level helper functions

template<class T>
VolumeOctreeNode* VolumeOctree::createTreeNodeFromTexture(const tgt::svec3& llf, const tgt::svec3& urb,
    const std::vector<const void*>& textureBuffers, const tgt::svec3& textureDim,
    bool octreeOptimization, uint16_t homogeneityThreshold,
    uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues,
    std::vector< std::vector<uint64_t> >& histograms) {
    tgtAssert(textureBuffers.size() == getNumChannels(), "number of texture buffers does not match channel count");
    tgtAssert(getNumChannels() <= MAX_CHANNELS, "more than max channels");
    tgtAssert(textureBuffers.front(), "null pointer passed");
    tgtAssert(brickPoolManager_, "brick pool manager");
    tgtAssert(avgValues && minValues && maxValues, "null pointer passed as avg/min/max value array");

    tgtAssert(tgt::hand(tgt::lessThan(llf, getOctreeDim())), "invalid llf");
    tgtAssert(tgt::hand(tgt::lessThanEqual(urb, getOctreeDim())), "invalid urb");
    tgtAssert(tgt::hand(tgt::lessThan(llf, urb)), "llf larger than or equal urb");

    // highest level (full resolution) has been reached => create brick from volume texture and terminate recursion
    const uint64_t virtualBrickAddress = brickPoolManager_->allocateBrick();
    uint16_t* brickBuffer = brickPoolManager_->getWritableBrick(virtualBrickAddress);
    switch(getNumChannels()) {
        case 1: extractBrickFromTexture<T, 1>(textureBuffers, textureDim, brickBuffer, getBrickDim(), llf,
                        avgValues, minValues, maxValues, histograms);
                break;
        case 2: extractBrickFromTexture<T, 2>(textureBuffers, textureDim, brickBuffer, getBrickDim(), llf,
                        avgValues, minValues, maxValues, histograms);
                break;
        case 3: extractBrickFromTexture<T, 3>(textureBuffers, textureDim, brickBuffer, getBrickDim(), llf,
                        avgValues, minValues, maxValues, histograms);
                break;
        case 4: extractBrickFromTexture<T, 4>(textureBuffers, textureDim, brickBuffer, getBrickDim(), llf,
                        avgValues, minValues, maxValues, histograms);
                break;
        default:
            tgtAssert(false, "Invalid number of channels");
    }
    brickPoolManager_->releaseBrick(virtualBrickAddress, OctreeBrickPoolManagerBase::WRITE);

    // determine whether node is homogeneous (in all channels)
    bool homogeneous = true;
    for (size_t i=0; i<getNumChannels(); i++) {
        tgtAssert(minValues[i] <= avgValues[i] && avgValues[i] <= maxValues[i], "invalid avg/min/max values");
        homogeneous &= (maxValues[i] - minValues[i] <= homogeneityThreshold);
    }

    VolumeOctreeNode* node = 0;

    // node not homogeneous => create leaf node with brick and shift virtual memory offset
    if (!homogeneous || !octreeOptimization) {
        node = VolumeOctreeBase::createNode(getNumChannels(), avgValues, minValues, maxValues, virtualBrickAddress);
    }
    else { // node homogeneous => store only avg value (without brick)
        brickPoolManager_->deleteBrick(virtualBrickAddress);
        node = VolumeOctreeBase::createNode(getNumChannels(), avgValues, minValues, maxValues);
    }

    tgtAssert(node, "node not created");

    return node;
}

template<class T, size_t numChannels>
void VolumeOctree::extractBrickFromTexture(const std::vector<const void*>& textures, const svec3& textureDim,
    uint16_t* brickBuffer, const svec3& brickDim, const svec3& brickOffsetInTexture,
    uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues,
    std::vector< std::vector<uint64_t> >& histograms) const
{
    tgtAssert(textures.size() == getNumChannels(), "number of channel textures does not match channel count");
    tgtAssert(textures.front(), "null pointer passed");
    tgtAssert(brickBuffer, "null pointer passed");
    tgtAssert(avgValues && minValues && maxValues, "null pointer passed as avg/min/max value array");
    tgtAssert(histograms.size() == getNumChannels(), "invalid number of channel histograms");
    tgtAssert(histograms.front().size() == NUM_HISTOGRAM_BUCKETS, "invalid histogram buffer size");

    tgtAssert(tgt::hand(tgt::lessThanEqual(brickDim, getOctreeDim())), "brick dimensions greater than octree dimensions");
    tgtAssert(tgt::hand(tgt::lessThanEqual(brickOffsetInTexture+brickDim, getOctreeDim())), "brick (partially) outside octree dimensions");

    tgtAssert(getNumChannels() <= MAX_CHANNELS, "more than max channels");

    tgtAssert(getNumChannels() == numChannels, "Invalid number of channels");

    // Use pointer to first element of texture/histogram vector in order to avoid repeated vector lookups
    const T* const * textureArray = reinterpret_cast<const T* const *>(&textures.front());

    size_t brickBufferSize = tgt::hmul(getBrickDim())*numChannels;
    size_t textureBufferSize = tgt::hmul(textureDim);

    // Copy voxels from channel textures to brick
    tgt::svec3 inTextureSize = tgt::min(brickOffsetInTexture + brickDim, textureDim) - brickOffsetInTexture;
    uint64_t numSignificantBrickVoxels = tgt::hmul(inTextureSize); //< number of brick voxels lying inside texture (NPOT)

    for (size_t channel=0; channel<numChannels; channel++) {

        uint64_t sum = 0;
        size_t minValue = 65535;
        size_t maxValue = 0;
        std::vector<uint64_t>& histogram = histograms[channel];
        const T* texture = textureArray[channel];

        for (size_t z = 0; z < inTextureSize.z; ++z) {
            for (size_t y = 0; y < inTextureSize.y; ++y) {
                const tgt::svec3 textureCoord = tgt::svec3(0,y,z)+brickOffsetInTexture;
                size_t textureLinearCoord = cubicCoordToLinear(textureCoord, textureDim);
                size_t textureLinearCoordEnd = textureLinearCoord+inTextureSize.x;
                size_t brickLinearCoord = numChannels * brickDim.x * (brickDim.y * z + y) + channel;

                for (; textureLinearCoord < textureLinearCoordEnd; ++textureLinearCoord) {
                    tgtAssert(brickLinearCoord < brickBufferSize, "invalid brick linear coord");

                    tgtAssert(textureLinearCoord < textureBufferSize, "invalid texture linear coord");

                    size_t value = convertVoxelValueToUInt16<T>(texture[textureLinearCoord]);

                    brickBuffer[brickLinearCoord] = value;
                    sum += value;
                    minValue = std::min(minValue, value);
                    maxValue = std::max(maxValue, value);

                    // 16 bit values, but only NUM_HISTOGRAM_BUCKETS_BITS bits for histogram
                    size_t rightShiftAmt = 16 - NUM_HISTOGRAM_BUCKETS_BITS;
                    size_t histValue = value >> rightShiftAmt;
                    histogram[histValue]++;

                    brickLinearCoord += numChannels;
                }
            }
        }

        avgValues[channel] = static_cast<uint16_t>(sum / numSignificantBrickVoxels);
        maxValues[channel] = maxValue;
        minValues[channel] = minValue;
        tgtAssert(minValues[channel] <= maxValues[channel], "min value is larger than max value");
        tgtAssert(minValues[channel] <= avgValues[channel] && avgValues[channel] <= maxValues[channel], "avg value outside min-max value range");
    }

    // Fill brick outside of texture with zeros
    tgt::svec3 end = brickDim;
    for(int dim=2; dim >= 0; --dim) {
        tgt::svec3 start(0);
        start[dim] = inTextureSize[dim];
        VRN_FOR_EACH_VOXEL(brickVoxel, start, end) {
            for (size_t channel=0; channel<numChannels; channel++) {
                size_t brickLinearCoord = cubicCoordToLinear(brickVoxel, brickDim)*numChannels;
                brickBuffer[brickLinearCoord + channel] = 0;
            }
        }
        end[dim] = start[dim];
    }

    // Initialize metadata if brick if empty
    if (numSignificantBrickVoxels == 0) {
        for (size_t channel=0; channel<numChannels; channel++) {
            avgValues[channel] = 0;
            minValues[channel] = 0;
            maxValues[channel] = 0;
        }
    }
}
VolumeOctreeNode* VolumeOctree::createParentNode(VolumeOctreeNode* children[8], bool octreeOptimization, uint16_t homogeneityThreshold,
    const tgt::svec3& brickUrb, uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues, HalfSampleAggregateFunction halfSampleFn) {
    switch(halfSampleFn) {
    case MEAN: return createParentNodeWithHalfsampling<HalfSampleMean>(children, octreeOptimization, homogeneityThreshold,
                brickUrb, avgValues, minValues, maxValues);
    case MAX: return createParentNodeWithHalfsampling<HalfSampleMax>(children, octreeOptimization, homogeneityThreshold,
                brickUrb, avgValues, minValues, maxValues);
    case MIN: return createParentNodeWithHalfsampling<HalfSampleMin>(children, octreeOptimization, homogeneityThreshold,
                brickUrb, avgValues, minValues, maxValues);
    case MEDIAN: return createParentNodeWithHalfsampling<HalfSampleMedian>(children, octreeOptimization, homogeneityThreshold,
                brickUrb, avgValues, minValues, maxValues);
    default:
        tgtAssert(false, "Invalid half sample mode");
        return nullptr;
    }
}

template<typename HalfSample>
VolumeOctreeNode* VolumeOctree::createParentNodeWithHalfsampling(VolumeOctreeNode* children[8], bool octreeOptimization, uint16_t homogeneityThreshold,
    const tgt::svec3& brickUrb, uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues) {
    switch(getNumChannels()) {
    case 1: return createParentNodeConstChannels<1, HalfSample>(children, octreeOptimization, homogeneityThreshold,
                brickUrb, avgValues, minValues, maxValues);
    case 2: return createParentNodeConstChannels<2, HalfSample>(children, octreeOptimization, homogeneityThreshold,
                brickUrb, avgValues, minValues, maxValues);
    case 3: return createParentNodeConstChannels<3, HalfSample>(children, octreeOptimization, homogeneityThreshold,
                brickUrb, avgValues, minValues, maxValues);
    case 4: return createParentNodeConstChannels<4, HalfSample>(children, octreeOptimization, homogeneityThreshold,
                brickUrb, avgValues, minValues, maxValues);
    default:
        tgtAssert(false, "more than 4 channels");
        return nullptr;
    }
}


template<size_t numChannels, typename HalfSample>
VolumeOctreeNode* VolumeOctree::createParentNodeConstChannels(VolumeOctreeNode* children[8], bool octreeOptimization, uint16_t homogeneityThreshold,
    const tgt::svec3& brickUrb, uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues) {
    tgtAssert(brickPoolManager_, "no brick pool manager");

    tgtAssert(getNumChannels() <= MAX_CHANNELS, "more than max channels");
    tgtAssert(getNumChannels() == numChannels, "Invalid number of channels");

    // compute parent avg/min/max values from children and determine whether parent is homogeneous
    uint64_t avgValuesUInt64[MAX_CHANNELS];
    for (size_t ch=0; ch<getNumChannels(); ch++) {
        avgValuesUInt64[ch] = 0;
        minValues[ch] = std::numeric_limits<uint16_t>::max();
        maxValues[ch] = 0;
    }
    size_t numNonEmptyChildren = 0; //< number of child nodes not completely outside volume
    for (size_t childID=0; childID<8; childID++) {
        VolumeOctreeNode* child = children[childID];
        tgtAssert(child, "null pointer");
        // ignore empty child node (completely outside volume) for avg/min/max calculation
        if (!child->inVolume()) {
            continue;
        }
        else {
            numNonEmptyChildren++;
            for (size_t ch=0; ch<getNumChannels(); ch++) {
                avgValuesUInt64[ch] += child->getAvgValue(ch);
                minValues[ch] = std::min(minValues[ch], child->getMinValue(ch));
                maxValues[ch] = std::max(maxValues[ch], child->getMaxValue(ch));
            }
        }
    }
    bool homogeneous = true;
    if (numNonEmptyChildren > 0) {
        for (size_t ch=0; ch<getNumChannels(); ch++) {
            avgValues[ch] = static_cast<uint16_t>(avgValuesUInt64[ch] / numNonEmptyChildren);
            tgtAssert(minValues[ch] <= avgValues[ch] && avgValues[ch] <= maxValues[ch], "invalid avg/min/max values");
            homogeneous &= ((maxValues[ch] - minValues[ch]) <= homogeneityThreshold);
        }
    }
    else {
        for (size_t ch=0; ch<getNumChannels(); ch++) {
            minValues[ch] = 0;
            maxValues[ch] = 0;
            avgValues[ch] = 0;
        }
    }

    if (homogeneous && octreeOptimization || numNonEmptyChildren == 0) { // node is homogeneous => create leaf node without brick
        VolumeOctreeNode* parent = 0;
        if (numNonEmptyChildren > 0)
            parent = VolumeOctreeBase::createNode(getNumChannels(), avgValues, minValues, maxValues);
        else
            parent = VolumeOctreeBase::createNode(getNumChannels());

        // delete child nodes
        for (size_t i=0; i<8; i++)
            deleteSubTree(children[i]);

        return parent;

    }
    else { // node is not homogeneous => create inner node with brick by merging child nodes

        // half sample child bricks
        svec3 brickDim = getBrickDim();
        svec3 halfBrickDim = getBrickDim() / svec3(2);

        // copy halfsampled bricks to dest buffer (child nodes/halfsampled bricks are expected to be zyx ordered)
        uint64_t brickVirtualMemoryAddress = brickPoolManager_->allocateBrick();

        VRN_FOR_EACH_VOXEL(childPos, svec3::zero, svec3::two) {
            // Reaquire brick buffer for every child in order to work around deadlock due to bad
            // LRU-caching behaviour of BrickPoolManagerDisk:
            uint16_t* brickBuffer = brickPoolManager_->getWritableBrick(brickVirtualMemoryAddress);

            for(size_t channel = 0; channel < numChannels; ++channel) {
                const size_t childID = cubicCoordToLinear(childPos, svec3::two);
                VolumeOctreeNode* child = children[childID];
                tgtAssert(child, "null pointer");

                tgt::svec3 inParentOffset = childPos * halfBrickDim;

                if (child->hasBrick()) { // node brick present => halfsample into buffer
                    tgt::svec3 llf = tgt::min(childPos * brickDim, brickUrb);
                    tgt::svec3 urb = tgt::min(llf + brickDim, brickUrb);
                    tgt::svec3 childUrb = urb - llf;
                    size_t childBrickAddr = child->getBrickAddress();

                    const uint16_t* childBrick = brickPoolManager_->getBrick(childBrickAddr);

                    tgt::svec3 maxPos = childUrb-tgt::svec3::one;
                    for (size_t z = 0; z < brickDim.z; z+=2) {
                        size_t zl = std::min(z  , maxPos.z) * numChannels;
                        size_t zh = std::min(z+1, maxPos.z) * numChannels;

                        size_t zp = (inParentOffset.z + z/2) * numChannels;
                        for (size_t y = 0; y < brickDim.y; y+=2) {
                            size_t yl = std::min(y  , maxPos.y) * numChannels;
                            size_t yh = std::min(y+1, maxPos.y) * numChannels;

                            size_t yp = (inParentOffset.y + y/2) * numChannels;

                            size_t yzll = brickDim.x * (yl + (brickDim.y * zl)) + channel;
                            size_t yzlh = brickDim.x * (yl + (brickDim.y * zh)) + channel;
                            size_t yzhl = brickDim.x * (yh + (brickDim.y * zl)) + channel;
                            size_t yzhh = brickDim.x * (yh + (brickDim.y * zh)) + channel;

                            size_t yzp = brickDim.x * (yp + (brickDim.y * zp)) + channel;

                            size_t x = 0;

                            auto halfSampleAt = [&] (size_t xl, size_t xh) {
                                uint16_t halfValue = HalfSample::halfsample(
                                    childBrick[xl+yzll],
                                    childBrick[xh+yzll],
                                    childBrick[xl+yzhl],
                                    childBrick[xh+yzhl],
                                    childBrick[xl+yzlh],
                                    childBrick[xh+yzlh],
                                    childBrick[xl+yzhh],
                                    childBrick[xh+yzhh]
                                );

                                size_t xp = (inParentOffset.x + x/2) * numChannels;

                                size_t parentBrickLinearPos = xp + yzp;
                                brickBuffer[parentBrickLinearPos] = halfValue;
                            };
                            for (; x < maxPos.x; x+=2) {
                                size_t xl =  x    * numChannels;
                                size_t xh = (x+1) * numChannels;

                                halfSampleAt(xl, xh);
                            }
                            for (; x < brickDim.x; x+=2) {
                                size_t xl = maxPos.x * numChannels;
                                size_t xh = maxPos.x * numChannels;

                                halfSampleAt(xl, xh);
                            }
                        }
                    }

                    brickPoolManager_->releaseBrick(childBrickAddr);
                }
                else { // no brick present => use node's avg values
                    uint16_t childAvgValue = child->getAvgValue(channel);
                    VRN_FOR_EACH_VOXEL(halfPos, svec3::zero, halfBrickDim) {
                        tgt::svec3 parentPos = inParentOffset + halfPos;
                        size_t parentBrickLinearPos = cubicCoordToLinear(parentPos, brickDim)*numChannels + channel;
                        brickBuffer[parentBrickLinearPos] = childAvgValue;
                    }
                }
            }
            brickPoolManager_->releaseBrick(brickVirtualMemoryAddress, OctreeBrickPoolManagerBase::WRITE);
        }

        // create parent node
        VolumeOctreeNode* parent = VolumeOctreeBase::createNode(getNumChannels(), avgValues, minValues, maxValues,
            brickVirtualMemoryAddress, children);

        return parent;
    }
}

void VolumeOctree::copyBrickToTexture(const uint16_t* brick, const tgt::svec3& brickDim,
    uint16_t* texture, const tgt::svec3& textureDim, const tgt::svec3& brickOffsetInTexture) const
{
    tgtAssert(brick, "null pointer passed");
    //tgtAssert(channel < getNumChannels(), "invalid channel");
    tgtAssert(texture, "null pointer passed");
    tgtAssert(tgt::hand(tgt::lessThanEqual(brickDim, getOctreeDim())), "brick dimensions greater than texture dimensions");
    tgtAssert(tgt::hand(tgt::lessThanEqual(brickOffsetInTexture+brickDim, getOctreeDim())), "brick (partially) outside texture");

    const size_t numChannels = getNumChannels();
    const size_t brickChannelSize = tgt::hmul(brickDim);
    const size_t brickBufferSize = brickChannelSize*getNumChannels();
    const size_t textureBufferSize = tgt::hmul(textureDim)*getNumChannels();

    VRN_FOR_EACH_VOXEL(brickVoxel, svec3(0,0,0), brickDim) {
        size_t brickLinearCoord = cubicCoordToLinear(brickVoxel, brickDim)*numChannels;
        tgtAssert(brickLinearCoord+numChannels-1 < brickBufferSize, "invalid brick linear coord");
        const tgt::svec3 textureCoord = brickVoxel+brickOffsetInTexture;
        if (tgt::hand(tgt::lessThan(textureCoord, textureDim))) {
            size_t textureLinearCoord = cubicCoordToLinear(textureCoord, textureDim)*numChannels;
            tgtAssert(textureLinearCoord+numChannels-1 < textureBufferSize, "invalid texture linear coord");
            std::copy(brick+brickLinearCoord, brick+brickLinearCoord+numChannels, texture+textureLinearCoord);
        }
    }
}

void VolumeOctree::serialize(Serializer& s) const {
    if (!brickPoolManager_)
        throw SerializationException("Unable to serialize octree: no brick pool manager assigned");

    // determine output path for node buffer
    const std::string octreeFile = s.getDocumentPath();
    const std::string octreePath = tgt::FileSystem::dirName(octreeFile);
    if (octreeFile.empty() || !tgt::FileSystem::dirExists(octreePath))
        throw SerializationException("Octree path does not exist: " + octreePath);

    // serialize base octree
    VolumeOctreeBase::serialize(s);

    // serialize tree to binary buffer
    char* nodeBuffer = 0;
    size_t bufferSize = 0;
    serializeNodeBuffer(nodeBuffer, bufferSize);
    tgtAssert(nodeBuffer, "null pointer returned");
    tgtAssert(bufferSize, "invalid buffer size returned");

    // write node buffer to file
    const std::string bufferFile = tgt::FileSystem::cleanupPath(octreePath + "/" + NODE_BUFFER_FILE_NAME);
    std::fstream fileStream(bufferFile.c_str(), std::ios_base::out | std::ios_base::binary);
    if (fileStream.fail()) {
        delete[] nodeBuffer;
        throw SerializationException("Failed to open file '" + bufferFile + "' for writing");
    }
    try {
        fileStream.write(nodeBuffer, bufferSize);
    }
    catch (std::exception& e) {
        delete[] nodeBuffer;
        fileStream.close();
        throw SerializationException("Failed to write node buffer to file '" + bufferFile + "': " + std::string(e.what()));
    }
    fileStream.close();
    delete[] nodeBuffer;
    nodeBuffer = 0;

    // serialize node buffer meta information
    s.serialize("nodeBufferName", std::string(NODE_BUFFER_FILE_NAME));
    s.serialize("nodeCount", rootNode_->getNodeCount());
    s.serialize("nodeBufferSize", bufferSize);

    // serialize brick pool manager
    s.serialize("brickPoolManager", brickPoolManager_);

    // serialize histograms
    s.serialize("histograms", histograms_);

}

void VolumeOctree::deserialize(Deserializer& s) {
    // determine output path for node buffer
    const std::string octreeFile = s.getDocumentPath();
    const std::string octreePath = tgt::FileSystem::dirName(octreeFile);
    if (octreeFile.empty() || !tgt::FileSystem::dirExists(octreePath))
        throw SerializationException("Octree path does not exist: " + octreePath);

    // deserialize base octree properties
    VolumeOctreeBase::deserialize(s);

    // deserialize node buffer meta information
    size_t nodeCount, nodeBufferSize;
    s.deserialize("nodeCount", nodeCount);
    if (nodeCount == 0)
        throw SerializationException("Invalid node count: " + itos(nodeCount));
    s.deserialize("nodeBufferSize", nodeBufferSize);
    if (nodeBufferSize == 0 || nodeBufferSize < nodeCount)
        throw SerializationException("Invalid node buffer size: " + itos(nodeBufferSize));
    std::string nodeBufferName;
    s.optionalDeserialize<std::string>("nodeBufferName", nodeBufferName, NODE_BUFFER_FILE_NAME);
    // load binary node buffer from file
    const std::string bufferFile = tgt::FileSystem::cleanupPath(octreePath + "/" + nodeBufferName);
    std::fstream fileStream(bufferFile.c_str(), std::ios_base::in | std::ios_base::binary);
    if (fileStream.fail())
        throw SerializationException("Failed to open node buffer file '" + bufferFile + "' for reading");
    char* nodeBuffer = new char[nodeBufferSize];
    try {
        fileStream.read(nodeBuffer, nodeBufferSize);
    }
    catch (std::exception& e) {
        delete[] nodeBuffer;
        fileStream.close();
        throw SerializationException("Failed to read node buffer from file '" + bufferFile + "': " + std::string(e.what()));
    }
    fileStream.close();

    // construct tree nodes from node buffer
    if (rootNode_)
        deleteSubTree(rootNode_);
    rootNode_ = 0;
    try {
        tgt::Stopwatch watch;
        watch.start();
        rootNode_ = deserializeNodeBuffer(nodeBuffer, nodeCount, nodeBufferSize);
        LDEBUG("Node buffer deserialization time: " << watch.getRuntime() << " msec");
    }
    catch (std::exception& e) {
        delete[] nodeBuffer;
        throw SerializationException("Failed to deserialize binary node buffer '" + bufferFile + "': " + std::string(e.what()));
    }
    tgtAssert(rootNode_, "no root node"); //< exception expected from deserializeNodeBuffer
    updateTreeMetaDataCache();
    delete[] nodeBuffer;
    nodeBuffer = 0;

    // deserialize brick pool manager
    delete brickPoolManager_;
    brickPoolManager_ = 0;
    tgt::Stopwatch watch;
    watch.start();
    s.deserialize("brickPoolManager", brickPoolManager_);
    LDEBUG("Brick pool manager deserialization time: " << watch.getRuntime() << " msec");
    if (!brickPoolManager_) {
        deleteSubTree(rootNode_);
        rootNode_ = 0;
        throw SerializationException("Brick pool manager not deserialized");
    }

    // deserialize histograms
    s.deserialize("histograms", histograms_);
    if (histograms_.size() != getNumChannels()) {
        deleteSubTree(rootNode_);
        rootNode_ = 0;
        throw SerializationException("Number of deserialized histograms does not match channel count [" +
            itos(histograms_.size()) + " != " + itos(getNumChannels()));
    }

}

void VolumeOctree::serializeNodeBuffer(char*& binaryBuffer, size_t& bufferSize) const {
    tgtAssert(rootNode_, "no root node");

    const size_t NODE_CONTENT_SIZE = rootNode_->getContentSize();
    const size_t NODE_SIZE = NODE_CONTENT_SIZE + sizeof(uint64_t); // content size + child group offset
    const size_t NODE_COUNT = rootNode_->getNodeCount();

    bufferSize = NODE_COUNT*NODE_SIZE;
    binaryBuffer = new char[bufferSize];

    // pair consisting of a octree node whose children still have to be added to the buffer, and the node's buffer offset
    typedef std::pair<const VolumeOctreeNode*, size_t> QuededNode;
    std::stack<QuededNode> workQueue;

    // start with root node, put encountered nodes into fifo queue, iterate until queue is empty
    rootNode_->serializeContentToBinaryBuffer(binaryBuffer);
    *reinterpret_cast<uint64_t*>(binaryBuffer + NODE_CONTENT_SIZE) = std::numeric_limits<uint64_t>::max(); //add max to chid group
    workQueue.push(QuededNode(rootNode_, 0));
    size_t curBufferOffset = 1; //< next after root node
    while (!workQueue.empty()) {
        // retrieve next node to process
        const VolumeOctreeNode* curNode = workQueue.top().first;
        size_t nodeOffset = workQueue.top().second;
        workQueue.pop();

        // no children => nothing to do
        if (!curNode->children_[0])
            continue;

        tgtAssert(nodeOffset < curBufferOffset,  "node offset not less than current offset");
        tgtAssert(curBufferOffset+8 <= NODE_COUNT, "invalid current buffer offset");

        // set curNode's child group pointer to current offset
        size_t nodeByteOffset = nodeOffset*NODE_SIZE;
        tgtAssert(nodeByteOffset < bufferSize, "invalid byte offset");
        *reinterpret_cast<uint64_t*>(binaryBuffer + nodeByteOffset + NODE_CONTENT_SIZE) = curBufferOffset;

        // create eight adjacent buffer entries at curOffset for children, and add them to work queue
        for (size_t childID = 0; childID < 8; childID++) {
            const VolumeOctreeNode* child = curNode->children_[childID];
            tgtAssert(child, "missing child");
            size_t childOffset = curBufferOffset+childID;

            size_t childByteOffset = childOffset*NODE_SIZE;
            tgtAssert(childByteOffset < bufferSize, "invalid child byte offset");
            child->serializeContentToBinaryBuffer(binaryBuffer + childByteOffset);
            // assign max uint64_t as child group offset to nodes without children
            *reinterpret_cast<uint64_t*>(binaryBuffer + childByteOffset + NODE_CONTENT_SIZE) = std::numeric_limits<uint64_t>::max();

            workQueue.push(QuededNode(child, childOffset));
        }
        curBufferOffset += 8;
    }
    tgtAssert(curBufferOffset == NODE_COUNT, "buffer offset does not equal number of tree nodes");

    // validate result
    /*LINFO("Validating node buffer against octree...");
    watch.reset();
    watch.start();
    try {
        compareNodeToBuffer(rootNode, 0);
    }
    catch (VoreenException& e) {
        LERROR(e.what());
    }
    LINFO("Validation time: " << watch.getRuntime() << " ms");*/
}

VolumeOctreeNode* VolumeOctree::deserializeNodeBuffer(const char* binaryBuffer, const size_t nodeCount, const size_t bufferSize) {
    tgtAssert(binaryBuffer, "null pointer passed");
    tgtAssert(nodeCount > 0, "invalid node count");
    tgtAssert(bufferSize > 0, "invalid buffer size");

    VolumeOctreeNode* rootNode = VolumeOctreeBase::createNode(getNumChannels());
    const size_t NODE_CONTENT_SIZE = rootNode->getContentSize();
    const size_t NODE_SIZE = NODE_CONTENT_SIZE + sizeof(uint64_t); // content size + child group offset

    // check buffer size against node count
    if (bufferSize != nodeCount*NODE_SIZE) {
        delete rootNode;
        throw SerializationException("Node buffer byte size does not match nodeCount*numBytesPerNode [" +
                                      itos(bufferSize) + " != " + itos(nodeCount) + "*" + itos(NODE_SIZE) + "]");
    }

    // pair consisting of a octree node whose children still have to be read from the buffer, and the children group's buffer offset
    typedef std::pair<VolumeOctreeNode*, uint64_t> QuededNode;
    std::stack<QuededNode> workQueue;

    // create root node from first buffer entry
    rootNode->deserializeContentFromBinaryBuffer(binaryBuffer);
    uint64_t childGroupOffset = *reinterpret_cast<const uint64_t*>(binaryBuffer + NODE_CONTENT_SIZE);
    if (childGroupOffset < std::numeric_limits<uint64_t>::max()) {
        if (childGroupOffset >= bufferSize-8*NODE_SIZE) {
            delete rootNode;
            throw SerializationException("Invalid child group offset: " + itos((size_t)childGroupOffset));
        }
        workQueue.push(std::pair<VolumeOctreeNode*, uint64_t>(rootNode, childGroupOffset));
    }

    // process work queue
    while (!workQueue.empty()) {
        // retrieve next node/child group offset to process
        VolumeOctreeNode* curNode = workQueue.top().first;
        uint64_t childGroupOffset = workQueue.top().second;
        tgtAssert(childGroupOffset < bufferSize-8*NODE_SIZE, "invalid child group offset");
        workQueue.pop();

        // iterate over child group and create nodes from buffer entries
        for (size_t i=0; i<8; i++) {
            // create child node
            VolumeOctreeNode* childNode = VolumeOctreeBase::createNode(getNumChannels());
            const char* childBuffer = binaryBuffer + (childGroupOffset+i)*NODE_SIZE;
            childNode->deserializeContentFromBinaryBuffer(childBuffer);

            // retrieve child group offset and add to work queue, if children present
            uint64_t grandChildGroupOffset = *reinterpret_cast<const uint64_t*>(childBuffer + NODE_CONTENT_SIZE);
            if (grandChildGroupOffset < std::numeric_limits<uint64_t>::max()) {
                if (grandChildGroupOffset >= bufferSize-8*NODE_SIZE) {
                    deleteSubTree(rootNode);
                    throw SerializationException("Invalid child group offset: " + itos((size_t)grandChildGroupOffset));
                }
                workQueue.push(std::pair<VolumeOctreeNode*, uint64_t>(childNode, grandChildGroupOffset));
            }

            curNode->children_[i] = childNode;
        }
    }
    tgtAssert(rootNode, "no root node");
    if (rootNode->getNodeCount() != nodeCount) {
        deleteSubTree(rootNode);
        throw SerializationException("Node count of deserialized octree does not match specified node count [" +
                                      itos(rootNode->getNodeCount()) + " != " + itos(nodeCount) + "]");
    }

    return rootNode;
}

const OctreeBrickPoolManagerBase* VolumeOctree::getBrickPoolManager() const {
    return brickPoolManager_;
}
OctreeBrickPoolManagerBase* VolumeOctree::getBrickPoolManager() {
    return brickPoolManager_;
}

uint16_t* VolumeOctree::acquireTempBrickBuffer() {
    return new uint16_t[getBrickMemorySize()/2];
}

void VolumeOctree::releaseTempBrickBuffer(uint16_t* buffer) {
    tgtAssert(buffer, "null pointer passed");

    delete[] buffer;
}



} // namespace
