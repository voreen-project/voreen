/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEOCTREE_H
#define VRN_VOLUMEOCTREE_H

#include "volumeoctreebase.h"
#include "octreebrickpoolmanager.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/io/progressreporter.h"
#include "voreen/core/utils/exception.h"

#include "tgt/stopwatch.h"

#include <vector>
#include <string>

namespace voreen {

/**
 * Basic multi-channel octree implementation that creates the octree from one or multiple channel volumes.
 *
 * The octree can be optimized: Bricks or subtrees that represent a volume region whose
 * max-min difference does not exceed a specified threshold are discarded.
 */
class VRN_CORE_API VolumeOctree : public VolumeOctreeBase {

public:

    /**
     * Standard constructor for the construction of a multi-channel octree.
     *
     * @param channelVolumes Input volumes the octree is constructed from. Up to four channels are supported.
     *        All channel volumes must have the same dimensions and format.
     * @param brickDim dimensions of one brick (number of voxels per dimension). Must be power-of-two.
     * @param homogeneityThreshold All regions within the volume with a value range (i.e., max-min)
     *        less or equal this threshold are considered homogeneous and will be discarded in the octree.
     *        The threshold is normalized, i.e., 1.0 represents the data type's full value range.
     *        A negative threshold disables octree optimization, resulting in a complete tree.
     * @param brickPoolManager Mandatory helper class that organizes the bricks in RAM/disk memory.
     *        The octree takes ownership of the passed manager and deletes it on its own destruction.
     * @param numThreads Number of threads to use during octree construction. Note: multi-threading requires the OpenMP module.
     * @param progressReporter Optional progress reporter that is updated during tree construction.
     *
     * @throws std::exception If the octree construction fails.
     */
    VolumeOctree(const std::vector<const VolumeBase*>& channelVolumes, size_t brickDim, float homogeneityThreshold = 0.01f,
            OctreeBrickPoolManagerBase* brickPoolManager = new OctreeBrickPoolManagerRAM(),
            size_t numThreads = 1, ProgressReporter* progessReporter = 0);

    /**
     * Convenience constructor for a single-channel octree.
     */
    VolumeOctree(const VolumeBase* volume, size_t brickDim, float homogeneityThreshold = 0.01f,
            OctreeBrickPoolManagerBase* brickPoolManager = new OctreeBrickPoolManagerRAM(),
            size_t numThreads = 1, ProgressReporter* progessReporter = 0);

    /**
     * Construct a VolumeOctree from preprocessed parts, i.e., an existing hierarchy of nodes whose bricks are stored
     * in the passed brickPoolManager.
     */
    VolumeOctree(VolumeOctreeNode* root, OctreeBrickPoolManagerBase* brickPoolManager, const tgt::svec3& brickDim,
            const tgt::svec3& volumeDim, size_t numChannels);

    /**
     * Decompose the VolumeOctree into its components. The caller takes ownership of the components.
     */
    std::pair<OctreeBrickPoolManagerBase*, VolumeOctreeNode*> decompose() &&;

    /// Default constructor for serialization only.
    VolumeOctree();
    VolumeOctree(VolumeOctree&& other);
    VolumeOctree& operator=(VolumeOctree&& other);

    virtual ~VolumeOctree();
    virtual VolumeOctree* create() const;

    virtual std::string getClassName() const { return "VolumeOctree"; }

    virtual size_t getActualTreeDepth() const;

    virtual size_t getNumNodes() const;

    virtual size_t getNumBricks() const;

    virtual uint64_t getBrickPoolMemoryAllocated() const;

    virtual uint64_t getBrickPoolMemoryUsed() const;

    virtual std::string getDescription() const;

    virtual const Histogram1D* getHistogram(size_t channel = 0) const;

    const OctreeBrickPoolManagerBase* getBrickPoolManager() const;


    virtual uint16_t getVoxel(const tgt::svec3& pos, size_t channel = 0, size_t level = 0) const;

    virtual const VolumeOctreeNode* getRootNode() const;

    virtual const VolumeOctreeNode* getNode(const tgt::vec3& point, size_t& level,
        tgt::svec3& voxelLLF, tgt::svec3& voxelURB, tgt::vec3& normLLF, tgt::vec3& normURB) const;

    virtual const uint16_t* getNodeBrick(const VolumeOctreeNode* node) const;

    virtual void releaseNodeBrick(const VolumeOctreeNode* node) const;

    virtual VolumeRAM* createVolume(size_t level = 0, clock_t timeLimit = 0, bool* complete = 0) const;

    virtual VolumeRAM* createSlice(SliceAlignment sliceAlignment, size_t sliceIndex, size_t level = 0,
        clock_t timeLimit = 0, bool* complete = 0, tgt::svec3 begin = tgt::svec3::zero, tgt::svec3 end = tgt::svec3(-1)) const;


    virtual void serialize(Serializer& s) const;

    virtual void deserialize(Deserializer& s);

protected:
    static const std::string loggerCat_;

private:
    void buildOctreeIteratively(const std::vector<const VolumeBase*>& volumes,
        bool octreeOptimization, uint16_t homogeneityThreshold, size_t numThreads, ProgressReporter* progessReporter);

    const VolumeOctreeNode* getNodeAtVoxel(const tgt::svec3& voxel, const size_t curLevel, const size_t targetLevel,
        const VolumeOctreeNode* node, const tgt::svec3& nodeLlf, const tgt::svec3& nodeUrb,
        size_t& resultLevel, tgt::svec3& resultLlf, tgt::svec3& resultUrb) const;

    void composeNodeTexture(const VolumeOctreeNode* node, const tgt::svec3& nodeOffset, size_t curLevel, size_t targetLevel,
        uint16_t* textureBuffer, const tgt::svec3& textureDim, clock_t timeLimit, tgt::Stopwatch& runtimeWatch, bool& complete) const;

    void composeNodeSliceTexture(SliceAlignment sliceAlignment, const VolumeOctreeNode* node,
        const tgt::svec3& nodeOffsetInTexture, size_t sliceIndexInNode, size_t curLevel, size_t targetLevel,
        uint16_t* textureBuffer, const tgt::svec3& textureDim, clock_t timeLimit, tgt::Stopwatch& runtimeWatch,
        bool& complete, const tgt::svec3 begin, const tgt::svec3 end) const;

    // low-level helper functions
    template<class T>
    VolumeOctreeNode* createTreeNodeFromTexture(const tgt::svec3& llf, const tgt::svec3& urb,
        const std::vector<const void*>& textureBuffers, const tgt::svec3& textureDim,
        bool octreeOptimization, uint16_t homogeneityThreshold,
        uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues,
        std::vector< std::vector<uint64_t> >& histograms);

    template<class T>
    void extractBrickFromTexture(const std::vector<const void*>& textures, const tgt::svec3& textureDim,
        uint16_t* brickBuffer, const tgt::svec3& brickDim, const tgt::svec3& brickOffsetInTexture,
        uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues,
        std::vector< std::vector<uint64_t> >& histograms) const;

    VolumeOctreeNode* createParentNode(VolumeOctreeNode* children[8], bool octreeOptimization, uint16_t homogeneityThreshold,
        uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues);

    void copyBrickToTexture(const uint16_t* brick, const tgt::svec3& brickDim,
        uint16_t* texture, const tgt::svec3& textureDim,
        const tgt::svec3& brickOffsetInTexture) const;

    void halfSampleBrick(const uint16_t* brick, const tgt::svec3& brickDim, uint16_t* halfSampledBrick) const;

    void deleteSubTree(VolumeOctreeNode* root) const;

    void serializeNodeBuffer(char*& binaryBuffer, size_t& bufferSize) const;

    VolumeOctreeNode* deserializeNodeBuffer(const char* binaryBuffer, const size_t nodeCount, const size_t bufferSize);

    uint16_t* acquireTempBrickBuffer();
    void releaseTempBrickBuffer(uint16_t* buffer);

    VolumeOctreeNode* rootNode_;

    OctreeBrickPoolManagerBase* brickPoolManager_;

    std::vector<Histogram1D*> histograms_;
};

} // namespace

#endif
