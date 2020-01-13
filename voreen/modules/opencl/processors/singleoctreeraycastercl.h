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

#ifndef VRN_SINGLEOCTREERAYCASTERCL_H
#define VRN_SINGLEOCTREERAYCASTERCL_H

#include "voreen/core/processors/renderprocessor.h"

#include "../utils/clwrapper.h"
#include "../properties/openclproperty.h"
#include "modules/opencl/processors/openclprocessor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/collectivesettingsproperty.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/utils/voreenqualitymode.h"

#include "tgt/texture.h"
#include "tgt/stopwatch.h"

#include <deque>

namespace voreen {

class VolumeOctreeBase;
class VolumeOctreeNode;
class MeshListGeometry;
class MeshGeometry;
class TransFunc1DKeys;

class VRN_CORE_API SingleOctreeRaycasterCL : public cl::OpenCLProcessor<RenderProcessor>, public QualityModeObserver {
public:
    SingleOctreeRaycasterCL();
    virtual ~SingleOctreeRaycasterCL();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "SingleOctreeRaycasterCL"; }
    virtual std::string getCategory() const   { return "Raycasting";              }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;        }

    virtual bool isReady() const;

    virtual void initializeCL();
    virtual void deinitializeCL();

protected:
    virtual void setDescriptions() {
        setDescription("Performs an out-of-core volume rendering using OpenCL. "
                       "The implementation is based on the GigaVoxel technique [Crassin et al., SIGGRAPH 2009]. "
                       "<p><strong>NOTE:</strong> This renderer requires an input volume with an octree representation, which "
                                                  "has to created using the OctreeCreator.</p>");
    }

    virtual void process();
    virtual void beforeProcess();
    virtual void afterProcess();

    virtual void initialize();
    virtual void deinitialize();

    /// Invalidates the processor when it switches out of interaction mode.
    virtual void qualityModeChanged();

protected:
    enum DisplayMode {
        FULL_FRAME = 0, ///< Current frame is completely rendered, using higher level bricks if necessary (interaction mode)
        REFINEMENT = 1  ///< Rendering is iteratively refined at each frame, using full resolution bricks
    };

    /// Consecutive states of the render process.
    enum RenderState {
        RENDER_STATE_FIRST_INTERACTION_FRAME = 0,   ///< first interactive frame
        RENDER_STATE_INTERACTION_FRAME = 1,         ///< interaction is in progress
        RENDER_STATE_FIRST_FULL_FRAME = 2,          ///< first full-quality frame after switching out of interaction mode
        RENDER_STATE_FULL_FRAME = 3,                ///< full-frame refinement, but first full-frame already rendered
        RENDER_STATE_FIRST_ITERATIVE_FRAME = 4,     ///< first iterative frame (starts after full-frame refinement)
        RENDER_STATE_ITERATIVE_FRAME = 5,           ///< iterative refinement, but first iterative frame already rendered
        RENDER_STATE_FINISHED = 6                   ///< current frame is finished (last iterative or full-frame has been rendered)
    };

    void rebuildKernels(const VolumeOctreeBase* octree);

    /**
     * Copies over the content of the srcPort to the destPort,
     * thereby implicitly rescaling the image to the dest dimensions.
     * To be used by subclasses for implementing coarseness (i.e., rendering with reduced dimensions in interaction mode).
     *
     * The entry/exit point ports are used for transforming the srcPortDepth values that are in ray coordinates (linear [0;1]).
     */
    void rescaleRendering(RenderPort& srcPort, RenderPort& srcPortDepth,
        RenderPort& entryPoints, RenderPort& exitPoints, RenderPort& destPort);

    virtual void initializeNodeBuffer();
    void initializeBrickBuffer();

    void initializeRayBuffer();

    void updateBrickBuffer(int keepLevel, size_t& numUsed, size_t& numRequested, size_t& numUploaded);

    void analyzeBrickFlagBuffer(size_t& numTraversed, size_t& numUsed, size_t& numRequested) const;
    void clearBrickFlagBuffer(size_t numNodes, bool clearMIPUsed);

    void renderFrame(DisplayMode displayMode, const VolumeBase* volume, const VolumeOctreeBase* octree, float samplingStepSize, size_t nodeLevel, bool useBricks);

    void generateTransfuncTexture();

    void invalidateNodeBuffer();
    void invalidateBrickBuffer();
    void invalidateTransFuncTex();
    void resetRenderState();

    void updatePropertyConfiguration();

    void resetChannelShift();

    void compareNodeToBuffer(const VolumeOctreeNode* node, size_t nodeOffset) const;

    // ports
    VolumePort volumeInport_;
    RenderPort entryPointsInport_;
    RenderPort exitPointsInport_;

    RenderPort renderOutport_;
    RenderPort internalRenderPort_;         ///< used for coarseness (rendering with reduced resolution in interaction mode).
    RenderPort internalRenderPortDepth_;    ///< necessary because OpenGL-OpenCL sharing does not work with depth textures

    // properties
    StringOptionProperty compositingMode_;
    CollectiveSettingsProperty renderingQuality_;
    CollectiveSettingsProperty interactionQuality_;

    TransFunc1DKeysProperty transferFunc0_;
    TransFunc1DKeysProperty transferFunc1_;
    TransFunc1DKeysProperty transferFunc2_;
    TransFunc1DKeysProperty transferFunc3_;

    FloatProperty samplingRate_;
    BoolProperty adaptiveSampling_;
    BoolProperty useIterativeRefinement_;
    IntProperty nodeLevelReduction_;
    IntProperty screenResolutionReduction_;

    FloatProperty interactionSamplingRate_;
    IntProperty interactionNodeLevelReduction_;
    IntProperty interactionScreenResolutionReduction_;

    BoolProperty useBricks_;
    StringOptionProperty textureFilterMode_;
    IntProperty brickBufferSizeMB_;
    IntProperty brickUploadBudgetMB_;
    IntProperty brickUploadTimeLimit_;
    BoolProperty keepInteractionBricks_;
    BoolProperty uploadBricksDuringInteraction_;

    BoolProperty applyChannelShift_;
    FloatVec3Property channelShift0_;
    FloatVec3Property channelShift1_;
    FloatVec3Property channelShift2_;
    FloatVec3Property channelShift3_;

    ButtonProperty resetChannelShift_;

    OpenCLProperty openclKernel_;
    CameraProperty cameraProperty_;

    BoolProperty enableDebugOutput_;

    // further members
    uint64_t* nodeBuffer_;
    size_t nodeBufferSize_;
    uint8_t* brickFlagBuffer_;
    size_t brickBufferSize_;            //< size of brick buffer in uint16_t elements (2 bytes)
    size_t numBrickBufferSlots_;        //< number of bricks that can be stored in brick buffer
    size_t brickUploadBufferSize_;      //< size of brick upload buffer in uint16_t elements (2 bytes)
    size_t numBrickUploadBufferSlots_;  //< number of bricks that can be uploaded per frame

    cl::Buffer* nodeBufferCL_;
    cl::Buffer* brickBufferCL_;
    cl::Buffer* brickFlagBufferCL_;

    cl::Buffer* rayBufferCL_;
    tgt::svec2 rayBufferDim_;

    tgt::Texture* transfuncTexture_;
    cl::ImageObject2D* transfuncTextureCL_;

    RenderState renderState_;
    tgt::Stopwatch refinementWatch_;

    /// Maps from node buffer index to a struct holding the actual octree node and additional information
    struct NodeInfo {
        NodeInfo() : nodeExists_(false), node_(0), parentNodeIndex_(0), level_(1000) {}
        bool nodeExists_;               ///< the IncrementalSingleOctreeRaycaster can have empty node entries
        const VolumeOctreeNode* node_;  ///< pointer to the actual octree node.
        size_t parentNodeIndex_;        ///< index of parent node.
        size_t level_;                  ///< the node's tree level (0 equals the level with the maximum resolution)
    };
    NodeInfo* nodeInfos_;

    // Multi-level LRU list of used bricks:
    // * for each octree level a separate list is maintained
    // * the least recently used brick is stored at the tail
    struct LRUEntry {
        LRUEntry(size_t nodeID) : nodeID_(nodeID), lastUsed_(0) {}
        size_t nodeID_;
        size_t lastUsed_; //< number of frames since the entry has been used (0 corresponds to current frame)
    };
    std::vector< std::list<LRUEntry> > leastRecentlyUsedBricks_;

    bool nodeBufferInvalid_;
    bool brickBufferInvalid_;
    bool transFuncTexInvalid_;

    tgt::Shader* rescaleShader_; ///< Shader used by the rescaleRendering() method

    tgt::Stopwatch frameClock_;

    static const std::string loggerCat_;

    // helper functions for querying/manipulating the node and flag buffers
    inline uint64_t createNodeBufferEntry(bool homogeneous, bool inPool, size_t childPointer, size_t brickAddress, const uint16_t* avgValues, size_t numChannels);
    inline void extractNodeValues(const uint64_t nodeEntry, bool& inPool, bool& homogeneous, size_t& childPointer, uint16_t& avgValue, size_t& brickPointer) const;
    inline void updateNodeChildPointer(uint64_t& nodeEntry, size_t childPointer);
    inline size_t getNodeChildPointer(uint64_t& nodeEntry) const;
    inline void updateNodeBrickPointer(uint64_t& nodeEntry, bool hasBrick, size_t brickPointer, const uint16_t* avgValues, size_t numChannels);
    inline size_t getNodeBrickPointer(uint64_t& nodeEntry) const;
    inline bool hasNodeBrick(uint64_t& nodeEntry) const;
    inline bool hasBrickBeenUsed(const uint8_t flagEntry) const;
    inline uint8_t setBrickUsed(const uint8_t flagEntry, const bool used);
    inline bool hasNodeBeenUsedInMIP(const uint8_t flagEntry, size_t numChannels) const;
    inline uint8_t setNodeBeenUsedInMIP(const uint8_t flagEntry, const bool used, size_t numChannels);
    inline bool isBrickRequested(const uint8_t flagEntry) const;
    inline uint8_t setBrickRequested(const uint8_t flagEntry, const bool requested);
};

} //namespace

#endif
