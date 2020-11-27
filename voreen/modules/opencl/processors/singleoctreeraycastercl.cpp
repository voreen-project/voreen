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

#include "singleoctreeraycastercl.h"

#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/octree/octreeutils.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/geometry/meshgeometry.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresizepoweroftwo.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/utils/stringutils.h"

#include "voreen/core/voreenapplication.h"

#include "modules/opencl/openclmodule.h"

#include "tgt/vector.h"
#include "tgt/tgt_math.h"
#include "tgt/stopwatch.h"
#include "tgt/shadermanager.h"
#include "tgt/textureunit.h"

#include <queue>
#include <stack>

using tgt::ivec2;
using tgt::ivec3;
using tgt::ivec4;
using tgt::vec2;
using tgt::vec3;
using tgt::vec4;
using tgt::svec3;

namespace {

const tgt::svec2 RENDER_WORKGROUP_SIZE(8, 8);

// 64-bit masks for node buffer entries
const uint64_t MASK_HOMOGENEOUS = 0x8000000000000000;  //< 10000000 00000000 00000000 00000000 00000000 00000000 00000000 00000000
const size_t   MASK_HOMOGENEOUS_SHIFT  = 63;
const size_t   MASK_HOMOGENEOUS_NUMBITS = 1;

const uint64_t MASK_INBRICKPOOL = 0x4000000000000000;  //< 01000000 00000000 00000000 00000000 00000000 00000000 00000000 00000000
const size_t   MASK_INBRICKPOOL_SHIFT  = 62;
const size_t   MASK_INBRICKPOOL_NUMBITS = 1;

// This is space for an index of the first of the 8 (consecutively stored)
// children of the current node in the node buffer. We can index children with
// a maximum offset of (1 >> MASK_CHILD_NUMBITS). For a tree depth of n we
// require 8^(n-1) (as in "power of") nodes for the leaf layer. Another 8^(n-1)
// nodes are sufficient for all non-leaf nodes. Thus, for 2*8^(n-1) =
// 2*(2^(3*(n-1)) = 2^(3n-2) we need 3n+1 bits.
const uint64_t MASK_CHILD =       0x3FFFFFE000000000;  //< 00111111 11111111 11111111 11100000 00000000 00000000 00000000 00000000
const size_t   MASK_CHILD_SHIFT  =  37;
const size_t   MASK_CHILD_NUMBITS = 25;
const size_t   MAX_ADDRESSABLE_NUM_NODES = 1ul << MASK_CHILD_NUMBITS;

// This is space for an index of the brick of the current node in the brick
// buffer. As all nodes can (potentially) have a brick, we support brick
// buffers with 1 >> (MASK_BRICK_NUMBITS) slots in total.
const uint64_t MASK_BRICK =       0x0000000001FFFFFF;  //< 00000000 00000000 00000000 00000000 00000001 11111111 11111111 11111111
const size_t   MASK_BRICK_SHIFT  =  0;
const size_t   MASK_BRICK_NUMBITS = 25;
const size_t   MAX_ADDRESSABLE_NUM_BRICKS = 1ul << MASK_BRICK_NUMBITS;

typedef struct {
    uint64_t MASK;
    size_t   SHIFT;
    size_t   NUMBITS;
} AVG_MASK;

// returns the node entry mask for avg values: the remaining 40 bits (lower five bytes, without brick pointer)
// are divided among the used channels
AVG_MASK getAvgMask(size_t channel, size_t numChannels) {
    tgtAssert(numChannels > 0 && numChannels <= 4, "invalid channel count");
    tgtAssert(channel >= 0 && channel < 4, "invalid channel");

    AVG_MASK mask;

    if (numChannels == 1) {
        mask.NUMBITS = 16;
    }
    else if (numChannels == 2) {
        mask.NUMBITS = 16;
    }
    else if (numChannels == 3) {
        mask.NUMBITS = 12;
    }
    else if (numChannels == 4) {
        mask.NUMBITS = 9;
    }
    else {
        tgtAssert(false, "should not get here (invalid channel count)");
        mask.NUMBITS = 0;
    }
    mask.SHIFT = (numChannels - channel - 1) * mask.NUMBITS;
    mask.MASK = ((1ul << mask.NUMBITS) - 1) << mask.SHIFT;

    tgtAssert(mask.NUMBITS * numChannels < 37, "Too many bits allocated for avg mask");

    return mask;
}

// 8-bit masks for node flag buffer entries
const uint8_t MASK_BRICK_INUSE =      1;  //< 00000001 (brick currently in use and should be kept on the GPU)
const uint8_t MASK_BRICK_REQUESTED =  2;  //< 00000010 (brick should be uploaded to the GPU)
const uint8_t MASK_MIP_0_USED      =  4;  //< 00000100
const uint8_t MASK_MIP_1_USED      =  8;  //< 00001000
const uint8_t MASK_MIP_2_USED      = 16;  //< 00010000
const uint8_t MASK_MIP_3_USED      = 32;  //< 00100000

// opencl literal conversion
std::string toOpenCLLiteral(const tgt::vec3& vec) {
    return "(float3)(" + voreen::ftos(vec.x) + "," + voreen::ftos(vec.y) + "," + voreen::ftos(vec.z) + ")";
}
std::string toOpenCLLiteral(const tgt::svec3& vec) {
    return "(uint3)(" + voreen::itos(vec.x) + "," + voreen::itos(vec.y) + "," + voreen::itos(vec.z) + ")";
}
std::string toOpenCLLiteral(const tgt::mat4& mat) {
    // note: tgt matrices are stored in row order
    std::string result = "(float16)(";
    for (size_t i=0; i<15; i++)
        result += voreen::ftos(mat.elem[i]) + ",";
    result += voreen::ftos(mat.elem[15]) + ")";

    return result;
}

//------------------------------------

struct RealWorldMappingCL {
    RealWorldMappingCL(voreen::RealWorldMapping mapping)
        : scale_(mapping.getScale())
        , offset_(mapping.getOffset())
    {}
    float scale_;
    float offset_;
};

} // namespace anonymous

namespace voreen {

inline uint64_t SingleOctreeRaycasterCL::createNodeBufferEntry(bool homogeneous, bool inPool, size_t childPointer, size_t brickAddress,
                                                               const uint16_t* avgValues, size_t numChannels)
{
    uint64_t result = 0;
    result |= (homogeneous ? 1LL : 0LL) << MASK_HOMOGENEOUS_SHIFT;
    result |= (inPool ? 1LL : 0LL)      << MASK_INBRICKPOOL_SHIFT;
    result |= ((uint64_t)childPointer   << MASK_CHILD_SHIFT) & MASK_CHILD;

    if (inPool) { //< brick present => store brick pointer
        result |= ((uint64_t)brickAddress   << MASK_BRICK_SHIFT) & MASK_BRICK_SHIFT;
    }
    else { //< no brick => store avg values
        tgtAssert(avgValues, "null pointer passed");
        for (size_t c=0; c<numChannels; c++) {
            result &= ~getAvgMask(c, numChannels).MASK;
            uint16_t avgScaled = avgValues[c] >> (16-getAvgMask(c, numChannels).NUMBITS); //< less than 16 bit available for storing avg value
            result |= ((uint64_t)avgScaled  << getAvgMask(c, numChannels).SHIFT) & getAvgMask(c, numChannels).MASK;
        }
    }

    return result;
}

inline void SingleOctreeRaycasterCL::extractNodeValues(const uint64_t nodeEntry, bool& inPool, bool& homogeneous, size_t& childPointer, uint16_t& avgValue, size_t& brickPointer) const {
    homogeneous = ((nodeEntry & MASK_HOMOGENEOUS) >> MASK_HOMOGENEOUS_SHIFT) == 1;
    inPool =      ((nodeEntry & MASK_INBRICKPOOL) >> MASK_INBRICKPOOL_SHIFT) == 1;
    childPointer = static_cast<size_t>((nodeEntry & MASK_CHILD) >> MASK_CHILD_SHIFT);

    if (inPool) { //< brick present => extract brick address
        brickPointer = static_cast<size_t>((nodeEntry & MASK_BRICK) >> MASK_BRICK_SHIFT);
    }
    else { //< no brick => extract avg values (TODO)
        avgValue = static_cast<uint16_t>( (nodeEntry & getAvgMask(0, 1).MASK) >> getAvgMask(0, 1).SHIFT );
    }
}

inline void SingleOctreeRaycasterCL::updateNodeChildPointer(uint64_t& nodeEntry, size_t childPointer) {
    // clear old value
    nodeEntry &= ~MASK_CHILD;
    // set new one
    nodeEntry |= (((uint64_t)childPointer) << MASK_CHILD_SHIFT) & MASK_CHILD;
}

inline size_t SingleOctreeRaycasterCL::getNodeChildPointer(uint64_t& nodeEntry) const {
    return static_cast<size_t>((nodeEntry & MASK_CHILD) >> MASK_CHILD_SHIFT);
}

inline void SingleOctreeRaycasterCL::updateNodeBrickPointer(uint64_t& nodeEntry, bool hasBrick, size_t brickPointer, const uint16_t* avgValues, size_t numChannels) {
    // has brick
    nodeEntry &= ~MASK_INBRICKPOOL;
    nodeEntry |= (((uint64_t)(hasBrick ? 1 : 0)) << MASK_INBRICKPOOL_SHIFT) & MASK_INBRICKPOOL;

    if (hasBrick) { // brick pointer
        nodeEntry &= ~MASK_BRICK;
        nodeEntry |= (((uint64_t)brickPointer) << MASK_BRICK_SHIFT) & MASK_BRICK;
    }
    else { // avg values
        tgtAssert(avgValues, "null pointer passed");
        for (size_t c=0; c<numChannels; c++) {
            nodeEntry &= ~getAvgMask(c, numChannels).MASK;
            uint16_t avgScaled = avgValues[c] >> (16-getAvgMask(c, numChannels).NUMBITS); //< less than 16 bit available for storing avg value
            nodeEntry |= ((uint64_t)avgScaled << getAvgMask(c, numChannels).SHIFT) & getAvgMask(c, numChannels).MASK;
        }
    }
}

inline size_t SingleOctreeRaycasterCL::getNodeBrickPointer(uint64_t& nodeEntry) const {
    return static_cast<size_t>((nodeEntry & MASK_BRICK) >> MASK_BRICK_SHIFT);
}

inline bool SingleOctreeRaycasterCL::hasNodeBrick(uint64_t& nodeEntry) const {
    return (nodeEntry & MASK_INBRICKPOOL) > 0;
}

// flag buffer access functions

inline bool SingleOctreeRaycasterCL::hasBrickBeenUsed(const uint8_t flagEntry) const {
    return (flagEntry & MASK_BRICK_INUSE) > 0;
}
inline uint8_t SingleOctreeRaycasterCL::setBrickUsed(const uint8_t flagEntry, const bool used) {
    return (flagEntry & ~MASK_BRICK_INUSE) | (used ? MASK_BRICK_INUSE : 0);
}
inline bool SingleOctreeRaycasterCL::hasNodeBeenUsedInMIP(const uint8_t flagEntry, size_t numChannels) const {
    uint8_t mask = MASK_MIP_0_USED | (numChannels > 1 ? MASK_MIP_1_USED : 0) | (numChannels > 2 ? MASK_MIP_2_USED : 0) | (numChannels > 3 ? MASK_MIP_3_USED : 0);
    return (flagEntry & mask) > 0;
}
inline uint8_t SingleOctreeRaycasterCL::setNodeBeenUsedInMIP(const uint8_t flagEntry, const bool used, size_t numChannels) {
    uint8_t mask = MASK_MIP_0_USED | (numChannels > 1 ? MASK_MIP_1_USED : 0) | (numChannels > 2 ? MASK_MIP_2_USED : 0) | (numChannels > 3 ? MASK_MIP_3_USED : 0);
    return (flagEntry & ~mask) | (used ? mask : 0);
}
inline bool SingleOctreeRaycasterCL::isBrickRequested(const uint8_t flagEntry) const {
    return (flagEntry & MASK_BRICK_REQUESTED) > 0;
}
inline uint8_t SingleOctreeRaycasterCL::setBrickRequested(const uint8_t flagEntry, const bool requested) {
    return (flagEntry & ~MASK_BRICK_REQUESTED) | (requested ? MASK_BRICK_REQUESTED : 0);
}


const std::string SingleOctreeRaycasterCL::loggerCat_("voreen.SingleOctreeRaycasterCL");

const int INVALID_BUFFERS = 50;

SingleOctreeRaycasterCL::SingleOctreeRaycasterCL()
    : OpenCLProcessor()
    , volumeInport_(Port::INPORT, "volumeInput", "Octree Input")
    , entryPointsInport_(Port::INPORT, "entryPointsInport", "Entry-points")
    , exitPointsInport_(Port::INPORT, "exitPointsInport", "Exit-points")
    , renderOutport_(Port::OUTPORT, "renderOutport", "Image Output")
    , internalRenderPort_(Port::OUTPORT, "internalRenderPort", "Internal RenderPort")
    , internalRenderPortDepth_(Port::OUTPORT, "internalRenderPortDepth", "Internal RenderPort Depth"/*, false, INVALID_RESULT, RenderPort::RENDERSIZE_DEFAULT, GL_LUMINANCE32F_ARB*/)
    , compositingMode_("compositingMode", "Compositing", Processor::INVALID_PROGRAM)
    , renderingQuality_("renderingQuality", "Rendering Quality")
    , interactionQuality_("interactionQuality", "Interaction Quality")
    , transferFunc0_("transferFunction", "Transfer Function")
    , transferFunc1_("transferFunction1", "Transfer Function 2")
    , transferFunc2_("transferFunction2", "Transfer Function 3")
    , transferFunc3_("transferFunction3", "Transfer Function 4")
    , samplingRate_("samplingRate", "Sampling Rate", 1.f, 0.1f, 5.f)
    , useIterativeRefinement_("useIterativeRefinement", "Iterative Refinement", true)
    , adaptiveSampling_("adaptiveSampling", "Adaptive Sampling", true, Processor::INVALID_PROGRAM)
    , nodeLevelReduction_("nodeLevelReduction", "Level of Detail Reduction", 0, 0, 3)
    , screenResolutionReduction_("screenResolutionReduction", "Screen Reduction", 1, 1, 4)
    , interactionSamplingRate_("interactionSamplingRate", "Sampling Rate", 0.75f, 0.1f, 5.f, Processor::VALID)
    , interactionNodeLevelReduction_("interactionNodeLevelReduction", "Level of Detail Reduction", 1, 0, 3, Processor::VALID)
    , interactionScreenResolutionReduction_("interactionScreenResolutionReduction", "Screen Resolution Reduction", 2, 1, 4, Processor::VALID)
    , useBricks_("useBricks", "Use Bricks", true, Processor::INVALID_PROGRAM)
    , textureFilterMode_("textureFilterMode", "Texture Filtering", Processor::INVALID_PROGRAM)
    , brickBufferSizeMB_("brickBufferSizeMB", "Brick Buffer Size (MB)", 512, 1, 2048, INVALID_BUFFERS)
    , brickUploadBudgetMB_("brickUploadBudgetMD", "Brick Upload Budget (MB)", 32, 1, 1024)
    , brickUploadTimeLimit_("brickUploadTimeLimit", "Brick Upload Time Limit (ms)", 250, 0, 1000)
    , keepInteractionBricks_("keepInteractionBricks", "Keep Interaction Level Bricks", true, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , uploadBricksDuringInteraction_("uploadBricksDuringInteracion", "Upload Bricks During Interaction", true, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , applyChannelShift_("applyChannelShift", "Apply Channel Shift", false, Processor::INVALID_PROGRAM, Property::LOD_ADVANCED)
    , channelShift0_("channelShift0", "Channel Shift 1", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
    , channelShift1_("channelShift1", "Channel Shift 2", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
    , channelShift2_("channelShift2", "Channel Shift 3", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
    , channelShift3_("channelShift3", "Channel Shift 4", tgt::vec3(0.f), tgt::vec3(-50.f), tgt::vec3(50.f), Processor::INVALID_RESULT, NumericProperty<tgt::vec3>::STATIC, Property::LOD_ADVANCED)
    , resetChannelShift_("resetchannelshift", "Reset Channel Shift")
    , openclKernel_("renderProgram", "OpenCL Kernel", VoreenApplication::app()->getBasePath("modules/opencl/cl/singleoctreeraycaster.cl"), Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , cameraProperty_("cameraProperty", "Camera", tgt::Camera(), /*false*/true, true, 500.f, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , enableDebugOutput_("enableDebugOutput", "Debug Output", false, Processor::VALID, Property::LOD_DEBUG)
    , nodeBuffer_(0)
    , brickFlagBuffer_(0)
    , brickBufferSize_(0)
    , numBrickBufferSlots_(0)
    , brickUploadBufferSize_(0)
    , numBrickUploadBufferSlots_(0)
    , nodeBufferCL_(0)
    , brickBufferCL_(0)
    , brickFlagBufferCL_(0)
    , rayBufferCL_(0)
    , rayBufferDim_(tgt::svec2::zero)
    , transfuncTexture_(0)
    , transfuncTextureCL_(0)
    , renderState_(RENDER_STATE_FIRST_FULL_FRAME)
    , nodeInfos_(0)
    , nodeBufferInvalid_(true)
    , brickBufferInvalid_(true)
    , transFuncTexInvalid_(true)
{
    addPort(volumeInport_);
    addPort(entryPointsInport_);
    addPort(exitPointsInport_);

    addPort(renderOutport_);
    addPrivateRenderPort(internalRenderPort_);
    addPrivateRenderPort(internalRenderPortDepth_);

    compositingMode_.addOption("dvr", "Direct Volume Rendering (DVR)",  "COMPOSITING_MODE_DVR");
    compositingMode_.addOption("mip", "Max Intensity Projection (MIP)", "COMPOSITING_MODE_MIP");
    compositingMode_.addOption("mop", "Max Opacity Projection (MOP)",   "COMPOSITING_MODE_MOP");
    addProperty(compositingMode_);

    // The same as "Very High", but does always renders up to the maximum level of detail.
    renderingQuality_.addSetting("Full",    &samplingRate_,                 1.f);
    renderingQuality_.addSetting("Full",    &useIterativeRefinement_,       true);
    renderingQuality_.addSetting("Full",    &nodeLevelReduction_,           0); //< no reduction
    renderingQuality_.addSetting("Full",    &screenResolutionReduction_,    1); //< no reduction (factor of 1.f)

    renderingQuality_.addSetting("Very High",    &samplingRate_,                 1.f);
    renderingQuality_.addSetting("Very High",    &useIterativeRefinement_,       true);
    renderingQuality_.addSetting("Very High",    &nodeLevelReduction_,           0); //< no reduction
    renderingQuality_.addSetting("Very High",    &screenResolutionReduction_,    1); //< no reduction (factor of 1.f)

    renderingQuality_.addSetting("High",    &samplingRate_,                 1.f);
    renderingQuality_.addSetting("High",    &useIterativeRefinement_,       true);
    renderingQuality_.addSetting("High",    &nodeLevelReduction_,           1); //< reduction by one level
    renderingQuality_.addSetting("High",    &screenResolutionReduction_,    1); //< no reduction

    renderingQuality_.addSetting("Medium",  &samplingRate_,                 0.75f);
    renderingQuality_.addSetting("Medium",  &useIterativeRefinement_,       false);
    renderingQuality_.addSetting("Medium",  &nodeLevelReduction_,           1); //< reduction by one level
    renderingQuality_.addSetting("Medium",  &screenResolutionReduction_,    2); //< reduction by a factor of 2.f

    renderingQuality_.addSetting("Low",     &samplingRate_,                 0.5f);
    renderingQuality_.addSetting("Low",     &useIterativeRefinement_,       false);
    renderingQuality_.addSetting("Low",     &nodeLevelReduction_,           2); //< reduction by two levels
    renderingQuality_.addSetting("Low",     &screenResolutionReduction_,    3); //< reduction by a factor of 3.f

    renderingQuality_.selectMode("High");
    renderingQuality_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::updatePropertyConfiguration));
    addProperty(renderingQuality_);

    // quality settings are the same as for the rendering quality (and should be kept that way!)
    interactionQuality_.addSetting("High",    &interactionSamplingRate_,                 1.f);
    interactionQuality_.addSetting("High",    &interactionNodeLevelReduction_,           1);
    interactionQuality_.addSetting("High",    &interactionScreenResolutionReduction_,    1);

    interactionQuality_.addSetting("Medium",  &interactionSamplingRate_,                 0.75f);
    interactionQuality_.addSetting("Medium",  &interactionNodeLevelReduction_,           1);
    interactionQuality_.addSetting("Medium",  &interactionScreenResolutionReduction_,    2);

    interactionQuality_.addSetting("Low",     &interactionSamplingRate_,                 0.5f);
    interactionQuality_.addSetting("Low",     &interactionNodeLevelReduction_,           2);
    interactionQuality_.addSetting("Low",     &interactionScreenResolutionReduction_,    3);

    interactionQuality_.selectMode("Medium");
    interactionQuality_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::updatePropertyConfiguration));
    addProperty(interactionQuality_);

    addProperty(transferFunc0_);
    addProperty(transferFunc1_);
    transferFunc1_.setVisibleFlag(false);
    addProperty(transferFunc2_);
    transferFunc2_.setVisibleFlag(false);
    addProperty(transferFunc3_);
    transferFunc3_.setVisibleFlag(false);

    transferFunc0_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::invalidateTransFuncTex));
    transferFunc1_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::invalidateTransFuncTex));
    transferFunc2_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::invalidateTransFuncTex));
    transferFunc3_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::invalidateTransFuncTex));

    addProperty(samplingRate_);
    addProperty(adaptiveSampling_);
    addProperty(useIterativeRefinement_);
    addProperty(nodeLevelReduction_);
    addProperty(screenResolutionReduction_);
    samplingRate_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));
    adaptiveSampling_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));
    useIterativeRefinement_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));
    nodeLevelReduction_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));
    screenResolutionReduction_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));
    //renderingQuality_.setGroupID("rendering-quality");
    samplingRate_.setGroupID("rendering-quality");
    adaptiveSampling_.setGroupID("rendering-quality");
    useIterativeRefinement_.setGroupID("rendering-quality");
    nodeLevelReduction_.setGroupID("rendering-quality");
    screenResolutionReduction_.setGroupID("rendering-quality");
    setPropertyGroupGuiName("rendering-quality", "Rendering Quality");

    addProperty(interactionSamplingRate_);
    addProperty(interactionNodeLevelReduction_);
    addProperty(interactionScreenResolutionReduction_);
    interactionSamplingRate_.setGroupID("interaction");
    interactionScreenResolutionReduction_.setGroupID("interaction");
    interactionNodeLevelReduction_.setGroupID("interaction");
    setPropertyGroupGuiName("interaction", "Interaction Quality");

    addProperty(useBricks_);
    useBricks_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::invalidateBrickBuffer));
    textureFilterMode_.addOption("nearest", "Nearest");
    textureFilterMode_.addOption("linear",  "Linear");
    textureFilterMode_.selectByKey("linear");
    addProperty(textureFilterMode_);
    addProperty(brickBufferSizeMB_);
    brickBufferSizeMB_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::invalidateBrickBuffer));
    addProperty(brickUploadBudgetMB_);
    brickUploadBudgetMB_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::invalidateBrickBuffer));
    addProperty(brickUploadTimeLimit_);
    addProperty(keepInteractionBricks_);
    addProperty(uploadBricksDuringInteraction_);
    useBricks_.setGroupID("bricking");
    textureFilterMode_.setGroupID("bricking");
    brickBufferSizeMB_.setGroupID("bricking");
    brickUploadBudgetMB_.setGroupID("bricking");
    brickBufferSizeMB_.setGroupID("bricking");
    brickUploadTimeLimit_.setGroupID("bricking");
    uploadBricksDuringInteraction_.setGroupID("bricking");
    keepInteractionBricks_.setGroupID("bricking");
    setPropertyGroupGuiName("bricking", "Brick Handling");

    addProperty(applyChannelShift_);
    addProperty(channelShift0_);
    addProperty(channelShift1_);
    addProperty(channelShift2_);
    addProperty(channelShift3_);
    addProperty(resetChannelShift_);
    resetChannelShift_.onClick(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetChannelShift));
    applyChannelShift_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::updatePropertyConfiguration));
    channelShift0_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));
    channelShift1_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));
    channelShift2_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));
    channelShift3_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));
    applyChannelShift_.setGroupID("channel-shift");
    channelShift0_.setGroupID("channel-shift");
    channelShift1_.setGroupID("channel-shift");
    channelShift2_.setGroupID("channel-shift");
    channelShift3_.setGroupID("channel-shift");
    resetChannelShift_.setGroupID("channel-shift");
    setPropertyGroupGuiName("channel-shift", "Channel Shift Correction (in voxels)");

    addProperty(openclKernel_);
    addProperty(cameraProperty_);
    cameraProperty_.onChange(MemberFunctionCallback<SingleOctreeRaycasterCL>(this, &SingleOctreeRaycasterCL::resetRenderState));

    addProperty(enableDebugOutput_);

    // check that node buffer bitmasks are complementary
    tgtAssert((MASK_HOMOGENEOUS & MASK_INBRICKPOOL) == 0, "bitmasks not complementary");
    tgtAssert((MASK_HOMOGENEOUS & MASK_CHILD)       == 0, "bitmasks not complementary");
    tgtAssert((MASK_HOMOGENEOUS & MASK_BRICK)       == 0, "bitmasks not complementary");

    tgtAssert((MASK_INBRICKPOOL & MASK_CHILD) == 0, "bitmasks not complementary");
    tgtAssert((MASK_INBRICKPOOL & MASK_BRICK) == 0, "bitmasks not complementary");
    tgtAssert((MASK_CHILD & MASK_BRICK) == 0, "bitmasks not complementary");

    for (size_t numChannels=0; numChannels<4; numChannels++) {
        for (size_t channel=0; channel < numChannels; channel++) {
            tgtAssert((MASK_HOMOGENEOUS & getAvgMask(channel, numChannels).MASK) == 0, "bitmasks not complementary");
            tgtAssert((MASK_INBRICKPOOL & getAvgMask(channel, numChannels).MASK) == 0, "bitmasks not complementary");
            tgtAssert((MASK_CHILD & getAvgMask(channel, numChannels).MASK) == 0, "bitmasks not complementary");
            //tgtAssert((MASK_BRICK & getAvgMask(channel, numChannels).MASK) == 0, "bitmasks not complementary");

            for (size_t i=0; i<channel; i++)
                tgtAssert((getAvgMask(i, numChannels).MASK & getAvgMask(channel, numChannels).MASK) == 0, "bitmasks not complementary");
        }
    }

    // check flag buffer masks
    tgtAssert((MASK_MIP_0_USED & MASK_BRICK_INUSE) == 0, "bitmasks not complementary");
    tgtAssert((MASK_MIP_0_USED & MASK_BRICK_REQUESTED) == 0, "bitmasks not complementary");
    tgtAssert((MASK_MIP_1_USED & MASK_BRICK_INUSE) == 0, "bitmasks not complementary");
    tgtAssert((MASK_MIP_1_USED & MASK_BRICK_REQUESTED) == 0, "bitmasks not complementary");
    tgtAssert((MASK_MIP_2_USED & MASK_BRICK_INUSE) == 0, "bitmasks not complementary");
    tgtAssert((MASK_MIP_2_USED & MASK_BRICK_REQUESTED) == 0, "bitmasks not complementary");
    tgtAssert((MASK_MIP_3_USED & MASK_BRICK_INUSE) == 0, "bitmasks not complementary");
    tgtAssert((MASK_MIP_3_USED & MASK_BRICK_REQUESTED) == 0, "bitmasks not complementary");


    tgtAssert((MASK_BRICK_REQUESTED & MASK_BRICK_INUSE) == 0, "bitmasks not complementary");
}

SingleOctreeRaycasterCL::~SingleOctreeRaycasterCL() {
    delete[] nodeBuffer_;
    nodeBuffer_ = 0;
    delete[] brickFlagBuffer_;
    brickFlagBuffer_ = 0;
}

Processor* SingleOctreeRaycasterCL::create() const {
    return new SingleOctreeRaycasterCL();
}

void SingleOctreeRaycasterCL::initialize() {
    OpenCLProcessor::initialize();

    QualityMode.addObserver(this);

    resetRenderState();

    updatePropertyConfiguration();
}

void SingleOctreeRaycasterCL::deinitialize() {
    QualityMode.removeObserver(this);
    OpenCLProcessor::deinitialize();
}

void SingleOctreeRaycasterCL::initializeCL() {
    //update brick buffer size property
    cl_ulong mem = OpenCLModule::getInstance()->getCLDevice().getInfo<cl_ulong>(CL_DEVICE_MAX_MEM_ALLOC_SIZE);
    mem /= 1048576; // 1024*1024 mem size in MB
    brickBufferSizeMB_.setMaxValue(mem);
    LINFO("Restrict brick buffer size to " << mem << "MB");

    // load rescale shader program
    rescaleShader_ = ShdrMgr.loadSeparate("passthrough.vert", "copyimage_octreecoarseness.frag", generateHeader(), false);

    // build OpenCL kernel
    rebuildKernels(0);
}

void SingleOctreeRaycasterCL::deinitializeCL() {
    delete nodeBufferCL_;
    nodeBufferCL_ = 0;
    delete brickBufferCL_;
    brickBufferCL_ = 0;
    delete brickFlagBufferCL_;
    brickFlagBufferCL_ = 0;

    delete[] nodeInfos_;
    nodeInfos_ = 0;

    ShdrMgr.dispose(rescaleShader_);
    rescaleShader_ = 0;

    delete transfuncTextureCL_;
    transfuncTextureCL_ = 0;
    delete transfuncTexture_;
    transfuncTexture_ = 0;

    delete rayBufferCL_;
    rayBufferCL_ = 0;
}

bool SingleOctreeRaycasterCL::isReady() const {
    return OpenCLProcessor::isReady();
}

void SingleOctreeRaycasterCL::beforeProcess() {
    OpenCLProcessor::beforeProcess();

    const VolumeBase* sourceVolume = volumeInport_.getData();
    if (!sourceVolume || !sourceVolume->hasRepresentation<VolumeOctreeBase>())
        return;
    const VolumeOctreeBase* octree = sourceVolume->getRepresentation<VolumeOctreeBase>();

    bool rebuildKernel = (getInvalidationLevel() >= Processor::INVALID_PROGRAM);

    if (volumeInport_.hasChanged()) {
        // update properties
        updatePropertyConfiguration();

        // assign octree source volume to transfer functions for extraction of derived data
        transferFunc0_.setVolume(sourceVolume, 0);
        if (sourceVolume->getNumChannels() > 1)
            transferFunc1_.setVolume(sourceVolume, 1);
        if (sourceVolume->getNumChannels() > 2)
            transferFunc2_.setVolume(sourceVolume, 2);
        if (sourceVolume->getNumChannels() > 3)
            transferFunc3_.setVolume(sourceVolume, 3);

        // adapt scene to octree/volume bounding box
        cameraProperty_.adaptInteractionToScene(sourceVolume->getBoundingBox().getBoundingBox(), tgt::min(sourceVolume->getSpacing()));

        nodeBufferInvalid_ = true;
        brickBufferInvalid_ = true;
        transFuncTexInvalid_ = true;
        rebuildKernel = true;
        resetRenderState();
    }

    if (entryPointsInport_.hasChanged() || exitPointsInport_.hasChanged())
        resetRenderState();

    // initialize octree buffers, if necessary
    if (nodeBufferInvalid_) {
        initializeNodeBuffer();
        initializeBrickBuffer();

        //octreeBuffersInvalid_ = false;
        nodeBufferInvalid_ = false;
        brickBufferInvalid_ = false;
        resetRenderState();
    }

    if (brickBufferInvalid_) {
        initializeBrickBuffer();

        // clear brick pointers
        if (nodeBufferSize_ > 0) {
            tgtAssert(nodeBuffer_, "no node buffer");
            tgtAssert(nodeInfos_, "no node infos");
            for (size_t i=0; i<nodeBufferSize_; i++) {
                if(nodeInfos_[i].nodeExists_) {
                    const VolumeOctreeNode* node = nodeInfos_[i].node_;
                    tgtAssert(node, "no node");
                    updateNodeBrickPointer(nodeBuffer_[i], false, 0, node->getAvgValues(), node->getNumChannels());
                }
            }

            tgtAssert(nodeBufferCL_, "no CL node buffer");
            tgtAssert(OpenCLModule::getInstance()->getCLCommandQueue(), "no command queue");
            OpenCLModule::getInstance()->getCLCommandQueue()->enqueueWriteBuffer(nodeBufferCL_, nodeBuffer_, true);
        }

        brickBufferInvalid_ = false;
        resetRenderState();
    }

    // initialize transfunc buffer, if necessary
    if (!transfuncTextureCL_ || transFuncTexInvalid_) {
        generateTransfuncTexture();
        resetRenderState();
    }

    // rebuild OpenCL kernel
    if (rebuildKernel) {
        rebuildKernels(octree);
        resetRenderState();
    }

}

void SingleOctreeRaycasterCL::process() {
    if (!volumeInport_.getData()->hasRepresentation<VolumeOctreeBase>()) {
        LWARNING("Input volume has no octree representation. Use OctreeCreator processor!");
        return;
    }
    const VolumeBase* inputVolume = volumeInport_.getData();
    tgtAssert(inputVolume, "no input volume");
    const VolumeOctreeBase* inputOctree = inputVolume->getRepresentation<VolumeOctreeBase>();
    tgtAssert(inputOctree, "no input octree");

    if (inputOctree->getNumChannels() > 4) {
        LERROR("Multi-channel octrees with more than 4 channels not supported");
        renderOutport_.clear();
        return;
    }
    if (inputOctree->getNumNodes() > MAX_ADDRESSABLE_NUM_NODES) {
        LERROR("Octree with more than " << MAX_ADDRESSABLE_NUM_NODES << " nodes not supported");
        renderOutport_.clear();
        return;
    }
    if (!openclKernel_.getProgram()) {
        LERROR("OpenCL program missing (full frame)");
        renderOutport_.clear();
        return;
    }
    if (!nodeBufferCL_ || !brickBufferCL_ || !brickFlagBufferCL_) {
        LERROR("OpenCL buffer(s) not initialized");
        renderOutport_.clear();
        return;
    }

    // resize outport to inport dimensions
    tgt::svec2 outputRenderSize = entryPointsInport_.getSize();
    if (outputRenderSize != tgt::svec2(renderOutport_.getSize()))
        renderOutport_.resize(outputRenderSize);
    LGL_ERROR;

    // resize internal render port according to coarseness and clear it
    int coarsenessFactor = QualityMode.isInteractionMode() ? interactionScreenResolutionReduction_.get() : screenResolutionReduction_.get();
    const tgt::svec2 internalRenderSize = renderOutport_.getSize() / coarsenessFactor;
    if (internalRenderSize != tgt::svec2(internalRenderPort_.getSize()))
        internalRenderPort_.resize(internalRenderSize);
    if (internalRenderSize != tgt::svec2(internalRenderPortDepth_.getSize()))
        internalRenderPortDepth_.resize(internalRenderSize);
    LGL_ERROR;

    glFinish();

    // interrupt refinement on viewport change
    if (renderState_ >= RENDER_STATE_FIRST_ITERATIVE_FRAME && rayBufferDim_ != (tgt::svec2)internalRenderPort_.getSize())
        resetRenderState();

    if (enableDebugOutput_.get())
        LINFO("--- next render pass ---");

    frameClock_.reset();
    frameClock_.start();

    // use average dimension resolution for calculating the sampling step size
    float samplingRate = QualityMode.isInteractionMode() ? interactionSamplingRate_.get() : samplingRate_.get();
    int avgDim = tgt::hadd((tgt::ivec3)inputOctree->getDimensions()) / 3;
    float samplingStepSize = 1.f / (float)(avgDim * samplingRate);

    // compute desired node level of detail
    size_t nodeLevel;
    if(renderingQuality_.getValue() == "Full" && !QualityMode.isInteractionMode()) {
        // Always show all levels in full quality mode
        nodeLevel = inputOctree->getNumLevels()-1;
    } else {
        // Otherwise compute level of detail from screen size and other parameters
        nodeLevel = computeNodeLevel(tgt::vec3(0.5f), inputOctree->getDimensions(),
                inputVolume->getTextureToWorldMatrix(), inputOctree->getNumLevels(),
                cameraProperty_.get().getPosition(), cameraProperty_.get().getProjectionMatrix(outputRenderSize), outputRenderSize);
        int nodeLevelReduction = QualityMode.isInteractionMode() ? interactionNodeLevelReduction_.get() : nodeLevelReduction_.get();
        nodeLevel -= std::min((size_t)nodeLevelReduction, nodeLevel);
    }

    // update brick buffer (upload bricks)
    size_t numBricksUsed = 0;
    size_t numBricksRequested = 0;
    size_t numBricksUploaded = 0;
    if ((renderState_ > RENDER_STATE_FIRST_INTERACTION_FRAME) && //< do not upload bricks before first interaction frame
        (renderState_ > RENDER_STATE_INTERACTION_FRAME || uploadBricksDuringInteraction_.get()) )
    {
        // determine level at which bricks should be kept in the buffer, even if they are currently unused (interaction)
        int keepLevel = -1;
        if (keepInteractionBricks_.get()) {
            keepLevel = (int)inputOctree->getNumLevels()-1 - (int)(nodeLevel);

            // make sure that not more than half of the brick buffer is currently occupied by interaction level bricks
            while (keepLevel < (int)inputOctree->getNumLevels()) {
                size_t numBricksKept = 0;
                for (size_t level=keepLevel; level < inputOctree->getNumLevels(); level++) {
                    tgtAssert(level < leastRecentlyUsedBricks_.size(), "invalid state of LRU list");
                    numBricksKept += leastRecentlyUsedBricks_.at(level).size();
                }
                if (numBricksKept <= (numBrickBufferSlots_ / 2))
                    break;
                else
                    keepLevel++;
            }
            tgtAssert(keepLevel >= 0 && keepLevel <= (int)inputOctree->getNumLevels(), "invalid keep level");
        }

        updateBrickBuffer(keepLevel, numBricksUsed, numBricksRequested, numBricksUploaded);
    }

    // update render state and derive display mode
    DisplayMode displayMode = (DisplayMode)-1;
    RenderState nextRenderState = renderState_;
    if (renderState_ == RENDER_STATE_FIRST_INTERACTION_FRAME) {
        displayMode = FULL_FRAME;
        nextRenderState = RENDER_STATE_INTERACTION_FRAME;
        if (enableDebugOutput_.get())
            LINFO("Type: First interaction frame");
    }
    if (renderState_ == RENDER_STATE_INTERACTION_FRAME) {
        displayMode = FULL_FRAME;
        if (enableDebugOutput_.get())
            LINFO("Type: Interaction frame");
    }
    if (renderState_ == RENDER_STATE_FIRST_FULL_FRAME) {
        displayMode = FULL_FRAME;
        nextRenderState = RENDER_STATE_FULL_FRAME;
        if (enableDebugOutput_.get())
            LINFO("Type: First full-frame");
    }
    if (renderState_ == RENDER_STATE_FULL_FRAME) {
        if (numBricksUploaded == 0) { // no more bricks could be uploaded => switch to iterative refinement
            if (useIterativeRefinement_.get())
                renderState_ = RENDER_STATE_FIRST_ITERATIVE_FRAME;
            else
                renderState_ = RENDER_STATE_FINISHED;
        }
        else { // new bricks uploaded => continue full frame rendering
            displayMode = FULL_FRAME;
            if (enableDebugOutput_.get())
                LINFO("Type: Full-frame");
        }
    }
    if (renderState_ == RENDER_STATE_FIRST_ITERATIVE_FRAME) {
        displayMode = REFINEMENT;
        nextRenderState = RENDER_STATE_ITERATIVE_FRAME;
        refinementWatch_.reset();
        refinementWatch_.start();
        if (enableDebugOutput_.get())
            LINFO("Type: First refinement");
    }
    if (renderState_ == RENDER_STATE_ITERATIVE_FRAME) {
        if (numBricksUploaded == 0) { // no further bricks requested => frame is finished
            renderState_ = RENDER_STATE_FINISHED;
            if (enableDebugOutput_.get())
                LINFO("Refinement finished: " + formatTime(refinementWatch_.getRuntime()));
        }
        else { // new bricks uploaded => continue refinement
            displayMode = REFINEMENT;
            if (enableDebugOutput_.get())
                LINFO("Type: Iterative Refinement");
        }
    }
    if (renderState_ == RENDER_STATE_FINISHED) {
        if (enableDebugOutput_.get())
            LINFO("Frame finished.");
        return;
    }
    tgtAssert(displayMode == FULL_FRAME || displayMode == REFINEMENT, "invalid display mode");

    // make sure that refinement kernel is available
    if (displayMode == REFINEMENT && !openclKernel_.getProgram(REFINEMENT)) {
        LWARNING("OpenCL program missing (refinement)");
        displayMode = FULL_FRAME;
    }

    // initialize ray buffer
    if (!rayBufferCL_ || rayBufferDim_ != (tgt::svec2)internalRenderPort_.getSize())
        initializeRayBuffer();

   /* if(displayMode != FULL_FRAME){
        size_t a=0,b=0,c=0;
        analyzeBrickFlagBuffer(a,b,c);
        std::cout << "number of traversed bricks: " << a << std::endl;
    }*/

    // clear flag buffer: keep traversed flags in refinement mode
    clearBrickFlagBuffer(nodeBufferSize_, displayMode == FULL_FRAME);

    // perform actual raycasting into internal render port
    renderFrame(displayMode, inputVolume, inputOctree, samplingStepSize, nodeLevel, useBricks_.get());

    // copy rendering from internal color and depth buffers to outport (thereby resize it)
    rescaleRendering(internalRenderPort_, internalRenderPortDepth_, entryPointsInport_, exitPointsInport_, renderOutport_);

    if (enableDebugOutput_.get())
        LINFO("Time: " << formatTime(frameClock_.getRuntime()));

    renderState_ = nextRenderState;

    //force repaint in animation
    if(QualityMode.isHighQualityMode()) { // && renderState_ < RENDER_STATE_FINISHED) if state is finished, function returns before if clause
        process();
    }
}

void SingleOctreeRaycasterCL::afterProcess() {
    OpenCLProcessor::afterProcess();

    if(QualityMode.isInteractionMode())
        processedInInteraction_ = true;

    // refinement not finished
    if (renderState_ < RENDER_STATE_FINISHED && openclKernel_.getProgram()) {
        invalidate();
        renderOutport_.invalidatePort();
    }
}

void SingleOctreeRaycasterCL::qualityModeChanged() {
    resetRenderState();
    if (!QualityMode.isInteractionMode() && processedInInteraction_) {
        processedInInteraction_ = false;
        invalidate();
    }

}

// private functions

void SingleOctreeRaycasterCL::initializeNodeBuffer() {
    tgtAssert(volumeInport_.hasData() && volumeInport_.getData()->hasRepresentation<VolumeOctreeBase>(), "no input octree");
    const VolumeOctreeBase* octree = volumeInport_.getData()->getRepresentation<VolumeOctreeBase>();

    const size_t nodeCount = octree->getNumNodes();
    const size_t numChannels = octree->getNumChannels();
    const size_t treeDepth = octree->getActualTreeDepth();
    tgtAssert(numChannels > 0 && numChannels <= 4, "invalid channel count");
    tgtAssert(nodeCount <= MAX_ADDRESSABLE_NUM_NODES, "Too many nodes");

    nodeBufferSize_ = nodeCount;

    delete nodeBufferCL_;
    nodeBufferCL_ = 0;

    delete[] nodeBuffer_;
    nodeBuffer_ = new uint64_t[nodeCount];

    delete[] nodeInfos_;
    nodeInfos_ = new NodeInfo[nodeCount];

    tgt::Stopwatch watch;
    watch.start();

    // pair consisting of a octree node whose children still have to be added to the buffer,
    // and the node's buffer offset
    typedef std::pair<const VolumeOctreeNode*, size_t> QueuedNode;
    std::queue<QueuedNode> workQueue;
    //std::stack<QueuedNode> workQueue;

    uint16_t avgValues[4]; //< is passed to createNodeBufferEntry
    bool useMaxValue = false; // compositingMode_.isSelected("mip") || compositingMode_.isSelected("mop");

    // start with root node, put encountered nodes into fifo queue, iterate until queue is empty
    const VolumeOctreeNode* rootNode = octree->getRootNode();
    tgtAssert(rootNode, "no root node");
    for (size_t c=0; c<numChannels; c++)
        avgValues[c] = useMaxValue ? rootNode->getMaxValue(c) : rootNode->getAvgValue(c);
    nodeBuffer_[0] = createNodeBufferEntry(rootNode->isHomogeneous(), false, 0, 0, avgValues, numChannels);
    nodeInfos_[0].node_ = octree->getRootNode();
    nodeInfos_[0].level_ = treeDepth-1;
    nodeInfos_[0].nodeExists_ = true;
    workQueue.push(QueuedNode(octree->getRootNode(), 0));
    size_t curOffset = 1; //< next after root node
    while (!workQueue.empty()) {
        // retrieve next node to process
        const VolumeOctreeNode* curNode = workQueue.front().first;
        size_t curNodeOffset = workQueue.front().second;
        tgtAssert(curNodeOffset < nodeCount, "invalid nodeOffset");
        tgtAssert(nodeInfos_[curNodeOffset].node_ == curNode && nodeInfos_[curNodeOffset].level_ < treeDepth, "invalid NodeInfo entry");

        const size_t nodeLevel = nodeInfos_[curNodeOffset].level_;

        /*const VolumeOctreeNode* curNode = workQueue.top().first;
        size_t nodeOffset = workQueue.top().second; */
        workQueue.pop();

        // no children => nothing to do
        if (!curNode->children_[0])
            continue;

        tgtAssert(curNodeOffset < curOffset,  "node offset not less than current offset");
        tgtAssert(curOffset+8 <= nodeCount, "invalid current buffer offset");

        // set curNode's child pointer to current offset
        updateNodeChildPointer(nodeBuffer_[curNodeOffset], curOffset);
        tgtAssert(((nodeBuffer_[curNodeOffset] & MASK_CHILD) >> MASK_CHILD_SHIFT) == curOffset, "child pointer not properly set");

        // create eight adjacent buffer entries at curOffset for children, and add them to work queue
        for (size_t childID = 0; childID < 8; childID++) {
            const VolumeOctreeNode* child = curNode->children_[childID];
            tgtAssert(child, "missing child");
            size_t childOffset = curOffset+childID;
            tgtAssert(childOffset < nodeCount, "invalid child node offset");

            for (size_t c=0; c<numChannels; c++)
                avgValues[c] = useMaxValue ? child->getMaxValue(c) : child->getAvgValue(c);
            uint64_t childEntry = createNodeBufferEntry(child->isHomogeneous(), false, 0, 0, avgValues, numChannels);
            tgtAssert(((childEntry & MASK_INBRICKPOOL) >> MASK_INBRICKPOOL_SHIFT) == 0, "in-brick-pool value not properly set");
            tgtAssert(((childEntry & MASK_HOMOGENEOUS) >> MASK_HOMOGENEOUS_SHIFT) == (!child->hasBrick()), "homogeneous value not properly set");
            tgtAssert(((childEntry & MASK_CHILD) >> MASK_CHILD_SHIFT) == 0, "child pointer not properly set");
            tgtAssert(((childEntry & getAvgMask(0,numChannels).MASK) >> getAvgMask(0,numChannels).SHIFT) == (avgValues[0] >> (16-getAvgMask(0,numChannels).NUMBITS)),
                "node avg value not properly set");
            //tgtAssert(((childEntry & MASK_BRICK) >> MASK_BRICK_SHIFT) == 0, "node brick pointer not properly set");
            nodeBuffer_[childOffset] = childEntry;
            tgtAssert(nodeLevel > 0, "invalid node level");
            nodeInfos_[childOffset].nodeExists_ = true;
            nodeInfos_[childOffset].node_ = child;
            nodeInfos_[childOffset].parentNodeIndex_ = curNodeOffset;
            nodeInfos_[childOffset].level_ = nodeLevel-1;

            workQueue.push(QueuedNode(child, childOffset));
        }
        curOffset += 8;
    }
    tgtAssert(curOffset == nodeCount, "buffer offset does not equal number of tree nodes");
    //tgtAssert(nodeBufferToNodeMap_.size() == nodeCount, "nodeBufferToNodeMap size does not equal number of tree nodes");

    LINFO("Node buffer size: \t" << formatMemorySize(nodeCount*sizeof(uint64_t)));
    LDEBUG("Node buffer generation time: " << watch.getRuntime() << " ms");

    // cl buffers
    cl::Context* context = OpenCLModule::getInstance()->getCLContext();
    if (!context) {
        LERROR("No CL context");
        return;
    }
    nodeBufferCL_ = new cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, nodeCount*8, nodeBuffer_);


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

void SingleOctreeRaycasterCL::initializeBrickBuffer() {
    tgtAssert(volumeInport_.hasData() && volumeInport_.getData()->hasRepresentation<VolumeOctreeBase>(), "no input octree");
    tgtAssert(nodeBuffer_, "no node buffer");

    const VolumeOctreeBase* octree = volumeInport_.getData()->getRepresentation<VolumeOctreeBase>();
    const size_t nodeCount = nodeBufferSize_;
    //tgtAssert(nodeBufferToNodeMap_.size() == nodeCount, "nodeBufferToNodeMap has invalid size");

    // clear OpenCL buffers
    delete brickBufferCL_;
    brickBufferCL_ = 0;
    delete brickFlagBufferCL_;
    brickFlagBufferCL_ = 0;

    // clear CPU buffers
    brickBufferSize_ = 0;
    numBrickBufferSlots_ = 0;

    brickUploadBufferSize_ = 0;
    numBrickUploadBufferSlots_ = 0;

    delete[] brickFlagBuffer_;
    brickFlagBuffer_ = 0;

    leastRecentlyUsedBricks_ = std::vector< std::list<LRUEntry> >(octree->getNumLevels());

    // allocate CPU buffers
    numBrickBufferSlots_ = tgt::ifloor((float)(brickBufferSizeMB_.get() << 20) / (float)octree->getBrickMemorySize());
    numBrickBufferSlots_ = std::min(numBrickBufferSlots_, nodeCount);
    if(numBrickBufferSlots_ >= MAX_ADDRESSABLE_NUM_BRICKS) {
        LWARNING("Restricting number of brick buffer slots to maximum addressable " << MAX_ADDRESSABLE_NUM_BRICKS);
        numBrickBufferSlots_ = MAX_ADDRESSABLE_NUM_BRICKS;
    }
    if (numBrickBufferSlots_ == 0) {
        LERROR("Brick buffer size (" << formatMemorySize(brickBufferSizeMB_.get() << 20) <<
            ") smaller than memory size of a single brick (" << formatMemorySize(octree->getBrickMemorySize()) << ")");
        return;
    }
    brickBufferSize_ = numBrickBufferSlots_ * (octree->getBrickMemorySize() / sizeof(uint16_t)); // brick buffer size in uint16_t

    numBrickUploadBufferSlots_ = std::min((brickUploadBudgetMB_.get() << 20) / octree->getBrickMemorySize(), numBrickBufferSlots_);
    if (numBrickUploadBufferSlots_ == 0) {
        LERROR("Brick upload buffer size (" << formatMemorySize(brickUploadBudgetMB_.get() << 20) <<
            ") smaller than memory size of a single brick (" << formatMemorySize(octree->getBrickMemorySize()) << ")");

        brickBufferSize_ = 0;
        numBrickBufferSlots_ = 0;

        return;
    }
    brickUploadBufferSize_ = std::min((brickUploadBudgetMB_.get() << 20) / sizeof(uint16_t), brickBufferSize_);

    brickFlagBuffer_ = new uint8_t[nodeCount];
    memset(brickFlagBuffer_, 0, nodeCount);

    // create OpenCL buffers
    cl::Context* context = OpenCLModule::getInstance()->getCLContext();
    if (!context) {
        LERROR("No CL context");
        return;
    }
    brickBufferCL_ = new cl::Buffer(context, CL_MEM_READ_WRITE, brickBufferSize_ * sizeof(uint16_t));
    brickFlagBufferCL_ = new cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, nodeCount, brickFlagBuffer_);

    LINFO("Brick buffer size: \t" << formatMemorySize(brickBufferSize_*sizeof(uint16_t)) << " / " << numBrickBufferSlots_ << " slots");
    LINFO("Brick upload buffer: \t" << formatMemorySize(brickUploadBufferSize_*sizeof(uint16_t)) << " / " << numBrickUploadBufferSlots_ << " slots");
    LINFO("Brick flag buffer size: " << formatMemorySize(nodeCount));
}

void SingleOctreeRaycasterCL::initializeRayBuffer() {
    LDEBUG("Initializing ray buffer");

    // free resources
    delete rayBufferCL_;
    rayBufferCL_ = 0;
    rayBufferDim_ = tgt::svec2::zero;

    // allocate refinement CL buffer
    cl::Context* context = OpenCLModule::getInstance()->getCLContext();
    if (!context) {
        LERROR("No CL context");
        return;
    }
    rayBufferDim_ = internalRenderPort_.getSize();



    // Rayinfo is defined in shader (octree_config.cl)
    size_t RAYINFO_NUM_ELEMENTS = 23;

    size_t bufferByteSize = tgt::hmul(rayBufferDim_) * RAYINFO_NUM_ELEMENTS * sizeof(cl_float);
    tgtAssert(bufferByteSize > 0, "invalid ray buffer size");
    rayBufferCL_ = new cl::Buffer(context, CL_MEM_READ_WRITE, bufferByteSize);
}

void SingleOctreeRaycasterCL::updateBrickBuffer(int keepLevel, size_t& numUsedInLastPass, size_t& numRequested, size_t& numBricksUploaded) {
    tgtAssert(volumeInport_.hasData() && volumeInport_.getData()->hasRepresentation<VolumeOctreeBase>(), "no input octree");

    numUsedInLastPass = 0;
    numRequested = 0;
    numBricksUploaded = 0;

    if (!brickBufferCL_)
        return;

    tgtAssert(nodeBuffer_, "no node buffer");
    tgtAssert(brickBufferCL_, "no brick buffer");
    tgtAssert(brickBufferSize_ > 0, "brick buffer size is 0");
    tgtAssert(numBrickBufferSlots_ > 0, "num brick buffer slots is 0");
    tgtAssert(brickUploadBufferSize_ > 0, "brick upload buffer size is 0");
    tgtAssert(numBrickUploadBufferSlots_ > 0, "num brick upload buffer slots is 0");
    tgtAssert(numBrickUploadBufferSlots_ <= numBrickBufferSlots_, "number of brick update slots exceeds number of brick slots");
    tgtAssert(brickFlagBuffer_, "no flag buffer");
    tgtAssert(nodeInfos_, "no NodeInfo buffer");

    const VolumeOctreeBase* octree = volumeInport_.getData()->getRepresentation<VolumeOctreeBase>();
    const size_t treeDepth = octree->getActualTreeDepth();
    const size_t numNodes = nodeBufferSize_;
    const size_t numChannels = octree->getNumChannels();
    const size_t brickMemorySize = octree->getBrickMemorySize();
    //tgtAssert(leastRecentlyUsedBricks_.size() == treeDepth, "invalid LRU list");
    //tgtAssert(nodeBufferToNodeMap_.size() == numNodes, "nodeBufferToNodeMap has invalid size");

    tgt::Stopwatch watch;
    watch.start();

    //
    // 1. Update LRU list:
    // - move used bricks that have been used during last rendering pass to front
    // - determine number of bricks used in last pass
    tgt::Stopwatch lruWatch;
    numUsedInLastPass = 0;
    for (size_t level=0; level < treeDepth; level++) {
        std::list<LRUEntry>& levelLRU = leastRecentlyUsedBricks_.at(level);
        for (std::list<LRUEntry>::iterator it = levelLRU.begin(); it!=levelLRU.end(); ) {
            LRUEntry& entry = *it;
            tgtAssert(entry.nodeID_ < numNodes, "invalid brick id");
            tgtAssert(nodeInfos_[entry.nodeID_].level_ == level, "node in wrong LRU list (wrong level)");
            if (hasBrickBeenUsed(brickFlagBuffer_[entry.nodeID_])) {
                numUsedInLastPass++;
                LRUEntry frontEntry(entry);
                it = levelLRU.erase(it);
                frontEntry.lastUsed_ = 0;
                levelLRU.push_front(frontEntry);
            }
            else
            {
                entry.lastUsed_++;
                ++it;
            }
        }
    }

    //
    // 2. Determine requested bricks, separated by level
    //
    std::vector< std::list<size_t> > requestedBricks(octree->getNumLevels());
    numRequested = 0;
    for (size_t nodeID = 0; nodeID < numNodes; nodeID++) {
        if (isBrickRequested(brickFlagBuffer_[nodeID])) {
            numRequested++;
            requestedBricks.at(nodeInfos_[nodeID].level_).push_back(nodeID);
        }
    }

    LINFO("Num requested: " << numRequested);

    // no bricks requested => keep brick buffer as it is
    if (numRequested == 0) {
        return;
    }

    // limit number of requested bricks to upload budget
    size_t numBricksToUpload = std::min(numRequested, numBrickUploadBufferSlots_);

    //
    // 3. Determine brick slots that are in use (have a GPU node pointing to them)
    //
    bool* brickSlotUsed = new bool[numBrickBufferSlots_];
    memset(brickSlotUsed, 0, numBrickBufferSlots_);  //< initialize with false
    size_t numSlotsOccupied = 0;
    for (size_t nodeID=0; nodeID<numNodes; nodeID++) {
        if (hasNodeBrick(nodeBuffer_[nodeID])) {
            size_t brickID = getNodeBrickPointer(nodeBuffer_[nodeID]);
            tgtAssert(brickID < numBrickBufferSlots_, "invalid brick id");
            tgtAssert(!brickSlotUsed[brickID], "brick slot already marked as used");
            brickSlotUsed[brickID] = true;
            numSlotsOccupied++;
        }
    }
    size_t numLRU = 0;
    for (size_t i=0; i<leastRecentlyUsedBricks_.size(); i++)
        numLRU += leastRecentlyUsedBricks_.at(i).size();
    tgtAssert(numSlotsOccupied == numLRU, "number of occupied slots does not match size of LRU list");
    tgtAssert(numSlotsOccupied >= numUsedInLastPass, "invalid number of occupied slots");

    //
    // 4. If not enough free capacity, mark brick slots to be freed:
    // a) that have not been used during last rendering pass
    // b) bricks of used lower-level nodes, if higher level nodes have been requested (lower resolution => higher priority)
    // NOTE: leastRecentlyUsedBricks_ is not modified at this step
    size_t numFreeBrickSlots = numBrickBufferSlots_ - numSlotsOccupied;
    // list of <nodeID, brickSlotID> pairs to be removed, separated by level
    std::vector< std::list<std::pair<size_t, size_t> > > brickKickoutList(octree->getNumLevels());
    size_t numToKickout = 0;
    if (numBricksToUpload > numFreeBrickSlots) {
        std::vector<std::list<LRUEntry> > tempLRUList(leastRecentlyUsedBricks_);

        const size_t additionalSlotsNeeded = numBricksToUpload - numFreeBrickSlots;
        const size_t numUnusedInLastPass = numSlotsOccupied - numUsedInLastPass;


        // a) remove *unused* bricks from tail of LRU list, beginning with the lowest level (highest resolution)
        for (size_t level=0; numToKickout < additionalSlotsNeeded; ) {
            while ((level < treeDepth) && (tempLRUList.at(level).empty() || tempLRUList.at(level).back().lastUsed_ == 0) )
                level++;
            if (level == treeDepth)
                break;

            // do not throw out unused bricks with a level above or equal the keepLevel (if set)
            if (keepLevel >= 0 && (int)level >= keepLevel)
                break;

            tgtAssert(level < treeDepth, "invalid LRU level");
            tgtAssert(!tempLRUList.at(level).empty(), "LRU list is empty");
            const LRUEntry entry = tempLRUList.at(level).back();
            tempLRUList.at(level).pop_back();

            size_t brickID = getNodeBrickPointer(nodeBuffer_[entry.nodeID_]);
            tgtAssert(brickID < numBrickBufferSlots_, "invalid brick id");
            tgtAssert(brickSlotUsed[brickID], "brick slot not marked as used");

            brickKickoutList.at(level).push_back(std::make_pair(entry.nodeID_, brickID));
            numToKickout++;
        }
        tgtAssert(numToKickout <= additionalSlotsNeeded, "removed too many bricks from buffer");

        // b) force-remove *used* bricks which have a lower level than requested bricks (if necessary)
        if ((numToKickout < additionalSlotsNeeded) /*&& forceReplaceUsedLowerLevelBricks_.get()*/) {
            tgtAssert(numFreeBrickSlots + numToKickout < numBricksToUpload, "invalid brick counts");

            // first remove numFreeBrickSlots from requestedBricks structure, beginning with highest level
            int curRequestedLevel = static_cast<int>(octree->getNumLevels()-1);
            for (size_t i=0; i<numFreeBrickSlots+numToKickout; i++) {
                while (requestedBricks.at(curRequestedLevel).empty()) {
                    tgtAssert(curRequestedLevel > 0, "requestedBrick vector in invalid state");
                    curRequestedLevel--;
                }
                requestedBricks.at(curRequestedLevel).pop_back();
            }

            // force-remove lower level bricks
            size_t forceRemoved = 0;
            size_t currentLRULevel = 0;
            while (numToKickout < additionalSlotsNeeded) {
                tgtAssert(numFreeBrickSlots + numToKickout < numBricksToUpload, "invalid brick counts");

                // find highest level at which bricks are used
                while (currentLRULevel < treeDepth && tempLRUList.at(currentLRULevel).empty()) {
                    //tgtAssert(currentLRULevel < treeDepth, "invalid LRU structure");
                    currentLRULevel++;
                }

                // find lowest level at which bricks are requested
                while (requestedBricks.at(curRequestedLevel).empty()) {
                    tgtAssert(curRequestedLevel > 0, "invalid requested bricks structure");
                    curRequestedLevel--;
                }

                // break, if no used bricks with a level higher than the remaining requested bricks are left
                if (curRequestedLevel <= (int)currentLRULevel)
                    break;

                // do not throw out bricks with a level above or equal the keepLevel (if set)
                if (keepLevel >= 0 && (int)currentLRULevel >= keepLevel)
                    break;

                // fetch least-recently used brick at current level
                tgtAssert(!tempLRUList.at(currentLRULevel).empty(), "LRU list is empty");
                const LRUEntry entry = tempLRUList.at(currentLRULevel).back();
                tempLRUList.at(currentLRULevel).pop_back();

                // only remove lower-level brick, if ancestor brick is present on the GPU
                /*size_t parentIndex = nodeInfos_[entry.nodeID_].parentNodeIndex_;
                tgtAssert(parentIndex < numNodes, "invalid parent index");
                size_t grandParentIndex = nodeInfos_[parentIndex].parentNodeIndex_;
                tgtAssert(grandParentIndex < numNodes, "invalid grandparent index");
                if (!hasNodeBrick(nodeBuffer_[parentIndex]) && !hasNodeBrick(nodeBuffer_[grandParentIndex]))
                    continue; */

                // add current entry to kick-out list
                size_t brickID = getNodeBrickPointer(nodeBuffer_[entry.nodeID_]);
                tgtAssert(brickID < numBrickBufferSlots_, "invalid brick id");
                tgtAssert(brickSlotUsed[brickID], "brick slot not marked as used");

                brickKickoutList.at(currentLRULevel).push_back(std::make_pair(entry.nodeID_, brickID));
                numToKickout++;

                requestedBricks.at(curRequestedLevel).pop_back();
            }
        } // force-removal of lower-level bricks
    }

    // actually remove marked bricks
    bool useMaxValue = false; // compositingMode_.isSelected("mip") || compositingMode_.isSelected("mop");
    uint16_t avgValues[4];
    for (size_t level=0; level<treeDepth; level++) {
        std::list<std::pair<size_t, size_t> >& currentLevelList = brickKickoutList.at(level);

        std::list<LRUEntry>& currentLRUList = leastRecentlyUsedBricks_.at(level);
        std::list<LRUEntry>::reverse_iterator lruIter =
            currentLRUList.rbegin(); //< currentLevelList contains bricks to be removed in reverse order

        while (!currentLevelList.empty()) {
            size_t nodeID = currentLevelList.front().first;
            size_t brickID = currentLevelList.front().second;
            currentLevelList.pop_front();

            tgtAssert(lruIter != currentLRUList.rend(), "inconsistent LRU lists");
            while (lruIter != currentLRUList.rend() && lruIter->nodeID_ != nodeID) {
                lruIter++;
                tgtAssert(lruIter != currentLRUList.rend(), "inconsistent LRU lists");
            }

            tgtAssert(brickID < numBrickBufferSlots_, "invalid brick id");
            tgtAssert(brickSlotUsed[brickID], "brick slot not marked as used");
            const VolumeOctreeNode* node = nodeInfos_[nodeID].node_;
            for (size_t c=0; c<numChannels; c++)
                avgValues[c] = useMaxValue ? node->getMaxValue(c) : node->getAvgValue(c);
            updateNodeBrickPointer(nodeBuffer_[nodeID], false, 0, avgValues, numChannels);
            brickSlotUsed[brickID] = false;
            numFreeBrickSlots++;
            numSlotsOccupied--;

            lruIter = std::list<LRUEntry>::reverse_iterator(currentLRUList.erase(--lruIter.base()));
        }
    }
    tgtAssert(numFreeBrickSlots+numSlotsOccupied == numBrickBufferSlots_, "invalid brick slot occupation count");


    // still no free brick slots => cannot upload any requested bricks
    if (numFreeBrickSlots == 0) {
        numBricksToUpload = 0;
        delete[] brickSlotUsed;
        /*if (enableDebugOutput_.get() && watch.getRuntime() > 0)
            LINFO("Update brick buffer time: " << watch.getRuntime() << " ms "); */
        return;
    }


    //
    // 5. Determine nodes whose bricks are to be uploaded to the GPU
    //
    // The nodes are implicitly ordered by inverse node level, giving higher-level nodes priority (i.e. lower res nodes)!
    // (see initializeNodeBuffer())
    //
    numBricksToUpload = std::min(numBricksToUpload, numFreeBrickSlots);
    std::vector<size_t> brickUploadNodeIds;
    size_t numInserted = 0;
    for (size_t nodeID=0; nodeID<numNodes && numInserted<numBricksToUpload; nodeID++) { //< higher levels come first
        const uint8_t brickFlag = brickFlagBuffer_[nodeID];
        if (isBrickRequested(brickFlag)) {
            tgtAssert(nodeInfos_[nodeID].node_ && nodeInfos_[nodeID].level_ < treeDepth, "invalid NodeInfo entry");
            const VolumeOctreeNode* node = nodeInfos_[nodeID].node_;
            const size_t nodeLevel = nodeInfos_[nodeID].level_;
            if (node->hasBrick()) {
                brickUploadNodeIds.push_back(nodeID);
                numInserted++;
                tgtAssert(brickUploadNodeIds.size() == numInserted, "brickUploadNodeIds has invalid size");
            }
            else {
                LWARNING("Non existing brick has been requested by GPU: nodeID=" << nodeID);
                //updateNodeBrickPointer(nodeBuffer_[nodeID], false, 0);
            }
        }
    }
    tgtAssert(brickUploadNodeIds.size() == numBricksToUpload, "brickUploadNodeIds has invalid size");


    //
    // 6. Write requested bricks to an temporary update buffer
    //
    tgt::Stopwatch updateBufferCreationWatch;
    updateBufferCreationWatch.start();
    uint16_t* brickUpdateBuffer = new uint16_t[numBricksToUpload * (brickMemorySize/sizeof(uint16_t))]; //< stores the bricks to be uploaded
    uint32_t* brickUpdateAddressBuffer = new uint32_t[numBricksToUpload]; //< stores the bricks' main buffer addresses

    // copy each requested brick to update buffer and assign next free brick slot as brick address
    size_t updateBufferPos = 0; //< current brick position in the temporary update buffer
    size_t nextFreeSlot = 0;    //< next free slot in the main brick buffer
    const int updateTimeLimit = brickUploadTimeLimit_.get(); //< time limit for brick buffer updates
    numBricksUploaded = 0;

    bool end_upload = false;
    //size_t lastRuntime = 0;
#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for (long i=0; i<static_cast<long>(numBricksToUpload); ++i) {
        if(end_upload) {
            continue;
        }
        const size_t nodeID = brickUploadNodeIds[i];
        tgtAssert(nodeID < numNodes, "invalid node id");
        const NodeInfo nodeInfo = nodeInfos_[nodeID];
        const VolumeOctreeNode* node = nodeInfo.node_;
        tgtAssert(node, "no node");
        tgtAssert(nodeInfo.level_ < treeDepth, "invalid level");
        size_t thisNextFreeSlot;
        size_t thisUpdateBufferPos;

        // find next free slot in main buffer
#ifdef VRN_MODULE_OPENMP
#pragma omp critical
#endif
        {
            while (brickSlotUsed[nextFreeSlot]) {
                nextFreeSlot++;
            }
            thisNextFreeSlot = nextFreeSlot;
            brickSlotUsed[thisNextFreeSlot] = true;

            thisUpdateBufferPos = updateBufferPos;
            updateBufferPos++;

        }
        tgtAssert(thisUpdateBufferPos < numBricksToUpload, "invalid update buffer pos");
        tgtAssert(thisNextFreeSlot < numBrickBufferSlots_, "invalid next free slot");

        // retrieve brick
        const uint16_t* brick = 0;
        try {
            brick = octree->getNodeBrick(nodeInfo.node_);
        }
        catch (VoreenException& e) {
            LWARNING("Failed to retrieve node brick " << nodeInfo.node_->getBrickAddress() << ": " + std::string(e.what()));
            continue;
        }
        tgtAssert(brick, "no brick returned (exception expected)");

        // copy brick to update buffer and store correspond brick address
        memcpy(brickUpdateBuffer + thisUpdateBufferPos*(brickMemorySize/sizeof(uint16_t)), brick, brickMemorySize);

        octree->releaseNodeBrick(node);

#ifdef VRN_MODULE_OPENMP
#pragma omp critical
#endif
        {
            brickUpdateAddressBuffer[thisUpdateBufferPos] = static_cast<uint32_t>(thisNextFreeSlot);

            // update node brick pointer
            updateNodeBrickPointer(nodeBuffer_[nodeID], true, thisNextFreeSlot, 0, numChannels);
            tgtAssert(getNodeBrickPointer(nodeBuffer_[nodeID]) == thisNextFreeSlot, "node brick pointer not properly stored");

            // add brick to LRU list
            leastRecentlyUsedBricks_.at(nodeInfo.level_).push_front(LRUEntry(nodeID));
            numSlotsOccupied++;
            numBricksUploaded++;

            if (updateTimeLimit > 0 && watch.getRuntime() > updateTimeLimit) {
                //LWARNING("Break: " << formatTime(watch.getRuntime()) << " (former: " << formatTime(lastRuntime) << ")");
                end_upload = true;
            }
        }
        //lastRuntime = watch.getRuntime();
    }

    tgtAssert(numBricksUploaded <= numBricksToUpload, "more bricks uploaded than expected");
    tgtAssert(numSlotsOccupied <= numBrickBufferSlots_, "more bricks occupied than slots available");
    updateBufferCreationWatch.stop();

    //
    // 7. Use OpenCL kernel to copy bricks from update buffer to main brick buffer on the GPU
    //
    cl::Context* context = OpenCLModule::getInstance()->getCLContext();
    cl::CommandQueue* commandQueue = OpenCLModule::getInstance()->getCLCommandQueue();
    tgtAssert(context, "No OpenCL context");
    tgtAssert(commandQueue, "No OpenCL command queue");
    tgtAssert(openclKernel_.getProgram(), "no OpenCL program");
    cl::Kernel* updateBrickKernel = openclKernel_.getProgram()->getKernel("updateBrickBuffer");
    if (!updateBrickKernel) {
        LERROR("Found no kernel 'updateBrickBuffer' in CL program");
        delete[] brickUpdateBuffer;
        delete[] brickUpdateAddressBuffer;
        delete[] brickSlotUsed;
        return;
    }

    tgt::Stopwatch gpuWatch;
    gpuWatch.start();
    if (numBricksUploaded > 0) {
        tgtAssert(numBricksUploaded <= numBrickUploadBufferSlots_, "number of bricks to upload exceeds upload budget");

        // create CL buffers for temporary update buffers (reducing their size to the actually used size!)
        cl::Buffer brickUpdateBufferCL(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
            numBricksUploaded*brickMemorySize, brickUpdateBuffer);
        cl::Buffer brickUpdateAddressBufferCL(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
            numBricksUploaded*4, brickUpdateAddressBuffer);

        updateBrickKernel->setArg(0, brickUpdateBufferCL);
        updateBrickKernel->setArg(1, brickUpdateAddressBufferCL);
        updateBrickKernel->setArg(2, static_cast<uint32_t>(numBricksUploaded));
        updateBrickKernel->setArg(3, brickBufferCL_);
        commandQueue->enqueue(updateBrickKernel, numBricksUploaded*16);
        commandQueue->finish();
    }
    gpuWatch.stop();

    // update node buffer
    commandQueue->enqueueWriteBuffer(nodeBufferCL_, nodeBuffer_, true);

    // free temporary CPU buffers
    delete[] brickUpdateBuffer;
    delete[] brickUpdateAddressBuffer;
    delete[] brickSlotUsed;

    if (enableDebugOutput_.get()) {
        LINFO("Update brick buffer time: " << watch.getRuntime() << " ms "
            << "(buffer creation: " << updateBufferCreationWatch.getRuntime() << " ms, "
            << "GPU upload: " << gpuWatch.getRuntime() << " ms)");

        LINFO("GPU bricks occupied/used/uploaded: " <<
            numSlotsOccupied  << " (" << formatMemorySize(numSlotsOccupied*brickMemorySize)  << ") / " <<
            numUsedInLastPass << " (" << formatMemorySize(numUsedInLastPass*brickMemorySize) << ") / " <<
            numBricksUploaded << " (" << formatMemorySize(numBricksUploaded*brickMemorySize) << ") "   );

        LINFO("Brick buffer: ");
        for (size_t level=0; level<treeDepth; level++)
            LINFO("- level " << level << ": " << leastRecentlyUsedBricks_.at(level).size());
    }

}

void SingleOctreeRaycasterCL::analyzeBrickFlagBuffer(size_t& numTraversed, size_t& numUsed, size_t& numRequested) const {
    tgtAssert(volumeInport_.hasData() && volumeInport_.getData()->hasRepresentation<VolumeOctreeBase>(), "no input octree");
    tgtAssert(brickFlagBuffer_, "no flag buffer");

    const VolumeOctreeBase* octree = volumeInport_.getData()->getRepresentation<VolumeOctreeBase>();
    const size_t numNodes = nodeBufferSize_;

    numTraversed = 0;
    numUsed = 0;
    numRequested = 0;
    for (size_t nodeID = 0; nodeID < numNodes; nodeID++) {
        uint8_t nodeFlag = brickFlagBuffer_[nodeID];
        if (hasNodeBeenUsedInMIP(nodeFlag, octree->getNumChannels()))
            numTraversed++;
        if (hasBrickBeenUsed(nodeFlag))
            numUsed++;
        if (isBrickRequested(nodeFlag))
            numRequested++;
    }
}

void SingleOctreeRaycasterCL::clearBrickFlagBuffer(size_t numNodes, bool clearUsedInMIP) {
    tgtAssert(brickFlagBuffer_ && brickFlagBufferCL_, "no flag buffer");
    tgtAssert(OpenCLModule::getInstance()->getCLCommandQueue(), "no command queue");

    if (clearUsedInMIP) {
        // reset entire buffer
        memset(brickFlagBuffer_, 0, numNodes);
    }
    else {
        // keep traversed flags
        for (size_t i=0; i<numNodes; i++)
            brickFlagBuffer_[i] &= ~(MASK_BRICK_INUSE | MASK_BRICK_REQUESTED);
    }
    OpenCLModule::getInstance()->getCLCommandQueue()->enqueueWriteBuffer(brickFlagBufferCL_, brickFlagBuffer_, true);
}

void SingleOctreeRaycasterCL::compareNodeToBuffer(const VolumeOctreeNode* node, size_t nodeOffset) const {
    tgtAssert(volumeInport_.hasData() && volumeInport_.getData()->hasRepresentation<VolumeOctreeBase>(), "no input octree");
    tgtAssert(node, "null pointer passed");
    tgtAssert(nodeBuffer_, "no node buffer");

    const VolumeOctreeBase* octree = volumeInport_.getData()->getRepresentation<VolumeOctreeBase>();
    if (nodeOffset >= nodeBufferSize_)
        throw VoreenException("Invalid node offset: " + itos(nodeOffset));

    bool inPool, homogeneous;
    size_t childPointer = 0;
    size_t brickPointer = 0;

    uint16_t avgValue;
    extractNodeValues(nodeBuffer_[nodeOffset], inPool, homogeneous, childPointer, avgValue, brickPointer);

    // compare node values
    if (inPool != 0)
        throw VoreenException("node at " + itos(nodeOffset) + ": mismatching inPool value (0 expected)");
    if (homogeneous == (!node->hasBrick()))
        throw VoreenException("node at " + itos(nodeOffset) + ": mismatching homogeneous value");
    uint16_t avgExpected = node->getAvgValue() >> (16-getAvgMask(0,1).NUMBITS);
    if (avgValue != avgExpected)
        throw VoreenException("node at " + itos(nodeOffset) + ": mismatching avg value [" + itos(avgValue) + "!=" + itos(avgExpected) + "]");
    if (brickPointer != 0)
        throw VoreenException("node at " + itos(nodeOffset) + ": mismatching brick pointer value [" + itos(avgValue) + "!= 0]");

    // recursively compare children
    if (childPointer == 0 && node->children_[0] != 0)
        throw VoreenException("node at " + itos(nodeOffset) + ": missing children (child pointer == 0)");
    else if (childPointer != 0 && node->children_[0] == 0)
        throw VoreenException("node at " + itos(nodeOffset) + ": unexpected children (child pointer != 0)");
    else if (node->children_[0] != 0) {
        tgtAssert(childPointer > 0, "invalid child pointer");
        for (int childID = 0; childID<8; childID++) {
            const VolumeOctreeNode* childNode = node->children_[childID];
            tgtAssert(childNode, "child is 0");
            size_t childOffset = childPointer + childID;
            compareNodeToBuffer(childNode, childOffset);
        }
    }
}


void SingleOctreeRaycasterCL::renderFrame(DisplayMode displayMode, const VolumeBase* volume, const VolumeOctreeBase* octree,
    float samplingStepSize, size_t nodeLevel, bool useBricks)
{
    tgtAssert(nodeBufferCL_ && brickFlagBufferCL_ && brickFlagBufferCL_ && rayBufferCL_, "missing CL buffer");
    tgtAssert(rayBufferDim_ == tgt::svec2(internalRenderPort_.getSize()), "invalid ray buffer size");

    const uint32_t coarsenessFactor = tgt::iround((float)renderOutport_.getSize().x / (float)internalRenderPort_.getSize().x);
    const bool sharing = OpenCLModule::getInstance()->getGLSharing();

    // retrieve CL kernel
    cl::Program* clProgram = openclKernel_.getProgram(displayMode);
    tgtAssert(clProgram, "no CL program");
    cl::Kernel* renderKernel = clProgram->getKernel("render");
    if (!renderKernel) {
        LERROR("Found no kernel 'render' in CL program");
        return;
    }
    tgtAssert(renderKernel, "no render kernel");

    // retrieve command queue
    const cl::Device& clDevice = OpenCLModule::getInstance()->getCLDevice();
    cl::Context* context = OpenCLModule::getInstance()->getCLContext();
    cl::CommandQueue* commandQueue = OpenCLModule::getInstance()->getCLCommandQueue();
    tgtAssert(context, "No OpenCL context");
    tgtAssert(commandQueue, "No OpenCL command queue");

    std::unique_ptr<cl::MemoryObject> entryTex;
    std::unique_ptr<cl::MemoryObject> exitTex;

    std::unique_ptr<cl::MemoryObject> outTex;
    std::unique_ptr<cl::MemoryObject> outTexDepth;
    if(sharing) {
        // entry/exit points
        LDEBUG("Binding entry/exit point textures");
        entryTex.reset(new cl::SharedTexture(context, CL_MEM_READ_ONLY, entryPointsInport_.getColorTexture()));
        exitTex.reset(new cl::SharedTexture(context, CL_MEM_READ_ONLY, exitPointsInport_.getColorTexture()));

        // output texture
        LDEBUG("Binding output textures");
        outTex.reset(new cl::SharedTexture(context, CL_MEM_WRITE_ONLY, internalRenderPort_.getColorTexture()));
        outTexDepth.reset(new cl::SharedTexture(context, CL_MEM_WRITE_ONLY, internalRenderPortDepth_.getColorTexture()));
    } else {
        // entry/exit points
        LDEBUG("Binding entry/exit point textures");
        entryPointsInport_.getColorTexture()->downloadTexture();
        exitPointsInport_.getColorTexture()->downloadTexture();
        entryTex.reset(new cl::ImageObject2D(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, entryPointsInport_.getColorTexture()));
        exitTex.reset(new cl::ImageObject2D(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, exitPointsInport_.getColorTexture()));

        // output texture
        LDEBUG("Binding output textures");
        outTex.reset(new cl::ImageObject2D(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, internalRenderPort_.getColorTexture()));
        outTexDepth.reset(new cl::ImageObject2D(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, internalRenderPortDepth_.getColorTexture()));
    }

    // transfer function
    tgtAssert(transfuncTexture_ && transfuncTextureCL_, "no transfunc tex");
    tgtAssert(transfuncTexture_->getHeight() == static_cast<int>(octree->getNumChannels()*3), "invalid transfunc texture height");
    if (!transfuncTextureCL_)
        generateTransfuncTexture();
    tgt::vec2 transFuncDomains[4];
    transFuncDomains[0] = transferFunc0_.get()->getDomain();
    transFuncDomains[1] = transferFunc1_.get()->getDomain();
    transFuncDomains[2] = transferFunc2_.get()->getDomain();
    transFuncDomains[3] = transferFunc3_.get()->getDomain();

    const RealWorldMapping rwm = volume->getRealWorldMapping();

    // configure kernel
    cl_uint paramID = 0;
    // entry/exit textures already set above...
    renderKernel->setArg(paramID++, entryTex.get());
    renderKernel->setArg(paramID++, exitTex.get());
    renderKernel->setArg(paramID++, internalRenderPort_.getSize());
    renderKernel->setArg(paramID++, coarsenessFactor);

    renderKernel->setArg(paramID++, samplingStepSize);
    tgt::ivec4 volumeDim((int)volume->getDimensions().x, (int)volume->getDimensions().y, (int)volume->getDimensions().z, 1);
    renderKernel->setArg(paramID++, volumeDim);

    renderKernel->setArg(paramID++, nodeBufferCL_);
    renderKernel->setArg(paramID++, brickBufferCL_);
    renderKernel->setArg(paramID++, brickFlagBufferCL_);
    renderKernel->setArg(paramID++, static_cast<const uint32_t>(nodeLevel));

    renderKernel->setArgTemplate<RealWorldMappingCL>(paramID++, RealWorldMappingCL(rwm));
    renderKernel->setArg(paramID++, transfuncTextureCL_);
    renderKernel->setArg(paramID++, 4*sizeof(tgt::vec2), &transFuncDomains);

    renderKernel->setArg(paramID++, outTex.get());
    renderKernel->setArg(paramID++, outTexDepth.get());

    renderKernel->setArg(paramID++, rayBufferCL_);
    renderKernel->setArg(paramID++, static_cast<uint32_t>(renderState_ == RENDER_STATE_FIRST_ITERATIVE_FRAME));

    if (applyChannelShift_.get()) {
        // channel shift is set in voxel coordinates
        renderKernel->setArg(paramID++, tgt::vec4(channelShift0_.get(), 1.f));
        renderKernel->setArg(paramID++, tgt::vec4(channelShift1_.get(), 1.f));
        renderKernel->setArg(paramID++, tgt::vec4(channelShift2_.get(), 1.f));
        renderKernel->setArg(paramID++, tgt::vec4(channelShift3_.get(), 1.f));
    }

    // determine global work size: next multiple of work group size that is larger than viewport dim
    tgt::svec2 globalWorkSize;
    globalWorkSize.x = tgt::iceil(internalRenderPort_.getSize().x / (float)RENDER_WORKGROUP_SIZE.x) * RENDER_WORKGROUP_SIZE.x;
    globalWorkSize.y = tgt::iceil(internalRenderPort_.getSize().y / (float)RENDER_WORKGROUP_SIZE.y) * RENDER_WORKGROUP_SIZE.y;

    // execute kernel
    if(sharing) {
        commandQueue->enqueueAcquireGLObject(entryTex.get());
        commandQueue->enqueueAcquireGLObject(exitTex.get());
        commandQueue->enqueueAcquireGLObject(outTex.get());
        commandQueue->enqueueAcquireGLObject(outTexDepth.get());
    }
    commandQueue->enqueue(renderKernel, globalWorkSize, RENDER_WORKGROUP_SIZE);
    if(sharing) {
        commandQueue->enqueueReleaseGLObject(entryTex.get());
        commandQueue->enqueueReleaseGLObject(exitTex.get());
        commandQueue->enqueueReleaseGLObject(outTex.get());
        commandQueue->enqueueReleaseGLObject(outTexDepth.get());
    } else {
        cl::ImageObject2D* outTexPtr = dynamic_cast<cl::ImageObject2D*>(outTex.get());
        cl::ImageObject2D* outTexDepthPtr = dynamic_cast<cl::ImageObject2D*>(outTexDepth.get());
        tgtAssert(outTexPtr, "Expected image ImageObject2D");
        tgtAssert(outTexDepthPtr, "Expected image ImageObject2D");
        commandQueue->enqueueReadImage(*outTexPtr, internalRenderPort_.getColorTexture(), false);
        commandQueue->enqueueReadImage(*outTexDepthPtr, internalRenderPortDepth_.getColorTexture(), true);
        internalRenderPort_.getColorTexture()->uploadTexture();
        internalRenderPortDepth_.getColorTexture()->uploadTexture();
    }

    // read back flag buffer
    commandQueue->enqueueReadBuffer(brickFlagBufferCL_, brickFlagBuffer_, true);

    commandQueue->finish();
}

void SingleOctreeRaycasterCL::generateTransfuncTexture() {
    delete transfuncTextureCL_;
    transfuncTextureCL_ = 0;
    delete transfuncTexture_;
    transfuncTexture_ = 0;

    LDEBUG("Generating transfunc texture");

    tgtAssert(transferFunc0_.get()->getTexture(), "no transfunc texture");
    const size_t numChannels = std::min<size_t>(volumeInport_.hasData() ? volumeInport_.getData()->getNumChannels() : 4, 4);
    const size_t transfuncDim = transferFunc0_.get()->getTexture()->getWidth();
    const size_t bytesPerTexel = 4*sizeof(uint8_t);         //< RGBA/uint8
    tgt::svec3 textureDim(transfuncDim, 3*numChannels, 1);  //< each transfunc is stored 3 times (for linear filtering)

    // combined texture buffer for all transfer functions
    GLubyte* textureData = new GLubyte[tgt::hmul(textureDim)*bytesPerTexel];
    memset(textureData, 0, tgt::hmul(textureDim)*bytesPerTexel);

    // 0. channel
    size_t channelOffset = 0;
    GLubyte* channelData = transferFunc0_.get()->getTexture()->downloadTextureToBuffer(GL_RGBA, GL_UNSIGNED_BYTE);
    for (size_t i=0; i<3; i++)
        memcpy(textureData + channelOffset + i*transfuncDim*bytesPerTexel, channelData, transfuncDim*bytesPerTexel);
    delete[] channelData;

    // 1. channel
    if (numChannels > 1) {
        const size_t channelOffset = 3*transfuncDim*bytesPerTexel;
        const size_t transfuncDim1 = transferFunc1_.get()->getDimensions().x;
        if (transfuncDim1 != transfuncDim) {
            LWARNING("Texture width of transfer function 1 differs from transfer function 0 [" << transfuncDim1 << " != " << transfuncDim << "]."
                << " Unable to use transfer function for channel 1.");
        }
        else {
            GLubyte* channelData = transferFunc1_.get()->getTexture()->downloadTextureToBuffer(GL_RGBA, GL_UNSIGNED_BYTE);
            for (size_t i=0; i<3; i++)
                memcpy(textureData + channelOffset + i*transfuncDim*bytesPerTexel, channelData, transfuncDim*bytesPerTexel);
            delete[] channelData;
        }
    }

    // 2. channel
    if (numChannels > 2) {
        const size_t channelOffset = 6*transfuncDim*bytesPerTexel;
        const size_t transfuncDim2 = transferFunc2_.get()->getDimensions().x;
        if (transfuncDim2 != transfuncDim) {
            LWARNING("Texture width of transfer function 2 differs from transfer function 0 [" << transfuncDim2 << " != " << transfuncDim << "]."
                << " Unable to use transfer function for channel 2.");
        }
        else {
            GLubyte* channelData = transferFunc2_.get()->getTexture()->downloadTextureToBuffer(GL_RGBA, GL_UNSIGNED_BYTE);
            for (size_t i=0; i<3; i++)
                memcpy(textureData + channelOffset + i*transfuncDim*bytesPerTexel, channelData, transfuncDim*bytesPerTexel);
            delete[] channelData;
        }
    }

    // 3. channel
    if (numChannels > 3) {
        const size_t channelOffset = 9*transfuncDim*bytesPerTexel;
        const size_t transfuncDim3 = transferFunc3_.get()->getDimensions().x;
        if (transfuncDim3 != transfuncDim) {
            LWARNING("Texture width of transfer function 3 differs from transfer function 0 [" << transfuncDim3 << " != " << transfuncDim << "]."
                << " Unable to use transfer function for channel 3.");
        }
        else {
            GLubyte* channelData = transferFunc3_.get()->getTexture()->downloadTextureToBuffer(GL_RGBA, GL_UNSIGNED_BYTE);
            for (size_t i=0; i<3; i++)
                memcpy(textureData + channelOffset + i*transfuncDim*bytesPerTexel, channelData, transfuncDim*bytesPerTexel);
            delete[] channelData;
        }
    }

    // create final texture
    transfuncTexture_ = new tgt::Texture(textureDim, GL_RGBA,GL_RGBA,GL_UNSIGNED_BYTE,
                                         tgt::Texture::LINEAR, tgt::Texture::REPEAT, textureData,true);
    LGL_ERROR;

    // create CL wrapper
    cl::Context* context = OpenCLModule::getInstance()->getCLContext();
    tgtAssert(context, "No OpenCL context");
    cl::ImageFormat imgf(CL_RGBA, CL_UNORM_INT8);
    transfuncTextureCL_ = new cl::ImageObject2D(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, imgf,
        transfuncTexture_->getWidth(), transfuncTexture_->getHeight(), 0, transfuncTexture_->getCpuTextureData());

    transFuncTexInvalid_ = false;
}

void SingleOctreeRaycasterCL::rebuildKernels(const VolumeOctreeBase* octree) {
    const size_t numChannels = octree ? octree->getNumChannels() : 1;

    std::ostringstream clDefinesCommon;
    clDefinesCommon << " -cl-fast-relaxed-math -cl-mad-enable";
    // set include path for modules (in "...")
#ifdef WIN32
    clDefinesCommon << " -I\"" << tgt::FileSystem::cleanupPath(VoreenApplication::app()->getBasePath() + "/modules/opencl/cl/") << "\"";
#else
    // TODO: no file paths with spaces supported
    std::string includePathString = tgt::FileSystem::cleanupPath(VoreenApplication::app()->getBasePath() + "/modules/opencl/cl/");
    if (includePathString.find(" ") != std::string::npos)
       LWARNING("CL include path contains spaces, might lead to problems: " << includePathString);
    //clDefinesCommon << " -I" << tgt::FileSystem::cleanupPath(VoreenApplication::app()->getBasePath() + "/modules/opencl/cl/");
    clDefinesCommon << " -I" << includePathString;
#endif

    // octree defines
    clDefinesCommon << " -DOCTREE_DIMENSIONS_DEF="    << (octree ? octree->getOctreeDim().x : 1);
    clDefinesCommon << " -DOCTREE_BRICKDIM_DEF="      << (octree ? octree->getBrickDim().x  : 1);
    clDefinesCommon << " -DOCTREE_DEPTH_DEF="         << (octree ? octree->getNumLevels()   : 1);
    clDefinesCommon << " -DOCTREE_NUMCHANNELS_DEF="   << numChannels;

    // node masks
    clDefinesCommon << " -DMASK_HOMOGENEOUS_DEF="         << MASK_HOMOGENEOUS;
    clDefinesCommon << " -DMASK_HOMOGENEOUS_SHIFT_DEF="   << MASK_HOMOGENEOUS_SHIFT;
    clDefinesCommon << " -DMASK_HOMOGENEOUS_NUMBITS_DEF=" << MASK_HOMOGENEOUS_NUMBITS;

    clDefinesCommon << " -DMASK_INBRICKPOOL_DEF="         << MASK_INBRICKPOOL;
    clDefinesCommon << " -DMASK_INBRICKPOOL_SHIFT_DEF="   << MASK_INBRICKPOOL_SHIFT;
    clDefinesCommon << " -DMASK_INBRICKPOOL_NUMBITS_DEF=" << MASK_INBRICKPOOL_NUMBITS;

    clDefinesCommon << " -DMASK_CHILD_DEF="               << MASK_CHILD;
    clDefinesCommon << " -DMASK_CHILD_SHIFT_DEF="         << MASK_CHILD_SHIFT;
    clDefinesCommon << " -DMASK_CHILD_NUMBITS_DEF="       << MASK_CHILD_NUMBITS;

    clDefinesCommon << " -DMASK_BRICK_DEF="               << MASK_BRICK;
    clDefinesCommon << " -DMASK_BRICK_SHIFT_DEF="         << MASK_BRICK_SHIFT;
    clDefinesCommon << " -DMASK_BRICK_NUMBITS_DEF="       << MASK_BRICK_NUMBITS;

    clDefinesCommon << " -DMASK_AVG_0_DEF="               << getAvgMask(0, numChannels).MASK;
    clDefinesCommon << " -DMASK_AVG_0_SHIFT_DEF="         << getAvgMask(0, numChannels).SHIFT;
    clDefinesCommon << " -DMASK_AVG_0_NUMBITS_DEF="       << getAvgMask(0, numChannels).NUMBITS;

    clDefinesCommon << " -DMASK_AVG_1_DEF="               << getAvgMask(1, numChannels).MASK;
    clDefinesCommon << " -DMASK_AVG_1_SHIFT_DEF="         << getAvgMask(1, numChannels).SHIFT;
    clDefinesCommon << " -DMASK_AVG_1_NUMBITS_DEF="       << getAvgMask(1, numChannels).NUMBITS;

    clDefinesCommon << " -DMASK_AVG_2_DEF="               << getAvgMask(2, numChannels).MASK;
    clDefinesCommon << " -DMASK_AVG_2_SHIFT_DEF="         << getAvgMask(2, numChannels).SHIFT;
    clDefinesCommon << " -DMASK_AVG_2_NUMBITS_DEF="       << getAvgMask(2, numChannels).NUMBITS;

    clDefinesCommon << " -DMASK_AVG_3_DEF="               << getAvgMask(3, numChannels).MASK;
    clDefinesCommon << " -DMASK_AVG_3_SHIFT_DEF="         << getAvgMask(3, numChannels).SHIFT;
    clDefinesCommon << " -DMASK_AVG_3_NUMBITS_DEF="       << getAvgMask(3, numChannels).NUMBITS;

    // brick flag buffer masks
    clDefinesCommon << " -DMASK_BRICK_INUSE_DEF="       << (uint32_t)MASK_BRICK_INUSE;
    clDefinesCommon << " -DMASK_BRICK_REQUESTED_DEF="   << (uint32_t)MASK_BRICK_REQUESTED;
    clDefinesCommon << " -DMASK_MIP_0_USED_DEF="    << (uint32_t)MASK_MIP_0_USED;
    clDefinesCommon << " -DMASK_MIP_1_USED_DEF="    << (uint32_t)MASK_MIP_1_USED;
    clDefinesCommon << " -DMASK_MIP_2_USED_DEF="    << (uint32_t)MASK_MIP_2_USED;
    clDefinesCommon << " -DMASK_MIP_3_USED_DEF="    << (uint32_t)MASK_MIP_3_USED;

    // mode defines
    clDefinesCommon << " -D" << compositingMode_.getValue();
    if (adaptiveSampling_.get())
        clDefinesCommon << " -DADAPTIVE_SAMPLING";
    if (textureFilterMode_.isSelected("linear"))
        clDefinesCommon << " -DTEXTURE_FILTER_LINEAR";
    if (useBricks_.get())
        clDefinesCommon << " -DUSE_BRICKS";
    if (applyChannelShift_.get())
        clDefinesCommon << " -DAPPLY_CHANNEL_SHIFT";

    // defines for full frame kernel (interaction);
    std::string clDefinesFullFrame = clDefinesCommon.str();
    clDefinesFullFrame += " -DUSE_ANCESTOR_NODES";
    openclKernel_.setDefines(clDefinesFullFrame, FULL_FRAME);

    // defines for refinement kernel
    std::string clDefinesRefinement = clDefinesCommon.str();
    clDefinesRefinement += " -DDISPLAY_MODE_REFINEMENT";
    openclKernel_.setDefines(clDefinesRefinement, REFINEMENT);

    //LINFO("kernel defines: " << clDefinesCommon.str());
    bool success = openclKernel_.rebuild();
    if (!success) {
        LERROR("Failed to build CL program");
        openclKernel_.clearProgram();
    }
    LGL_ERROR;
}

void SingleOctreeRaycasterCL::rescaleRendering(RenderPort& srcPort, RenderPort& srcPortDepth,
    RenderPort& entryPoints, RenderPort& exitPoints, RenderPort& destPort)
{
    // activate and clear output render target
    destPort.activateTarget();
    destPort.clearTarget();


    // activate shader and set uniforms
    tgtAssert(rescaleShader_, "bypass shader not loaded");
    rescaleShader_->activate();

    setGlobalShaderParameters(rescaleShader_, 0, destPort.getSize());

    rescaleShader_->setIgnoreUniformLocationError(true);
    // provide values needed for correct depth value calculation
    float n = cameraProperty_.get().getNearDist();
    float f = cameraProperty_.get().getFarDist();
    rescaleShader_->setUniform("const_to_z_e_1", 0.5f + 0.5f*((f+n)/(f-n)));
    rescaleShader_->setUniform("const_to_z_e_2", ((f-n)/(f*n)));
    rescaleShader_->setUniform("const_to_z_w_1", ((f*n)/(f-n)));
    rescaleShader_->setUniform("const_to_z_w_2", 0.5f*((f+n)/(f-n))+0.5f);
    rescaleShader_->setIgnoreUniformLocationError(false);

    // bind input rendering to texture units
    tgt::TextureUnit colorUnit, depthUnit;
    srcPort.bindColorTexture(colorUnit.getEnum(), GL_LINEAR);
    srcPortDepth.bindColorTexture(depthUnit.getEnum(), GL_LINEAR);
    srcPort.setTextureParameters(rescaleShader_, "texParams_");
    rescaleShader_->setUniform("colorTex_", colorUnit.getUnitNumber());
    rescaleShader_->setUniform("depthTex_", depthUnit.getUnitNumber());

    tgt::TextureUnit depthUnitEntry, depthUnitExit;
    entryPoints.bindDepthTexture(depthUnitEntry, GL_LINEAR);
    exitPoints.bindDepthTexture(depthUnitExit, GL_LINEAR);
    entryPoints.setTextureParameters(rescaleShader_, "texParamsEEP_");
    rescaleShader_->setUniform("depthTexEntry_", depthUnitEntry.getUnitNumber());
    rescaleShader_->setUniform("depthTexExit_", depthUnitExit.getUnitNumber());

    // render screen aligned quad
    renderQuad();

    // cleanup
    rescaleShader_->deactivate();
    destPort.deactivateTarget();
    tgt::TextureUnit::setZeroUnit();

    // check for OpenGL errors
    LGL_ERROR;
}

void SingleOctreeRaycasterCL::invalidateNodeBuffer() {
    nodeBufferInvalid_ = true;
}

void SingleOctreeRaycasterCL::invalidateBrickBuffer() {
    brickBufferInvalid_ = true;
}

void SingleOctreeRaycasterCL::invalidateTransFuncTex() {
    transFuncTexInvalid_ = true;
}

void SingleOctreeRaycasterCL::resetRenderState() {
    // if not in interaction mode, reset to first full-quality state (first full-frame)
    // otherwise, reset to first interaction frame
    renderState_ = (QualityMode.isInteractionMode() ? RENDER_STATE_FIRST_INTERACTION_FRAME : RENDER_STATE_FIRST_FULL_FRAME);
}

void SingleOctreeRaycasterCL::updatePropertyConfiguration() {
    const VolumeOctreeBase* octree = 0;
    if (volumeInport_.hasData())
        octree = volumeInport_.getData()->getRepresentation<VolumeOctreeBase>();

    size_t numChannels = 1;
    if (octree)
        numChannels = octree->getNumChannels();

    transferFunc1_.setVisibleFlag(numChannels >= 2);
    transferFunc2_.setVisibleFlag(numChannels >= 3);
    transferFunc3_.setVisibleFlag(numChannels >= 4);

    channelShift1_.setVisibleFlag(numChannels >= 2);
    channelShift2_.setVisibleFlag(numChannels >= 3);
    channelShift3_.setVisibleFlag(numChannels >= 4);

    bool channelShift = applyChannelShift_.get();
    channelShift0_.setReadOnlyFlag(!channelShift);
    channelShift1_.setReadOnlyFlag(!channelShift);
    channelShift2_.setReadOnlyFlag(!channelShift);
    channelShift3_.setReadOnlyFlag(!channelShift);
    resetChannelShift_.setReadOnlyFlag(!channelShift);

    // limit interaction quality to render quality
    int renderingQualityIndex = renderingQuality_.getSelectedIndex();
    int interactionQualityIndex = interactionQuality_.getSelectedIndex();
    if (interactionQualityIndex <= renderingQualityIndex-2) //< account for "full" and "very high" quality option not available in interaction quality
        interactionQuality_.selectByIndex(renderingQualityIndex-2);

    nodeLevelReduction_.setReadOnlyFlag(renderingQuality_.getValue() == "Full");
}

void SingleOctreeRaycasterCL::resetChannelShift() {
    channelShift0_.set(tgt::vec3::zero);
    channelShift1_.set(tgt::vec3::zero);
    channelShift2_.set(tgt::vec3::zero);
    channelShift3_.set(tgt::vec3::zero);
}

}   // namespace
