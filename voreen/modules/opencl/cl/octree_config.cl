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

//set guard to prevent double include
#ifndef OCTREE_CONFIG_FILE
#define OCTREE_CONFIG_FILE

// structs
typedef struct {
    float scale_;
    float offset_;
} RealWorldMapping;

typedef struct {
    ulong value_;
    ulong offset_;
    float3 llf_;
    float3 urb_;
} OctreeNode;

typedef struct {
    float intensity[4];
    float firsthit;
} RayResult;

typedef struct {
    float4 color;
    float param;
    RayResult current;
    RayResult pending;
} RayInfo;
__constant uint RAYINFO_NUM_ELEMENTS = 15; //< number of single floats stored in an RayInfo object

// octree properties (set by CPU)
__constant uint  OCTREE_DIMENSIONS  = OCTREE_DIMENSIONS_DEF;    //< voxel dimensions of the octree (cubic, power-of-two, >= volume dim)
__constant uint  OCTREE_BRICKDIM    = OCTREE_BRICKDIM_DEF;      //< brick dimensions (cubic, power-of-two)
__constant uint  OCTREE_DEPTH       = OCTREE_DEPTH_DEF;         //< number of levels of the octree
__constant uint  OCTREE_NUMCHANNELS = OCTREE_NUMCHANNELS_DEF;   //< number of channels of the octree [1..4]

//--------------------------------------------------
// 64-bit masks for node buffer entries (set by CPU)

// flag indicating whether a node is homogeneous
__constant ulong MASK_HOMOGENEOUS           = MASK_HOMOGENEOUS_DEF;
__constant uint  MASK_HOMOGENEOUS_SHIFT     = MASK_HOMOGENEOUS_SHIFT_DEF;
__constant uint  MASK_HOMOGENEOUS_NUMBITS   = MASK_HOMOGENEOUS_NUMBITS_DEF;

// flag indicating whether a node's brick is in the GPU buffer
__constant ulong MASK_INBRICKPOOL           = MASK_INBRICKPOOL_DEF;
__constant uint  MASK_INBRICKPOOL_SHIFT     = MASK_INBRICKPOOL_SHIFT_DEF;
__constant uint  MASK_INBRICKPOOL_NUMBITS   = MASK_INBRICKPOOL_NUMBITS_DEF;

// offset of a node's child group in the node buffer
__constant ulong MASK_CHILD                 = MASK_CHILD_DEF;
__constant uint  MASK_CHILD_SHIFT           = MASK_CHILD_SHIFT_DEF;
__constant uint  MASK_CHILD_NUMBITS         = MASK_CHILD_NUMBITS_DEF;

// address of a node's brick in the GPU buffer
__constant ulong MASK_BRICK                 = MASK_BRICK_DEF;
__constant uint  MASK_BRICK_SHIFT           = MASK_BRICK_SHIFT_DEF;
__constant uint  MASK_BRICK_NUMBITS         = MASK_BRICK_NUMBITS_DEF;

// node's avg values for all four channels
__constant ulong MASK_AVG_0                 = MASK_AVG_0_DEF;
__constant uint  MASK_AVG_0_SHIFT           = MASK_AVG_0_SHIFT_DEF;
__constant uint  MASK_AVG_0_NUMBITS         = MASK_AVG_0_NUMBITS_DEF;
__constant ulong MASK_AVG_1                 = MASK_AVG_1_DEF;
__constant uint  MASK_AVG_1_SHIFT           = MASK_AVG_1_SHIFT_DEF;
__constant uint  MASK_AVG_1_NUMBITS         = MASK_AVG_1_NUMBITS_DEF;
__constant ulong MASK_AVG_2                 = MASK_AVG_2_DEF;
__constant uint  MASK_AVG_2_SHIFT           = MASK_AVG_2_SHIFT_DEF;
__constant uint  MASK_AVG_2_NUMBITS         = MASK_AVG_2_NUMBITS_DEF;
__constant ulong MASK_AVG_3                 = MASK_AVG_3_DEF;
__constant uint  MASK_AVG_3_SHIFT           = MASK_AVG_3_SHIFT_DEF;
__constant uint  MASK_AVG_3_NUMBITS         = MASK_AVG_3_NUMBITS_DEF;

// 8-bit masks for node flag buffer entries
__constant uchar MASK_BRICK_INUSE     = MASK_BRICK_INUSE_DEF;
__constant uchar MASK_BRICK_REQUESTED = MASK_BRICK_REQUESTED_DEF;
__constant uchar MASK_MIP_0_USED      = MASK_MIP_0_USED_DEF;
__constant uchar MASK_MIP_1_USED      = MASK_MIP_1_USED_DEF;
__constant uchar MASK_MIP_2_USED      = MASK_MIP_2_USED_DEF;
__constant uchar MASK_MIP_3_USED      = MASK_MIP_3_USED_DEF;
__constant uchar MASK_MIP_A_USED      = MASK_MIP_0_USED_DEF
                                          #if OCTREE_NUMCHANNELS_DEF > 1
                                              | MASK_MIP_1_USED_DEF
                                          #endif
                                          #if OCTREE_NUMCHANNELS_DEF > 2
                                              | MASK_MIP_2_USED_DEF
                                          #endif
                                          #if OCTREE_NUMCHANNELS_DEF > 3
                                              | MASK_MIP_3_USED_DEF
                                          #endif
                                          ;
// ---------------------
// node access functions

bool isHomogeneous(const ulong node) {
    return (node & MASK_HOMOGENEOUS) > 0;
}

bool hasBrick(const ulong node) {
    return (node & MASK_INBRICKPOOL) > 0;
}

/// Returns the offset of the node's child group in the node buffer, i.e. the offset of its first child.
ulong getNodeChildGroupOffset(const ulong node) {
    return ((node & MASK_CHILD) >> MASK_CHILD_SHIFT);
}

/// Returns true, if the node has no children, i.e. if its child group pointer is 0.
bool isLeafNode(const ulong node) {
    return (getNodeChildGroupOffset(node) == 0);
}

/// Returns the address of the node's brick in the brick buffer.
/// If the node has no brick, its brick offset is undefined.
ulong getNodeBrickAddress(const ulong node) {
    return ((node & MASK_BRICK) >> MASK_BRICK_SHIFT);
}

__global const ushort* getNodeBrick(const ulong node, __global const ushort* brickBuffer) {
    ulong brickOffset = getNodeBrickAddress(node);
    return &brickBuffer[brickOffset * (OCTREE_BRICKDIM*OCTREE_BRICKDIM*OCTREE_BRICKDIM * OCTREE_NUMCHANNELS)];
}

/// Returns the node's avg value for a specific channel, converted to a normalized float.
float getNodeAvgValue0(const ulong node) {
    ulong avgInt = (node & MASK_AVG_0) >> MASK_AVG_0_SHIFT;
    float avgNorm = clamp((float)(avgInt) * 1.f/(float)((1 << MASK_AVG_0_NUMBITS)-1), 0.f, 1.f);
    //float avgNorm = clamp((float)(avgInt) * (1.f/4095.f), 0.f, 1.f); //< for 12 bit precision
    return avgNorm;
}
float getNodeAvgValue1(const ulong node) {
    ulong avgInt = (node & MASK_AVG_1) >> MASK_AVG_1_SHIFT;
    float avgNorm = clamp((float)(avgInt) * 1.f/(float)((1 << MASK_AVG_1_NUMBITS)-1), 0.f, 1.f);
    return avgNorm;
}
float getNodeAvgValue2(const ulong node) {
    ulong avgInt = (node & MASK_AVG_2) >> MASK_AVG_2_SHIFT;
    float avgNorm = clamp((float)(avgInt) * 1.f/(float)((1 << MASK_AVG_2_NUMBITS)-1), 0.f, 1.f);
    return avgNorm;
}
float getNodeAvgValue3(const ulong node) {
    ulong avgInt = (node & MASK_AVG_3) >> MASK_AVG_3_SHIFT;
    float avgNorm = clamp((float)(avgInt) * 1.f/(float)((1 << MASK_AVG_3_NUMBITS)-1), 0.f, 1.f);
    return avgNorm;
}

/// Returns the node's avg values for all channels, converted to a normalized float.
void getNodeAvgValues(const ulong node, float* channelAvgValues) {
    channelAvgValues[0] = getNodeAvgValue0(node);
#if OCTREE_NUMCHANNELS_DEF > 1
    channelAvgValues[1] = getNodeAvgValue1(node);
#endif
#if OCTREE_NUMCHANNELS_DEF > 2
    channelAvgValues[2] = getNodeAvgValue2(node);
#endif
#if OCTREE_NUMCHANNELS_DEF > 3
    channelAvgValues[3] = getNodeAvgValue3(node);
#endif
}

// ----------------------------
// flag buffer access functions

bool isBrickInUse(const uchar flagEntry) {
    return (flagEntry & MASK_BRICK_INUSE) > 0;
}
void setBrickUsed(global uchar* flagEntry, const bool used) {
    *flagEntry = (*flagEntry & ~MASK_BRICK_INUSE) | (used ? MASK_BRICK_INUSE : 0);
}
bool isBrickRequested(const uchar flagEntry) {
    return (flagEntry & MASK_BRICK_REQUESTED) > 0;
}
void setBrickRequested(global uchar* flagEntry, const bool requested) {
    *flagEntry = (*flagEntry & ~MASK_BRICK_REQUESTED) | (requested ? MASK_BRICK_REQUESTED : 0);
}
bool hasNodeBeenUsedInMIPAll(const uchar flagEntry) {
    return (flagEntry & MASK_MIP_A_USED) > 0;
}
bool hasNodeBeenUsedInMIP(const uchar flagEntry, const int ch) {
    if(ch == 0)
        return (flagEntry & MASK_MIP_0_USED) > 0;
    else if(ch == 1)
        return (flagEntry & MASK_MIP_1_USED) > 0;
    else if(ch == 2)
        return (flagEntry & MASK_MIP_2_USED) > 0;
    else if(ch == 3)
        return (flagEntry & MASK_MIP_3_USED) > 0;
    return false;
}
void setNodeUsedInMIPAll(global uchar* flagEntry, const bool used) {
    *flagEntry = (*flagEntry & ~MASK_MIP_A_USED) | (used ? MASK_MIP_A_USED : 0);
}
void setNodeUsedInMIP(global uchar* flagEntry, const bool used, const int ch) {
    if(ch == 0)
        *flagEntry = (*flagEntry & ~MASK_MIP_0_USED) | (used ? MASK_MIP_0_USED : 0);
    else if(ch == 1)
        *flagEntry = (*flagEntry & ~MASK_MIP_1_USED) | (used ? MASK_MIP_1_USED : 0);
    else if(ch == 2)
        *flagEntry = (*flagEntry & ~MASK_MIP_2_USED) | (used ? MASK_MIP_2_USED : 0);
    else if(ch == 3)
        *flagEntry = (*flagEntry & ~MASK_MIP_3_USED) | (used ? MASK_MIP_3_USED : 0);
}

//--------------------------------------
// RayInfo storage

void fetchRayFromBuffer(const global float* buffer, const int2 pos, const uint2 bufferDim, RayInfo* ray) {
    const uint bufferIndex = (pos.y*bufferDim.x + pos.x) * RAYINFO_NUM_ELEMENTS;
    uint offset = 0;

    ray->color.x  = buffer[bufferIndex + offset++];
    ray->color.y  = buffer[bufferIndex + offset++];
    ray->color.z  = buffer[bufferIndex + offset++];
    ray->color.w  = buffer[bufferIndex + offset++];

    ray->param    = buffer[bufferIndex + offset++];
    ray->current.firsthit = buffer[bufferIndex + offset++];
    ray->pending.firsthit = buffer[bufferIndex + offset++];

    for(int c=0; c<4; ++c) {
        ray->current.intensity[c] = buffer[bufferIndex + offset++];
        ray->pending.intensity[c] = buffer[bufferIndex + offset++];
    }
}

void writeRayToBuffer(global float* buffer, const int2 pos, const uint2 bufferDim, const RayInfo* ray) {
    const uint bufferIndex = (pos.y*bufferDim.x + pos.x) * RAYINFO_NUM_ELEMENTS;
    uint offset = 0;

    buffer[bufferIndex + offset++] = ray->color.x;
    buffer[bufferIndex + offset++] = ray->color.y;
    buffer[bufferIndex + offset++] = ray->color.z;
    buffer[bufferIndex + offset++] = ray->color.w;

    buffer[bufferIndex + offset++] = ray->param;

    buffer[bufferIndex + offset++] = ray->current.firsthit;
    buffer[bufferIndex + offset++] = ray->pending.firsthit;

    for(int c=0; c<4; ++c) {
        buffer[bufferIndex + offset++] = ray->current.intensity[c];
        buffer[bufferIndex + offset++] = ray->pending.intensity[c];
    }
}

//--------------------------------------
// low-level helper functions

uint cubicCoordToLinear(uint3 coord, const uint3 dim) {
    //coord = clamp(coord, (uint3)0, dim-(uint3)(1));
    return coord.z*dim.x*dim.y + coord.y*dim.x + coord.x;
}

uint3 linearCoordToCubic(uint coord, const uint3 dim) {
    uint3 result;
    result.z = coord / (dim.x*dim.y);
    coord = coord % (dim.x*dim.y);
    result.y = coord / dim.x;
    result.x = coord % dim.x;
    return result;
}

float minFloat3(float3 vec) {
    return min(vec.x, min(vec.y, vec.z));
}

float minFloat4(float4 vec) {
    return min(vec.x, min(vec.y, min(vec.z, vec.w)));
}

/// multiplication of a 4x4 row-matrix by a 4-component vector
float4 matVecMult4x4(const float16 mat, const float4 vec) {
    return (float4)(dot(mat.s0123, vec),
                    dot(mat.s4567, vec),
                    dot(mat.s89AB, vec),
                    dot(mat.sCDEF, vec));
}

/// multiplication of a 4x4 row-matrix by a 3-component vector (1.0 is added as w-component)
float3 matVecMult4x3(const float16 mat, const float3 vec) {
    float4 res4 = matVecMult4x4(mat, (float4)(vec, 1.f));
    return (res4 / res4.w).xyz;
}

//-------------------------------------
// brick filtering

float filterBrickNearest(float3 samplePosInBrick, const uint channel, __global const ushort* brick) {
    uint3 sampleCoordsInBrick;
    sampleCoordsInBrick.x = clamp((uint)round(samplePosInBrick.x*(float)(OCTREE_BRICKDIM)-0.5f), (uint)0, (uint)(OCTREE_BRICKDIM-1));
    sampleCoordsInBrick.y = clamp((uint)round(samplePosInBrick.y*(float)(OCTREE_BRICKDIM)-0.5f), (uint)0, (uint)(OCTREE_BRICKDIM-1));
    sampleCoordsInBrick.z = clamp((uint)round(samplePosInBrick.z*(float)(OCTREE_BRICKDIM)-0.5f), (uint)0, (uint)(OCTREE_BRICKDIM-1));

    uint inBrickOffset = cubicCoordToLinear(sampleCoordsInBrick, (uint3)(OCTREE_BRICKDIM))*OCTREE_NUMCHANNELS + channel;
    return brick[inBrickOffset] * (1.f/(float)(1 << 16)); // 65535.f
}

float filterBrickLinear(float3 samplePosInBrick, const uint channel, __global const ushort* brick) {
    const float3 sampleCoordsInBrick = samplePosInBrick * ((float3)(OCTREE_BRICKDIM)-1.f);

    // cubic coordinates of llf/urb corner voxels surrounding the sample pos
    uint3 llf = (uint3)(floor(sampleCoordsInBrick.x),
                        floor(sampleCoordsInBrick.y),
                        floor(sampleCoordsInBrick.z));
    uint3 urb = llf + (uint3)(1);
    llf = clamp(llf, (uint3)(0), (uint3)(OCTREE_BRICKDIM-1));
    urb = min(urb, (uint3)(OCTREE_BRICKDIM-1));

    // linear coodinates of the eight corner voxels
    uint llfOffset = cubicCoordToLinear((uint3)(llf.x, llf.y, llf.z), (uint3)(OCTREE_BRICKDIM))*OCTREE_NUMCHANNELS + channel;
    uint llbOffset = cubicCoordToLinear((uint3)(llf.x, llf.y, urb.z), (uint3)(OCTREE_BRICKDIM))*OCTREE_NUMCHANNELS + channel;
    uint lrfOffset = cubicCoordToLinear((uint3)(llf.x, urb.y, llf.z), (uint3)(OCTREE_BRICKDIM))*OCTREE_NUMCHANNELS + channel;
    uint lrbOffset = cubicCoordToLinear((uint3)(llf.x, urb.y, urb.z), (uint3)(OCTREE_BRICKDIM))*OCTREE_NUMCHANNELS + channel;
    uint ulfOffset = cubicCoordToLinear((uint3)(urb.x, llf.y, llf.z), (uint3)(OCTREE_BRICKDIM))*OCTREE_NUMCHANNELS + channel;
    uint ulbOffset = cubicCoordToLinear((uint3)(urb.x, llf.y, urb.z), (uint3)(OCTREE_BRICKDIM))*OCTREE_NUMCHANNELS + channel;
    uint urfOffset = cubicCoordToLinear((uint3)(urb.x, urb.y, llf.z), (uint3)(OCTREE_BRICKDIM))*OCTREE_NUMCHANNELS + channel;
    uint urbOffset = cubicCoordToLinear((uint3)(urb.x, urb.y, urb.z), (uint3)(OCTREE_BRICKDIM))*OCTREE_NUMCHANNELS + channel;

    // voxel intensities of the eight corner voxels
    float llfVoxel = brick[llfOffset];
    float llbVoxel = brick[llbOffset];
    float lrfVoxel = brick[lrfOffset];
    float lrbVoxel = brick[lrbOffset];
    float ulfVoxel = brick[ulfOffset];
    float ulbVoxel = brick[ulbOffset];
    float urfVoxel = brick[urfOffset];
    float urbVoxel = brick[urbOffset];

    // interpolate between corner voxels
    float3 p = sampleCoordsInBrick - (float3)(llf.x,llf.y,llf.z);
    float sample = llfVoxel * (1.f-p.x)*(1.f-p.y)*(1.f-p.z) +
                   llbVoxel * (1.f-p.x)*(1.f-p.y)*(    p.z) +
                   lrfVoxel * (1.f-p.x)*(    p.y)*(1.f-p.z) +
                   lrbVoxel * (1.f-p.x)*(    p.y)*(    p.z) +
                   ulfVoxel * (    p.x)*(1.f-p.y)*(1.f-p.z) +
                   ulbVoxel * (    p.x)*(1.f-p.y)*(    p.z) +
                   urfVoxel * (    p.x)*(    p.y)*(1.f-p.z) +
                   urbVoxel * (    p.x)*(    p.y)*(    p.z);

    sample *= (1.f/(float)(1 << 16)); // 65535.f;

    return sample;
}

void filterBrick(float3 samplePosInBrick, __global const ushort* const brick, float* channelIntensities) {
#ifdef TEXTURE_FILTER_LINEAR
    channelIntensities[0] = filterBrickLinear(samplePosInBrick, 0, brick);
    #if OCTREE_NUMCHANNELS_DEF > 1
        channelIntensities[1] = filterBrickLinear(samplePosInBrick, 1, brick);
    #endif
    #if OCTREE_NUMCHANNELS_DEF > 2
        channelIntensities[2] = filterBrickLinear(samplePosInBrick, 2, brick);
    #endif
    #if OCTREE_NUMCHANNELS_DEF > 3
        channelIntensities[3] = filterBrickLinear(samplePosInBrick, 3, brick);
    #endif
#else
    channelIntensities[0] = filterBrickNearest(samplePosInBrick, 0, brick);
    #if OCTREE_NUMCHANNELS_DEF > 1
        channelIntensities[1] = filterBrickNearest(samplePosInBrick, 1, brick);
    #endif
    #if OCTREE_NUMCHANNELS_DEF > 2
        channelIntensities[2] = filterBrickNearest(samplePosInBrick, 2, brick);
    #endif
    #if OCTREE_NUMCHANNELS_DEF > 3
        channelIntensities[3] = filterBrickNearest(samplePosInBrick, 3, brick);
    #endif
#endif
}


//-------------------------------------
// transfer function

__constant sampler_t transFuncSampler = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

float4 applyChannelTransFunc(const float normIntensity, read_only image2d_t transFuncTex,
    const float2 transFuncDomain, const RealWorldMapping realWorldMapping, const uint channel)
{
    // transform intensity value into TF domain
    float realWorldIntensity = (normIntensity * realWorldMapping.scale_) + realWorldMapping.offset_;
    float tfIntensity = (realWorldIntensity - transFuncDomain.x) / (transFuncDomain.y - transFuncDomain.x);
    float yCoord = (float)(channel*3 + 1) / (float)(3*OCTREE_NUMCHANNELS);

    // apply transfer function
    float4 color = read_imagef(transFuncTex, transFuncSampler, (float2)(tfIntensity, yCoord));
    return color;
}

void applyTransFuncs(const float* channelIntensities, read_only image2d_t transFuncTex,
    const float8 transFuncDomains, const RealWorldMapping realWorldMapping, float4* channelColors)
{
    channelColors[0] = applyChannelTransFunc(channelIntensities[0], transFuncTex,
        (float2)(transFuncDomains.s0, transFuncDomains.s1), realWorldMapping, 0);
#if OCTREE_NUMCHANNELS_DEF > 1
    channelColors[1] = applyChannelTransFunc(channelIntensities[1], transFuncTex,
        (float2)(transFuncDomains.s2, transFuncDomains.s3), realWorldMapping, 1);
#endif
#if OCTREE_NUMCHANNELS_DEF > 2
    channelColors[2] = applyChannelTransFunc(channelIntensities[2], transFuncTex,
        (float2)(transFuncDomains.s4, transFuncDomains.s5), realWorldMapping, 2);
#endif
#if OCTREE_NUMCHANNELS_DEF > 3
    channelColors[3] = applyChannelTransFunc(channelIntensities[3], transFuncTex,
        (float2)(transFuncDomains.s6, transFuncDomains.s7), realWorldMapping, 3);
#endif

}

//-------------------------------------

float3 computeNodeExitPoint(const float3 nodeLLF, const float3 nodeURB, const float3 nodeEntry, const float3 rayDir) {
    //tgtAssert(inRange(nodeEntry, nodeLLF, nodeURB), "node entry point outside node");

    float3 exitPlanes = (float3)(rayDir.x >= 0.f ? nodeURB.x : nodeLLF.x,
                                    rayDir.y >= 0.f ? nodeURB.y : nodeLLF.y,
                                    rayDir.z >= 0.f ? nodeURB.z : nodeLLF.z);
    float3 tNodeExit;
    tNodeExit.x = rayDir.x != 0.f ? ((exitPlanes.x - nodeEntry.x) / rayDir.x) : 1e6f;
    tNodeExit.y = rayDir.y != 0.f ? ((exitPlanes.y - nodeEntry.y) / rayDir.y) : 1e6f;
    tNodeExit.z = rayDir.z != 0.f ? ((exitPlanes.z - nodeEntry.z) / rayDir.z) : 1e6f;
    //tgtAssert(tgt::hand(tgt::greaterThanEqual(tNodeExit, vec3::zero)), "at least one negative node exit parameter");

    float tNodeExitMin = minFloat3(tNodeExit);
    //tgtAssert(inRange(tNodeExitMin, 0.f, 1.f), "minimum node exit parameter outside range [0.0;1.0]");

    float3 nodeExit = nodeEntry + (tNodeExitMin - 1e-6f)*rayDir;
    //tgtAssert(inRange(nodeExit, nodeLLF, nodeURB), "node exit point outside node");
    return nodeExit;
}

//-------------------------------------

// Get the node at the provided position and requestlevel (or first available
// parent node, if no node at the requested level exists)
OctreeNode getNodeAtSamplePos(const float3 samplePos, const uint requestedLevel, __global const ulong* const nodeBuffer,
    uint* returnedLevel)
{

    OctreeNode currentNode;
    currentNode.value_ = nodeBuffer[0];
    currentNode.offset_ = 0;
    currentNode.llf_ = (float3)(0.f, 0.f, 0.f);
    currentNode.urb_ = (float3)(1.f, 1.f, 1.f);

    uint currentLevel = 0;

    // iteratively descent to children
    while (true) {
        const ulong childGroupOffset = getNodeChildGroupOffset(currentNode.value_);
        if (currentLevel == requestedLevel || childGroupOffset == 0) { // current node level requested, or current node is leaf => stop descent
            *returnedLevel = currentLevel;
            return currentNode;
        }
        else { // descent to next level

            float3 nodeDim = currentNode.urb_ - currentNode.llf_;
            float3 nodeHalfDim = nodeDim * (float3)(0.5f);
            float3 nodeCenter = currentNode.llf_ + nodeHalfDim;

            // select child node
            uint3 childNodeGroupCoord; //< "coordinates" of the child node within its child group
            childNodeGroupCoord.x = (samplePos.x >= nodeCenter.x) ? 1 : 0;
            childNodeGroupCoord.y = (samplePos.y >= nodeCenter.y) ? 1 : 0;
            childNodeGroupCoord.z = (samplePos.z >= nodeCenter.z) ? 1 : 0;
            currentNode.offset_ = childGroupOffset + cubicCoordToLinear(childNodeGroupCoord, (uint3)(2));
            currentNode.value_ = nodeBuffer[currentNode.offset_];
            currentLevel++;

            // update LLF, URB
            currentNode.llf_.x += nodeHalfDim.x*(float)childNodeGroupCoord.x;
            currentNode.llf_.y += nodeHalfDim.y*(float)childNodeGroupCoord.y;
            currentNode.llf_.z += nodeHalfDim.z*(float)childNodeGroupCoord.z;
            currentNode.urb_ = currentNode.llf_ + nodeHalfDim;
        }
    }

    // should not get here
    OctreeNode dummy;
    return dummy;
}

// Get the best (i.e., highest resolution) node (up to requestedLevel) with a
// brick on the GPU at the provided sample position. Also request child nodes
// at the position if their bricks are not on the GPU currently.
OctreeNode getBestAvailableNodeAndRequest(const float3 samplePos, const uint requestedLevel, __global const ulong* const nodeBuffer,
        __global uchar* const brickFlagBuffer, uint* returnedLevel)
{

    // Storage for (up to MAX_REQUEST_LEVEL) offsets of nodes that have a
    // higher resolution brick than the highest resolution currently available
#define MAX_REQUEST_LEVEL 8 // log2(256) => then it is definitely better to use node averages
    ulong usefulNodeOffsets[MAX_REQUEST_LEVEL];
    uint nextUsefulNodeOffsetsEntry = 0;

    OctreeNode currentNode;
    currentNode.value_ = nodeBuffer[0];
    currentNode.offset_ = 0;
    currentNode.llf_ = (float3)(0.f, 0.f, 0.f);
    currentNode.urb_ = (float3)(1.f, 1.f, 1.f);

    uint currentLevel = 0;

    OctreeNode bestAvailableNode = currentNode;
    uint bestAvailableLevel = 0;

    // iteratively descent to children
    while (true) {
        const ulong childGroupOffset = getNodeChildGroupOffset(currentNode.value_);

        if(!isHomogeneous(currentNode.value_) && !hasBrick(currentNode.value_)) {
            // Node has an associated brick, but brick is not on GPU

            // It is always true that currentLevel <= requestedLevel, so the
            // following ensure that we store the lowest (highest resolution)
            // MAX_REQUEST_LEVEL levels at most
            if(currentLevel+MAX_REQUEST_LEVEL > requestedLevel) {
                usefulNodeOffsets[nextUsefulNodeOffsetsEntry] = currentNode.offset_;
                nextUsefulNodeOffsetsEntry += 1;
            }
        } else {
            bestAvailableNode = currentNode;
            bestAvailableLevel = currentLevel;

            // All nodes above this level are not useful anymore
            nextUsefulNodeOffsetsEntry = 0;
        }

        if (currentLevel == requestedLevel || childGroupOffset == 0) {
            // current node level requested, or current node is leaf => stop descent

            // Request bricks of all useful nodes (i.e., those that have a
            // higher than currently available resolution brick)
            for(uint i=0; i<nextUsefulNodeOffsetsEntry; ++i) {
                setBrickRequested(brickFlagBuffer + usefulNodeOffsets[i], true);
            }
            *returnedLevel = bestAvailableLevel;
            return bestAvailableNode;
        }
        else { // descent to next level
            float3 nodeDim = currentNode.urb_ - currentNode.llf_;
            float3 nodeHalfDim = nodeDim * (float3)(0.5f);
            float3 nodeCenter = currentNode.llf_ + nodeHalfDim;

            // select child node
            uint3 childNodeGroupCoord; //< "coordinates" of the child node within its child group
            childNodeGroupCoord.x = (samplePos.x >= nodeCenter.x) ? 1 : 0;
            childNodeGroupCoord.y = (samplePos.y >= nodeCenter.y) ? 1 : 0;
            childNodeGroupCoord.z = (samplePos.z >= nodeCenter.z) ? 1 : 0;
            currentNode.offset_ = childGroupOffset + cubicCoordToLinear(childNodeGroupCoord, (uint3)(2));
            currentNode.value_ = nodeBuffer[currentNode.offset_];
            currentLevel++;

            // update LLF, URB
            currentNode.llf_.x += nodeHalfDim.x*(float)childNodeGroupCoord.x;
            currentNode.llf_.y += nodeHalfDim.y*(float)childNodeGroupCoord.y;
            currentNode.llf_.z += nodeHalfDim.z*(float)childNodeGroupCoord.z;
            currentNode.urb_ = currentNode.llf_ + nodeHalfDim;
        }
    }

    // should not get here
    OctreeNode dummy;
    return dummy;
}

OctreeNode fetchNextRayNode(const float3 samplePos, const float rayParam, const float3 rayDir, const float samplingStepSize,
    const uint nodeLevel, __global const ulong* const nodeBuffer, __global const ushort* const brickBuffer, __global uchar* const brickFlagBuffer,
    float* exitParam, float* samplingStepSizeNode, bool* nodeHasBrick, __global const ushort** brick)
{

    uint resultNodeLevel = 0;
    *nodeHasBrick = false;
    *brick = 0;

#if defined(USE_BRICKS) && defined(USE_ANCESTOR_NODES) && !defined(DISPLAY_MODE_REFINEMENT)
    uint bestNodeLevel;
    OctreeNode resultNode = getBestAvailableNodeAndRequest(samplePos, nodeLevel, nodeBuffer, brickFlagBuffer, &resultNodeLevel);
    if (!isHomogeneous(resultNode.value_)) { // node not homogeneous => fetch brick
        if(hasBrick(resultNode.value_)) {
            *brick = getNodeBrick(resultNode.value_, brickBuffer);
            *nodeHasBrick = true;
        }
    }
#else
    // retrieve node
    OctreeNode resultNode = getNodeAtSamplePos(samplePos, nodeLevel, nodeBuffer, &resultNodeLevel);

    #ifdef USE_BRICKS
    if (!isHomogeneous(resultNode.value_)) { // node not homogeneous => fetch brick
        *nodeHasBrick = hasBrick(resultNode.value_);
        if (*nodeHasBrick) { //< brick is in GPU buffer => use it
            *brick = getNodeBrick(resultNode.value_, brickBuffer);
        } else { // brick is not in GPU buffer
            setBrickRequested(brickFlagBuffer + resultNode.offset_, true);
        }
    }
    #endif

#endif

    // compute node exit point and ray end parameter for the node
    if (length(rayDir) > 0.f) {
        float3 nodeExit = computeNodeExitPoint(resultNode.llf_, resultNode.urb_, samplePos, rayDir);
        // determine ray parameter of last sample before node exit point
        float tOffset = minFloat3((nodeExit - samplePos) / rayDir);
        tOffset = floor(tOffset / samplingStepSize) * samplingStepSize;
        *exitParam = rayParam + tOffset;
    }
    else {
        *exitParam = rayParam;
    }

#ifdef ADAPTIVE_SAMPLING
    // adapt sampling step size to current node level/resolution
    *samplingStepSizeNode = samplingStepSize * (float)(1<<(OCTREE_DEPTH-1-min(resultNodeLevel, OCTREE_DEPTH-1)));
#else
    *samplingStepSizeNode = samplingStepSize;
#endif
    return resultNode;
}

#endif
//end of guard
