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

//#pragma OPENCL EXTENSION cl_amd_printf : enable

//include basic functions
#include "octree_config.cl"

//include macros
#ifdef COMPOSITING_MODE_DVR
    #include "octree_dvr.cl"
#elif COMPOSITING_MODE_MIP
    #include "octree_mip.cl"
#elif COMPOSITING_MODE_MOP
    #include "octree_mop.cl"
#else
    ERROR: unknown compositing mode
#endif

// sampler configuration for 2D texture lookup
__constant sampler_t imageSampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

__constant float SAMPLING_BASE_INTERVAL_RCP = 200.f;


//-------------------------------------
// ray traversal functions

#ifdef APPLY_CHANNEL_SHIFT
    void traverseRayCS(const float3 entry, const float3 exit, const uint2 viewportSize,
                     __global const ulong* nodeBuffer,
                     __global const ushort* brickBuffer,
                     __global uchar* brickFlagBuffer,
                     const int nodeLevelOfDetail,

                     const RealWorldMapping realWorldMapping,
                     read_only image2d_t transFunc,
                     const float8 transFuncDomains,

                     const float samplingStepSize,

                     const float3 channelShifts[OCTREE_NUMCHANNELS_DEF],

                     RayInfo* ray)
    {
        // calculate ray parameters
        float3 direction = exit - entry;
        float tEnd = length(direction);
        if (tEnd > 0.f)
            direction = normalize(direction);

        // current node information
        OctreeNode currentNode[OCTREE_NUMCHANNELS_DEF];                     ///< current octree node (per channel)
        uint currentNodeLevel[OCTREE_NUMCHANNELS_DEF];                      ///< level of the current (fetched) node (per channel)
        float3 currentNodeDimRec[OCTREE_NUMCHANNELS_DEF];                   ///< reciproke dimension of current node in texture coordinates
        bool currentNodeHasBrick[OCTREE_NUMCHANNELS_DEF];                   ///< does current node have a brick?
        __global const ushort* currentBrick[OCTREE_NUMCHANNELS_DEF];        ///< brick of current node (uint16_t)
        // ray casting loop
        float samplingStepSizeNode = samplingStepSize;    ///< adapted sampling step size for current node
        float tEndCurrentNode[OCTREE_NUMCHANNELS_DEF];    ///< end ray parameter for current node, i.e. t value of last sample within the current node
        //only used in MIP/MOP
        bool sampleNode[OCTREE_NUMCHANNELS_DEF]; ///< only MIP/MOP
        bool rayFinished = true;                 ///< only MIP/MOP

        for (int ch=0; ch<OCTREE_NUMCHANNELS; ch++) {
            currentNodeHasBrick[ch] = false;
            currentBrick[ch] = 0;
            tEndCurrentNode[ch] = -1.f;
            sampleNode[ch] = true;
        }

        // initialize intensity values
#ifdef COMPOSITING_MODE_MOP
        float4 channelColors[OCTREE_NUMCHANNELS_DEF];
        applyTransFuncs(ray->channelIntensities, transFunc, transFuncDomains, realWorldMapping, channelColors);
#endif
        // main ray-casting loop
        while (ray->param <= tEnd) {

            float previousStepSize = samplingStepSizeNode; ///< only used DVR

            // determine sampling point for each channel + fetch next node, if necessary
            float3 sample[OCTREE_NUMCHANNELS_DEF];
            for (int ch=0; ch<OCTREE_NUMCHANNELS; ch++) {

                sample[ch] = entry + ray->param*direction + channelShifts[ch];

                // retrieve next node for ray, if ray has passed last sample within current node (for respective channel)
                if (ray->param > (tEndCurrentNode[ch]+1.e-6f)) {

                    float tSamplingStepSizeNode;
                    currentNode[ch] = fetchNextRayNode(sample[ch], ray->param, direction, samplingStepSize,
                                        nodeLevelOfDetail, nodeBuffer, brickBuffer, brickFlagBuffer,
                                        &tEndCurrentNode[ch], &tSamplingStepSizeNode, &currentNodeHasBrick[ch],
                                        &currentNodeLevel[ch], &currentBrick[ch]);
                    if (ch == 0) // first channel determines sampling rate
                        samplingStepSizeNode = tSamplingStepSizeNode;

            #ifdef DISPLAY_MODE_REFINEMENT
                    //macro define
                    getBrickInRefinementCS
            #else
                    // mark available brick as used (=> should be kept on GPU)
                    if (currentNodeHasBrick[ch])
                        setBrickUsed(brickFlagBuffer + currentNode[ch].offset_, true);
            #endif
                    currentNodeDimRec[ch] = 1.f / (currentNode[ch].urb_-currentNode[ch].llf_);

                } // node retrieval

            } // channel sampling points

            // retrieve sample intensity (for each channel)
            float channelIntensities[OCTREE_NUMCHANNELS_DEF];
            for (int ch=0; ch<OCTREE_NUMCHANNELS; ch++) {
                if (sampleNode[ch]) {
                    //set intensity zero, if we are outside the volume TODO: optimize
                    if(any(clamp(sample[ch],(float3)(0.f,0.f,0.f),(float3)(1.f,1.f,1.f)) != sample[ch])) {
                        channelIntensities[ch] = 0.f;
                    }
                    else {
                        if (currentNodeHasBrick[ch]) { // sample brick
                            float3 samplePosInBrick = (sample[ch] - currentNode[ch].llf_) * currentNodeDimRec[ch];
                            #ifdef TEXTURE_FILTER_LINEAR
                            channelIntensities[ch] = filterBrickLinear(samplePosInBrick, ch, currentBrick[ch]);
                            #else
                            channelIntensities[ch] = filterBrickNearest(samplePosInBrick, ch, currentBrick[ch]);
                            #endif
                        }
                        else { // node has no brick => use average value for DVR, 0 for MIP, resp.
                            float channelIntensitiesTemp[OCTREE_NUMCHANNELS_DEF];
                            getNodeAvgValues(currentNode[ch].value_, channelIntensitiesTemp);
                            channelIntensities[ch] = channelIntensitiesTemp[ch];
                        }
                    }
                } else {
                    channelIntensities[ch] = 0.f;
                }
            }

            //macro define
            applyTFandCombineColors

        //only needed in the adaptive sampling case
    #ifdef ADAPTIVE_SAMPLING
            // switch back to base sampling step size,
            // if node step size would yield a sampling point beyond node exit point (for any channel)
            for (int ch=0; ch<OCTREE_NUMCHANNELS; ch++) {
                if (ray->param+samplingStepSizeNode > (tEndCurrentNode[ch]+1.e-6f)) {
                    samplingStepSizeNode = samplingStepSize;
                    break;
                }
            }
    #endif

            // advance along ray
            ray->param += samplingStepSizeNode;

        } // ray-casting loop

        // Make sure firsthit is in normalized coordinates, i.e., an
        // interpolation parameter between entry and exit points instead of a
        // depth in texture coordinates.
        if(tEnd > 0.0f) {
            ray->firsthit /= tEnd;
        }

        //macro define
        postRaycastingLoop
    }


#else
    void traverseRay(const float3 entry, const float3 exit, const uint2 viewportSize,
                     __global const ulong* nodeBuffer,
                     __global const ushort* brickBuffer,
                     __global uchar* brickFlagBuffer,
                     const int nodeLevelOfDetail,

                     const RealWorldMapping realWorldMapping,
                     read_only image2d_t transFunc,
                     const float8 transFuncDomains,

                     const float samplingStepSize,

                     RayInfo* ray)
    {
        // calculate ray parameters
        float3 direction = exit - entry;
        float tEnd = length(direction);
        if (tEnd > 0.f)
            direction = normalize(direction);

        // current node information
        OctreeNode currentNode;                     ///< current octree node
        uint currentNodeLevel;                      ///< level of the current (fetched) octree node
        float3 currentNodeDimRec;                   ///< reciproke dimension of current node in texture coordinates
        bool currentNodeHasBrick = false;           ///< does current node have a brick?
        __global const ushort* currentBrick = 0;    ///< brick of current node (uint16_t)
        // ray casting loop
        float tEndCurrentNode = -1.f;                   ///< end ray parameter for current node, i.e. t value of last sample within the current node
        float samplingStepSizeNode = samplingStepSize;  ///< adapted sampling step size for current node
        bool sampleNode = true;                   ///< only MIP/MOP
        bool rayFinished = true;                 ///< only MIP/MOP

        // initialize intensity values
#ifdef COMPOSITING_MODE_MOP
        float4 channelColors[OCTREE_NUMCHANNELS_DEF];
        applyTransFuncs(ray->channelIntensities, transFunc, transFuncDomains, realWorldMapping, channelColors);
#endif

        // main ray-casting loop
        while (ray->param <= tEnd) {

            float previousStepSize = samplingStepSizeNode; ///< only used DVR

            float3 sample = entry + ray->param*direction;
            { // bracket to simulate for-loop

            // retrieve next node, if ray has passed last sample within current node
                if (ray->param > (tEndCurrentNode+1.e-6f)) {

                    currentNode = fetchNextRayNode(sample, ray->param, direction, samplingStepSize,
                                            nodeLevelOfDetail, nodeBuffer, brickBuffer, brickFlagBuffer,
                                            &tEndCurrentNode, &samplingStepSizeNode, &currentNodeHasBrick, &currentNodeLevel, &currentBrick);

            #ifdef DISPLAY_MODE_REFINEMENT
                    //macro define
                    getBrickInRefinement
            #else
                    // mark available brick as used (=> should be kept on GPU)
                    if (currentNodeHasBrick)
                        setBrickUsed(brickFlagBuffer + currentNode.offset_, true);
            #endif
                    currentNodeDimRec = 1.f / (currentNode.urb_-currentNode.llf_);
                } // node retrieval

            } // channel sampling points

            // retrieve sample intensity (for each channel)
            float channelIntensities[OCTREE_NUMCHANNELS_DEF];
            if (sampleNode) { //always true in dvr
                //set intensity zero, if we are outside the volume TODO: optimize
                if(any(clamp(sample,(float3)(0.f,0.f,0.f),(float3)(1.f,1.f,1.f)) != sample)) {
                    for (int ch=0; ch<OCTREE_NUMCHANNELS; ch++)
                        channelIntensities[ch] = 0.f;
                }
                else {
                    if (currentNodeHasBrick) { // sample brick
                        float3 samplePosInBrick = (sample - currentNode.llf_) * currentNodeDimRec;
                        filterBrick(samplePosInBrick, currentBrick, channelIntensities);
                    }
                    else { // node has no brick => use average value for DVR, 0 for MIP, resp.
                        getNodeAvgValues(currentNode.value_, channelIntensities);
                    }
                }

                //macro define
                applyTFandCombineColors
            }

        //only needed in the adaptive sampling case
    #ifdef ADAPTIVE_SAMPLING
            // switch back to base sampling step size,
            // if node step size would yield a sampling point beyond node exit point
            if (ray->param+samplingStepSizeNode > (tEndCurrentNode+1.e-6f))
                samplingStepSizeNode = samplingStepSize;
    #endif

            // advance along ray
            ray->param += samplingStepSizeNode;

        } // ray-casting loop

        //macro define
        postRaycastingLoop
    }
#endif

//-------------------------------------------------------------------------------------------
// main kernel

__kernel void render( read_only image2d_t entryTex
                    , read_only image2d_t exitTex
                    , const uint2 viewportSize
                    , const uint coarsenessFactor

                    , const float samplingStepSize
                    , const uint4 volumeDimensions

                    , __global const ulong* nodeBuffer
                    , __global const ushort* brickBuffer
                    , __global uchar* brickFlagBuffer
                    , const uint nodeLevelOfDetail

                    , const RealWorldMapping realWorldMapping
                    , read_only image2d_t transFunc
                    , const float8 transFuncDomains

                    , write_only image2d_t output
                    , write_only image2d_t outputDepth

                    , __global float* rayBuffer
                    , const uint firstRefinementFrame   //< bool flag

#ifdef APPLY_CHANNEL_SHIFT
                    , const float4 channelShift0
                    , const float4 channelShift1
                    , const float4 channelShift2
                    , const float4 channelShift3
#endif
                    )
{

    // determine fragment position: each workgroup is assigned to one tile of the image
    const uint2 localSize = (uint2)(get_local_size(0), get_local_size(1));
    const uint2 groupID = (uint2)(get_group_id(0), get_group_id(1));
    const uint2 localID = (uint2)(get_local_id(0), get_local_id(1));
    int2 fragPos = (int2)(groupID.x*localSize.x + localID.x, groupID.y*localSize.y + localID.y);

    if (fragPos.x >= viewportSize.x || fragPos.y >= viewportSize.y)
        return;

    // read entry/exit points
    float4 entry = read_imagef(entryTex, imageSampler, fragPos*(int2)(coarsenessFactor));
    float4 exit = read_imagef(exitTex, imageSampler, fragPos*(int2)(coarsenessFactor));

    // transform entry/exit points first from texture to voxel coordinates and then from NPOT to POT
    // (URB voxel of the volume differs from the URB of the octree for NPOT volumes)
    const float4 textureOffset = (float4) (0.5, 0.5, 0.5, 0.0);
    const float4 volDimF = (float4) ((float) volumeDimensions.x,
                                     (float) volumeDimensions.y,
                                     (float) volumeDimensions.z,
                                     1.f);
    const float4 volToPOTF = (float4) ((float) OCTREE_DIMENSIONS,
                                       (float) OCTREE_DIMENSIONS,
                                       (float) OCTREE_DIMENSIONS,
                                       1.f);
    entry = ((entry * volDimF) - textureOffset) / volToPOTF;
    exit = ((exit * volDimF) - textureOffset) / volToPOTF;


    // initialize ray
    RayInfo ray;
    ray.param = 0.f;             ///< ray parameter
    ray.color = (float4)(0.f);   ///< resulting color
    ray.firsthit = 1.0f;          ///< first hit point
    ray.channelIntensities[0] = 0.f;
    ray.channelIntensities[1] = 0.f;
    ray.channelIntensities[2] = 0.f;
    ray.channelIntensities[3] = 0.f;

    ray.level = -1.0f;

// copy channel shifts to array (channel shift is already in voxel coordinates and thus only has to be converted to POT)
#ifdef APPLY_CHANNEL_SHIFT
    float3 channelShifts[OCTREE_NUMCHANNELS_DEF];
    channelShifts[0] = channelShift0.xyz / volToPOTF.xyz;
    #if OCTREE_NUMCHANNELS_DEF > 1
    channelShifts[1] = channelShift1.xyz / volToPOTF.xyz;
    #endif
    #if OCTREE_NUMCHANNELS_DEF > 2
    channelShifts[2] = channelShift2.xyz / volToPOTF.xyz;
    #endif
    #if OCTREE_NUMCHANNELS_DEF > 3
    channelShifts[3] = channelShift3.xyz / volToPOTF.xyz;
    #endif
#endif

#ifdef DISPLAY_MODE_REFINEMENT
    // fetch intermediate ray result from refinement buffer
    fetchRayFromBuffer(rayBuffer, fragPos, viewportSize, &ray);

    if (firstRefinementFrame) {
        ray.param = 0.f;
        #ifdef COMPOSITING_MODE_DVR
        ray.color = (float4)(0.f);
        ray.firsthit = 1.f;
        #endif
    }

    if (ray.param >= 1.8f) // 1.8 \approx sqrt(3), which is the diagonal of the "texture space cube"
        return;
#endif

    // traverse ray
    if (entry.w > 0.f && (entry.x != exit.x || entry.y != exit.y || entry.z != exit.z)) {
    #ifdef APPLY_CHANNEL_SHIFT
        traverseRayCS(entry.xyz, exit.xyz, viewportSize, nodeBuffer, brickBuffer, brickFlagBuffer,
                    nodeLevelOfDetail, realWorldMapping,
                    transFunc, transFuncDomains, samplingStepSize,
                    channelShifts,
                                &ray);
    #else
    traverseRay(entry.xyz, exit.xyz, viewportSize, nodeBuffer, brickBuffer, brickFlagBuffer,
                    nodeLevelOfDetail, realWorldMapping,
                    transFunc, transFuncDomains, samplingStepSize,
                                &ray);
    #endif
    }

// store ray result in ray buffer
writeRayToBuffer(rayBuffer, fragPos, viewportSize, &ray);

// write fragment
#ifdef DISPLAY_MODE_REFINEMENT
    // DVR:     output finished rays only
    // MIP/MOP: always output ray result
    bool outputRay = true;
    #ifdef COMPOSITING_MODE_DVR
    if (ray.param < 1.f)
        outputRay = false;
    #endif

    if (outputRay) {
        write_imagef(output, fragPos, ray.color);
        write_imagef(outputDepth, fragPos, ray.firsthit);
    }

#else // full frame mode => always output ray
    write_imagef(output, fragPos, ray.color);
    write_imagef(outputDepth, fragPos, ray.firsthit);
#endif

}

//-------------------------------------------------------------------------------------------
// helper kernels

__kernel void updateBrickBuffer(__global const ushort* brickUpdateBuffer,
                                __global const uint* brickUpdateAddressBuffer,
                                const uint numUpdateBricks,
                                __global ushort* brickBuffer
                                )
{
    // each brick is divided among 16 work units
    const uint BRICK_NUM_BLOCKS = 16;

    const uint brickID =          get_global_id(0) / BRICK_NUM_BLOCKS;
    const uint brickBlockOffset = get_global_id(0) % BRICK_NUM_BLOCKS;
    if (brickID >= numUpdateBricks)
        return;

    const uint BRICK_NUMVOXELS = OCTREE_BRICKDIM*OCTREE_BRICKDIM*OCTREE_BRICKDIM*OCTREE_NUMCHANNELS;
    const uint srcOffset =   brickID*BRICK_NUMVOXELS + brickBlockOffset*BRICK_NUMVOXELS/BRICK_NUM_BLOCKS;
    const uint destOffset =  brickUpdateAddressBuffer[brickID]*BRICK_NUMVOXELS +
                             brickBlockOffset*BRICK_NUMVOXELS/BRICK_NUM_BLOCKS;
    for (int i=0; i<BRICK_NUMVOXELS/BRICK_NUM_BLOCKS; i++)
        brickBuffer[destOffset+i] = brickUpdateBuffer[srcOffset+i];
}

