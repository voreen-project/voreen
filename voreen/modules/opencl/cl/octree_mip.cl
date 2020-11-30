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

#include "octree_config.cl"

/**
 * skip node, if brick is missing (inhomogeneous node)
 * un-request traversed brick
 * mark available brick (or homogeneous node) as traversed
 */
#define getBrickInRefinement \
    if (!currentNodeHasBrick && !isHomogeneous(currentNode.value_)) {\
        sampleNode = false;\
        if (hasNodeBeenUsedInMIPAll(brickFlagBuffer[currentNode.offset_])) \
                setBrickRequested(brickFlagBuffer + currentNode.offset_, false);\
        else\
            rayFinished = false;\
    }\
    else {\
        sampleNode = true;\
        setNodeUsedInMIPAll(brickFlagBuffer + currentNode.offset_, true);\
    }

/**
 * skip node, if brick is missing (inhomogeneous node)
 * un-request traversed brick
 * mark available brick (or homogeneous node) as traversed
 * sampleNode = false; is missing in this case
 */
#define getBrickInRefinementCS \
    if (!currentNodeHasBrick && !isHomogeneous(currentNode.value_)) {\
        sampleNode = false;\
        if (hasNodeBeenUsedInMIPAll(brickFlagBuffer[currentNode.offset_])) \
                setBrickRequested(brickFlagBuffer + currentNode.offset_, false);\
        else\
            rayFinished = false;\
    }\
    else {\
        sampleNode = true;\
        setNodeUsedInMIPAll(brickFlagBuffer + currentNode.offset_, true);\
    }

/**
 * Update pending intensity and depth value
 * Apply pending results if new pending intensity is already higher than the previous.
 */
#define applyTFandCombineColors \
    for (int ch=0; ch<OCTREE_NUMCHANNELS; ch++) {\
        if (ray->pending.intensity[ch] < channelIntensities[ch]) {\
            ray->pending.intensity[ch] = channelIntensities[ch];\
            ray->pending.firsthit = min(ray->pending.firsthit, ray->param/tEnd);\
            if(ray->current.intensity[ch] < ray->pending.intensity[ch]) {\
                ray->current.intensity[ch] = ray->pending.intensity[ch];\
                ray->current.firsthit = ray->pending.firsthit;\
            }\
        }\
    }\

/**
 * Retrieve and composit results based on current intensities.
 * If the ray is finished, use the previously pending, now current results.
 */
#define postRaycastingLoop\
    if(rayFinished) {\
        for (int ch=0; ch<OCTREE_NUMCHANNELS; ch++) {\
            ray->current.intensity[ch] = ray->pending.intensity[ch];\
            ray->current.firsthit = ray->pending.firsthit;\
        }\
        ray->param = 2.f;\
    } else {\
        ray->param = 0.f;\
    }\
    float4 maxIntensityColor[OCTREE_NUMCHANNELS_DEF];\
    applyTransFuncs(ray->current.intensity, transFunc, transFuncDomains, realWorldMapping, maxIntensityColor);\
    ray->color = (float4)(0.f);\
    for (int ch=0; ch<OCTREE_NUMCHANNELS; ch++) {\
        float4 channelColor = maxIntensityColor[ch];\
        if (channelColor.w > 0.f) {\
            ray->color += channelColor;\
        }\
    }\
    ray->color = min(ray->color, (float4)(1.f));\

