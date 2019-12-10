/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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
 * skipNode = false; is missing in this case
 */
#define getBrickInRefinementCS \
                if (!currentNodeHasBrick[ch] && !isHomogeneous(currentNode[ch].value_)) {\
                    sampleNode[ch] = false;\
                    if (hasNodeBeenUsedInMIP(brickFlagBuffer[currentNode[ch].offset_],ch)) \
                        setBrickRequested(brickFlagBuffer + currentNode[ch].offset_, false);\
                    else\
                        rayFinished = false;\
                }\
                else {\
                    sampleNode[ch] = true;\
                    setNodeUsedInMIP(brickFlagBuffer + currentNode[ch].offset_, true, ch);\
                }

/**
 * retrieve channel color values
 * update channel color and depth values
 * close bracket opend in main loop
 */
#define applyTFandCombineColors \
            float4 sampleIntensityColor[OCTREE_NUMCHANNELS_DEF];\
            applyTransFuncs(channelIntensities, transFunc, transFuncDomains, realWorldMapping, sampleIntensityColor);\
\
            for (int ch=0; ch<OCTREE_NUMCHANNELS; ch++) {\
                if (sampleIntensityColor[ch].w > channelColors[ch].w) {\
                    channelColors[ch] = sampleIntensityColor[ch];\
                    ray->channelIntensities[ch] = channelIntensities[ch];\
                    ray->firsthit = min(ray->firsthit, ray->param);\
                }\
            }\

/**
 * Set final ray color and mark ray as finished
 */
#define postRaycastingLoop \
            ray->color = (float4)(0.f);\
            for (int i=0; i<OCTREE_NUMCHANNELS; i++) {\
                float4 channelColor = channelColors[i];\
                if (channelColor.w > 0.f) {\
                    ray->color += channelColor;\
                }\
            }\
            ray->color = min(ray->color, (float4)(1.f));\
\
            ray->param = rayFinished ? 2.f : 0.f;
