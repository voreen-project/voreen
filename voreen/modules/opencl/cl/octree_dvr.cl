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

#include "octree_config.cl"

/**
 * break ray traversal on missing brick
 * traversal flag not used in DVR
 */
#define getBrickInRefinement \
            if (!currentNodeHasBrick && !isHomogeneous(currentNode.value_)) {\
                return;\
            }

/**
 * break ray traversal on missing brick
 * traversal flag not used in DVR
 */
#define getBrickInRefinementCS \
            if (!currentNodeHasBrick[ch] && !isHomogeneous(currentNode[ch].value_)) {\
                return;\
            }

/**
 * apply channel transfer functions
 * perform compositing
 *   save first-hit point
 * early ray termination
 */
#define applyTFandCombineColors \
        float4 channelColors[OCTREE_NUMCHANNELS_DEF];\
        applyTransFuncs(channelIntensities, transFunc, transFuncDomains, realWorldMapping, channelColors);\
\
        float4 sampleColor = (float4)(0.f);\
        for (int i=0; i<OCTREE_NUMCHANNELS; i++) {\
            if (channelColors[i].w > 0.f) {\
                sampleColor.xyz += channelColors[i].xyz;\
                sampleColor.w = max(channelColors[i].w, sampleColor.w);\
            }\
        }\
        sampleColor.xyz = min(sampleColor.xyz, (float3)(1.f));\
\
        if (sampleColor.w > 0.f) {\
            sampleColor.w = 1.f - pow(1.f - sampleColor.w, previousStepSize * SAMPLING_BASE_INTERVAL_RCP);\
\
            ray->color.xyz = ray->color.xyz + (1.f - ray->color.w) * sampleColor.w * sampleColor.xyz;\
            ray->color.w   = ray->color.w +   (1.f - ray->color.w) * sampleColor.w;\
\
            ray->firsthit = min(ray->param, ray->firsthit);\
        }\
\
        if (ray->color.w >= 0.95f) {\
            ray->color.w = 1.f;\
            ray->param = tEnd + samplingStepSize;\
        }

/**
 * Sets the ray param to 2.f = finished
 */
#define postRaycastingLoop ray->param = 2.f;
