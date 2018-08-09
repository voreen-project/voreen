/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "modules/mod_sampler2d.frag"
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"

#include "modules/mod_depth.frag"

uniform float samplingStepSize_;
const float SAMPLING_BASE_INTERVAL_RCP = 200.0;

uniform sampler2D entryPoints_;          // ray entry points
uniform sampler2D entryPointsDepth_;     // ray entry points depth
uniform TextureParameters entryParameters_;
uniform sampler2D exitPoints_;           // ray exit points
uniform sampler2D exitPointsDepth_;      // ray exit points depth
uniform TextureParameters exitParameters_;

uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;               // texture lookup parameters for volume_

uniform TransFuncParameters transferFunc_;
uniform TF_SAMPLER_TYPE transferFuncTex_;

/***
 * Performs direct volume rendering and
 * returns the final fragment color.
 ***/
vec4 directRendering(in vec3 first, in vec3 last, vec2 p) {

    vec4 result = vec4(0.0);
    float depthT = -1.0;
    bool finished = false;

    // calculate ray parameters
    float stepIncr = samplingStepSize_;
    float t = 0.0;
    vec3 direction = last.rgb - first.rgb;
    float tend = length(direction);
    direction = normalize(direction);

    // raycasting loop
    for (int loop=0; !finished && loop<255*255; loop++) {

        vec3 sampleVal = first.rgb + t * direction;
        vec4 voxel = getVoxel(volume_, volumeStruct_, sampleVal);

        // no shading is applied
        vec4 color = applyTF(transferFunc_, transferFuncTex_, voxel);

        // perform compositing
        if (color.a > 0.0) {
            // accomodate for variable sampling rates (base interval defined by mod_compositing.frag)
            color.a = 1.0 - pow(1.0 - color.a, samplingStepSize_ * SAMPLING_BASE_INTERVAL_RCP);
            result.rgb = result.rgb + (1.0 - result.a) * color.a * color.rgb;
            result.a = result.a + (1.0 -result.a) * color.a;
        }

        // save first hit ray parameter for depth value calculation
        if (depthT < 0.0 && result.a > 0.0)
            depthT = t;

        // early ray termination
        if (result.a >= 1.0) {
            result.a = 1.0;
            finished = true;
        }

        t += stepIncr;
        finished = finished || (t > tend);
    }

    // calculate depth value from ray parameter
    gl_FragDepth = 1.0;
    if (depthT >= 0.0) {
        float depthEntry = textureLookup2Dnormalized(entryPointsDepth_, entryParameters_, p).x;
        float depthExit = textureLookup2Dnormalized(exitPointsDepth_, exitParameters_, p).x;
        gl_FragDepth = calculateDepthValue(depthT/tend, depthEntry, depthExit);
    }

    return result;
}

/***
 * The main method.
 ***/
void main() {
    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    vec3 frontPos = textureLookup2Dnormalized(entryPoints_, entryParameters_, p).rgb;
    vec3 backPos = textureLookup2Dnormalized(exitPoints_, exitParameters_, p).rgb;

    //determine whether the ray has to be casted
    if (frontPos == backPos) {
        //background need no raycasting
        discard;
    } else {
        //fragCoords are lying inside the boundingbox
        FragData0 = directRendering(frontPos, backPos, p);
    }
}
