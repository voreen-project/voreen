/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "modules/mod_sampler2d.frag"
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"

#include "modules/mod_depth.frag"

uniform sampler2D entryPoints_;                    // ray entry points
uniform sampler2D entryPointsDepth_;               // ray entry points depth
uniform TextureParameters entryParameters_;

uniform sampler2D exitPoints_;                        // ray exit points
uniform sampler2D exitPointsDepth_;                // ray exit points depth
uniform TextureParameters exitParameters_;

uniform float samplingStepSize_;

uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;            // texture lookup parameters for volume_

uniform TransFuncParameters transferFunc_;
uniform TF_SAMPLER_TYPE transferFuncTex_;

/***
 * Calculates texture for early ray termination
 * result.rgb = the position where the early ray termination was triggered
 * result.a = the depth of the last hit point seen from the front
 ***/
vec4 ertTexture(in vec3 first, in vec3 last, vec2 p) {

    vec4 result = vec4(0.0);
    float depthT = -1.0;
    bool finished = false;

    // calculate ray parameters
    float stepIncr = samplingStepSize_;
    float tend;
    float t = 0.0;
    vec3 direction = last.rgb - first.rgb;

    tend = length(direction);
    direction = normalize(direction);

    vec3 sample;
    vec4 voxel, color;

    //Perform back-to-front ray casting to determine furthest away structures
    while(!finished){
            sample = last.rgb - t * direction;
            voxel = getVoxel(volume_, volumeStruct_, sample);

            // no shading is applied
            color = applyTF(transferFunc_, transferFuncTex_, voxel);

            // encountered first hit point from the rear
            if (color.a > 0.0) {
                finished = true;
            }

            t += stepIncr;
            finished = finished || (t > tend);
    }

    //Recalculate triggered distance from the back to distance seen from the front
    float last_hit_distance = t - stepIncr;
    last_hit_distance = max(tend - last_hit_distance, 0.0);

    //Recalculate tend based on our determined last hit point, i.e. don't do front-to-back if we already traversed the entire ray.
    tend = length(sample - first.rgb);

    t = 0.0;
    finished = false;
    sample = first.rgb;

    //Perform front-to-back ray casting to determine when opacity has reached the ERT threshold
    while(!finished){
            sample = first.rgb + t * direction;
            voxel = getVoxel(volume_, volumeStruct_, sample);

            // no shading is applied
            color = applyTF(transferFunc_, transferFuncTex_, voxel);

            // perform compositing
            if (color.a > 0.0) {
                // accomodate for variable sampling rates (base interval defined by mod_compositing.frag)
                result.a = result.a + (1.0 -result.a) * color.a;
            }

            // early ray termination
            if (result.a >= 0.90) {
                result.a = 1.0;
                finished = true;
            }

            t += stepIncr;
            finished = finished || (t > tend);
    }

    result.rgb = sample;
    result.a = last_hit_distance;

    // calculate depth value from ray parameter
    gl_FragDepth = 1.0;
    if (p.x == 0.0)
        gl_FragDepth = calculateDepthValue(depthT/tend, textureLookup2Dnormalized(entryPointsDepth_, entryParameters_, p).z,
                                                        textureLookup2Dnormalized(exitPointsDepth_, exitParameters_, p).z);

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
        FragData0 = ertTexture(frontPos, backPos, p);
    }
}
