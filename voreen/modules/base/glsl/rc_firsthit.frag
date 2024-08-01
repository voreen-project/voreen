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

#include "modules/mod_raysetup.frag"

#include "modules/mod_depth.frag"
#include "modules/mod_compositing.frag"
#include "modules/mod_gradients.frag"
#include "modules/mod_shading.frag"

// declare entry and exit parameters
uniform sampler2D entryPoints_;            // ray entry points
uniform sampler2D entryPointsDepth_;       // ray entry points depth

uniform float samplingStepSize_;

// declare volume
uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;    // texture lookup parameters for volume_

uniform TransFuncParameters transferFunc_;
uniform TF_SAMPLER_TYPE transferFuncTex_;

/***
 * ...
 ***/
vec4 getFirstHitColor(in vec3 firstHitPos) {
    vec4 result = vec4(0.0);
    vec4 voxel = getVoxel(volume_, volumeStruct_, firstHitPos);

    // calculate gradients
    //if(t == 0.0)
        //voxel.xyz = fixClipBorderGradient(firstHitPos, rayDirection, entryPoints_, entryParameters_);
    //else
        //voxel.xyz = CALC_GRADIENT(vec3(0.0), entryPoints_, 0.0);

    // apply classification
    result = RC_APPLY_CLASSIFICATION(transferFunc_, transferFuncTex_, voxel, voxel.a);

    // apply shading
    result.rgb = APPLY_SHADING(voxel.xyz, texToPhysical(firstHitPos, volumeStruct_), volumeStruct_.lightPositionPhysical_, volumeStruct_.cameraPositionPhysical_, result.rgb, result.rgb, result.rgb);

    return result;
}

/***
 * The main method.
 ***/
void main() {

    vec3 firstHitPos = textureLookup2D(entryPoints_, gl_FragCoord.xy).rgb;

    if (firstHitPos == vec3(0.0))
        // background needs no rendering
        discard;
    else {
        FragData0 = getFirstHitColor(firstHitPos);
        // FIXME: somehow the depth values seem to be not correct (tr)
        //gl_FragDepth = calculateDepthValue(0.5, entryPointsDepth_, entryPointsDepth_);
        //gl_FragDepth = textureLookup2D(entryPointsDepth_, gl_FragCoord.xy).z;

        // obtain coordinates of proxy front face in world coordinates
        vec4 front = vec4((gl_FragCoord.x*screenDimRCP_.x*2.0)-1.0,
                          (gl_FragCoord.y*screenDimRCP_.y*2.0)-1.0,
                          (textureLookup2D(entryPointsDepth_, gl_FragCoord.xy).x*2.0)-1.0,
                          1.0);
        front = gl_ModelViewProjectionMatrixInverse * front;

        front = gl_ModelViewProjectionMatrix * front;
        gl_FragDepth = ((front.z/front.w)+1.0)/2.0;

    }
}
