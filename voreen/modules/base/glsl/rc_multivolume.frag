/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

/**
 * IMPORTANT: We are assuming each entry/exit texture has the same parameters.
 */

// variables for storing compositing results
vec4 result = vec4(0.0);
vec4 result1 = vec4(0.0);
vec4 result2 = vec4(0.0);

uniform float samplingStepSize_;

// declare entry and exit parameters
uniform sampler2D entryGlobalPoints_;            // global ray entry points
uniform sampler2D entryGlobalPointsDepth_;       // global ray entry points depth
uniform sampler2D exitGlobalPoints_;             // global ray exit points
uniform sampler2D exitGlobalPointsDepth_;        // global ray exit points depth

uniform TextureParameters entryExitParameters_;

// declare volumes
#ifdef VOLUME_1_ACTIVE
uniform VolumeParameters volumeStruct1_;
uniform sampler3D volume1_;    // volume dataset 1
uniform TransFuncParameters transferFunc1_;
uniform TF_SAMPLER_TYPE_1 transferFuncTex1_;
uniform sampler2D entryPoints1_;                 // entry points 1
uniform sampler2D exitPoints1_;                  // exit points 1
#endif

#ifdef VOLUME_2_ACTIVE
uniform VolumeParameters volumeStruct2_;
uniform sampler3D volume2_;    // volume dataset 2
uniform TransFuncParameters transferFunc2_;
uniform TF_SAMPLER_TYPE_2 transferFuncTex2_;
uniform sampler2D entryPoints2_;                 // entry points 2
uniform sampler2D exitPoints2_;                  // exit points 2
#endif

#ifdef VOLUME_3_ACTIVE
uniform VolumeParameters volumeStruct3_;
uniform sampler3D volume3_;    // volume dataset 3
uniform TransFuncParameters transferFunc3_;
uniform TF_SAMPLER_TYPE_3 transferFuncTex3_;
uniform sampler2D entryPoints3_;                 // entry points 3
uniform sampler2D exitPoints3_;                  // exit points 3
#endif

#ifdef VOLUME_4_ACTIVE
uniform VolumeParameters volumeStruct4_;
uniform sampler3D volume4_;    // volume dataset 4
uniform TransFuncParameters transferFunc4_;
uniform TF_SAMPLER_TYPE_4 transferFuncTex4_;
uniform sampler2D entryPoints4_;                 // entry points 4
uniform sampler2D exitPoints4_;                  // exit points 4
#endif

/***
 * Performs the ray traversal
 * returns the final fragment color.
 ***/
void rayTraversal(in vec3 first, in vec3 last, float entryDepth, float exitDepth
#ifdef VOLUME_1_ACTIVE
                 , in vec3 vol1Entry, in vec3 vol1Exit
#endif
#ifdef VOLUME_2_ACTIVE
                 , in vec3 vol2Entry, in vec3 vol2Exit
#endif
#ifdef VOLUME_3_ACTIVE
                 , in vec3 vol3Entry, in vec3 vol3Exit
#endif
#ifdef VOLUME_4_ACTIVE
                 , in vec3 vol4Entry, in vec3 vol4Exit
#endif
) {
    // calculate the required ray parameters
    float t     = 0.0;
    float tIncr = 0.0;
    float tEnd  = 1.0;
    vec3 rayDirection;
    raySetup(first, last, samplingStepSize_, rayDirection, tIncr, tEnd);

#ifdef VOLUME_1_ACTIVE
    vec3 vol1first = worldToTex(first, volumeStruct1_);
    vec3 vol1dir = worldToTex(last, volumeStruct1_) - vol1first;
    float vol1Begin = length(vol1Entry - first) / tEnd;
    float vol1End = length(vol1Exit - first) / tEnd;
    float lastIntensity1 = 0.0f; //used for pre-integrated transfer-functions
#endif
#ifdef VOLUME_2_ACTIVE
    vec3 vol2first = worldToTex(first, volumeStruct2_);
    vec3 vol2dir = worldToTex(last, volumeStruct2_) - vol2first;
    float vol2Begin = length(vol2Entry - first) / tEnd;
    float vol2End = length(vol2Exit - first) / tEnd;
    float lastIntensity2 = 0.0f; //used for pre-integrated transfer-functions
#endif
#ifdef VOLUME_3_ACTIVE
    vec3 vol3first = worldToTex(first, volumeStruct3_);
    vec3 vol3dir = worldToTex(last, volumeStruct3_) - vol3first;
    float vol3Begin = length(vol3Entry - first) / tEnd;
    float vol3End = length(vol3Exit - first) / tEnd;
    float lastIntensity3 = 0.0f; //used for pre-integrated transfer-functions
#endif
#ifdef VOLUME_4_ACTIVE
    vec3 vol4first = worldToTex(first, volumeStruct4_);
    vec3 vol4dir = worldToTex(last, volumeStruct4_) - vol4first;
    float vol4Begin = length(vol4Entry - first) / tEnd;
    float vol4End = length(vol4Exit - first) / tEnd;
    float lastIntensity4 = 0.0f; //used for pre-integrated transfer-functions
#endif
    float maxIntensity1 = 0;
    float maxIntensity2 = 0;
    float maxIntensity3 = 0;
    float realT;
    float tDepth = -1.0f;
    bool finished = false;
    WHILE(!finished) {
        realT = t / tEnd;

        vec3 worldSamplePos = first + t * rayDirection;

#ifdef VOLUME_1_ACTIVE
        vec3 samplePos1 = vol1first + realT * vol1dir;

        if(realT >= vol1Begin && realT <= vol1End) {
            vec4 voxel1 = getVoxel(volume1_, volumeStruct1_, samplePos1);

#ifdef CLASSIFICATION_REQUIRES_GRADIENT
        // calculate gradients
        voxel1.xyz = CALC_GRADIENT(volume1_, volumeStruct1_, samplePos1);
#endif

        // apply classification
        vec4 color = RC_APPLY_CLASSIFICATION(transferFunc1_, transferFuncTex1_, voxel1, lastIntensity1);

        // if opacity greater zero, apply compositing
        if (color.a > 0.0) {
#ifndef CLASSIFICATION_REQUIRES_GRADIENT
            voxel1.xyz = CALC_GRADIENT(volume1_, volumeStruct1_, samplePos1);
#endif

            // apply shading
            color.rgb = APPLY_SHADING_1(voxel1.xyz, texToPhysical(samplePos1, volumeStruct1_), volumeStruct1_.lightPositionPhysical_, volumeStruct1_.cameraPositionPhysical_, color.rgb, color.rgb, vec3(1.0,1.0,1.0));

            result = RC_APPLY_COMPOSITING1(result, voxel1.a, mvOpacityCorrection(color), worldSamplePos, voxel1.xyz, t, samplingStepSize_, tDepth, maxIntensity1, 0, 0);
            result1 = RC_APPLY_COMPOSITING2(result1, voxel1.a, mvOpacityCorrection(color), worldSamplePos, voxel1.xyz, t, samplingStepSize_, tDepth, maxIntensity2, 0, 0);
            result2 = RC_APPLY_COMPOSITING3(result2, voxel1.a, mvOpacityCorrection(color), worldSamplePos, voxel1.xyz, t, samplingStepSize_, tDepth, maxIntensity3, 0, 0);
        }
        lastIntensity1 = voxel1.a;
        }
#endif

#ifdef VOLUME_2_ACTIVE
        //second sample:
        vec3 samplePos2 = vol2first + realT * vol2dir;

        if(realT >= vol2Begin && realT <= vol2End) {
            vec4 voxel2 = getVoxel(volume2_, volumeStruct2_, samplePos2);

#ifdef CLASSIFICATION_REQUIRES_GRADIENT
            // calculate gradients
            voxel2.xyz = CALC_GRADIENT(volume2_, volumeStruct2_, samplePos2);
#endif

            // apply classification
            vec4 color2 = RC_APPLY_CLASSIFICATION2(transferFunc2_, transferFuncTex2_, voxel2, lastIntensity2);

            // if opacity greater zero, apply compositing
            if (color2.a > 0.0) {
#ifndef CLASSIFICATION_REQUIRES_GRADIENT
                voxel2.xyz = CALC_GRADIENT(volume2_, volumeStruct2_, samplePos2);
#endif

                // apply shading
                color2.rgb = APPLY_SHADING_2(voxel2.xyz, texToPhysical(samplePos2, volumeStruct2_), volumeStruct2_.lightPositionPhysical_, volumeStruct2_.cameraPositionPhysical_, color2.rgb, color2.rgb, vec3(1.0,1.0,1.0));

                result = RC_APPLY_COMPOSITING1(result, voxel2.a, mvOpacityCorrection(color2), worldSamplePos, voxel2.xyz, t, samplingStepSize_, tDepth, maxIntensity1, 0, 0);
                result1 = RC_APPLY_COMPOSITING2(result1, voxel2.a, mvOpacityCorrection(color2), worldSamplePos, voxel2.xyz, t, samplingStepSize_, tDepth, maxIntensity2, 0, 0);
                result2 = RC_APPLY_COMPOSITING3(result2, voxel2.a, mvOpacityCorrection(color2), worldSamplePos, voxel2.xyz, t, samplingStepSize_, tDepth, maxIntensity3, 0, 0);
            }
            lastIntensity2 = voxel2.a;
        }
#endif

#ifdef VOLUME_3_ACTIVE
        //second sample:
        vec3 samplePos3 = vol3first + realT * vol3dir;

        if(realT >= vol3Begin && realT <= vol3End) {
            vec4 voxel3 = getVoxel(volume3_, volumeStruct3_, samplePos3);

#ifdef CLASSIFICATION_REQUIRES_GRADIENT
            // calculate gradients
            voxel3.xyz = CALC_GRADIENT(volume3_, volumeStruct3_, samplePos3);
#endif

            // apply classification
            vec4 color3 = RC_APPLY_CLASSIFICATION3(transferFunc3_, transferFuncTex3_, voxel3, lastIntensity3);

            // if opacity greater zero, apply compositing
            if (color3.a > 0.0) {
#ifndef CLASSIFICATION_REQUIRES_GRADIENT
                voxel3.xyz = CALC_GRADIENT(volume3_, volumeStruct3_, samplePos3);
#endif

                // apply shading
                color3.rgb = APPLY_SHADING_3(voxel3.xyz, texToPhysical(samplePos3, volumeStruct3_), volumeStruct3_.lightPositionPhysical_, volumeStruct3_.cameraPositionPhysical_, color3.rgb, color3.rgb, vec3(1.0,1.0,1.0));

                result = RC_APPLY_COMPOSITING1(result, voxel3.a, mvOpacityCorrection(color3), worldSamplePos, voxel3.xyz, t, samplingStepSize_, tDepth, maxIntensity1, 0, 0);
                result1 = RC_APPLY_COMPOSITING2(result1, voxel3.a, mvOpacityCorrection(color3), worldSamplePos, voxel3.xyz, t, samplingStepSize_, tDepth, maxIntensity2, 0, 0);
                result2 = RC_APPLY_COMPOSITING3(result2, voxel3.a, mvOpacityCorrection(color3), worldSamplePos, voxel3.xyz, t, samplingStepSize_, tDepth, maxIntensity3, 0, 0);
            }
            lastIntensity3 = voxel3.a;
        }
#endif

#ifdef VOLUME_4_ACTIVE
        //second sample:
        vec3 samplePos4 = vol4first + realT * vol4dir;

        if(realT >= vol4Begin && realT <= vol4End) {
            vec4 voxel4 = getVoxel(volume4_, volumeStruct4_, samplePos4);

#ifdef CLASSIFICATION_REQUIRES_GRADIENT
            // calculate gradients
            voxel4.xyz = CALC_GRADIENT(volume4_, volumeStruct4_, samplePos4);
#endif

            // apply classification
            vec4 color4 = RC_APPLY_CLASSIFICATION4(transferFunc4_, transferFuncTex4_, voxel4, lastIntensity4);

            // if opacity greater zero, apply compositing
            if (color4.a > 0.0) {
#ifndef CLASSIFICATION_REQUIRES_GRADIENT
                voxel4.xyz = CALC_GRADIENT(volume4_, volumeStruct4_, samplePos4);
#endif

                // apply shading
                color4.rgb = APPLY_SHADING_4(voxel4.xyz, texToPhysical(samplePos4, volumeStruct4_), volumeStruct4_.lightPositionPhysical_, volumeStruct4_.cameraPositionPhysical_, color4.rgb, color4.rgb, vec3(1.0,1.0,1.0));

                result = RC_APPLY_COMPOSITING1(result, voxel4.a, mvOpacityCorrection(color4), worldSamplePos, voxel4.xyz, t, samplingStepSize_, tDepth, maxIntensity1, 0, 0);
                result1 = RC_APPLY_COMPOSITING2(result1, voxel4.a, mvOpacityCorrection(color4), worldSamplePos, voxel4.xyz, t, samplingStepSize_, tDepth, maxIntensity2, 0, 0);
                result2 = RC_APPLY_COMPOSITING3(result2, voxel4.a, mvOpacityCorrection(color4), worldSamplePos, voxel4.xyz, t, samplingStepSize_, tDepth, maxIntensity3, 0, 0);
            }
            lastIntensity4 = voxel4.a;
        }
#endif
#ifdef USE_EARLY_TERMINATION
        finished = earlyRayTermination(result.a, result1.a, result2.a, EARLY_RAY_TERMINATION_OPACITY);
#endif
        t += tIncr;
        finished = finished || (t > tEnd);
    } END_WHILE
    gl_FragDepth = getDepthValue(tDepth, tEnd, entryDepth, exitDepth);
}

/***
 * The main method.
 ***/
void main() {
    // fetch entry/exit points
    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    vec3 frontPosGlobal = textureLookup2Dnormalized(entryGlobalPoints_, entryExitParameters_, p).rgb;
    vec3 backPosGlobal = textureLookup2Dnormalized(exitGlobalPoints_, entryExitParameters_, p).rgb;

    // determine whether the ray has to be casted, background needs no raycasting
    if (frontPosGlobal == backPosGlobal)
        discard;

    // fragCoords are lying inside the bounding box
#ifdef VOLUME_1_ACTIVE
    vec3 volume1Entry = textureLookup2Dnormalized(entryPoints1_, entryExitParameters_, p).rgb;
    vec3 volume1Exit  = textureLookup2Dnormalized(exitPoints1_ , entryExitParameters_, p).rgb;
#endif
#ifdef VOLUME_2_ACTIVE
    vec3 volume2Entry = textureLookup2Dnormalized(entryPoints2_, entryExitParameters_, p).rgb;
    vec3 volume2Exit  = textureLookup2Dnormalized(exitPoints2_ , entryExitParameters_, p).rgb;
#endif
#ifdef VOLUME_3_ACTIVE
    vec3 volume3Entry = textureLookup2Dnormalized(entryPoints3_, entryExitParameters_, p).rgb;
    vec3 volume3Exit  = textureLookup2Dnormalized(exitPoints3_ , entryExitParameters_, p).rgb;
#endif
#ifdef VOLUME_4_ACTIVE
    vec3 volume4Entry = textureLookup2Dnormalized(entryPoints4_, entryExitParameters_, p).rgb;
    vec3 volume4Exit  = textureLookup2Dnormalized(exitPoints4_ , entryExitParameters_, p).rgb;
#endif
    float entryDepth = textureLookup2Dnormalized(entryGlobalPointsDepth_, entryExitParameters_, p).x;
    float exitDepth = textureLookup2Dnormalized(exitGlobalPointsDepth_, entryExitParameters_, p).x;

    //main raytraversel
    rayTraversal(frontPosGlobal, backPosGlobal, entryDepth, exitDepth
#ifdef VOLUME_1_ACTIVE
                 , volume1Entry, volume1Exit
#endif
#ifdef VOLUME_2_ACTIVE
                 , volume2Entry, volume2Exit
#endif
#ifdef VOLUME_3_ACTIVE
                 , volume3Entry, volume3Exit
#endif
#ifdef VOLUME_4_ACTIVE
                 , volume4Entry, volume4Exit
#endif
                );

    #ifdef OP0
        FragData0 = result;
    #endif
    #ifdef OP1
        FragData1 = result1;
    #endif
    #ifdef OP2
        FragData2 = result2;
    #endif
}
