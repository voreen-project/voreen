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

// variables for storing compositing results
vec4 result1 = vec4(0.0);
vec4 result2 = vec4(0.0);
vec4 result3 = vec4(0.0);

// declare entry and exit parameters
uniform sampler2D entryPoints_;            // ray entry points
uniform sampler2D entryPointsDepth_;       // ray entry points depth
uniform TextureParameters entryParameters_;

uniform sampler2D exitPoints_;             // ray exit points
uniform sampler2D exitPointsDepth_;        // ray exit points depth
uniform TextureParameters exitParameters_;

uniform float samplingStepSize_;

uniform TransFuncParameters transferFunc1_;
#if defined (TF_SAMPLER_TYPE)
uniform TF_SAMPLER_TYPE transferFuncTex1_;
#endif

#if NUM_CHANNELS >= 2
uniform TransFuncParameters transferFunc2_;
#if defined (TF_SAMPLER_TYPE)
uniform TF_SAMPLER_TYPE transferFuncTex2_;
#endif
#endif

#if NUM_CHANNELS >= 3
uniform TransFuncParameters transferFunc3_;
#if defined (TF_SAMPLER_TYPE)
uniform TF_SAMPLER_TYPE transferFuncTex3_;
#endif
#endif

#if NUM_CHANNELS >= 4
uniform TransFuncParameters transferFunc4_;
#if defined (TF_SAMPLER_TYPE)
uniform TF_SAMPLER_TYPE transferFuncTex4_;
#endif
#endif

uniform float isoValue_;

uniform float gammaValue1_;                      // used for MIDA raycasting
uniform float gammaValue2_;                     // used for MIDA raycasting
uniform float gammaValue3_;                     // used for MIDA raycasting

//pass radius values
uniform float radiusMin_;
uniform float radiusMax_;
uniform float shiftFactor_;

// declare volume
uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;    // volume data with parameters

uniform vec3 channelShift1_;
uniform vec3 channelShift2_;
uniform vec3 channelShift3_;
uniform vec3 channelShift4_;


#define PI 3.1415926538

#define MAXRADIUS 6371/shiftFactor_
#define MINRADIUS 3485/shiftFactor_

/***
 * converts cartesian coordinates to spherical coordinates.
 * returns a vec3, where x is the radius, y is theta and z is phi.
 ***/
vec3 convCCtoSC(vec3 cc) {
    cc = (cc - vec3(0.5))*2*MAXRADIUS;
    float r = sqrt(pow(cc.x, 2) + pow(cc.y, 2) + pow(cc.z, 2));
    float theta = acos(cc.z/r);
    float phi = atan(cc.y, cc.x);
    return vec3(phi, r, theta);
}

vec4 getVoxelChannelShift(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos){
#ifndef ENABLE_CHANNEL_SHIFT
    return getVoxel(volume, volumeStruct, samplePos);
#else
    vec4 result;
    result.x = getVoxel(volume, volumeStruct, samplePos+channelShift1_).x;
#if NUM_CHANNELS >= 2
    result.y = getVoxel(volume, volumeStruct, samplePos+channelShift2_).y;
#endif
#if NUM_CHANNELS >= 3
    result.z = getVoxel(volume, volumeStruct, samplePos+channelShift3_).z;
#endif
#if NUM_CHANNELS >= 4
    result.w = getVoxel(volume, volumeStruct, samplePos+channelShift4_).w;
#endif
    return result;
#endif
}

/***
 * Performs the ray traversal
 * returns the final fragment color.
 ***/
void rayTraversal(in vec3 first, in vec3 last, float entryDepth, float exitDepth) {
    // calculate the required ray parameters
    float t     = 0.0;
    float tIncr = 0.0;
    float tEnd  = 1.0;
    float tDepth = -1.0;
    vec4 lastIntensity = vec4(0.0f); //used for pre-integrated transfer-functions
    //TODO: transform to real-world value?

    vec4 maxIntensity1 = vec4(0.0);   // used for MIDA raycasting
    vec4 maxIntensity2 = vec4(0.0);  // used for MIDA raycasting
    vec4 maxIntensity3 = vec4(0.0);  // used for MIDA raycasting

    vec3 rayDirection;
    raySetup(first, last, samplingStepSize_, rayDirection, tIncr, tEnd);

#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT1
    vec4 result1c1=vec4(0);
#if NUM_CHANNELS >= 2
    vec4 result1c2=vec4(0);
#endif
#if NUM_CHANNELS >= 3
    vec4 result1c3=vec4(0);
#endif
#if NUM_CHANNELS >= 4
    vec4 result1c4=vec4(0);
#endif


#define RESULT1_1 result1c1
#define RESULT1_2 result1c2
#define RESULT1_3 result1c3
#define RESULT1_4 result1c4
#else
#define RESULT1_1 result1
#define RESULT1_2 result1
#define RESULT1_3 result1
#define RESULT1_4 result1
#endif

#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT2
    vec4 result2c1=vec4(0);
#if NUM_CHANNELS >= 2
    vec4 result2c2=vec4(0);
#endif
#if NUM_CHANNELS >= 3
    vec4 result2c3=vec4(0);
#endif
#if NUM_CHANNELS >= 4
    vec4 result2c4=vec4(0);
#endif
#define RESULT2_1 result2c1
#define RESULT2_2 result2c2
#define RESULT2_3 result2c3
#define RESULT2_4 result2c4
#else
#define RESULT2_1 result2
#define RESULT2_2 result2
#define RESULT2_3 result2
#define RESULT2_4 result2
#endif

#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT3
    vec4 result3c1=vec4(0);
#if NUM_CHANNELS >= 2
    vec4 result3c2=vec4(0);
#endif
#if NUM_CHANNELS >= 3
    vec4 result3c3=vec4(0);
#endif
#if NUM_CHANNELS >= 4
    vec4 result3c4=vec4(0);
#endif
#define RESULT3_1 result3c1
#define RESULT3_2 result3c2
#define RESULT3_3 result3c3
#define RESULT3_4 result3c4getVoxel
#else
#define RESULT3_1 result3
#define RESULT3_2 result3
#define RESULT3_3 result3
#define RESULT3_4 result3
#endif


    bool finished = false;
    WHILE(!finished) {
        //setup radius min max
        float rMin = radiusMin_/shiftFactor_;
        float rMax = radiusMax_/shiftFactor_;

        //sample position in spherical Coordinates
        vec3 samplePos = first + t * rayDirection;
        vec3 sphericalSamplePos = convCCtoSC(samplePos);

        //normalise
        float rNN = sphericalSamplePos.y;
        sphericalSamplePos.x = ((sphericalSamplePos.x + PI) / (2 * PI));    //phi (2PI, 360)
        sphericalSamplePos.y = (sphericalSamplePos.y - MINRADIUS)/(MAXRADIUS - MINRADIUS); //radius
        sphericalSamplePos.z = sphericalSamplePos.z / PI;                   //theta (PI, 180)

        vec4 voxel = getVoxelChannelShift(volume_, volumeStruct_, sphericalSamplePos);
        vec3 gradient;

        //check if radius is in bounds. Added two values, because of inaccuracy
        bool radiusCheck = true;
        if(rNN > rMax || rNN < rMin)
            radiusCheck = false;
#ifdef CLASSIFICATION_REQUIRES_GRADIENT
        // calculate gradients
        gradient = CALC_GRADIENT(volume_, volumeStruct_, samplePos);
#endif
        bool hasOpacity = false;
        // apply classification
        vec4 color1 = RC_APPLY_CLASSIFICATION(transferFunc1_, transferFuncTex1_, vec4(gradient, voxel.r), lastIntensity.r);
        hasOpacity = (hasOpacity || (color1.a > 0)) && radiusCheck;
#if NUM_CHANNELS >= 2
        vec4 color2 = RC_APPLY_CLASSIFICATION(transferFunc2_, transferFuncTex2_, vec4(gradient, voxel.g), lastIntensity.g);
        hasOpacity = (hasOpacity || (color1.a > 0)) && radiusCheck;
#endif
#if NUM_CHANNELS >= 3
        vec4 color3 = RC_APPLY_CLASSIFICATION(transferFunc3_, transferFuncTex3_, vec4(gradient, voxel.b), lastIntensity.b);
        hasOpacity = (hasOpacity || (color1.a > 0)) && radiusCheck;
#endif
#if NUM_CHANNELS >= 4
        vec4 color4 = RC_APPLY_CLASSIFICATION(transferFunc4_, transferFuncTex4_, vec4(gradient, voxel.a), lastIntensity.a);
        hasOpacity = (hasOpacity || (color1.a > 0)) && radiusCheck;
#endif
        // if opacity greater zero, apply compositing
        if (hasOpacity) {
            if(t == 0.0) // the gradient fix is only need for shading purposes and will mess with 2D-TFs
                gradient = fixClipBorderGradient(samplePos, rayDirection, entryPoints_, entryParameters_);
#ifndef CLASSIFICATION_REQUIRES_GRADIENT
            else
                gradient = CALC_GRADIENT(volume_, volumeStruct_, samplePos);
#endif

            // apply shading
            color1.rgb = APPLY_SHADING(gradient, texToPhysical(samplePos, volumeStruct_), volumeStruct_.lightPositionPhysical_, volumeStruct_.cameraPositionPhysical_, color1.rgb, color1.rgb, vec3(1.0,1.0,1.0));
            RESULT1_1 = RC_APPLY_COMPOSITING1(RESULT1_1, voxel.r, color1, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity1.r, gammaValue1_, isoValue_);
            RESULT2_1 = RC_APPLY_COMPOSITING2(RESULT2_1, voxel.r, color1, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity2.r, gammaValue2_, isoValue_);
            RESULT3_1 = RC_APPLY_COMPOSITING3(RESULT3_1, voxel.r, color1, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity3.r, gammaValue3_, isoValue_);
#if NUM_CHANNELS >= 2
            color2.rgb = APPLY_SHADING(gradient, texToPhysical(samplePos, volumeStruct_), volumeStruct_.lightPositionPhysical_, volumeStruct_.cameraPositionPhysical_, color2.rgb, color2.rgb, vec3(1.0,1.0,1.0));
            RESULT1_2 = RC_APPLY_COMPOSITING1(RESULT1_2, voxel.g, color2, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity1.g, gammaValue1_, isoValue_);
            RESULT2_2 = RC_APPLY_COMPOSITING2(RESULT2_2, voxel.g, color2, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity2.g, gammaValue2_, isoValue_);
            RESULT3_2 = RC_APPLY_COMPOSITING3(RESULT3_2, voxel.g, color2, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity3.g, gammaValue3_, isoValue_);
#endif
#if NUM_CHANNELS >= 3
            color3.rgb = APPLY_SHADING(gradient, texToPhysical(samplePos, volumeStruct_), volumeStruct_.lightPositionPhysical_, volumeStruct_.cameraPositionPhysical_, color3.rgb, color3.rgb, vec3(1.0,1.0,1.0));
            RESULT1_3 = RC_APPLY_COMPOSITING1(RESULT1_3, voxel.b, color3, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity1.b, gammaValue1_, isoValue_);
            RESULT2_3 = RC_APPLY_COMPOSITING2(RESULT2_3, voxel.b, color3, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity2.b, gammaValue2_, isoValue_);
            RESULT3_3 = RC_APPLY_COMPOSITING3(RESULT3_3, voxel.b, color3, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity3.b, gammaValue3_, isoValue_);
#endif
#if NUM_CHANNELS >= 4
            color4.rgb = APPLY_SHADING(gradient, texToPhysical(samplePos, volumeStruct_), volumeStruct_.lightPositionPhysical_, volumeStruct_.cameraPositionPhysical_, color4.rgb, color4.rgb, vec3(1.0,1.0,1.0));
            RESULT1_4 = RC_APPLY_COMPOSITING1(RESULT1_4, voxel.a, color4, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity1.a, gammaValue1_, isoValue_);
            RESULT2_4 = RC_APPLY_COMPOSITING2(RESULT2_4, voxel.a, color4, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity2.a, gammaValue2_, isoValue_);
            RESULT3_4 = RC_APPLY_COMPOSITING3(RESULT3_4, voxel.a, color4, samplePos, gradient, t, samplingStepSize_, tDepth, maxIntensity3.a, gammaValue3_, isoValue_);
#endif
        }
        lastIntensity = voxel;
#ifdef USE_EARLY_TERMINATION
/**

 ===================================
               FIX EARLY RAY TERMINATION!!!!!!!!!!!!!!!!!!!!!!!!
==========================
*/
        finished = earlyRayTermination(RESULT1_1.a, RESULT2_1.a, RESULT3_1.a, EARLY_RAY_TERMINATION_OPACITY);
#endif
        t += tIncr;
        finished = finished || (t > tEnd);
    } END_WHILE


    result1 = RESULT1_1;
    result2 = RESULT2_1;
    result3 = RESULT3_1;
#if NUM_CHANNELS >= 2
#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT1
    result1 += result1c2;
#endif
#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT2
    result2 += result2c2;
#endif
#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT3
    result3 += result3c2;
#endif
#endif

#if NUM_CHANNELS >= 3
#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT1
    result1 += result1c3;
#endif
#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT2
    result2 += result2c3;
#endif
#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT3
    result3 += result3c3;
#endif
#endif

#if NUM_CHANNELS >= 4
#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT1
    result1 += result1c4;
#endif
#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT2
    result2 += result2c4;
#endif
#ifdef SEPERATE_CHANNEL_COMPOSITION_CHANNEL_OUTPUT3
    result3 += result3c4;
#endif
#endif

    gl_FragDepth = getDepthValue(tDepth, tEnd, entryDepth, exitDepth);
}

void main() {
    // fetch entry/exit points
    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    vec3 frontPos = textureLookup2Dnormalized(entryPoints_, entryParameters_, p).rgb;
    vec3 backPos = textureLookup2Dnormalized(exitPoints_, exitParameters_, p).rgb;
    float entryDepth = textureLookup2Dnormalized(entryPointsDepth_, entryParameters_, p).x;
    float exitDepth = textureLookup2Dnormalized(exitPointsDepth_, exitParameters_, p).x;

    // determine whether the ray has to be casted
    if (frontPos == backPos)
        // background needs no raycasting
        discard;
    else {
        // fragCoords are lying inside the bounding box
        rayTraversal(frontPos, backPos, entryDepth, exitDepth);
    }

    #ifdef OP0
        FragData0 = result1;
    #endif
    #ifdef OP1
        FragData1 = result2;
    #endif
    #ifdef OP2
        FragData2 = result3;
    #endif
}
