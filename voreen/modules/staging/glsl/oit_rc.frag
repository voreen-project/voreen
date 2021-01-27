/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#extension GL_ARB_gpu_shader5 : require

#include "modules/mod_sampler2d.frag"
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"

#include "modules/mod_raysetup.frag"

#include "modules/mod_depth.frag"
#include "modules/mod_compositing.frag"
#include "modules/mod_gradients.frag"
#include "modules/mod_shading.frag"

#include "mod_oit.frag"

uniform mat4 viewMatrixInverse_;
uniform mat4 projectionMatrixInverse_;
uniform mat4 viewMatrix_;
uniform mat4 projectionMatrix_;
uniform vec3 cameraPosition_;
uniform float near_;
uniform float far_;

uniform float samplingStepSize_;
uniform float clippingGradientDepth_; // fix gradients this far into the volume
vec3 worldPos_;
vec3 worldDir_;
float curDepth_;
float firstDepth_;

struct RayState {
    bool active_;
    vec3 pos_; // current sampling position
    vec3 dir_; // length(dir_) = samplingStepSize_
    float lastIntensity_; //used for pre-integrated transfer-functions
    vec3 cameraPositionPhysical_;
    vec3 lightPositionPhysical_;
    vec3 pgNormal_; // normal of the proxy-geometry (used for correct shading on clipping surfaces)
    float entryDepth_; // used to determine if clipping surface-based shading has to be used
    float skinExitDepth_; // used to determine if clipping surface-based shading has to be used
};

// declare volumes
#ifdef VOLUME_1_ACTIVE
uniform VolumeParameters volumeStruct1_;
uniform sampler3D volume1_;    // volume dataset 1
uniform TransFuncParameters transferFunc1_;
uniform TF_SAMPLER_TYPE_1 transferFuncTex1_;
RayState ray1;
#endif

#ifdef VOLUME_2_ACTIVE
uniform VolumeParameters volumeStruct2_;
uniform sampler3D volume2_;    // volume dataset 2
uniform TransFuncParameters transferFunc2_;
uniform TF_SAMPLER_TYPE_2 transferFuncTex2_;
RayState ray2;
#endif

#ifdef VOLUME_3_ACTIVE
uniform VolumeParameters volumeStruct3_;
uniform sampler3D volume3_;    // volume dataset 3
uniform TransFuncParameters transferFunc3_;
uniform TF_SAMPLER_TYPE_3 transferFuncTex3_;
RayState ray3;
#endif

#ifdef VOLUME_4_ACTIVE
uniform VolumeParameters volumeStruct4_;
uniform sampler3D volume4_;    // volume dataset 4
uniform TransFuncParameters transferFunc4_;
uniform TF_SAMPLER_TYPE_4 transferFuncTex4_;
RayState ray4;
#endif

/***
 * Performs the ray traversal
 * returns the final fragment color.
 ***/
vec4 rayTraversal(in float exitDepth, in vec4 inColor) {
    // variables for storing compositing results
    vec4 result = inColor;
    vec4 result1 = vec4(0.0);
    vec4 result2 = vec4(0.0);

    float tDepth = -1.0f;
    bool finished = false;
    WHILE(!finished) {
        //realT = curDepth_ / tEnd;

        //--------------VOLUME 1---------------------
#ifdef VOLUME_1_ACTIVE
        if(ray1.active_) {
            vec4 voxel1 = getVoxel(volume1_, volumeStruct1_, ray1.pos_);

#ifdef CLASSIFICATION_REQUIRES_GRADIENT
        // calculate gradients
        voxel1.xyz = CALC_GRADIENT(volume1_, volumeStruct1_, ray1.pos_);
#endif
            // apply classification
            vec4 color = RC_APPLY_CLASSIFICATION(transferFunc1_, transferFuncTex1_, voxel1, ray1.lastIntensity_);

            // if opacity greater zero, apply compositing
            if (color.a > 0.0) {
                // calculate gradients
                if(curDepth_ <= ray1.skinExitDepth_)
                    voxel1.xyz = ray1.pgNormal_;
#ifndef CLASSIFICATION_REQUIRES_GRADIENT
                else
                    voxel1.xyz = CALC_GRADIENT(volume1_, volumeStruct1_, ray1.pos_);
#endif

                // apply shading
                color.rgb = APPLY_SHADING_1(voxel1.xyz, texToPhysical(ray1.pos_, volumeStruct1_), ray1.lightPositionPhysical_, ray1.cameraPositionPhysical_, color.rgb, color.rgb, vec3(1.0,1.0,1.0));

                result = RC_APPLY_COMPOSITING1(result, voxel.x, color, worldPos_, voxel1.xyz, curDepth_, samplingStepSize_, tDepth, 0, 0, 0);
                //result1 = RC_APPLY_COMPOSITING2(result1, color, worldPos_, voxel1.xyz, curDepth_, samplingStepSize_, tDepth);
                //result2 = RC_APPLY_COMPOSITING3(result2, color, worldPos_, voxel1.xyz, curDepth_, samplingStepSize_, tDepth);
            }
            ray1.lastIntensity_ = voxel1.a;
            ray1.pos_ += ray1.dir_;
        }
#endif

#ifdef VOLUME_2_ACTIVE
        if(ray2.active_) {
            vec4 voxel2 = getVoxel(volume2_, volumeStruct2_, ray2.pos_);

#ifdef CLASSIFICATION_REQUIRES_GRADIENT
        // calculate gradients
        voxel2.xyz = CALC_GRADIENT(volume2_, volumeStruct2_, ray2.pos_);
#endif

            // apply classification
            vec4 color2 = RC_APPLY_CLASSIFICATION2(transferFunc2_, transferFuncTex2_, voxel2, ray2.lastIntensity_);

            // if opacity greater zero, apply compositing
            if (color2.a > 0.0) {
                // calculate gradients
                if(curDepth_ <= ray2.skinExitDepth_)
                    voxel2.xyz = ray2.pgNormal_;
#ifndef CLASSIFICATION_REQUIRES_GRADIENT
                else
                    voxel2.xyz = CALC_GRADIENT(volume2_, volumeStruct2_, ray2.pos_);
#endif

                // apply shading
                color2.rgb = APPLY_SHADING_2(voxel2.xyz, texToPhysical(ray2.pos_, volumeStruct2_), ray2.lightPositionPhysical_, ray2.cameraPositionPhysical_, color2.rgb, color2.rgb, vec3(1.0,1.0,1.0));

                result = RC_APPLY_COMPOSITING1(result, voxel2.a, color2, worldPos_, voxel2.xyz, curDepth_, samplingStepSize_, tDepth, 0, 0, 0);
                //result1 = RC_APPLY_COMPOSITING2(result1, color2, worldPos_, voxel2.xyz, curDepth_, samplingStepSize_, tDepth);
                //result2 = RC_APPLY_COMPOSITING3(result2, color2, worldPos_, voxel2.xyz, curDepth_, samplingStepSize_, tDepth);
            }
            ray2.lastIntensity_ = voxel2.a;
            ray2.pos_ += ray2.dir_;
        }
#endif

#ifdef VOLUME_3_ACTIVE
        if(ray3.active_) {
            vec4 voxel3 = getVoxel(volume3_, volumeStruct3_, ray3.pos_);

#ifdef CLASSIFICATION_REQUIRES_GRADIENT
        // calculate gradients
        voxel3.xyz = CALC_GRADIENT(volume3_, volumeStruct3_, ray3.pos_);
#endif

            // apply classification
            vec4 color3 = RC_APPLY_CLASSIFICATION3(transferFunc3_, transferFuncTex3_, voxel3, ray3.lastIntensity_);

            // if opacity greater zero, apply compositing
            if (color3.a > 0.0) {
                // calculate gradients
                if(curDepth_ <= ray3.skinExitDepth_)
                    voxel3.xyz = ray3.pgNormal_;
#ifndef CLASSIFICATION_REQUIRES_GRADIENT
                else
                    voxel3.xyz = CALC_GRADIENT(volume3_, volumeStruct3_, ray3.pos_);
#endif

                // apply shading
                color3.rgb = APPLY_SHADING_3(voxel3.xyz, texToPhysical(ray3.pos_, volumeStruct3_), ray3.lightPositionPhysical_, ray3.cameraPositionPhysical_, color3.rgb, color3.rgb, vec3(1.0,1.0,1.0));

                result = RC_APPLY_COMPOSITING1(result, voxel3.a, color3, worldPos_, voxel3.xyz, curDepth_, samplingStepSize_, tDepth, 0, 0, 0);
                //result1 = RC_APPLY_COMPOSITING2(result1, color3, worldPos_, voxel3.xyz, curDepth_, samplingStepSize_, tDepth);
                //result2 = RC_APPLY_COMPOSITING3(result2, color3, worldPos_, voxel3.xyz, curDepth_, samplingStepSize_, tDepth);
            }
            ray3.lastIntensity_ = voxel3.a;
            ray3.pos_ += ray3.dir_;
        }
#endif

#ifdef VOLUME_4_ACTIVE
        if(ray4.active_) {
            vec4 voxel4 = getVoxel(volume4_, volumeStruct4_, ray4.pos_);

#ifdef CLASSIFICATION_REQUIRES_GRADIENT
            // calculate gradients
            voxel4.xyz = CALC_GRADIENT(volume4_, volumeStruct4_, ray4.pos_);
#endif

            // apply classification
            vec4 color4 = RC_APPLY_CLASSIFICATION4(transferFunc4_, transferFuncTex4_, voxel4, ray4.lastIntensity_);

            // if opacity greater zero, apply compositing
            if (color4.a > 0.0) {
                // calculate gradients
                if(curDepth_ <= ray4.skinExitDepth_)
                    voxel4.xyz = ray4.pgNormal_;
#ifndef CLASSIFICATION_REQUIRES_GRADIENT
                else
                    voxel4.xyz = CALC_GRADIENT(volume4_, volumeStruct4_, ray4.pos_);
#endif

                // apply shading
                color4.rgb = APPLY_SHADING_4(voxel4.xyz, texToPhysical(ray4.pos_, volumeStruct4_), ray4.lightPositionPhysical_, ray4.cameraPositionPhysical_, color4.rgb, color4.rgb, vec3(1.0,1.0,1.0));

                result = RC_APPLY_COMPOSITING1(result, voxel4.a, color4, worldPos_, voxel4.xyz, curDepth_, samplingStepSize_, tDepth, 0, 0, 0);
                //result1 = RC_APPLY_COMPOSITING2(result1, color4, worldPos_, voxel4.xyz, curDepth_, samplingStepSize_, tDepth);
                //result2 = RC_APPLY_COMPOSITING3(result2, color4, worldPos_, voxel4.xyz, curDepth_, samplingStepSize_, tDepth);
            }
            ray4.lastIntensity_ = voxel4.a;
            ray4.pos_ += ray4.dir_;
        }
#endif

        if(firstDepth_ < 0.0 && result.a > 0.0) {
            vec2 tmp = (projectionMatrix_ * viewMatrix_ * vec4(worldPos_, 1.0)).zw;
            firstDepth_ = tmp.x / tmp.y;
        }
        worldPos_ += worldDir_;

        finished = earlyRayTermination(result.a, EARLY_RAY_TERMINATION_OPACITY);
        curDepth_ += samplingStepSize_;
        finished = finished || (curDepth_ > exitDepth);
    } END_WHILE
    return result;
}

vec4 blend(vec4 prev, vec4 cur) {
    cur.rgb = cur.rgb * cur.w;

    return (prev + cur * (1.0f - prev.a));
}

//vec4 blend(vec4 prev, vec4 cur) {
    ////cur.rgb = cur.rgb * cur.w;

    ////return (prev + cur * (1.0f - prev.a));

    //vec4 result = prev;

    //// apply opacity correction to accomodate for variable sampling intervals
    ////cur.a = 1.0 - pow(1.0 - cur.a, 4.0* SAMPLING_BASE_INTERVAL_RCP);

    //result.rgb = result.rgb + (1.0 - result.a) * cur.a * cur.rgb;
    //result.a = result.a + (1.0 -result.a) * cur.a;
    //return result;
//}

void main() {
    vec2 coord = gl_FragCoord.xy - vec2(0.5);

    uint head = imageLoad(headOffset, ivec2(coord)).x;

    if(head == 0U)
        discard;
    else {
        int numVolumesActive = 0;

#ifdef VOLUME_1_ACTIVE
        ray1.active_ = false;
        ray1.lastIntensity_ = 0.0;
#endif

#ifdef VOLUME_2_ACTIVE
        ray2.active_ = false;
        ray2.lastIntensity_ = 0.0;
#endif

#ifdef VOLUME_3_ACTIVE
        ray3.active_ = false;
        ray3.lastIntensity_ = 0.0;
#endif

#ifdef VOLUME_4_ACTIVE
        ray4.active_ = false;
        ray4.lastIntensity_ = 0.0;
#endif
        // Calculate near-plane position in world space:
        vec2 p = gl_FragCoord.xy * screenDimRCP_;
        vec4 wGeom = vec4(2.0*p - 1.0, -1.0, 1.0);

        // reverse perspective division by w (which is equal to the camera-space z)
        float origZG = 2.0*far_*near_ / ((far_ + near_) - wGeom.z*(far_ - near_));
        wGeom *= origZG;

        // unproject and reverse camera-transform
        worldPos_ = (viewMatrixInverse_ * (projectionMatrixInverse_ * wGeom)).xyz;
        worldDir_ = normalize(worldPos_ - cameraPosition_) * samplingStepSize_;

        // composite:
        vec4 color  = vec4(0.0);
        curDepth_ = distance(worldPos_, cameraPosition_);
        firstDepth_ = -1.f;

        while(head != 0U) {
            LinkedListStruct cur = linkedList_[head];
            head = cur.next_;

            if(numVolumesActive > 0)
                color = rayTraversal(cur.depth_, color);
            else {
                // empty segment - update curDepth_ and worldPos_
                float numSteps = ceil((cur.depth_ - curDepth_) / samplingStepSize_);
                curDepth_ += samplingStepSize_ * numSteps;
                worldPos_ += worldDir_ * numSteps;
            }

            if(firstDepth_ < 0.0 && cur.proxyGeometryId_ == 0) {
                vec2 tmp = (projectionMatrix_ * viewMatrix_ * vec4(worldPos_, 1.0)).zw;
                firstDepth_ = tmp.x / tmp.y;
            }

            if(cur.proxyGeometryId_ != 0) {
                ProxyGeometryStruct pg = proxyGeometries_[cur.proxyGeometryId_];

                // proxy-geometry fragment, toggle volumeXActive:
                switch(pg.volumeId_) {
#ifdef VOLUME_1_ACTIVE
                    case 1:
                        ray1.active_ = !ray1.active_;

                        if(ray1.active_) {
                            ++numVolumesActive;

                            ray1.pos_ = (pg.worldToTextureMatrix_ * vec4(worldPos_, 1.0)).xyz;
                            ray1.dir_ = (pg.worldToTextureMatrix_ * vec4(worldDir_, 0.0)).xyz;
                            ray1.cameraPositionPhysical_ = pg.cameraPositionPhysical_;
                            ray1.lightPositionPhysical_ = pg.lightPositionPhysical_;
                            ray1.lastIntensity_ = 0.0;

                            ray1.pgNormal_ = getNormal(cur);
                            ray1.pgNormal_ *= (volumeStruct1_.datasetDimensionsRCP_ * volumeStruct1_.datasetSpacingRCP_); // transform to physical space
                            ray1.pgNormal_ = normalize(ray1.pgNormal_);

                            ray1.entryDepth_ = curDepth_;
                            vec3 dirPhysical = ray1.dir_ * volumeStruct1_.datasetDimensions_ * volumeStruct1_.datasetSpacing_;
                            ray1.skinExitDepth_ = curDepth_ + (clippingGradientDepth_ / abs(dot(normalize(dirPhysical), ray1.pgNormal_)));
                        }
                        else {
                            --numVolumesActive;
                        }
                    break;
#endif
#ifdef VOLUME_2_ACTIVE
                    case 2:
                        ray2.active_ = !ray2.active_;

                        if(ray2.active_) {
                            ++numVolumesActive;

                            ray2.pos_ = (pg.worldToTextureMatrix_ * vec4(worldPos_, 1.0)).xyz;
                            ray2.dir_ = (pg.worldToTextureMatrix_ * vec4(worldDir_, 0.0)).xyz;
                            ray2.cameraPositionPhysical_ = pg.cameraPositionPhysical_;
                            ray2.lightPositionPhysical_ = pg.lightPositionPhysical_;
                            ray2.lastIntensity_ = 0.0;

                            ray2.pgNormal_ = getNormal(cur);
                            ray2.pgNormal_ *= (volumeStruct2_.datasetDimensionsRCP_ * volumeStruct2_.datasetSpacingRCP_); // transform to physical space
                            ray2.pgNormal_ = normalize(ray2.pgNormal_);

                            ray2.entryDepth_ = curDepth_;
                            vec3 dirPhysical = ray2.dir_ * volumeStruct2_.datasetDimensions_ * volumeStruct2_.datasetSpacing_;
                            ray2.skinExitDepth_ = curDepth_ + (clippingGradientDepth_ / abs(dot(normalize(dirPhysical), ray2.pgNormal_)));
                        }
                        else
                            --numVolumesActive;
                    break;
#endif
#ifdef VOLUME_3_ACTIVE
                    case 3:
                        ray3.active_ = !ray3.active_;

                        if(ray3.active_) {
                            ++numVolumesActive;

                            ray3.pos_ = (pg.worldToTextureMatrix_ * vec4(worldPos_, 1.0)).xyz;
                            ray3.dir_ = (pg.worldToTextureMatrix_ * vec4(worldDir_, 0.0)).xyz;
                            ray3.cameraPositionPhysical_ = pg.cameraPositionPhysical_;
                            ray3.lightPositionPhysical_ = pg.lightPositionPhysical_;
                            ray3.lastIntensity_ = 0.0;

                            ray3.pgNormal_ = getNormal(cur);
                            ray3.pgNormal_ *= (volumeStruct3_.datasetDimensionsRCP_ * volumeStruct3_.datasetSpacingRCP_); // transform to physical space
                            ray3.pgNormal_ = normalize(ray3.pgNormal_);

                            ray3.entryDepth_ = curDepth_;
                            vec3 dirPhysical = ray3.dir_ * volumeStruct3_.datasetDimensions_ * volumeStruct3_.datasetSpacing_;
                            ray3.skinExitDepth_ = curDepth_ + (clippingGradientDepth_ / abs(dot(normalize(dirPhysical), ray3.pgNormal_)));
                        }
                        else
                            --numVolumesActive;
                    break;
#endif
#ifdef VOLUME_4_ACTIVE
                    case 4:
                        ray4.active_ = !ray4.active_;

                        if(ray4.active_) {
                            ++numVolumesActive;

                            ray4.pos_ = (pg.worldToTextureMatrix_ * vec4(worldPos_, 1.0)).xyz;
                            ray4.dir_ = (pg.worldToTextureMatrix_ * vec4(worldDir_, 0.0)).xyz;
                            ray4.cameraPositionPhysical_ = pg.cameraPositionPhysical_;
                            ray4.lightPositionPhysical_ = pg.lightPositionPhysical_;
                            ray4.lastIntensity_ = 0.0;

                            ray4.pgNormal_ = getNormal(cur);
                            ray4.pgNormal_ *= (volumeStruct4_.datasetDimensionsRCP_ * volumeStruct4_.datasetSpacingRCP_); // transform to physical space
                            ray4.pgNormal_ = normalize(ray4.pgNormal_);

                            ray4.entryDepth_ = curDepth_;
                            vec3 dirPhysical = ray4.dir_ * volumeStruct4_.datasetDimensions_ * volumeStruct4_.datasetSpacing_;
                            ray4.skinExitDepth_ = curDepth_ + (clippingGradientDepth_ / abs(dot(normalize(dirPhysical), ray4.pgNormal_)));
                        }
                        else
                            --numVolumesActive;
                    break;
#endif
                }
            }
            else
                color = blend(color, getColor(cur)); // geometry fragment

            if(earlyRayTermination(color.a, EARLY_RAY_TERMINATION_OPACITY))
                break;
        }
        FragData0 = color;

        // take first non-proxy-geometry fragment depth as depth value
        if(firstDepth_ >= 0.0)
            gl_FragDepth = firstDepth_;
        else
            gl_FragDepth = 1.0;
    }

    //#ifdef OP0
        //FragData0 = result;
    //#endif
    //#ifdef OP1
        //FragData1 = result1;
    //#endif
    //#ifdef OP2
        //FragData2 = result2;
    //#endif
}

