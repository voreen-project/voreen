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

#ifdef GLSL_VERSION_400
#extension GL_ARB_gpu_shader5 : require
#extension GL_EXT_shader_image_load_store : require
#extension GL_EXT_bindable_uniform : require
#endif

#include "modules/mod_sampler2d.frag"
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"

#include "modules/mod_raysetup.frag"
#include "modules/mod_depth.frag"
#include "modules/mod_compositing.frag"
#include "modules/mod_gradients.frag"
#include "modules/mod_shading.frag"

uniform sampler2D entryPoints_;                      // ray entry points
uniform sampler2D entryPointsDepth_;              // ray entry points depth
uniform TextureParameters entryParameters_;

uniform sampler2D exitPoints_;                        // ray exit points
uniform sampler2D exitPointsDepth_;                // ray exit points depth
uniform TextureParameters exitParameters_;

uniform float samplingStepSize_;

uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;        // volume dataset

uniform TransFuncParameters transferFunc_;
uniform TF_SAMPLER_TYPE transferFuncTex_;

#if defined(ZONE_RENDERING) && defined(GLSL_VERSION_400)
uniform sampler2D zoneMap_;
uniform TextureParameters zoneMapParameters_;
uniform sampler2D renderTarget_;
uniform TextureParameters renderTargetParameters_;
uniform int currentRenderedZone_;
uniform int edgeFactor_;
uniform ivec2 posDirZoneVec_;
uniform ivec2 negDirZoneVec_;
uniform ivec2 projLightPos_;
#endif

#if defined(USE_EARLY_RAY_TERMINATION_SWEEP) && defined(GLSL_VERSION_400)
uniform sampler2D ertPoints_;                      // ray ert points
uniform TextureParameters ertParameters_;      // ray ert  parameters
#endif

// shading uniforms
uniform sampler1D albedoFunc_;
uniform sampler1D surfacenessFunc_;
uniform float ambientIntensity_;

#ifdef GLSL_VERSION_400
//Light map storage
uniform layout(size4x32) image2D lightMaps_[2];
uniform layout(size4x32) image2D lightPosMaps_[2];
uniform TextureParameters lightMapParameters_;
uniform ivec2 lightMapReadXWriteY;

uniform int lineWidth_;
uniform int firstLine_;
uniform float shadowBrightness_;
uniform ivec2 propagationDirection_;     //Propagation direction in screen coordinates, (1, 0), (-1,0), (0, 1) or (0, -1);
uniform mat4 objToLightMat_;
uniform mat4 lightViewMat_;
uniform vec2 lightMapSize_;

//Return true/false if valid/invalid coordinates
bool checkValidCoordinates(in ivec2 coords){
    return (coords.x >= 0 && coords.y >= 0 && coords.x < screenDim_.x && coords.y < screenDim_.y);
}

//Bilinear interpolation
vec4 imageLoadBilinear(in layout(size4x32) image2D map, in vec2 fcoord)
{
    ivec4 coord;
    coord.xy = ivec2(floor(fcoord) + 0.5);
    coord.zw = coord.xy + 1;

    vec4 tex11 = imageLoad(map, coord.xy);
    vec4 tex21 = imageLoad(map, coord.zy);
    vec4 tex12 = imageLoad(map, coord.xw);
    vec4 tex22 = imageLoad(map, coord.zw);

    vec2 factor = fcoord - vec2(coord.xy);

    return mix(mix(tex11, tex21, factor.x), mix(tex12, tex22, factor.x), factor.y);
}

//Loads the visibility from the current image using bilinear interpolation
vec4 loadVisibilityBI(in vec2 coords){
    return imageLoadBilinear(lightMaps_[lightMapReadXWriteY.x], coords);
}

//Loads the visibility from the current image using bilinear interpolation
vec3 loadPositionBI(in vec2 coords){
    return imageLoadBilinear(lightPosMaps_[lightMapReadXWriteY.x], coords).rgb;
}

//Loads the visibility from the current image
vec4 loadVisibility(in ivec2 coords){
    return imageLoad(lightMaps_[lightMapReadXWriteY.x], coords);
}

//Loads the latest saved voxel position from the position image
vec3 loadPosition(in ivec2 coords){
    return imageLoad(lightPosMaps_[lightMapReadXWriteY.x], coords).rgb;
}

//Stores the visibility from the light source (currently just the color from that slice)
void storeVisibility(in vec4 val, in ivec2 coords){
    imageStore(lightMaps_[lightMapReadXWriteY.y], coords, val);
}

//Stores the visibility position
void storePosition(in vec3 pos, in ivec2 coords){
    imageStore(lightPosMaps_[lightMapReadXWriteY.y], coords, vec4(pos, 1.0));
}

//Convert from normalized coordinates [0, 1] to real coordinates [0, tex_dim]
void convertFromNormToRealCoords(inout vec2 pos, in vec2 dim){
    pos *= dim;
}

//Get light texture coord using projective texturing.
vec2 getLightTexCoord(in vec3 pos, in mat4 transform){
    vec3 worldCoord = pos - 0.5;

    vec4 light_pos = transform * vec4(worldCoord, 1.0);

    light_pos /= light_pos.w;
    light_pos.xy *= 0.5;
    light_pos.xy += 0.5;

    convertFromNormToRealCoords(light_pos.xy, lightMapSize_);

    return light_pos.xy;
}
#endif

// RGB to YUV color conversion
vec3 rgb2yuv(vec3 colorRGB) {
    vec3 colorYUV = vec3(0.0);
    // r=y, g=u, b=v
    colorYUV.r = 0.299 * colorRGB.r + 0.587 * colorRGB.g + 0.114 * colorRGB.b;
    colorYUV.g = 0.436 * (colorRGB.b - colorYUV.r) / (1.0 - 0.114);
    colorYUV.b = 0.615 * (colorRGB.r - colorYUV.r) / (1.0 - 0.299);
    return colorYUV;
}

// YUV to RGB color conversion
vec3 yuv2rgb(vec3 colorYUV) {
    vec3 colorRGB = vec3(0.0);
    // r=y, g=u, b=v
    colorRGB.r = colorYUV.r + 1.13983 * colorYUV.b;
    colorRGB.g = colorYUV.r - 0.39465 * colorYUV.g - 0.58060 * colorYUV.b;
    colorRGB.b = colorYUV.r + 2.03211 * colorYUV.g;
    return colorRGB;
}

vec3 calcGradientFDUnscaled(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos, float v) {
    vec3 gradient = vec3(0.0);
    vec3 offset = volumeStruct.datasetDimensionsRCP_;
    float v0 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(offset.x, 0.0, 0.0)).a;
    float v1 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, offset.y, 0)).a;
    float v2 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, offset.z)).a;
    gradient = vec3(v - v0, v - v1, v - v2);
    return gradient;
}

//Get voxel color
vec4 getColor(in vec3 pos, in float intensity){
    vec4 color = applyTF(transferFunc_, transferFuncTex_, intensity);
    vec3 colorYUV = rgb2yuv(color.rgb);

    // get light source distance for attenuation and normalize light vector
    vec3 vpos = texToPhysical(pos, volumeStruct_);
    vec3 L = volumeStruct_.lightPositionPhysical_-vpos;
    float d = length(L);
    L /= d;
    colorYUV.r *= getAttenuation(d);

    // apply ambient
    colorYUV.r += ambientIntensity_;

    if (applyTF(surfacenessFunc_, intensity).a == 1.0) {
        // apply local illumination
        vec3 grad = calcGradientFDUnscaled(volume_, volumeStruct_, pos, intensity);
        float gradLength = length(grad);
        vec3 N = grad / gradLength;
        vec3 V = volumeStruct_.cameraPositionPhysical_-vpos;

        // apply diffuse
        float NdotL = max(dot(N, L), 0.0);
        colorYUV.r *= NdotL;

        // apply specular
        //vec3 H = normalize(V + L);
        //float NdotH = pow(max(dot(N, H), 0.0), 220.0);
        //color += NdotH;
    }

    return vec4(yuv2rgb(colorYUV), color.a);
}

#ifdef GLSL_VERSION_400
//Add visibillity contribution to vis and color buffer
void addContribution(inout vec4 buf, inout vec4 color, in vec4 vis, in vec4 vis_org){

    //For nice color effect... why before?
    color.rgb *= vis.rgb;

#if defined(BILINEAR_INTERPOLATION) && !defined(PIECEWISE_INTEGRATION)
    vis = vis_org;
#endif

#if defined(COLOR_CONTRIBUTION)
    buf.rgb = (1.0 - color.a) * vis.rgb + color.a * color.rgb * shadowBrightness_;
    buf.a = (1.0 - color.a) * vis.a + color.a;
#else
    buf = vec4((1.0 - color.a) * vis.a);
#endif
}

//Piecewise integration towards previous slice, return integrated visibillity
vec4 castLightRaySeg(in vec3 pos, in vec3 stored_pos, in vec4 vis, in int steps){
    vec3 direction = pos - stored_pos;
    float tend = length(direction);
    direction = normalize(direction);
    vec4 vis_buf = vec4(vis);

    float step_inc = tend/float(steps-1);
    float t = 0.0;
    tend -= step_inc;

    vec3 samplePos;
    vec4 color;

    while(t < tend){
        samplePos = stored_pos + t * direction;
        float intensity = textureLookup3DUnnormalized(volume_, volumeStruct_, samplePos).a;

        //get color, no shading applied
        color = getColor(samplePos, intensity);

        addContribution(vis_buf, color, vis_buf, vis);

        t += step_inc;
    }

    return vis_buf;
}
#endif

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
    float tend;
    float t = 0.0;
    vec3 direction = last.rgb - first.rgb;

#if defined(USE_EARLY_RAY_TERMINATION_SWEEP) && defined(GLSL_VERSION_400)
    tend = length(textureLookup2Dnormalized(ertPoints_, ertParameters_, p).rgb - first.rgb);
#else
    tend = length(direction);
#endif
    direction = normalize(direction);

    vec3 samplePos = first.rgb;
    vec4 voxel; //Voxel value
    vec4 color; //Color of voxel from TF

#ifdef GLSL_VERSION_400
    //mat4 objTransform = gl_ProjectionMatrix * lightViewMat_ * volumeStruct_.physicalToWorldMatrix_;
    mat4 objTransform = objToLightMat_;

    //Visibility buffer
    vec4 vis_buf = vec4(1.0);

    //Check if valid line to save visibility on (last one in one draw call).
    ivec2 abs_prop = abs(propagationDirection_);
    int currentLine = int(dot(gl_FragCoord.xy, vec2(abs_prop)));
    float lineInOrder = mod(float(abs(firstLine_ - currentLine)), float(lineWidth_));
    lineInOrder++;
    bool isSaveLine = (lineInOrder == lineWidth_);
    if(lineInOrder == 1) lineInOrder = lineWidth_;

    ivec2 prevCoord = ivec2(floor(getLightTexCoord(samplePos, objTransform)));

    //Current visibility value
    vec4 vis = loadVisibility(prevCoord);
    vec4 visPI = vis;

    vec2 lightCoord; //Current coordinate in light texture
#endif

#if defined(PIECEWISE_INTEGRATION) && defined(GLSL_VERSION_400)
    //Current visibility position
    vec3 prevPos = first.rgb;
    vec3 vis_pos = loadPosition(prevCoord);
#endif

    // 2 nested loops allow for more than 255 iterations,
    // but is slower than while (t < tend)
    for (int loop0=0; !finished && loop0<255; loop0++) {
        for (int loop1=0; !finished && loop1<255; loop1++) {

            samplePos = first.rgb + t * direction;
            float intensity = textureLookup3DUnnormalized(volume_, volumeStruct_, samplePos).a;

#ifdef GLSL_VERSION_400
            lightCoord = vec2(getLightTexCoord(samplePos, objTransform));
#endif

            color = getColor(samplePos, intensity);

            // back to front compositing
#if defined(BILINEAR_INTERPOLATION) && defined(PIECEWISE_INTEGRATION) && defined(GLSL_VERSION_400)
            vis_pos = loadPositionBI(lightCoord);
            visPI = castLightRaySeg(samplePos, vis_pos, vis, int(lineInOrder));
#elif defined(BILINEAR_INTERPOLATION) && defined(GLSL_VERSION_400)
            visPI = loadVisibilityBI(lightCoord);
#elif !defined(PIECEWISE_INTEGRATION) && defined(GLSL_VERSION_400)
            visPI = vis;
#endif

#if defined(GLSL_VERSION_400)
            // simulate scattering based on albedo
            float albedo = applyTF(albedoFunc_, intensity).a;
            vis = albedo*visPI + (1.0-albedo)*vis;
#endif

 #ifdef GLSL_VERSION_400
            // add contribution
            addContribution(vis_buf, color, visPI, vis);
#endif

            //opacity correction
            color.a = 1.0 - pow(1.0 - color.a, samplingStepSize_ * SAMPLING_BASE_INTERVAL_RCP);

            result.rgb = result.rgb + (1.0 - result.a) * color.a *color.rgb;
            result.a = result.a + (1.0 -result.a) * color.a;

            // save first hit ray parameter for depth value calculation
            if (depthT < 0.0 && result.a > 0.0){
                depthT = t;
            }

#ifndef GLSL_VERSION_400
            // early ray termination
           if (result.a >= 1.0) {
                result.a = 1.0;
                finished = true;
            }
#endif

#ifdef GLSL_VERSION_400
            if(prevCoord != ivec2(floor(lightCoord))){
                if(isSaveLine){
#endif
#if defined(VIS_ACCUMULATED_RAY) && defined(GLSL_VERSION_400)
                    storeVisibility(vec4(1.0-result.a), prevCoord);
#elif defined(GLSL_VERSION_400)
                    storeVisibility(vis_buf, prevCoord);
#endif
#if defined(PIECEWISE_INTEGRATION) && defined(GLSL_VERSION_400)
                    storePosition(prevPos+((samplePos-prevPos)/2.0), prevCoord);
                    prevPos = samplePos;
#endif
#ifdef GLSL_VERSION_400
                }
                prevCoord = ivec2(lightCoord);
                vis = loadVisibility(prevCoord);
#endif
#if !defined(BILINEAR_INTERPOLATION) && defined(PIECEWISE_INTEGRATION) && defined(GLSL_VERSION_400)
                vis_pos = loadPosition(prevCoord);
                visPI = castLightRaySeg(samplePos, vis_pos, vis, int(lineInOrder));
#endif
#ifdef GLSL_VERSION_400
            }
#endif
            t += stepIncr;
            finished = finished || (t > tend);
        }
    }

    // calculate depth value from ray parameter
    gl_FragDepth = 1.0;
    if (depthT >= 0.0)
        gl_FragDepth = calculateDepthValue(depthT/tend, textureLookup2Dnormalized(entryPointsDepth_, entryParameters_, p).z,
                                                        textureLookup2Dnormalized(exitPointsDepth_, exitParameters_, p).z);

#if defined(ZONE_RENDERING) && defined(GLSL_VERSION_400)
    vec2 zone = textureLookup2Dnormalized(zoneMap_, zoneMapParameters_, p).rg;
    if((currentRenderedZone_ == 2 && zone != vec2(0.0, 1.0)) ||
       (currentRenderedZone_ == 3 && zone != vec2(0.0, 0.0)) ||
       (currentRenderedZone_ == 4 && zone != vec2(1.0, 1.0))){
        vec4 currentColor = textureLookup2Dnormalized(renderTarget_, renderTargetParameters_, p);
        if(currentColor.a > 0.0){ //Perform edge blending
            vec2 closestToArea = gl_FragCoord.xy;
            vec2 closestToNeighbor = closestToArea;
            vec2 propInv = abs(vec2(propagationDirection_.yx));
            vec2 traverseVec = (dot(closestToArea, propInv) < dot(vec2(projLightPos_), propInv)) ?  vec2(negDirZoneVec_) : vec2(posDirZoneVec_);
            closestToArea = closestToNeighbor = traverseVec * abs(dot(closestToArea, propInv.yx) - dot(vec2(projLightPos_), propInv.yx));
            closestToNeighbor = edgeFactor_ * closestToArea;
            closestToArea += vec2(projLightPos_);
            closestToNeighbor += vec2(projLightPos_);

            vec2 blendFactor = smoothstep(closestToArea, closestToNeighbor, gl_FragCoord.xy);
            result.rgb = mix(result.rgb, currentColor.rgb, vec3(dot(blendFactor, propInv)));
        }
    }
#endif

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
