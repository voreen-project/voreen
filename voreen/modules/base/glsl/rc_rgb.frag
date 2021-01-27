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

#include "modules/mod_sampler2d.frag"
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"

#include "modules/mod_raysetup.frag"

#include "modules/mod_depth.frag"
#include "modules/mod_compositing.frag"
#include "modules/mod_gradients.frag"
#include "modules/mod_shading.frag"

uniform sampler2D entryPoints_;                    // ray entry points
uniform sampler2D entryPointsDepth_;               // ray entry points depth
uniform TextureParameters entryParameters_;

uniform sampler2D exitPoints_;                     // ray exit points
uniform sampler2D exitPointsDepth_;                // ray exit points depth
uniform TextureParameters exitParameters_;

uniform float samplingStepSize_;

uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;            // volume dataset

uniform TransFuncParameters transferFunc_;
uniform TF_SAMPLER_TYPE transferFuncTex_;

uniform bool applyColorModulation_;

vec3 rgb2hsv(vec3 colorRGB) {
    float minComponent = min(colorRGB.r, min(colorRGB.g, colorRGB.b));
    float maxComponent = max(colorRGB.r, max(colorRGB.g, colorRGB.b));
    float delta = maxComponent - minComponent;

    vec3 result;
    result.b = maxComponent;
    if (maxComponent != 0.0) result.g = delta / maxComponent;
    else result.g = 0.0;
    if (result.g == 0.0) result.r = 0.0; // no hue
    else {
        if (colorRGB.r == maxComponent) result.r = (colorRGB.g - colorRGB.b) / delta;
        else if (colorRGB.g == maxComponent) result.r = 2 + (colorRGB.b - colorRGB.r) / delta;
        else if (colorRGB.b == maxComponent) result.r = 4 + (colorRGB.r - colorRGB.g) / delta;
        result.r *= 60.0;
        if (result.r < 0.0) result.r += 360.0;
        result.r /= 360.0;
    }
    return result;
}

/***
 * Performs direct volume rendering and
 * returns the final fragment color.
 ***/
vec4 directRendering(vec2 p, in vec3 first, in vec3 last, float entryDepth, float exitDepth) {

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

    // 2 nested loops allow for more than 255 iterations,
    // but is slower than while (t < tend)
    for (int loop0=0; !finished && loop0<255; loop0++) {
        for (int loop1=0; !finished && loop1<255; loop1++) {

            vec3 sampleVal = first.rgb + t * direction;
            vec4 voxel = texture(volume_, sampleVal);

            // no shading is applied
            vec4 color = voxel;
            vec4 tfColor = applyTF(transferFunc_, transferFuncTex_, rgb2hsv(color.rgb).r);
            if (applyColorModulation_)
                color.rgb *= tfColor.rgb;
            color.a = tfColor.a;

            // perform compositing
            if (color.a > 0.0) {
                // to accomodate for variable slice spacing
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
    }

    // calculate depth value from ray parameter
    gl_FragDepth = 1.0;
    if (depthT >= 0.0)
        gl_FragDepth = getDepthValue(depthT, tend, entryDepth, exitDepth);


    return result;
}

/***
 * The main method.
 ***/
void main() {
    // fetch entry/exit points
    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    vec3 frontPos = textureLookup2Dnormalized(entryPoints_, entryParameters_, p).rgb;
    vec3 backPos = textureLookup2Dnormalized(exitPoints_, exitParameters_, p).rgb;
    float entryDepth = textureLookup2Dnormalized(entryPointsDepth_, entryParameters_, p).x;
    float exitDepth = textureLookup2Dnormalized(exitPointsDepth_, exitParameters_, p).x;

    //determine whether the ray has to be casted
    if (frontPos == backPos) {
        //background need no raycasting
        discard;
    } else {
        //fragCoords are lying inside the boundingbox
        FragData0 = directRendering(p, frontPos, backPos, entryDepth, exitDepth);
    }
}
