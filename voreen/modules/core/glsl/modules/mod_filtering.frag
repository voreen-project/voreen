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
  * RGB-to-Gray conversion, the alpha-channel is not modified.
  */
vec4 rgbToGrayScale(in vec4 color) {
    float brightness = (0.30 * color.r) + (0.59 * color.g) + (0.11 * color.b);
    return vec4(brightness, brightness, brightness, color.a);
}

/**
 * RGB-to-Gray conversion, modulated by an additional saturation factor.
 */
 vec4 rgbToGrayScaleSaturated(in vec4 color, in float saturation) {
    vec4 grayscale = rgbToGrayScale(color);
    return mix(grayscale, color, saturation);
}

/***
 * Max value of halfKernelDim is 3
 */
vec4 median(in sampler2D texture, in TextureParameters texParams, in vec2 fragCoord, in int halfKernelDim) {
    vec4[49] values;
    int curValID = 0;
    for (int x=-halfKernelDim; x<=halfKernelDim; x++)
        for (int y=-halfKernelDim; y<=halfKernelDim; y++)
            values[curValID++] = textureLookup2Dscreen(texture, texParams,
                                                       vec2(fragCoord.x+float(x), fragCoord.y+float(y)));
    // sort fetched values
    for (int i=curValID-2; i>=0; i--)
        for (int j=0; j<i; j++) {
            vec4 tmp = min(values[j], values[j+1]);
            values[j+1] = values[j] + values[j+1] - tmp;
            values[j] = tmp;
        }
    return values[int(floor(curValID/2.0))];
}


vec4 mean(in sampler2D texture, in TextureParameters texParams, in vec2 fragCoord, in int halfKernelDim) {
    vec4 result = vec4(0.0);
    for (int x=-halfKernelDim; x<=halfKernelDim; x++)
        for (int y=-halfKernelDim; y<=halfKernelDim; y++)
            result += textureLookup2Dscreen(texture, texParams,
                                            vec2(fragCoord.x+float(x), fragCoord.y+float(y)));
    result /= pow((halfKernelDim*2.0)+1.0, 2.0);
    return result;
}

vec4 erosion(in sampler2D texture, in TextureParameters texParams, in vec2 fragCoord, in int kernelRadius) {
    vec4 result = vec4(1.0);
    for (int x=-kernelRadius; x<=kernelRadius; x++)
        for (int y=-kernelRadius; y<=kernelRadius; y++) {
            vec4 curValue = textureLookup2Dscreen(texture, texParams,
                                                  vec2(fragCoord.x+float(x), fragCoord.y+float(y)));
            result = min(curValue, result);
        }
    return result;
}

vec4 dilation(in sampler2D texture, in TextureParameters texParams, in vec2 fragCoord, in int kernelRadius) {
    vec4 result = vec4(0.0);
    for (int x=-kernelRadius; x<=kernelRadius; x++)
        for (int y=-kernelRadius; y<=kernelRadius; y++) {
            vec4 curValue = textureLookup2Dscreen(texture, texParams,
                                                  vec2(fragCoord.x+float(x), fragCoord.y+float(y)));
            result = max(curValue, result);
        }
    return result;
}

vec4 erosion_circle(in sampler2D texture, in TextureParameters texParams, in vec2 fragCoord, in int kernelRadius) {
    vec4 result = vec4(1.0);
    for (int x=-kernelRadius; x<=kernelRadius; x++)
        for (int y=-kernelRadius; y<=kernelRadius; y++) {
            if(length(vec2(x, y)) > kernelRadius + 0.5)
                continue;
            vec4 curValue = textureLookup2Dscreen(texture, texParams, vec2(fragCoord.x+float(x), fragCoord.y+float(y)));
            result = min(curValue, result);
        }
    return result;
}

vec4 dilation_circle(in sampler2D texture, in TextureParameters texParams, in vec2 fragCoord, in int kernelRadius) {
    vec4 result = vec4(0.0);
    for (int x=-kernelRadius; x<=kernelRadius; x++)
        for (int y=-kernelRadius; y<=kernelRadius; y++) {
            if(length(vec2(x, y)) > kernelRadius + 0.5)
                continue;
            vec4 curValue = textureLookup2Dscreen(texture, texParams, vec2(fragCoord.x+float(x), fragCoord.y+float(y)));
            result = max(curValue, result);
        }
    return result;
}



