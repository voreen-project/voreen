/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

/*
 * This module contains compositing functions used for multimodal sliceviewing.
 */

vec4 maximumAlpha(vec4 previous, vec4 current, float factor) {
    if(previous.a > current.a)
        return previous;
    else
        return current;
}

vec4 minimumAlpha(vec4 previous, vec4 current, float factor) {
    if(previous.a < current.a)
        return previous;
    else
        return current;
}

vec4 replace(vec4 previous, vec4 current, float factor) {
    return current;
}

vec4 deactivate(vec4 previous, vec4 current, float factor) {
    return previous;
}

vec4 add(vec4 previous, vec4 current, float factor) {
    return previous + current;
}

vec4 subtract(vec4 previous, vec4 current, float factor) {
    return previous - current;
}

vec4 multiply(vec4 previous, vec4 current, float factor) {
    return previous * current;
}

vec4 divide(vec4 previous, vec4 current, float factor) {
    return previous / current;
}

vec4 difference(vec4 previous, vec4 current, float factor) {
    return abs(previous - current);
}

vec4 decal(vec4 previous, vec4 current, float factor) {
    vec4 ret = vec4(vec3(previous.rgb + (1.0 - previous.a) * current.a * current.rgb), 1.0);
    ret.a = previous.a + (1.0 - previous.a) * current.a;
    return ret;
}

vec4 over(vec4 previous, vec4 current, float factor) {
    float alpha0 = 1.0 - (1.0 - current.a)*(1.0 - previous.a);
    return (current.a / alpha0) * current + (1.0 - (current.a / alpha0)) * previous;
}

vec3 rgbToHsv(vec3 rgb) {
    float max = max(rgb.r, max(rgb.g, rgb.b));
    float min = min(rgb.r, min(rgb.g, rgb.b));

    float h = 0.0;
    if (max == min)
        h = 0.0;
    else if (max == rgb.r)
        h = 60*((rgb.g - rgb.b)/(max - min));
    else if (max == rgb.g)
        h = 60*(2 + (rgb.b - rgb.r)/(max - min));
    else
        h = 60*(4 + (rgb.r - rgb.g)/(max - min));

    if (h < 0.0)
        h += 360.0;
    float s = (max == 0.0) ? 0.0 : (max - min)/max;
    float v = max;

    return vec3(h, s, v);
}

vec3 hsvToRgb(vec3 hsv) {
    float hi = floor(hsv.x/60);
    float f = (hsv.x/60) - hi;

    float p = hsv.z * (1-hsv.y);
    float q = hsv.z * (1-(hsv.y*f));
    float t = hsv.z * (1-(hsv.y*(1-f)));

    if (hi == 1)
        return vec3(q, hsv.z, p);
    else if (hi == 2)
        return vec3(p, hsv.z, t);
    else if (hi == 3)
        return vec3(p, q, hsv.z);
    else if (hi == 4)
        return vec3(t, p, hsv.z);
    else if (hi == 5)
        return vec3(hsv.z, p, q);
    else
        return vec3(hsv.z, t, p);
}

vec4 setHue(vec4 previous, vec4 current, float factor) {
    vec3 hsvCur = rgbToHsv(current.rgb);
    vec3 hsvPrev = rgbToHsv(previous.rgb);
    vec3 hsv = vec3(hsvCur.x, hsvCur.y, hsvPrev.z);
    return vec4(hsvToRgb(hsv), previous.a);
}

vec4 blend(vec4 previous, vec4 current, float factor) {
    return mix(previous, current, factor);
}

//-------------------

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
