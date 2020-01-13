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

#ifndef MOD_COLORTABLES_FRAG
#define MOD_COLORTABLES_FRAG

#define COLOR_TABLE_RAINBOW     0
#define COLOR_TABLE_HOT_METAL   1
#define COLOR_TABLE_TEXTURE     2

#define COLOR_MODE_MAGNITUDE  0
#define COLOR_MODE_DIRECTION  1
#define COLOR_MODE_MONOCHROME 2

uniform float minValue_;
uniform float maxValue_;
uniform float maxMagnitude_;
uniform vec2 thresholds_;
#if COLOR_MODE == COLOR_MODE_MONOCHROME
uniform vec4 color_;
#endif
bool useThresholds_ = false;

#if COLOR_MODE == COLOR_MODE_MAGNITUDE
#if COLOR_TABLE == COLOR_TABLE_TEXTURE
uniform sampler1D colorTexture_;
#endif
#endif

#ifdef COLOR_MODE_DIRECTION
vec4 dirToColor(in vec3 dir) {
    vec3 tmp = dir/vec3(2.0)+vec3(0.5);
    return vec4(tmp, 1.0);
}
#endif

#if COLOR_MODE == COLOR_MODE_MAGNITUDE
vec4 getColorFromFlowMagnitude(const float m) {
    #if COLOR_TABLE == COLOR_TABLE_TEXTURE
    return texture(colorTexture_, m);
    #else
        #if COLOR_TABLE == COLOR_TABLE_HOT_METAL
        const int colorTableSize = 4;
        vec3 colorTable[colorTableSize];
        colorTable[0] = vec3(0.0, 0.0, 0.0);    // black
        colorTable[1] = vec3(1.0, 0.0, 0.0);    // red
        colorTable[2] = vec3(1.0, 1.0, 0.0);    // yellow
        colorTable[3] = vec3(1.0, 1.0, 1.0);    // white
        #elif COLOR_TABLE == COLOR_TABLE_RAINBOW
        const int colorTableSize = 6;
        vec3 colorTable[colorTableSize];
        colorTable[0] = vec3(0.0, 0.0, 0.0);    // black
        colorTable[1] = vec3(0.0, 0.0, 1.0);    // blue
        colorTable[2] = vec3(0.0, 1.0, 1.0);    // cyan
        colorTable[3] = vec3(0.0, 1.0, 0.0);    // green
        colorTable[4] = vec3(1.0, 1.0, 0.0);    // yellow
        colorTable[5] = vec3(1.0, 0.0, 0.0);    // red
        #endif

        float numColors = float(colorTableSize - 1);
        float v = clamp(m * numColors, 0.0, numColors);
        ivec2 limits = clamp(ivec2(int(v), int(ceil(v))), 0, colorTableSize);
        vec3 color = mix(colorTable[limits.x], colorTable[limits.y], fract(v));
        return vec4(color, 1.0);
    #endif
}
#endif // magnitude

vec4 getFlowColor(const vec3 direction) {
    float magnitude = length(direction);
    float maxMagnitude = maxMagnitude_;
    if (maxMagnitude <= 0.0) {
        float maxComponentValue = max(abs(maxValue_), abs(minValue_));
        maxMagnitude = sqrt(3.0 * (maxComponentValue * maxComponentValue));
    }
    useThresholds_ = (thresholds_ != vec2(0.0));
    if (useThresholds_ == true) {
        maxMagnitude = thresholds_.y;
        if ((magnitude < thresholds_.x) || (magnitude > thresholds_.y))
            return vec4(0.0, 0.0, 0.0, 1.0);
    }
    if  (magnitude > maxMagnitude)
        return vec4(1.0, 1.0, 1.0, 1.0);
    else
#if COLOR_MODE == COLOR_MODE_MAGNITUDE
    if ((useThresholds_ == true) && (thresholds_.y > thresholds_.x))
        return getColorFromFlowMagnitude((magnitude - thresholds_.x) / (thresholds_.y - thresholds_.x));
    else
        return getColorFromFlowMagnitude(magnitude / maxMagnitude);
#endif
#if COLOR_MODE == COLOR_MODE_DIRECTION
        return dirToColor(normalize(direction));
#endif
#if COLOR_MODE == COLOR_MODE_MONOCHROME
        return color_;
#endif
}

#endif
