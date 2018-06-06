/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

// Adapted from:
// http://horde3d.org/wiki/index.php5?title=Shading_Technique_-_FXAA

#include "modules/mod_sampler2d.frag"

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform TextureParameters texParams_;

uniform float spanMax_;
uniform float reduceMul_;

void main() {
    const float FXAA_REDUCE_MIN = (1.0/128.0);

    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    vec2 texcoordOffset = screenDimRCP_;

    vec3 rgbNW = texture(colorTex_, p + (vec2(-1.0, -1.0) * texcoordOffset)).xyz;
    vec3 rgbNE = texture(colorTex_, p + (vec2(+1.0, -1.0) * texcoordOffset)).xyz;
    vec3 rgbSW = texture(colorTex_, p + (vec2(-1.0, +1.0) * texcoordOffset)).xyz;
    vec3 rgbSE = texture(colorTex_, p + (vec2(+1.0, +1.0) * texcoordOffset)).xyz;
    vec3 rgbM = texture(colorTex_, p).xyz;

    vec3 luma = vec3(0.299, 0.587, 0.114);
    float lumaNW = dot(rgbNW, luma);
    float lumaNE = dot(rgbNE, luma);
    float lumaSW = dot(rgbSW, luma);
    float lumaSE = dot(rgbSE, luma);
    float lumaM = dot( rgbM, luma);

    float lumaMin = min(lumaM, min(min(lumaNW, lumaNE), min(lumaSW, lumaSE)));
    float lumaMax = max(lumaM, max(max(lumaNW, lumaNE), max(lumaSW, lumaSE)));

    vec2 dir;
    dir.x = -((lumaNW + lumaNE) - (lumaSW + lumaSE));
    dir.y = ((lumaNW + lumaSW) - (lumaNE + lumaSE));

    float dirReduce = max((lumaNW + lumaNE + lumaSW + lumaSE) * (0.25 / reduceMul_), FXAA_REDUCE_MIN);

    float rcpDirMin = 1.0/(min(abs(dir.x), abs(dir.y)) + dirReduce);

    dir = min(vec2(spanMax_, spanMax_),
    max(vec2(-spanMax_, -spanMax_), dir * rcpDirMin)) * texcoordOffset;

    vec3 rgbA = (1.0/2.0) * (
        texture(colorTex_, p + dir * (1.0/3.0 - 0.5)).xyz +
        texture(colorTex_, p + dir * (2.0/3.0 - 0.5)).xyz);
    vec3 rgbB = rgbA * (1.0/2.0) + (1.0/4.0) * (
        texture(colorTex_, p + dir * (0.0/3.0 - 0.5)).xyz +
        texture(colorTex_, p + dir * (3.0/3.0 - 0.5)).xyz);
    float lumaB = dot(rgbB, luma);

    if((lumaB < lumaMin) || (lumaB > lumaMax)){
        FragData0.xyz = rgbA;
    } else {
        FragData0.xyz = rgbB;
    }
    FragData0.a = 1.0;

    gl_FragDepth = textureLookup2Dnormalized(depthTex_, texParams_, p).x;
}
