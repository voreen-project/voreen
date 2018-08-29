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

uniform sampler2D firstHitPointsNormalEstimation_;

#ifdef VRN_GL_NVIDIA_VERSION
#if (VRN_GL_NVIDIA_VERSION == 9755 || VRN_GL_NVIDIA_VERSION == 1001419)
#define VRN_WORKAROUND_NORMALEST_UNROLL
#endif
#endif

// Image-based normal estimation based on depth buffer values
vec3 calcFirstHitNormal(in vec2 pos /*in vec3 direction*/) {

    const int N = 5;
    vec4 v[4];

// Workaround for compiler bug on NVIDIA 97.55/Linux (cgc 1.6.0011) and cgc 1.5.0019 (Feb 22
// 2007), where the compiler would suck up all memory.
#ifndef VRN_WORKAROUND_NORMALEST_UNROLL
    for (int i=0; i < 4; i++)
        v[i] = vec4(0);
#else
    // manual unrolling
    v[0] = vec4(0.0);
    v[1] = vec4(0.0);
    v[2] = vec4(0.0);
    v[3] = vec4(0.0);
#endif

    vec2 offset[4];
    offset[0] = vec2(1, 0);
    offset[1] = vec2(0, 1);
    offset[2] = vec2(-1, 0);
    offset[3] = vec2(0, -1);

    // TODO it is very likely that the index of gl_TextureMatrix must be set with
    // the proper texture unit used. Since this must be a compile time constant
    // on almost all cards a define must be used. (roland)
    vec2 p = (gl_TextureMatrix[0] * gl_FragCoord).xy;
    for (int i=1; i <= N; i++) {
#ifndef VRN_WORKAROUND_NORMALEST_UNROLL
        for (int j=0; j < 4; j++) {

            vec3 point = textureLookup2D(firstHitPointsNormalEstimation_, p + offset[j] * float(i)).xyz;
            point = (point - 0.5) * 1.0;

            //v[j].w += sign(length(point));
            if (point != vec3(0.0)) v[j].w++;
            v[j].xyz += point;
        }
#else
        // manual unrolling
        vec3 point = textureLookup2D(firstHitPointsNormalEstimation_, p + offset[0] * float(i)).xyz;
        point = (point - 0.5) * 1.0;
        //v[0].w += sign(length(point));
        if (point != vec3(0.0)) v[0].w++;
        v[0].xyz += point;

        point = textureLookup2D(firstHitPointsNormalEstimation_, p + offset[1] * float(i)).xyz;
        point = (point - 0.5) * 1.0;
        //v[1].w += sign(length(point));
        if (point != vec3(0.0)) v[1].w++;
        v[1].xyz += point;

        point = textureLookup2D(firstHitPointsNormalEstimation_, p + offset[2] * float(i)).xyz;
        point = (point - 0.5) * 1.0;
        //v[2].w += sign(length(point));
        if (point != vec3(0.0)) v[2].w++;
        v[2].xyz += point;

        point = textureLookup2D(firstHitPointsNormalEstimation_, p + offset[3] * float(i)).xyz;
        point = (point - 0.5) * 1.0;
        //v[3].w += sign(length(point));
        if (point != vec3(0.0)) v[3].w++;
        v[3].xyz += point;
#endif
    }

    vec3 c0 = textureLookup2D(firstHitPointsNormalEstimation_, p).rgb;
    c0 = (c0 - 0.5) * 1.0;
    vec4 n = vec4(0);
#ifndef VRN_WORKAROUND_NORMALEST_UNROLL
    for (int i=0; i < 2; i++) {
        if (v[i*2].w > 0.0 && v[i*2+1].w > 0.0) {
            n += vec4(cross(normalize(v[i*2].xyz / v[i*2].w - c0),
                            normalize(v[i*2+1].xyz / v[i*2+1].w - c0)), 1.0);
        }
    }
#else
    // manual unrolling
    if (v[0].w > 0.0 && v[1].w > 0.0)
        n += vec4(cross(normalize(v[0].xyz / v[0].w - c0),
                        normalize(v[1].xyz / v[1].w - c0)), 1.0);
    if (v[2].w > 0.0 && v[3].w > 0.0)
        n += vec4(cross(normalize(v[2].xyz / v[2].w - c0),
                        normalize(v[3].xyz / v[3].w - c0)), 1.0);
#endif
    n.xyz = normalize(n.xyz / n.w);
//    n.xyz = faceforward(n.xyz, direction, n.xyz);
    return normalize(n.xyz) * 0.5 + vec3(0.5);
}

void main() {

    FragData0 = vec4(calcFirstHitNormal(gl_FragCoord.xy), 1.0);
}
