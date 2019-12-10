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

#include "modules/mod_sampler2d.frag"

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform TextureParameters colorParams_;

// TODO
//uniform sampler2D randomTex_;

uniform mat4 projectionMatrixInverse_;
uniform mat4 viewMatrix_;
uniform mat4 viewMatrixInverse_;

uniform float occFac_;
uniform float colFac_;
uniform int numIts_;

vec2 fragCoord;

float rand(vec2 co) {
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

// unproject eye position from screen coords
vec3 posFromDepth(vec2 coord){
    float d = texture2D(depthTex_, coord).x;

    // use NDCs from -1 to 1, unproject to view coords
    vec4 ndc = vec4(2.0*coord - 1.0, 2.0*d - 1.0, 1.0);
    ndc = (projectionMatrixInverse_ * ndc);
    ndc /= ndc.w;
    return ndc.xyz;
}

// estimate normal from depth
vec3 readNormal(in vec2 coord) {
    vec2 offset1 = vec2(0.0, colorParams_.dimensionsRCP_.y);
    vec2 offset2 = vec2(colorParams_.dimensionsRCP_.x, 0.0);

    float depth  = texture2D(depthTex_, fragCoord).x;
    float depth1 = texture2D(depthTex_, fragCoord + offset1).x;
    float depth2 = texture2D(depthTex_, fragCoord + offset2).x;

    vec3 p1 = vec3(offset1, depth1 - depth);
    vec3 p2 = vec3(offset2, depth2 - depth);

    vec3 normal = cross(p2, p1);

    return normalize(normal);
}

//Ambient Occlusion form factor:
float occFF(in vec3 ddiff,in vec3 cnorm, in float c1, in float c2){
    float rd = length(ddiff);
    if(rd < 0.03)
        return 0.0;
    vec3 vv = normalize(ddiff);
    return clamp(dot(readNormal(fragCoord + vec2(c1,c2)), vv), 0.0, 1.0) * clamp(dot(cnorm, vv), 0.0, 1.0) * (1.0 - 1.0 / sqrt(1.0 / (rd*rd) + 1.0));
}

void main() {
    fragCoord = gl_FragCoord.xy * colorParams_.dimensionsRCP_;
    vec4 current = texture2D(colorTex_, fragCoord);
    if(current.a == 0.0)
        discard;

    //read current normal,position and color.
    vec3 n = readNormal(fragCoord);
    vec3 p = posFromDepth(fragCoord);
    vec3 col = current.xyz;

    vec2 random = vec2(rand(fragCoord), rand(fragCoord.yx));
    random = 2.0 * random - 1.0;

    //initialize variables:
    vec4 ao = vec4(0.0);

    float incx = colorParams_.dimensionsRCP_.x;
    float incy = colorParams_.dimensionsRCP_.y;
    float pw = incx;
    float ph = incy;

    float cdepth = 0.1 * abs(p.z);
    float sumWeight = 0.f;

    for(int i = -numIts_; i < numIts_; ++i) {
        for(int j = -numIts_; j < numIts_; ++j) {
            if((length(vec2(i, j)) > float(numIts_) + 0.5) || (i == 0 && j == 0))
                continue;

            float npw = (i * pw + 0.5 * colorParams_.dimensionsRCP_.x * random.x) / cdepth;
            float nph = (j * ph + 0.5 * colorParams_.dimensionsRCP_.y * random.y) / cdepth;

            random = 0.5 * random + 0.5;
            random = vec2(rand(random), rand(random.yx));
            random = 2.0 * random - 1.0;

            vec3 ddiff  = posFromDepth(fragCoord + vec2(npw,nph)) - p;

            float ff = (1.0 / length(vec2(i, j))) * occFF(ddiff,n,npw,nph);
            ao.a   += ff;
            ao.xyz += ff * texture2D(colorTex_, fragCoord+vec2(npw,nph)).xyz;

            sumWeight += 1.0 / length(vec2(i, j));
        }
    }

    ao /= sumWeight;

    gl_FragColor = vec4(col - occFac_ * vec3(ao.a) + colFac_ * ao.xyz, current.a);
    gl_FragDepth = texture2D(depthTex_, fragCoord).x;
}
