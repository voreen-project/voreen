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

#include "mod_colorcoding.frag"

#define PLANE_SAGITTAL  0
#define PLANE_CORONAL   2
#define PLANE_AXIAL     1

uniform sampler3D vectorField3D_;
uniform int sliceNo_;
uniform bool useBilinearInterpolation_;
uniform vec3 flowDimensions_;

#ifdef PLANE
#if PLANE == PLANE_AXIAL
int i = 0, j = 2, k = 1;
#elif PLANE == PLANE_CORONAL
int i = 0, j = 1, k = 2;
#else   // default is permutation for sagittal plane section through volume
int i = 2, j = 1, k = 0;
#endif
#else
int i = 2, j = 1, k = 0;
#endif

/*
vec2 lookupVectorField2D(const vec2 pos) {
    vec2 p = vec2((pos * 2.0) - vec2(1.0));            // map to [-1.0, 1.0]
    float radius = length(p);

    if ((0.0 <= radius) && (radius <= 1.0))
        return vec2(-p.y, p.x);

    return vec2(0.0, 0.0);
}

vec3 lookupVectorField(const vec3 pos) {
    vec2 v = lookupVectorField2D(vec2(pos[i], pos[j]));
    vec3 res = vec3(0.0);
    res[i] = v.x;
    res[j] = v.y;
    return res;
}
*/
/*
float randomSeed = 0.1415926538979;

void seedRandom(const float seed) {
    randomSeed *= fract(3.14519265358979 * seed);
}

float rand() {
    float a = 7271.2651132 * (randomSeed + 0.12345678945687);
    randomSeed = (mod(a, 671.0) + 1.0) / 671.0;
    return randomSeed;
}
*/

float getFlowMagnitude(const vec3 r) {
    vec3 v = texture3D(vectorField3D_, r).xyz;
    if (v != vec3(0.0)) {
        v = (v * (maxValue_ - minValue_)) + minValue_;
        return length(v);
    }
    return 0.0;
}

vec4 bilinearInterpolation(const vec3 r0) {
    vec3 voxelCoord = r0 * flowDimensions_;

    // fill in the components at the right positions. the latter
    // depends on the current slice orientation
    //
    vec3 ll = floor(voxelCoord) / flowDimensions_;   // lower-left
    ll[k] = r0[k];

    vec3 ur = ceil(voxelCoord) / flowDimensions_;   // upper-right
    ur[k] = r0[k];

    vec3 lr = vec3(0.0);    // lower-right
    lr[i] = ur[i];
    lr[j] = ll[j];
    lr[k] = r0[k];

    vec3 ul = vec3(0.0);    // upper-left
    ul[i] = ll[i];
    ul[j] = ur[j];
    ul[k] = r0[k];

    vec2 v = fract(vec2(voxelCoord[i], voxelCoord[j]));
    float mll = getFlowMagnitude(ll);
    float mlr = getFlowMagnitude(lr);
    float mur = getFlowMagnitude(ur);
    float mul = getFlowMagnitude(ul);

    vec4 c1 = mix(getColorFromFlowMagnitude(mll), getColorFromFlowMagnitude(mlr), v.x);
    vec4 c2 = mix(getColorFromFlowMagnitude(mul), getColorFromFlowMagnitude(mur), v.x);
    return mix(c1, c2, v.y);
}

vec4 getColorFromFlowDirection(const vec3 r) {
    vec3 v = texture3D(vectorField3D_, r).xyz;
    return vec4(abs(v), 1.0);
}

void main() {
    // get the fragments position within the vector field
    //
    vec3 r0 = vec3(0.0);
    r0[i] = gl_TexCoord[0].s;
    r0[j] = gl_TexCoord[0].t;
    //r0[k] = 0.5;
    r0[k] = float(sliceNo_) / flowDimensions_[k];

    // look up the flow at the position and calculate the color
    //
#if COLOR_MODE == COLOR_MODE_DIRECTION
    FragData0 = getColorFromFlowDirection(r0);
#elif COLOR_MODE == COLOR_MODE_MONOCHROME
    FragData0 = clamp(color_, vec4(0.0), vec4(1.0));
#else
    if (useBilinearInterpolation_ == true)
        FragData0 = bilinearInterpolation(r0);
    else {
        vec3 v = texture3D(vectorField3D_, r).xyz;
        FragData0 = getFlowColor(v);
    }
#endif
}
