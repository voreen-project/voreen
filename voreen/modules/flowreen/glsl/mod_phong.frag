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

#ifndef MOD_PHONG_FRAG
#define MOD_PHONG_FRAG

varying vec3 normalVec;
varying vec3 lightVec;
varying vec3 viewVec;

uniform bool useForLines_;
uniform bool usePhongShading_;
uniform float ka_;
uniform float kd_;
uniform float ks_;
uniform float shininess_;

vec3 phongShading() {
    // normalize the vectors again which are interpolated per pixel and are
    // therefore different for each fragment when transfered from vertex to
    // fragment shader!
    // lightVec does not need to be re-normalized, as it remains the same for
    // all fragments.
    //
    vec3 n = normalize(normalVec);  // normalize again as the normal is interpolated per pixel!
    vec3 v = normalize(viewVec);    // the same for the view vector: it is interpolated per pixel!

    float NdotL = max(dot(lightVec, n), 0.0);
    float specular = 0.0;
    if (NdotL > 0.0)    // omit calculations if no specular highlights are visible
    {
        vec3 reflectVec = normalize((2.0 * NdotL * n) - lightVec);
        float RdotV = dot(reflectVec, v);
        specular = max(pow(RdotV, shininess_), 0.0);
    }

    return vec3(1.0, NdotL, specular);
}

/**
 * Calculate Phong terms for simple OpenGL line primitives according to
 * the work of Zöckler, Stalling and Hege in "Interactive Visualization Of
 * 3D-Vector Fields Using Illuminated Stream Lines", from 1996.
 */
vec3 phongShadingForLines() {
    // normalize the vectors again which are interpolated per pixel and are
    // therefore different for each fragment when transfered from vertex to
    // fragment shader!
    // lightVec does not need to be re-normalized, as it remains the same for
    // all fragments.
    //
    vec3 t = normalize(normalVec);  // normalize again as the normal is interpolated per pixel!
    vec3 v = normalize(viewVec);    // the same for the view vector: it is interpolated per pixel!
    float LdotT = clamp(dot(lightVec, t), -1.0, 1.0);
    float NdotL = max(sqrt(1.0 - (LdotT * LdotT)), 0.0);

    float VdotT = clamp(dot(v, t), -1.0, 1.0);
    float VdotN = max(sqrt(1.0 - (VdotT * VdotT)), 0.0);
    float RdotV = (LdotT * VdotT) - (NdotL * VdotN);
    float specular = max(pow(RdotV, shininess_), 0.0);
    const float p = 4.8;
    return vec3(1.0, pow(NdotL, p), specular);
}

#endif
