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

float fetchCurvature(in sampler3D volume, in VolumeParameters curvatureVolumeParameters, in vec3 samplePos) {
    float curvature = getVoxel(volume, curvatureVolumeParameters, samplePos).a;
    curvature -= 0.5;
    curvature *= 2.0;
    return curvature;
}

mat3 computeHessian(in sampler3D volume, in VolumeParameters gradientVolumeParameters, in vec3 samplePos) {
    vec3 offset = gradientVolumeParameters.datasetDimensionsRCP_;
    vec3 gradientr = textureLookup3DUnnormalized(volume, gradientVolumeParameters, samplePos + vec3(offset.x,0.0,0.0)).xyz;
    vec3 gradientl = textureLookup3DUnnormalized(volume, gradientVolumeParameters, samplePos + vec3(-offset.x,0.0,0.0)).xyz;
    vec3 gradientu = textureLookup3DUnnormalized(volume, gradientVolumeParameters, samplePos + vec3(0.0,offset.y,0.0)).xyz;
    vec3 gradientd = textureLookup3DUnnormalized(volume, gradientVolumeParameters, samplePos + vec3(0.0,-offset.y,0.0)).xyz;
    vec3 gradientf = textureLookup3DUnnormalized(volume, gradientVolumeParameters, samplePos + vec3(0.0,0.0,offset.z)).xyz;
    vec3 gradientb = textureLookup3DUnnormalized(volume, gradientVolumeParameters, samplePos + vec3(0.0,0.0,-offset.z)).xyz;
    mat3 H;
    H[0][0] = gradientl.x - gradientr.x;
    H[0][1] = gradientl.y - gradientr.y;
    H[0][2] = gradientl.z - gradientr.z;
    H[1][0] = gradientl.y - gradientr.y;
    H[1][1] = gradientd.y - gradientu.y;
    H[1][2] = gradientd.z - gradientu.z;
    H[2][0] = gradientl.z - gradientr.z;
    H[2][1] = gradientd.z - gradientu.z;
    H[2][2] = gradientb.z - gradientf.z;
    return H;
}

vec2 computeCurvature(in sampler3D volume, in VolumeParameters gradientVolumeParameters, in vec3 samplePos) {
    mat3 H = computeHessian(volume, gradientVolumeParameters, samplePos);
    vec3 gradient = (textureLookup3DUnnormalized(volume, gradientVolumeParameters, samplePos).xyz-vec3(0.5))*2.0;

    vec3 n = -gradient / length(gradient);
    mat3 nnT;
    nnT[0][0] = n[0]*n[0];
    nnT[0][1] = n[0]*n[1];
    nnT[0][2] = n[0]*n[2];
    nnT[1][0] = n[1]*n[0];
    nnT[1][1] = n[1]*n[1];
    nnT[1][2] = n[1]*n[2];
    nnT[2][0] = n[2]*n[0];
    nnT[2][1] = n[2]*n[1];
    nnT[2][2] = n[2]*n[2];
    mat3 P = mat3(1.0) - nnT;

    mat3 G = -P*H*P / length(gradient);

    // compute trace of G
    float trace = G[0][0] + G[1][1] + G[2][2];

    // compute Frobenius norm of G
    float F = 0.0;
    F += pow(abs(G[0][0]), 2.0);
    F += pow(abs(G[0][1]), 2.0);
    F += pow(abs(G[0][2]), 2.0);
    F += pow(abs(G[1][0]), 2.0);
    F += pow(abs(G[1][1]), 2.0);
    F += pow(abs(G[1][2]), 2.0);
    F += pow(abs(G[2][0]), 2.0);
    F += pow(abs(G[2][1]), 2.0);
    F += pow(abs(G[2][2]), 2.0);
    F = sqrt(F);

    vec2 curvature;
    curvature.x = (trace + sqrt(2.0 * (F*F) - (trace*trace))) / 2.0; // kappa1
    curvature.y = (trace - sqrt(2.0 * (F*F) - (trace*trace))) / 2.0; // kappa2
    return curvature;
}

float computeViewCurvature(in sampler3D volume, in VolumeParameters gradientVolumeParameters, in vec3 samplePos, in vec3 v, in float T) {
    mat3 H = computeHessian(volume, gradientVolumeParameters, samplePos);
    vec3 gradient = (textureLookup3DUnnormalized(volume, gradientVolumeParameters, samplePos).xyz-vec3(0.5))*2.0;

    vec3 n = -gradient / length(gradient);
    mat3 nnT;
    nnT[0][0] = n[0]*n[0];
    nnT[0][1] = n[0]*n[1];
    nnT[0][2] = n[0]*n[2];
    nnT[1][0] = n[1]*n[0];
    nnT[1][1] = n[1]*n[1];
    nnT[1][2] = n[1]*n[2];
    nnT[2][0] = n[2]*n[0];
    nnT[2][1] = n[2]*n[1];
    nnT[2][2] = n[2]*n[2];
    mat3 P = mat3(1.0) - nnT;

    mat3 G = -P*H*P / length(gradient);

    float kv = dot(v, G * v) / dot(v, P * v);
    float result = sqrt(T * kv * (2.0 - (T * kv)));
    return result;
}
