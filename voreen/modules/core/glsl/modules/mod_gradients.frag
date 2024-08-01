/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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
 * This module contains all functions which can be used for calculating
 * gradients on the fly within a raycaster.
 * The functions below are referenced by CALC_GRADIENT which is
 * used in the raycaster fragment shaders.
 */

/**
 * Calculates a correct gradient for clipping and volume boundaries.
 *
 * @param samplePos the sample's position in texture space
 * @param rayDirection the ray direction
 */
vec3 fixClipBorderGradient(vec3 samplePos, vec3 rayDirection, sampler2D entryPoints, TextureParameters entryParameters) {
    vec3 v0 = normalize(textureLookup2Dnormalized(entryPoints, entryParameters, vec2(gl_FragCoord.x+2.0, gl_FragCoord.y)*screenDimRCP_ ).rgb - samplePos);
    vec3 v1 = normalize(textureLookup2Dnormalized(entryPoints, entryParameters, vec2(gl_FragCoord.x, gl_FragCoord.y+2.0)*screenDimRCP_ ).rgb - samplePos);
    //FIXME: error handling if v0 or v1 is (0,0,0)
    vec3 gradient = cross(v0, v1);
    gradient = faceforward(gradient, rayDirection, gradient);
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space based on the red
 * channel using forward differences.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientRFD(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 offset = volumeStruct.datasetDimensionsRCP_;
    float v = textureLookup3DUnnormalized(volume, volumeStruct, samplePos).r;
    float v0 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(offset.x, 0.0, 0.0)).r;
    float v1 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, offset.y, 0)).r;
    float v2 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, offset.z)).r;
    vec3 gradient = vec3(v - v0, v - v1, v - v2);
    gradient *= volumeStruct.datasetSpacingRCP_;
    gradient *= volumeStruct.rwmScale_;
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space based on the green
 * channel using forward differences.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientGFD(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 offset = volumeStruct.datasetDimensionsRCP_;
    float v = textureLookup3DUnnormalized(volume, volumeStruct, samplePos).g;
    float v0 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(offset.x, 0.0, 0.0)).g;
    float v1 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, offset.y, 0)).g;
    float v2 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, offset.z)).g;
    vec3 gradient = vec3(v - v0, v - v1, v - v2);
    gradient *= volumeStruct.datasetSpacingRCP_;
    gradient *= volumeStruct.rwmScale_;
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space based on the blue
 * channel using forward differences.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientBFD(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 offset = volumeStruct.datasetDimensionsRCP_;
    float v = textureLookup3DUnnormalized(volume, volumeStruct, samplePos).b;
    float v0 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(offset.x, 0.0, 0.0)).b;
    float v1 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, offset.y, 0)).b;
    float v2 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, offset.z)).b;
    vec3 gradient = vec3(v - v0, v - v1, v - v2);
    gradient *= volumeStruct.datasetSpacingRCP_;
    gradient *= volumeStruct.rwmScale_;
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space based on the alpha
 * channel using forward differences.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientAFD(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 offset = volumeStruct.datasetDimensionsRCP_;
    float v = textureLookup3DUnnormalized(volume, volumeStruct, samplePos).a;
    float v0 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(offset.x, 0.0, 0.0)).a;
    float v1 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, offset.y, 0)).a;
    float v2 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, offset.z)).a;
    vec3 gradient = vec3(v - v0, v - v1, v - v2);
    gradient *= volumeStruct.datasetSpacingRCP_;
    gradient *= volumeStruct.rwmScale_;
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space based on the red
 * channel using central differences.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientR(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 offset = volumeStruct.datasetDimensionsRCP_;
    float v0 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(offset.x, 0.0, 0.0)).r;
    float v1 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, offset.y, 0)).r;
    float v2 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, offset.z)).r;
    float v3 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x, 0, 0)).r;
    float v4 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, -offset.y, 0)).r;
    float v5 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, -offset.z)).r;
    vec3 gradient = vec3(v3 - v0, v4 - v1, v5 - v2);
    gradient *= volumeStruct.datasetSpacingRCP_ * 0.5;
    gradient *= volumeStruct.rwmScale_;
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space based on the green
 * channel using central differences.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientG(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 offset = volumeStruct.datasetDimensionsRCP_;
    float v0 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(offset.x, 0.0, 0.0)).g;
    float v1 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, offset.y, 0)).g;
    float v2 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, offset.z)).g;
    float v3 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x, 0, 0)).g;
    float v4 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, -offset.y, 0)).g;
    float v5 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, -offset.z)).g;
    vec3 gradient = vec3(v3 - v0, v4 - v1, v5 - v2);
    gradient *= volumeStruct.datasetSpacingRCP_ * 0.5;
    gradient *= volumeStruct.rwmScale_;
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space based on the blue
 * channel using central differences.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientB(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 offset = volumeStruct.datasetDimensionsRCP_;
    float v0 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(offset.x, 0.0, 0.0)).b;
    float v1 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, offset.y, 0)).b;
    float v2 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, offset.z)).b;
    float v3 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x, 0, 0)).b;
    float v4 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, -offset.y, 0)).b;
    float v5 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, -offset.z)).b;
    vec3 gradient = vec3(v3 - v0, v4 - v1, v5 - v2);
    gradient *= volumeStruct.datasetSpacingRCP_ * 0.5;
    gradient *= volumeStruct.rwmScale_;
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space based on the alpha
 * channel using central differences.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientA(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 offset = volumeStruct.datasetDimensionsRCP_;
    float v0 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(offset.x, 0.0, 0.0)).a;
    float v1 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, offset.y, 0)).a;
    float v2 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, offset.z)).a;
    float v3 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x, 0, 0)).a;
    float v4 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, -offset.y, 0)).a;
    float v5 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(0, 0, -offset.z)).a;
    vec3 gradient = vec3(v3 - v0, v4 - v1, v5 - v2);
    gradient *= volumeStruct.datasetSpacingRCP_ * 0.5;
    gradient *= volumeStruct.rwmScale_;
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space based on the alpha
 * channel using the sobel operator.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientSobel(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 offset = volumeStruct.datasetDimensionsRCP_;

    /*
     * [ 1 2 3 ] [ 10 11 12 ] [ 19 20 21 ]
     * [ 4 5 6 ] [ 13 14 15 ] [ 22 23 24 ]
     * [ 7 8 9 ] [ 16 17 18 ] [ 25 26 27 ]
     */

    float v01 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3( offset.x, -offset.y,  offset.z)).a;
    float v02 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3( offset.x,       0.0,  offset.z)).a;
    float v03 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3( offset.x,  offset.y,  offset.z)).a;
    float v04 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3( offset.x, -offset.y,       0.0)).a;
    float v05 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3( offset.x,       0.0,       0.0)).a;
    float v06 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3( offset.x,  offset.y,       0.0)).a;
    float v07 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3( offset.x, -offset.y, -offset.z)).a;
    float v08 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3( offset.x,       0.0, -offset.z)).a;
    float v09 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3( offset.x,  offset.y, -offset.z)).a;

    float v10 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(      0.0, -offset.y,  offset.z)).a;
    float v11 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(      0.0,       0.0,  offset.z)).a;
    float v12 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(      0.0,  offset.y,  offset.z)).a;
    float v13 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(      0.0, -offset.y,       0.0)).a;
    //float v14 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(      0.0,       0.0,       0.0)).a;
    float v15 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(      0.0,  offset.y,       0.0)).a;
    float v16 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(      0.0, -offset.y, -offset.z)).a;
    float v17 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(      0.0,       0.0, -offset.z)).a;
    float v18 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(      0.0,  offset.y, -offset.z)).a;

    float v19 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x, -offset.y,  offset.z)).a;
    float v20 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x,       0.0,  offset.z)).a;
    float v21 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x,  offset.y,  offset.z)).a;
    float v22 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x, -offset.y,       0.0)).a;
    float v23 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x,       0.0,       0.0)).a;
    float v24 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x,  offset.y,       0.0)).a;
    float v25 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x, -offset.y, -offset.z)).a;
    float v26 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x,       0.0, -offset.z)).a;
    float v27 = textureLookup3DUnnormalized(volume, volumeStruct, samplePos + vec3(-offset.x,  offset.y, -offset.z)).a;

    float gx =      v01 + 2.0*v02 +     v03    /* 0 0 0 */   -     v19 - 2.0*v20 -     v21
               +2.0*v04 + 4.0*v05 + 2.0*v06    /* 0 0 0 */   - 2.0*v22 - 4.0*v23 - 2.0*v24
                   +v07 + 2.0*v08 +     v09    /* 0 0 0 */   -     v25 - 2.0*v26 -     v27;

    float gy =     -v01 /* 0 */ +     v03    - 2.0*v10 /* 0 */ + 2.0*v12    -     v19 /* 0 */ +     v21
               -2.0*v04 /* 0 */ + 2.0*v06    - 4.0*v13 /* 0 */ + 4.0*v15    - 2.0*v22 /* 0 */ + 2.0*v24
                   -v07 /* 0 */ +     v09    - 2.0*v16 /* 0 */ + 2.0*v18    -     v25 /* 0 */ +     v27;

    float gz = v01 + 2.0*v02 + v03    + 2.0*v10 + 4.0*v11 + 2.0*v12    + v19 + 2.0*v20 + v21
                   /* 0 0 0 */                 /* 0 0 0 */                    /* 0 0 0 */
              -v07 - 2.0*v08 - v09    - 2.0*v16 - 4.0*v17 - 2.0*v18    - v25 - 2.0*v26 - v27;

    vec3 gradient = -vec3(gx, gy, gz);
    gradient *= volumeStruct.datasetSpacingRCP_ * 0.5;
    gradient /=  16.0; // sum of all positive elements
    gradient *= volumeStruct.rwmScale_;
    return gradient;
}

/**
 * Calculates a voxel's gradient in volume object space by calculating the
 * gradients for all eight neighbours based on the alpha channel using
 * central differences and filtering between the reuslts.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradientFiltered(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    vec3 delta = volumeStruct.datasetDimensionsRCP_;
    vec3 g0 = calcGradientA(volume, volumeStruct, samplePos);
    vec3 g1 = calcGradientA(volume, volumeStruct, samplePos+vec3(-delta.x, -delta.y, -delta.z));
    vec3 g2 = calcGradientA(volume, volumeStruct, samplePos+vec3( delta.x,  delta.y,  delta.z));
    vec3 g3 = calcGradientA(volume, volumeStruct, samplePos+vec3(-delta.x,  delta.y, -delta.z));
    vec3 g4 = calcGradientA(volume, volumeStruct, samplePos+vec3( delta.x, -delta.y,  delta.z));
    vec3 g5 = calcGradientA(volume, volumeStruct, samplePos+vec3(-delta.x, -delta.y,  delta.z));
    vec3 g6 = calcGradientA(volume, volumeStruct, samplePos+vec3( delta.x,  delta.y, -delta.z));
    vec3 g7 = calcGradientA(volume, volumeStruct, samplePos+vec3(-delta.x,  delta.y,  delta.z));
    vec3 g8 = calcGradientA(volume, volumeStruct, samplePos+vec3( delta.x, -delta.y, -delta.z));
    vec3 mix0 = mix(mix(g1, g2, 0.5), mix(g3, g4, 0.5), 0.5);
    vec3 mix1 = mix(mix(g5, g6, 0.5), mix(g7, g8, 0.5), 0.5);
    return mix(g0, mix(mix0, mix1, 0.5), 0.75);
}

/**
 * Calculates a voxel's gradient in volume object space based on the alpha
 * channel using central differences.
 *
 * @param volume the voxel's volume
 * @param volumeStruct additional information about the passed volume
 * @param samplePos the sample's position in texture space
 */
vec3 calcGradient(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    return calcGradientA(volume, volumeStruct, samplePos);
}

