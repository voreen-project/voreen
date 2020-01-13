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

/**
 * This module contains all functions which can be used for compositing
 * voxels within a raycaster.
 * The functions are referenced by RC_APPLY_COMPOSITING as used in the
 * raycaster fragment shaders.
 */

/**
 * Reciprocal value of the reference sampling interval used
 * for the opacity correction that is necessary to compensate
 * variable sampling rates.
 *
 * See Engel et. al.: "Real-Time Volume Graphics" - Ch 9.1.3
 */
const float SAMPLING_BASE_INTERVAL_RCP = 200.0;

/**
 * Performs regular DVR compositing. Expects the current voxels color
 * and the intermediate result. Returns the result after compositing.
 *
 */
vec4 compositeDVR(in vec4 curResult, in vec4 color, in float t, in float samplingStepSize, inout float tDepth) {
    vec4 result = curResult;

    // apply opacity correction to accomodate for variable sampling intervals
    color.a = 1.0 - pow(1.0 - color.a, samplingStepSize * SAMPLING_BASE_INTERVAL_RCP);

    result.rgb = result.rgb + (1.0 - result.a) * color.a * color.rgb;
    result.a = result.a + (1.0 -result.a) * color.a;
    // save first hit ray parameter for depth value calculation
    if (tDepth < 0.0)
        tDepth = t;
    return result;
}

/**
 * Performs regular DVR compositing, but without applying the opacity correction (e.g. for pre-integrated TF).
 * Expects the current voxels color and the intermediate result. Returns the result after compositing.
 *
 */
vec4 compositeDVRNoOpacityCorrection(in vec4 curResult, in vec4 color, in float t, inout float tDepth) {
    vec4 result = curResult;

    // apply opacity correction to accomodate for variable sampling intervals
    //color.a = 1.0 - pow(1.0 - color.a, samplingStepSize * SAMPLING_BASE_INTERVAL_RCP);

    result.rgb = result.rgb + (1.0 - result.a) * color.a * color.rgb;
    result.a = result.a + (1.0 -result.a) * color.a;
    // save first hit ray parameter for depth value calculation
    if (tDepth < 0.0)
        tDepth = t;
    return result;
}

/**
 * Performs MIP (maximum intensity projection) compositing. Expects the current
 * voxels color and the intermediate result. Returns the result after compositing.
 *
 */
vec4 compositeMIP(in vec4 curResult, float voxel, in vec4 color, in float t, inout float tDepth, inout float maxIntensity) {
    vec4 result;
    if (voxel > maxIntensity) {
        result = color;
        maxIntensity = voxel;
        // save ray parameter for depth value calculation
        tDepth = t;
    }
    else result = curResult;
    return result;
}

vec4 compositeMOP(in vec4 curResult, in vec4 color, in float t, inout float tDepth) {
    vec4 result;
    if (color.a > curResult.a) {
        result = color;
        // save ray parameter for depth value calculation
        tDepth = t;
    }
    else result = curResult;
    return result;
}

/// Acutally performs the MIDA raycasting compositing
vec4 compositeMIDAhelper(in vec4 curResult, in float voxel, in vec4 color, inout float maxIntensity, in float t, inout float tDepth, in float gamma) {
    vec4 result = curResult;
    float delta_i = 0.0;
    if (voxel > maxIntensity) {
        delta_i = voxel - maxIntensity;
    }
    float beta_i = 1.0 - delta_i * (1.0 + gamma);

    result.rgb = beta_i * result.rgb + (1.0 - beta_i * result.a) * color.a * color.rgb;
    result.a = beta_i * result.a + (1.0 - beta_i * result.a) * color.a;

    if (tDepth < 0.0)
        tDepth = t;

    return result;

}

/**
 * Performs MIDA (maximum intensity difference accumulation) compositing as described in
 * the paper "Instant Volume Visualization using Maximum Intensity Difference Accumulation"
 * by Bruckner et al. published in Eurographics 2009
 */
vec4 compositeMIDA(in vec4 curResult, in float voxel, in vec4 color, inout float maxIntensity, in float t, in float samplingStepSize, inout float tDepth, in float gamma) {
    // apply opacity correction to accomodate for variable sampling intervals
    color.a = 1.0 - pow(1.0 - color.a, samplingStepSize * SAMPLING_BASE_INTERVAL_RCP);

    vec4 result = curResult;

    if (gamma <= 0.0) {
        result = compositeMIDAhelper(result, voxel, color, maxIntensity, t, tDepth, gamma);
    }
    else {
        vec4 mipResult = compositeMIP(result, voxel, color, t, tDepth, maxIntensity);
        vec4 midaResult = compositeMIDAhelper(result, voxel, color, maxIntensity, t, tDepth, 0.0);

        result = gamma * mipResult + (1.0 - gamma) * midaResult;
    }

    return result;
}

/**
 * Performs isosurface rendering compositing. Expects the current voxels color
 * and the intermediate result. Returns the result after compositing.
 *
 */
vec4 compositeISO(in vec4 curResult, in vec4 color, in float t, inout float tDepth, in float isoValue) {
    vec4 result = curResult;
    float epsilon = 0.02;
    if (color.a >= isoValue-epsilon && color.a <= isoValue+epsilon) {
        result = color;
        result.a = 1.0;
        // save ray parameter for depth value calculation
        tDepth = t;
    }
    return result;
}

/**
 * Performs first hit point compositing.
 */
vec4 compositeFHP(in vec3 samplePos, in vec4 curResult, in float t, inout float tDepth) {
    vec4 result = curResult;
    // save first hit point
    if (result.xyz == vec3(0.0)) {
        result.xyz = samplePos;
        result.a = 1.0;
        // save first hit ray parameter for depth value calculation
        if (tDepth < 0.0)
            tDepth = t;
    }
    return result;
}

/**
 * Performs first hit normals (gradients) compositing.
 */
vec4 compositeFHN(in vec3 gradient, in vec4 curResult, in float t, inout float tDepth) {
    vec4 result = curResult;
    // save first hit normal
    if (result.xyz == vec3(0.0)) {
        result.xyz = normalize(gradient) / 2.0 + 0.5;
        result.a = 1.0;
        // save first hit ray parameter for depth value calculation
        if (tDepth < 0.0)
            tDepth = t;
    }
    return result;
}

/**
 * Performs first hit texture coordinate compositing.
 */
vec4 compositeFHT(in vec3 texCoords, in vec4 curResult, in float t, inout float tDepth) {
    vec4 result = curResult;
    // save first hit normal
    if (result.xyz == vec3(0.0)) {
        result.xyz = texCoords;
        result.a = 1.0;
        // save first hit ray parameter for depth value calculation
        if (tDepth < 0.0)
            tDepth = t;
    }
    return result;
}
