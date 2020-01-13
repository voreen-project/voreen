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
__constant float SAMPLING_BASE_INTERVAL_RCP = 200.0;

/**
 * Performs regular DVR compositing. Expects the current voxels color
 * and the intermediate result. Returns the result after compositing.
 *
 */
float4 compositeDVR(float4 curResult, float4 color, float t, float samplingStepSize, float tDepth) {
    float4 result = curResult;

    // apply opacity correction to accomodate for variable sampling intervals
    color.w = 1.0 - pow(1.0f - color.w, samplingStepSize * SAMPLING_BASE_INTERVAL_RCP);

    result.xyz = result.xyz + (1.0 - result.w) * color.w * color.xyz;
    result.w = result.w + (1.0 -result.w) * color.w;
    // save first hit ray parameter for depth value calculation
    if (tDepth < 0.0)
        tDepth = t;
    return result;
}

/**
 * Performs regular DVR compositing without opacity correction (e.g. for pre-integrated TF).
 * Expects the current voxels color and the intermediate result.
 * Returns the result after compositing.
 */
float4 compositeDVRNoOpacityCorrection(float4 curResult, float4 color, float t, float tDepth) {
    float4 result = curResult;

    // apply opacity correction to accomodate for variable sampling intervals
    //color.w = 1.0 - pow(1.0f - color.w, samplingStepSize * SAMPLING_BASE_INTERVAL_RCP);

    result.xyz = result.xyz + (1.0 - result.w) * color.w * color.xyz;
    result.w = result.w + (1.0 -result.w) * color.w;
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
float4 compositeMIP(float4 curResult, float4 color, float t, float tDepth) {
    float4 result;
    if (color.w > curResult.w) {
        result = color;
        // save ray parameter for depth value calculation
        tDepth = t;
    }
    else result = curResult;
    return result;
}

/// Acutally performs the MIDA raycasting compositing
float4 compositeMIDAhelper(float4 curResult, float4 voxel, float4 color, float f_max_i, float t, float tDepth, float gamma) {
    float4 result = curResult;
    float delta_i = 0.0;
    if (voxel.w > f_max_i) {
        delta_i = voxel.w - f_max_i;
    }
    float beta_i = 1.0 - delta_i * (1.0 + gamma);

    result.xyz = beta_i * result.xyz + (1.0 - beta_i * result.w) * color.w * color.xyz;
    result.w = beta_i * result.w + (1.0 - beta_i * result.w) * color.w;

    if (tDepth < 0.0)
        tDepth = t;

    return result;

}

/**
 * Performs MIDA (maximum intensity difference accumulation) compositing as described in
 * the paper "Instant Volume Visualization using Maximum Intensity Difference Accumulation"
 * by Bruckner et al. published in Eurographics 2009
 */
float4 compositeMIDA(float4 curResult, float4 voxel, float4 color, float f_max_i, float t, float samplingStepSize, float tDepth, float gamma) {
    // apply opacity correction to accomodate for variable sampling intervals
    color.w = 1.0 - pow(1.0f - color.w, samplingStepSize * SAMPLING_BASE_INTERVAL_RCP);

    float4 result = curResult;

    if (gamma <= 0.0) {
        result = compositeMIDAhelper(result, voxel, color, f_max_i, t, tDepth, gamma);
    }
    else {
        float4 mipResult = compositeMIP(result, color, t, tDepth);
        float4 midaResult = compositeMIDAhelper(result, voxel, color, f_max_i, t, tDepth, 0.0);

        result = gamma * mipResult + (1.0 - gamma) * midaResult;
    }

    return result;
}

/**
 * Performs isosurface rendering compositing. Expects the current voxels color
 * and the intermediate result. Returns the result after compositing.
 *
 */
float4 compositeISO(float4 curResult, float4 color, float t, float tDepth, float isoValue) {
    float4 result = curResult;
    float epsilon = 0.02;
    if (color.w >= isoValue-epsilon && color.w <= isoValue+epsilon) {
        result = color;
        result.w = 1.0;
        // save ray parameter for depth value calculation
        tDepth = t;
    }
    return result;
}

/**
 * Performs first hit point compositing.
 */
float4 compositeFHP(float4 samplePos, float4 curResult, float t, float tDepth) {
    float4 result = curResult;
    // save first hit point
    if (result.x == 0.0f && result.y == 0.0f && result.z == 0.0f) {
        result.xyz = samplePos.xyz;
        result.w = 1.0f;
        // save first hit ray parameter for depth value calculation
        if (tDepth < 0.0f)
            tDepth = t;
    }
    return result;
}

/**
 * Performs first hit normals (gradients) compositing.
 */
float4 compositeFHN(float3 gradient, float4 curResult, float t, float tDepth) {
    float4 result = curResult;
    // save first hit normal
    if (result.x == 0.0f && result.y == 0.0f && result.z == 0.0f) {
        result.xyz = normalize(gradient) / 2.0f + 0.5f;
        result.w = 1.0f;
        // save first hit ray parameter for depth value calculation
        if (tDepth < 0.0f)
            tDepth = t;
    }
    return result;
}

/**
 * Performs first hit texture coordinate compositing.
 */
float4 compositeFHT(float3 texCoords, float4 curResult, float t, float tDepth) {
    float4 result = curResult;
    // save first hit normal
    if (result.x == 0.0f && result.y == 0.0f && result.z == 0.0f) {
        result.xyz = texCoords;
        result.w = 1.0f;
        // save first hit ray parameter for depth value calculation
        if (tDepth < 0.0f)
            tDepth = t;
    }
    return result;
}
