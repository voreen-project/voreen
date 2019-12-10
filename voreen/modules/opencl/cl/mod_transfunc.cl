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

/**
 * This module contains all functions which can be used for performing
 * classifications of a voxel within a raycaster.
 * The functions below are referenced by RC_APPLY_CLASSIFICATION which
 * is used in the raycaster fragment shaders.
 */

//unused
struct TransFuncParameters {
    float3 domainLower_;
    float3 domainUpper_;
};

//normal TF:
float4 applyTF(image2d_t transfunc, sampler_t sampler, float intensity) {
    return read_imagef(transfunc, sampler, (float2)(intensity, 0.5f));
}

//pre-integrated TF:
float4 applyTFpi(image2d_t transfunc, sampler_t sampler, float intensity, float prev) {
    return read_imagef(transfunc, sampler, (float2)(intensity, prev));
}
