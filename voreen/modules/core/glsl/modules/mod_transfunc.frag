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

/**
 * This module contains all functions which can be used for performing
 * classifications of a voxel within a raycaster.
 * The functions below are referenced by RC_APPLY_CLASSIFICATION which
 * is used in the raycaster fragment shaders.
 */

struct TransFuncParameters {
    vec3 domainLower_;
    vec3 domainUpper_;
};

float realWorldToTexture(TransFuncParameters tf, float v) {
    if(v <= tf.domainLower_.x)
        return 0.0;
    else if(v >= tf.domainUpper_.x)
        return 1.0;
    else if(tf.domainUpper_.x == tf.domainLower_.x)
        return 0.0;
    else
        return (v - tf.domainLower_.x) / (tf.domainUpper_.x - tf.domainLower_.x);
}

vec2 realWorldToTexture(TransFuncParameters tf, vec2 v) {
    float x = realWorldToTexture(tf, v.x);

    if(v.y <= tf.domainLower_.y)
        return vec2(x, 0.0);
    else if(v.y >= tf.domainUpper_.y)
        return vec2(x, 1.0);
    else if(tf.domainUpper_.y == tf.domainLower_.y)
        return vec2(x, 0.0);
    else
        return vec2(x, (v.y - tf.domainLower_.y) / (tf.domainUpper_.y - tf.domainLower_.y));
}

vec4 applyTF(TransFuncParameters transfunc, sampler1D tex, float intensity) {
    intensity = realWorldToTexture(transfunc, intensity);
    #if __VERSION__ > 130
        return texture(tex, intensity);
    #else
        return texture1D(tex, intensity);
    #endif
}

vec4 applyTF(TransFuncParameters transfunc, sampler1D tex, vec4 intensity) {
    intensity.a = realWorldToTexture(transfunc, intensity.a);
    #if __VERSION__ > 130
        return texture(tex, intensity.a);
    #else
        return texture1D(tex, intensity.a);
    #endif
}

vec4 applyTF(TransFuncParameters transfunc, sampler2D tex, float intensity, float gradientMagnitude) {
    #if __VERSION__ > 130
        return texture(tex, realWorldToTexture(transfunc, vec2(intensity, gradientMagnitude)));
    #else
        return texture2D(tex, realWorldToTexture(transfunc, vec2(intensity, gradientMagnitude)));
    #endif
}

vec4 applyTF(TransFuncParameters transfunc, sampler2D tex, vec4 intensityGradient) {
    #if __VERSION__ > 130
        return texture(tex, realWorldToTexture(transfunc, vec2(intensityGradient.a, length(intensityGradient.rgb))));
    #else
        return texture2D(tex, realWorldToTexture(transfunc, vec2(intensityGradient.a, length(intensityGradient.rgb))));
    #endif
}

//pre-integrated TF:
vec4 applyTFpi(TransFuncParameters transfunc, sampler2D tex, float intensity, float prev) {
    intensity = realWorldToTexture(transfunc, intensity);
    prev = realWorldToTexture(transfunc, prev);
    #if __VERSION__ > 130
        //return texture(tex, vec2(prev, intensity));
        return texture(tex, vec2(intensity, prev));
    #else
        return texture2D(tex, vec2(intensity, prev));
    #endif
}
