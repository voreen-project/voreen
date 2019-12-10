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

#version 330

#include "mod_geometryoutput.frag"

// structs for light source and materials
struct LightSource {
    vec4 position_;        // light position eye space (w=0 -> directional)
    vec3 ambientColor_;    // ambient color (r,g,b)
    vec3 diffuseColor_;    // diffuse color (r,g,b)
    vec3 specularColor_;   // specular color (r,g,b)
    //vec3 attenuation_;     // attenuation (constant, linear, quadratic)
};

struct Material {
    vec4 ambientColor_;
    vec4 diffuseColor_;
    vec4 specularColor_;
    float shininess_;
};


// input varyings
in vec3 frag_EyePosition;
in vec3 frag_EyeNormal;

// uniforms
uniform vec4 color_;

// lighting uniforms
uniform bool lightingEnabled_;
uniform LightSource lightSource_;
uniform Material material_;


// ---------- functions for phong lighting ---------------------------- //

/**
 * Returns attenuation based on the currently set opengl values.
 * Incorporates constant, linear and quadratic attenuation.
 *
 * @param d Distance to the light source.
 */
/*float getAttenuation(in float d) {
    float att = 1.0 / (lightSource_.attenuation_.x +
                       lightSource_.attenuation_.y * d +
                       lightSource_.attenuation_.z * d * d);
    return min(att, 1.0);
}*/


/**
 * Returns the ambient term, considering the user defined lighting
 * parameters.
 *
 * @param ka The ambient color to be used, which is fetched from the
 * transfer function.
 */
vec4 getAmbientTerm(in vec4 ka) {
    return ka * vec4(lightSource_.ambientColor_,1.0);
}


/**
 * Returns the diffuse term, considering the user defined lighting
 * parameters.
 *
 * @param kd The diffuse color to be used, which is fetched from the
 * transfer function.
 * @param N The surface normal used for lambert shading.
 * @param L The normalized light vector used for lambert shading.
 */
vec4 getDiffuseTerm(in vec4 kd, in vec3 N, in vec3 L) {
    float NdotL = max(dot(N, L), 0.0);
    return kd * vec4(lightSource_.diffuseColor_,1.0) * NdotL;
}

/**
 * Returns the specular term, considering the user defined lighting
 * parameters.
 *
 * @param ks The specular material color to be used.
 * @param N The surface normal used.
 * @param L The normalized light vector used.
 * @param V The viewing vector used.
 * @param alpha The shininess coefficient used.
 */
vec4 getSpecularTerm(in vec4 ks, in vec3 N, in vec3 L, in vec3 V, in float alpha) {
    vec3 H = normalize(V + L);
    float NdotH = pow(max(dot(N, H), 0.0), alpha);
    return ks * vec4(lightSource_.specularColor_,1.0) * NdotH;
}

/**
 * Calculates phong shading in eye space by considering the uniforms.
 * Attenuation is applied, if the w component of the light position is 1.
 *
 * @param n Normal (does not need to be normalized) in eye space.
 * @param pos Position in eye space
 * @param ma Material ambient color
 * @param md Material diffuse color
 * @param ms Material specular color
 */
vec4 phongShading(in vec3 n, in vec3 pos, in vec4 ma, in vec4 md, in vec4 ms) {

    vec3 N = normalize(n);
    vec3 L = lightSource_.position_.xyz - pos;

    // we are in eye space
    vec3 V = normalize(-pos);

    // get light source distance for attenuation and normalize light vector
    float d = length(L);
    L /= d;

    vec4 shadedColor = vec4(0.0);
    shadedColor += getAmbientTerm(ma);
    shadedColor += getDiffuseTerm(md, N, L);
    shadedColor += getSpecularTerm(ms, N, L, V, material_.shininess_);
    //shadedColor *= getAttenuation(d);

    return shadedColor;
}

void main() {
    vec4 color;
    if (!lightingEnabled_) {
        color = color_;
    }
    else {
        vec4 ma = material_.ambientColor_;
        vec4 md = material_.diffuseColor_;
        vec4 ms = material_.specularColor_;

        color = phongShading(frag_EyeNormal, frag_EyePosition, ma, md, ms);
    }
    outputFragment(color);
}

