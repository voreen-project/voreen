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
 * This module contains all functions which can be used for shading in OpenCL.
 * The functions below, which implement a full shading model as for
 * instance Phong or Toon shading, are referenced by APPLY_SHADING
 * which is used in the openCL program.
 */

/// This struct contains all information about a light source.
typedef struct _LightSource {
    float4 position_;        // light position in world space
    float4 ambientColor_;    // ambient color (r,g,b)
    float4 diffuseColor_;    // diffuse color (r,g,b)
    float4 specularColor_;   // specular color (r,g,b)
    float4 attenuation_;     // attenuation (constant, linear, quadratic)
    float shininess_;
    float p1,p2,p3;
} LightSource;

/**
 * Returns attenuation based on the currently set opengl values.
 * Incorporates constant, linear and quadratic attenuation.
 *
 * @param d Distance to the light source.
 */
float getAttenuation(float d, __constant LightSource* lightSource_) {
    float att = 1.0f / ((*lightSource_).attenuation_.x +
                       (*lightSource_).attenuation_.y * d +
                       (*lightSource_).attenuation_.z * d * d);
    return min(att, 1.0f);
}


/**
 * Returns the ambient term, considering the user defined lighting
 * parameters.
 *
 * @param ka The ambient color to be used, which is fetched from the
 * transfer function.
 */
float3 getAmbientTerm(float3 ka, __constant LightSource* lightSource_) {
    return ka * (*lightSource_).ambientColor_.xyz;
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
float3 getDiffuseTerm(float3 kd, float3 N, float3 L, __constant LightSource* lightSource_) {
    float NdotL = max(dot(N, L), 0.0f);
    return kd * (*lightSource_).diffuseColor_.xyz * NdotL;
}


/**
 * This function implements the soft lighting technique described by
 * Josip Basic in the technote 'A cheap soft lighting for real-time 3D
 * environments.
 *
 * @param kd The diffuse color to be used, which is fetched from the
 * transfer function.
 * @param N The surface normal used for lambert shading.
 * @param L The normalized light vector used for lambert shading.
 */
float3 getLerpDiffuseTerm(float3 kd, float3 N, float3 L, __constant LightSource* lightSource_) {
    float alpha = 0.5f;
    float3 NV = mix(N, L, alpha);
    float NVdotL = max(dot(NV, L), 0.0f);
    return kd * (*lightSource_).diffuseColor_.xyz * NVdotL;
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
float3 getSpecularTerm(float3 ks, float3 N, float3 L, float3 V, __constant LightSource* lightSource_) {
    float3 H = normalize(V + L);
    float NdotH = pow(max(dot(N, H), 0.0f), (*lightSource_).shininess_);
    return ks * (*lightSource_).specularColor_.xyz * NdotH;
}

/**
 * Calculates phong shading by considering the user defined lighting parameters.
 * The front material's shininess parameter is used in the calculation of the specular term.
 * Attenuation is applied, if the symbol PHONG_APPLY_ATTENUATION is defined.
 *
 * @param n Normal (does not need to be normalized).
 * @param pos Position
 * @param lPos Light Position
 * @param cPos Camera Position
 * @param ka The ambient material color to be used.
 * @param kd The diffuse material color to be used.
 * @param ks The specular material color to be used.
 */
float3 phongShading(float3 n, float3 pos, float3 lPos, float3 cPos, float3 ka, float3 kd, float3 ks, __constant LightSource* lightSource_) {
    float3 N = fast_normalize(n);
    float3 L = lPos - pos;          // using light position in volume physical space
    float3 V = fast_normalize(cPos - pos);  // using camera position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = length(L);
    L /= d;

    float3 shadedColor = (float3)(0.0f);
    shadedColor += getAmbientTerm(ka, lightSource_);
    shadedColor += getDiffuseTerm(kd, N, L, lightSource_);
    shadedColor += getSpecularTerm(ks, N, L, V, lightSource_);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d, lightSource_);
    #endif
    return shadedColor;
}


/**
 * Calculates phong shading without considering the ambient term.
 * The front material's shininess parameter is used in the calculation of the specular term.
 * Attenuation is applied, if the symbol PHONG_APPLY_ATTENUATION is defined.
 *
 * @param n Normal (does not need to be normalized).
 * @param pos Position
 * @param lPos Light Position
 * @param cPos Camera Position
 * @param kd The diffuse material color to be used.
 * @param ks The specular material color to be used.
 */
float3 phongShadingDS(float3 n, float3 pos, float3 lPos, float3 cPos, float3 kd, float3 ks, __constant LightSource* lightSource_) {
    float3 N = fast_normalize(n);
    float3 L = lPos - pos;              // using light position in volume physical space
    float3 V = fast_normalize(cPos - pos);  // using camera position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = length(L);
    L /= d;

    float3 shadedColor = (float3)(0.0f);
    shadedColor += getDiffuseTerm(kd, N, L, lightSource_);
    shadedColor += getSpecularTerm(ks, N, L, V, lightSource_);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d, lightSource_);
    #endif
    return shadedColor;
}


/**
 * Calculates phong shading with only considering the specular term.
 * The front material's shininess parameter is used in the calculation of the specular term.
 * Attenuation is applied, if the symbol PHONG_APPLY_ATTENUATION is defined.
 *
 * @param n Normal (does not need to be normalized).
 * @param pos Position
 * @param lPos Light Position
 * @param cPos Camera Position
 * @param ks The specular material color to be used.
 */
float3 phongShadingS(float3 n, float3 pos, float3 lPos, float3 cPos, float3 ks, __constant LightSource* lightSource_) {
    float3 N = fast_normalize(n);
    float3 L = lPos - pos;              // using light position in volume physical space
    float3 V = fast_normalize(cPos - pos);  // using camera position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = length(L);
    L /= d;

    float3 shadedColor = (float3)(0.0f);
    shadedColor += getSpecularTerm(ks, N, L, V, lightSource_);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d, lightSource_);
    #endif
    return shadedColor;
}


/**
 * Calculates phong shading without considering the specular term.
 * Attenuation is applied, if the symbol PHONG_APPLY_ATTENUATION is defined.
 *
 * @param n Normal (does not need to be normalized).
 * @param pos Position
 * @param lPos Light Position
 * @param cPos Camera Position
 * @param kd The diffuse material color to be used.
 * @param ka The ambient material color to be used.
 */
float3 phongShadingDA(float3 n, float3 pos, float3 lPos, float3 cPos, float3 kd, float3 ka, __constant LightSource* lightSource_) {
    float3 N = fast_normalize(n);
    float3 L = lPos - pos;    // using light position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = fast_length(L);
    L /= d;

    float3 shadedColor = (float3)(0.0f);
    shadedColor += getAmbientTerm(ka, lightSource_);
    shadedColor += getDiffuseTerm(kd, N, L, lightSource_);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d, lightSource_);
    #endif
    return shadedColor;
}


/**
 * Calculates lambertian shading.
 * Attenuation is applied, if the symbol PHONG_APPLY_ATTENUATION is defined.
 *
 * @param n Normal (does not need to be normalized).
 * @param pos Position
 * @param lPos Light Position
 * @param cPos Camera Position
 * @param kd The diffuse material color to be used.
 */
float3 phongShadingD(float3 n, float3 pos, float3 lPos, float3 cPos, float3 kd, __constant LightSource* lightSource_) {
    float3 N = fast_normalize(n);
    float3 L = lPos - pos;    // using light position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = fast_length(L);
    L /= d;

    float3 shadedColor = (float3)(0.0f);
    shadedColor += getDiffuseTerm(kd, N, L, lightSource_);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d, lightSource_);
    #endif
    return shadedColor;
}

/**
 * Calculates a quantized toon shading.
 *
 * @param n Normal (does not need to be normalized).
 * @param pos Position
 * @param lPos Light Position
 * @param cPos Camera Position
 * @param numShades The number of different shadings.
 */
float3 toonShading(float3 n, float3 pos, float3 lPos, float3 cPos, float3 kd, int numShades, __constant LightSource* lightSource_) {
    float3 N = fast_normalize(n);
    float3 L = lPos - pos;    // using light position in volume physical space
    float NdotL = max(dot(N,L),0.0f);

    kd *= (*lightSource_).diffuseColor_.xyz;

    for (int i=1; i <= numShades; i++) {
        if (NdotL <= (float)(i) / (float)(numShades)) {
            NdotL = (float)(i) / (float)(numShades);
            break;
        }
    }
    return (float3)(kd.xyz * NdotL);
}

/**
 * Calculates Cook-Torrance shading by considering the user defined lighting
 * parameters.
 * The front material's shininess parameter is used as parameter for the material roughness.
 *
 * @param n Normal (does not need to be normalized).
 * @param pos Position
 * @param lPos Light Position
 * @param cPos Camera Position
 * @param ka The ambient material color to be used.
 * @param kd The diffuse material color to be used.
 * @param ks The specular material color to be used.
 */
float3 cookTorranceShading(float3 n, float3 pos, float3 lPos, float3 cPos, float3 ka, float3 kd, float3 ks, __constant LightSource* lightSource_) {
    //'Shininess' is mapped as combined parameter for roughness and specular reflection

    // Material roughness [0, 1]
    float m = 1.0f - ((*lightSource_).shininess_ / 128.0f);

    //Specular and diffuse [0, 1]
    //s + d = 1
    float s = ((*lightSource_).shininess_ / 128.0f);
    float d = 1.0f - s;

    //to achieve a similar intensity as the phong-model
    //a correction factor is used

    float correctionFactor = 2.0f;

    float f = 1.0f;     //reflection coefficient at vertical angle of entry
    float dw = 1.0f;    //solid angle;

    float3 N = fast_normalize(n);
    float3 L = fast_normalize(lPos - pos);
    float3 V = fast_normalize(cPos - pos);
    float3 H = fast_normalize(V+L);

    //solid angle depending on light distance
    float distance = fast_length(L);
    dw = distance * 1.0f;

    float NH = dot(N,H);
    float NV = dot(N,V);
    float NL = max(0.0000001f, dot(N,L));
    float VH = dot(V,H);

    //Fresnel-term: Schlick-approximation
    float F = f + (1.0f - f) * pow(1.0f - NV, 5.0f);

    //geometric attenuation
    float G1 = (2.0f * NH * NL) / VH;
    float G2 = (2.0f * NH * NV) / VH;
    float G = min(1.0f, max(0.0f, min(G1, G2)));

    //Beckmann microfacette distribution-function
    float mq = m*m + 1e-5;       // Add small value
    float NHq = NH*NH + 1e-5;    // to avoid division by zero

    float D1 = 1.0f / (4.0f * mq * NHq * NHq);
    float D2 = exp(-(1.0f - NHq) / (mq * NHq));
    float D = D1 * D2;

    //diffuse reflection
    float3 Rd = kd * (*lightSource_).diffuseColor_.xyz * correctionFactor * NL;

    //specular reflection
    float3 Rs = F * D * G /(3.1415926f * NV * NL) * ks * (*lightSource_).specularColor_.xyz;

    //final shaded color
    float3 shadedColor = (float3)(0.0f);
    shadedColor += getAmbientTerm(ka, lightSource_);
    shadedColor += (NL * dw * (s*Rs + d*Rd));
    return shadedColor;
}

/**
 * Calculates Oren-Nayar shading by considering the user defined lighting parameters.
 * The front material's shininess parameter is used as parameter for the material roughness.
 *
 * @param n Normal (does not need to be normalized).
 * @param pos Position
 * @param lPos Light Position
 * @param cPos Camera Position
 * @param ka The ambient material color to be used.
 * @param kd The diffuse material color to be used.
 */
float3 orenNayarShading(float3 n, float3 pos, float3 lPos, float3 cPos, float3 ka, float3 kd, __constant LightSource* lightSource_) {
    // Material Roughness
    float m = (*lightSource_).shininess_ / 128.0f;

    float3 N = fast_normalize(n);
    float3 L = fast_normalize(lPos - pos);
    float3 V = fast_normalize(cPos - pos);
    float3 H = fast_normalize(V + L);

    float VdotN = dot( V, N );
    float LdotN = dot( L, N );
    float NdotL = max(0.0000001f, dot( N, L ));

    float cosPhi = dot(V - N * VdotN, L - N * LdotN);
    float alpha = max( acos( VdotN ), acos( LdotN ) );
    float beta  = min( acos( VdotN ), acos( LdotN ) );

    float mq = m*m;

    float A = 1.0f - 0.5f * (mq / (mq + 0.57f));
    float B = 0.45f * (mq / (mq + 0.09f));

    float3 result = (float3)(0.0f);
    result += (A + B * max( 0.0f, cosPhi ) * sin(alpha) * tan(beta));

    float3 shadedColor = (float3)(0.0f);
    shadedColor += getAmbientTerm(ka, lightSource_);
    shadedColor += kd * (*lightSource_).diffuseColor_.xyz * NdotL * result;
    return shadedColor;
}

/**
 * Calculates Ward Isotropic shading by considering the user defined lighting parameters.
 * The front material's shininess parameter is used as parameter for the isotropic material roughness.
 *
 * @param n Normal (does not need to be normalized).
 * @param pos Position
 * @param lPos Light Position
 * @param cPos Camera Position
 * @param ka The ambient material color to be used.
 * @param kd The diffuse material color to be used.
 * @param ks The specular material color to be used.
 */
float3 wardShading(float3 n, float3 pos, float3 lPos, float3 cPos, float3 ka, float3 kd, float3 ks, __constant LightSource* lightSource_) {
    // Material roughness (isotropic)
    float alphaX = (*lightSource_).shininess_ / 128.0f;
    //float alphaY = alphaX * 0.5;
    float alphaY = alphaX * 2.0f;

    float3 N = fast_normalize(n);
    float3 L = fast_normalize(lPos - pos);
    float3 V = fast_normalize(cPos - pos);
    float3 H = fast_normalize(V + L);


    // Coordinateframe
    float3 epsilon = (float3)( 1.0f, 0.0f, 0.0f );
    float3 X = fast_normalize( cross( N, epsilon ) );  // Tangent
    float3 Y = fast_normalize( cross( N, X ));         // Orthogonal

    float VdotN = dot( V, N ); // Angle theta_r
    float LdotN = dot( L, N ); // Angle theta_i
    float HdotN = dot( H, N ); // Angle delta
    float HdotL = dot( H, L );
    float HdotX = dot( H, X );
    float HdotY = dot( H, Y );
    float NdotL = max(0.0000001f, dot( N, L));

    alphaX += 0.0000001f; // Add small value
    alphaY += 0.0000001f; // to avoid division by zero

    // Diffuse reflection
    float3 Rd = kd * (*lightSource_).diffuseColor_.xyz;

    // Specular Exponent
    float A = HdotX / alphaX;
          A *= A;

    float B = HdotY / alphaY;
          B *= B;

    float exponent = -2.0f * ( ( A + B ) / ( 1.0f + HdotN ) );

    // Specular Nenner
    float denominator  = 4.0f * 3.14159f;
          denominator  *= alphaX;
          denominator  *= alphaY;
          denominator  *= sqrt( LdotN * VdotN );

    if (denominator < 0.0f)
        denominator = min(0.0000001f, denominator);
    else
        denominator = max(0.0000001f, denominator);

    // Specular reflection
    float3 Rs = ks * (*lightSource_).specularColor_.xyz * ( exp( exponent ) / denominator );

    float3 shadedColor = (float3)(0.0f);
    shadedColor += getAmbientTerm(ka, lightSource_);
    shadedColor += (NdotL * (Rs + Rd));
    return shadedColor;
}
