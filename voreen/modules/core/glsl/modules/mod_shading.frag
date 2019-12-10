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
 * This module contains all functions which can be used for shading.
 * The functions below, which implement a full shading model as for
 * instance Phong or Toon shading, are referenced by APPLY_SHADING
 * which is used in the raycaster fragment shaders.
 */

/// This struct contains all information about a light source.
struct LightSource {
    vec3 position_;        // light position in world space
    vec3 ambientColor_;    // ambient color (r,g,b)
    vec3 diffuseColor_;    // diffuse color (r,g,b)
    vec3 specularColor_;   // specular color (r,g,b)
    vec3 attenuation_;     // attenuation (constant, linear, quadratic)
};

// uniforms needed for shading
uniform float shininess_;       // material shininess parameter
uniform LightSource lightSource_;

/**
 * Returns attenuation based on the currently set opengl values.
 * Incorporates constant, linear and quadratic attenuation.
 *
 * @param d Distance to the light source.
 */
float getAttenuation(in float d) {
    float att = 1.0 / (lightSource_.attenuation_.x +
                       lightSource_.attenuation_.y * d +
                       lightSource_.attenuation_.z * d * d);
    return min(att, 1.0);
}


/**
 * Returns the ambient term, considering the user defined lighting
 * parameters.
 *
 * @param ka The ambient color to be used, which is fetched from the
 * transfer function.
 */
vec3 getAmbientTerm(in vec3 ka) {
    return ka * lightSource_.ambientColor_;
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
vec3 getDiffuseTerm(in vec3 kd, in vec3 N, in vec3 L) {
    float NdotL = max(dot(N, L), 0.0);
    return kd * lightSource_.diffuseColor_ * NdotL;
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
vec3 getLerpDiffuseTerm(in vec3 kd, in vec3 N, in vec3 L) {
    float alpha = 0.5;
    vec3 NV = mix(N, L, alpha);
    float NVdotL = max(dot(NV, L), 0.0);
    return kd * lightSource_.diffuseColor_ * NVdotL;
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
vec3 getSpecularTerm(in vec3 ks, in vec3 N, in vec3 L, in vec3 V, in float alpha) {
    vec3 H = normalize(V + L);
    float NdotH = pow(max(dot(N, H), 0.0), alpha);
    return ks * lightSource_.specularColor_ * NdotH;
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
vec3 phongShading(in vec3 n, in vec3 pos, in vec3 lPos, in vec3 cPos, in vec3 ka, in vec3 kd, in vec3 ks) {
    vec3 N = normalize(n);
    vec3 L = lPos - pos;              // using light position in volume physical space
    vec3 V = normalize(cPos - pos);  // using camera position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = length(L);
    L /= d;

    vec3 shadedColor = vec3(0.0);
    shadedColor += getAmbientTerm(ka);
    shadedColor += getDiffuseTerm(kd, N, L);
    shadedColor += getSpecularTerm(ks, N, L, V, shininess_);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d);
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
vec3 phongShadingDS(in vec3 n, in vec3 pos, in vec3 lPos, in vec3 cPos, in vec3 kd, in vec3 ks) {
    vec3 N = normalize(n);
    vec3 L = lPos - pos;              // using light position in volume physical space
    vec3 V = normalize(cPos - pos);  // using camera position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = length(L);
    L /= d;

    vec3 shadedColor = vec3(0.0);
    shadedColor += getDiffuseTerm(kd, N, L);
    shadedColor += getSpecularTerm(ks, N, L, V, shininess_);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d);
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
vec3 phongShadingS(in vec3 n, in vec3 pos, in vec3 lPos, in vec3 cPos, in vec3 ks) {
    vec3 N = normalize(n);
    vec3 L = lPos - pos;              // using light position in volume physical space
    vec3 V = normalize(cPos - pos);  // using camera position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = length(L);
    L /= d;

    vec3 shadedColor = vec3(0.0);
    shadedColor += getSpecularTerm(ks, N, L, V, shininess_);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d);
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
vec3 phongShadingDA(in vec3 n, in vec3 pos, in vec3 lPos, in vec3 cPos, in vec3 kd, in vec3 ka) {
    vec3 N = normalize(n);
    vec3 L = lPos - pos;    // using light position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = length(L);
    L /= d;

    vec3 shadedColor = vec3(0.0);
    shadedColor += getAmbientTerm(ka);
    shadedColor += getDiffuseTerm(kd, N, L);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d);
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
vec3 phongShadingD(in vec3 n, in vec3 pos, in vec3 lPos, in vec3 cPos, in vec3 kd) {
    vec3 N = normalize(n);
    vec3 L = lPos - pos;    // using light position in volume physical space

    // get light source distance for attenuation and normalize light vector
    float d = length(L);
    L /= d;

    vec3 shadedColor = vec3(0.0);
    shadedColor += getDiffuseTerm(kd, N, L);
    #ifdef PHONG_APPLY_ATTENUATION
        shadedColor *= getAttenuation(d);
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
vec3 toonShading(in vec3 n, in vec3 pos, in vec3 lPos, in vec3 cPos, in vec3 kd, in int numShades) {
    vec3 N = normalize(n);
    vec3 L = lPos - pos;    // using light position in volume physical space
    float NdotL = max(dot(N,L),0.0);

    kd *= lightSource_.diffuseColor_;

    for (int i=1; i <= numShades; i++) {
        if (NdotL <= float(i) / float(numShades)) {
            NdotL = float(i) / float(numShades);
            break;
        }
    }
    return vec3(kd.rgb * NdotL);
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
vec3 cookTorranceShading(in vec3 n, in vec3 pos, in vec3 lPos, in vec3 cPos, in vec3 ka, in vec3 kd, in vec3 ks) {
    //'Shininess' is mapped as combined parameter for roughness and specular reflection

    // Material roughness [0, 1]
    float m = 1.0 - (shininess_ / 128.0);

    //Specular and diffuse [0, 1]
    //s + d = 1
    float s = (shininess_ / 128.0);
    float d = 1.0 - s;

    //Um eine aehnliche Helligkeit, wie das Phong-Modell
    //zu erreichen, wird ueber den Korrekturfaktor die Helligkeit
    //der Lichtquelle angepasst.
    float korrekturFaktor = 2.0;

    float f = 1.0;     //Refelxionskoefficient bei senkrechter Einfallrichtung
    float dw = 1.0;    //Raumwinkel;

    vec3 N = normalize(n);
    vec3 L = normalize(lPos - pos);
    vec3 V = normalize(cPos - pos);
    vec3 H = normalize(V+L);

    //Raumwinkel abhaengig von der Lichtentfernung
    float distance = length(L);
    dw = distance * 1.0;

    float NH = dot(N,H);
    float NV = dot(N,V);
    float NL = max(0.0000001, dot(N,L));
    float VH = dot(V,H);

    //Fresnel-Term: Schlick-Approximation
    float F = f + (1.0 - f) * pow(1.0 - NV, 5.0);

    //Geometrische Abschwaechung
    float G1 = (2.0 * NH * NL) / VH;
    float G2 = (2.0 * NH * NV) / VH;
    float G = min(1.0, max(0.0, min(G1, G2)));

    //Beckmann microfacette distribution-function
    float mq = m*m + 1e-5;       // Add small value
    float NHq = NH*NH + 1e-5;    // to avoid division by zero

    float D1 = 1.0 / (4.0 * mq * NHq * NHq);
    float D2 = exp(-(1.0 - NHq) / (mq * NHq));
    float D = D1 * D2;

    //Diffuse reflection
    vec3 Rd = kd * lightSource_.diffuseColor_ * korrekturFaktor * NL;

    //Specular reflection
    vec3 Rs = F * D * G /(3.1415926 * NV * NL) * ks * lightSource_.specularColor_;

    //Final shaded color
    vec3 shadedColor = vec3(0.0);
    shadedColor += getAmbientTerm(ka);
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
vec3 orenNayarShading(in vec3 n, in vec3 pos, in vec3 lPos, in vec3 cPos, in vec3 ka, in vec3 kd) {
    // Material Roughness
    float m = shininess_ / 128.0;

    vec3 N = normalize(n);
    vec3 L = normalize(lPos - pos);
    vec3 V = normalize(cPos - pos);
    vec3 H = normalize(V + L);

    float VdotN = dot( V, N );
    float LdotN = dot( L, N );
    float NdotL = max(0.0000001, dot( N, L ));

    float cosPhi   = dot(V - N * VdotN, L - N * LdotN);
    float alpha = max( acos( VdotN ), acos( LdotN ) );
    float beta  = min( acos( VdotN ), acos( LdotN ) );

    float mq = m*m;

    float A = 1.0 - 0.5 * (mq / (mq + 0.57));
    float B = 0.45 * (mq / (mq + 0.09));

    vec3 result = vec3(0.0);
    result += (A + B * max( 0.0, cosPhi ) * sin(alpha) * tan(beta));

    vec3 shadedColor = vec3(0.0);
    shadedColor += getAmbientTerm(ka);
    shadedColor += kd * lightSource_.diffuseColor_ * NdotL * result;
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
vec3 wardShading(in vec3 n, in vec3 pos, in vec3 lPos, in vec3 cPos, in vec3 ka, in vec3 kd, in vec3 ks) {
    // Material roughness (isotropic)
    float alphaX = shininess_ / 128.0;
    //float alphaY = alphaX * 0.5;
    float alphaY = alphaX * 2.0;

    vec3 N = normalize(n);
    vec3 L = normalize(lPos - pos);
    vec3 V = normalize(cPos - pos);
    vec3 H = normalize(V + L);


    // Coordinateframe
    vec3 epsilon = vec3( 1.0, 0.0, 0.0 );
    vec3 X = normalize( cross( N, epsilon ) );  // Tangent
    vec3 Y = normalize( cross( N, X ));         // Orthogonal

    float VdotN = dot( V, N ); // Angle theta_r
    float LdotN = dot( L, N ); // Angle theta_i
    float HdotN = dot( H, N ); // Angle delta
    float HdotL = dot( H, L );
    float HdotX = dot( H, X );
    float HdotY = dot( H, Y );
    float NdotL = max(0.0000001, dot( N, L));

    alphaX += 0.0000001; // Add small value
    alphaY += 0.0000001; // to avoid division by zero

    // Diffuse reflection
    vec3 Rd = kd * lightSource_.diffuseColor_;

    // Specular Exponent
    float A = HdotX / alphaX;
          A *= A;

    float B = HdotY / alphaY;
          B *= B;

    float exponent = -2.0 * ( ( A + B ) / ( 1.0 + HdotN ) );

    // Specular Nenner
    float nenner  = 4.0 * 3.14159;
          nenner  *= alphaX;
          nenner  *= alphaY;
          nenner  *= sqrt( LdotN * VdotN );

    if (nenner < 0.0)
        nenner = min(0.0000001, nenner);
    else
        nenner = max(0.0000001, nenner);

    // Specular reflection
    vec3 Rs = ks * lightSource_.specularColor_ * ( exp( exponent ) / nenner );

    vec3 shadedColor = vec3(0.0);
    shadedColor += getAmbientTerm(ka);
    shadedColor += (NdotL * (Rs + Rd));
    return shadedColor;
}

