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

#define PLANE_SAGITTAL  0
#define PLANE_CORONAL   1
#define PLANE_AXIAL     2

uniform vec3 flowDimensions_;
uniform vec2 textureDimensions_;
uniform int sliceNo_;
uniform vec2 textureScaling_; /** Scaling for texture coordinates (usually (1.0, 1.0), but smaller for viewports...*/
uniform vec2 randPos_;
uniform float spotScale_;     /** size of spot relative to texture width (e.g. 0.05 = 1/20) */
uniform float minValue_;
uniform float maxValue_;
uniform float maxMagnitude_;
uniform sampler3D vectorField3D_;
uniform bool useSpotBending_;

varying float discard_;

#if PLANE == PLANE_AXIAL
int i = 0, j = 1, k = 2;
#elif PLANE == PLANE_CORONAL
int i = 0, j = 2, k = 1;
#else   // default is permutation for sagittal plane section through volume
int i = 2, j = 1, k = 0;
#endif

/*
vec2 lookupVectorField2D(const vec2 pos)
{
    vec2 p = (pos / flowDimensions_.xz);    // map to [0.0, 1.0]
    p = ((p * 2.0) - vec2(1.0));            // map to [-1.0, 1.0]
    float radius = length(p);
    p = normalize(p);

    if ((0.0 <= radius) && (radius <= 1.0))
        return (vec2(-p.y, p.x) * (radius * 30.0));    // within [-30.0, 30.0]

    return vec2(0.0, 0.0);
}

vec3 lookupVectorField(const vec3 pos)
{
    vec2 v = lookupVectorField2D(pos.xz);
    return vec3(v.x, 0.0, v.y);
}
*/
vec3 lookupVectorField(const vec3 r)
{
    vec3 texCoord = clamp(r / flowDimensions_, 0.0, 1.0);
    vec3 v = texture3D(vectorField3D_, texCoord).xyz;
    if (v != vec3(0.0)) {
        float range = maxValue_ - minValue_;
        v = ((v * range) + minValue_);
    }
    return v;
}

vec3 viewportPosToFlowPos(const vec2 vp) {
    vec3 p = vec3(0.0);
    vec2 c = mod(vp, textureDimensions_);
    // map to [0, 1]
    //
    p[i] = (c.x / textureDimensions_.x);
    p[j] = (c.y / textureDimensions_.y);
    p[k] = float(sliceNo_) / flowDimensions_[k];
    return (p * flowDimensions_);
}

vec2 flowPosToViewportPos(const vec3 fp) {
    vec3 aux = (fp / flowDimensions_);  // map to [0, 1]
    return clamp((vec2(aux[i], aux[j]) * textureDimensions_), vec2(0.0), textureDimensions_);
}

bool isInsideFlowVolume(const vec3 r) {
    if (any(lessThan(r, vec3(0.0))) == true)
        return false;
    if (any(greaterThanEqual(r, flowDimensions_)) == true)
        return false;
    return true;
}

vec2 safeNormalize(const vec2 input) {
    if (input == vec2(0.0))
        return vec2(0.0);
    return normalize(input);
}

vec3 safeNormalize(const vec3 input) {
    if (input == vec3(0.0))
        return vec3(0.0);
    return normalize(input);
}

/**
 * Integrate a streamline using a 4th order Runge-Kutta method and
 * return the value at the given distance starting from position r0
 * using stepwidth h.
 */
vec3 computeStreamline3DRungeKutta(const vec3 r0, const float distance, const float h)
{
    if ((distance == 0.0) || (abs(distance) < h) || (h <= 0.0))
        return r0;

    int N = int(ceil(abs(distance) / h)); // number of partitions
    vec3 r = r0;
    vec3 k1 = vec3(0.0);
    vec3 k2 = vec3(0.0);
    vec3 k3 = vec3(0.0);
    vec3 k4 = vec3(0.0);

    vec3 prevR = r;
    bool lookupPos = (distance > 0.0);  // integrate along the streamline in positive direction?

    for (int i = 0; i < N; ++i) {
        vec3 v = lookupVectorField(r);
        if ((v == vec3(0.0)) || (isInsideFlowVolume(r) == false)) {
            discard_ = 1.0;
            return r;
        }

        k1 = safeNormalize(v) * h;
        if (lookupPos == true) {
            k2 = safeNormalize(lookupVectorField(r + (k1 / 2.0))) * h;
            k3 = safeNormalize(lookupVectorField(r + (k2 / 2.0))) * h;
            k4 = safeNormalize(lookupVectorField(r + k3)) * h;
            r += ((k1 / 6.0) + (k2 / 3.0) + (k3 / 3.0) + (k4 / 6.0));
        } else {
            k2 = safeNormalize(lookupVectorField(r - (k1 / 2.0))) * h;
            k3 = safeNormalize(lookupVectorField(r - (k2 / 2.0))) * h;
            k4 = safeNormalize(lookupVectorField(r - k3)) * h;
            r -= ((k1 / 6.0) + (k2 / 3.0) + (k3 / 3.0) + (k4 / 6.0));
        }

        if ((r == prevR) || (isInsideFlowVolume(r) == false)) {
            discard_ = 1.0;
            return prevR;
        } else
            prevR = r;
    }   // for

    return r;
}

vec3 computeStreamline(const vec3 r0, const float distance, const float h) {
    return computeStreamline3DRungeKutta(r0, distance, h);
}

void main() {
    float meshWidth = textureDimensions_.x * spotScale_;   // scaling factor for the mesh

    // Center the vertices around the center of the spot mesh and
    // map to range within [-0.5, 0.5] x [-0.5, 0.5]
    //
    vec2 vertexPos = (gl_Vertex.xy - vec2(0.5));
    vec2 p = randPos_ * textureDimensions_;    // ensure valid random position
    vec2 r = (vertexPos * meshWidth) + p;

    vec3 r0 = viewportPosToFlowPos(p);
    vec3 v3 = lookupVectorField(r0);
    vec2 v = vec2(v3[i], v3[j]);
    float magnitude = length(v);
    v = safeNormalize(v);

    if (v == vec2(0.0)) {
        discard_ = 1.0;
    } else {
        vec2 rake = vec2(-v.y, v.x);

        float maxWidth = (1.0 / spotScale_) - 1.0;
        float w = (1.0 + (maxWidth * (magnitude / maxMagnitude_))) * meshWidth;
        //float w = meshWidth * (1.0 + magnitude);
        float h = (meshWidth * meshWidth) / w;

        float len = vertexPos.x * w;
        float height = vertexPos.y * h;

        if (useSpotBending_ == true){
            vec3 convertedLength = viewportPosToFlowPos(vec2(w, 0.0) * vertexPos);
            vec3 rRake = viewportPosToFlowPos(p + (rake * height));

            if (lookupVectorField(rRake) != vec3(0.0)) {
                vec3 r3 = computeStreamline(rRake, convertedLength[i], 0.5);
                r = flowPosToViewportPos(r3);
            } else {
                discard_ = 1.0;
                r = (p + (v * meshWidth * vertexPos.x) + (rake * meshWidth * vertexPos.y));
            }
        } else
            r = (p + (v * len) + (rake * height));
        discard_ = 0.0;
    }

    // finally ransform the mesh to unified screen coordinates [-1.0, 1.0]
    //
    gl_Position = (gl_ModelViewProjectionMatrix * vec4(r, 0.0, 1.0));

    // texture coordinate for output texture
    //
    gl_TexCoord[1] = (gl_Position + vec4(1.0)) * 0.5;   // map to [0.0, 1.0]
    gl_TexCoord[1].st *= clamp(textureScaling_, 0.0, 1.0);

    // texture coordinates for the spot texture itself. gl_Vertex is within [0.0, 1.0]
    //
    gl_TexCoord[2] = gl_Vertex;
}
