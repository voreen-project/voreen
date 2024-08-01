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

in vec4 frag_color;
in vec2 frag_normCoords;
in vec2 frag_halfLineDimensions;
in float frag_radiusLeftSquared;
in float frag_radiusRightSquared;
in vec2 frag_texPosSeed;

out vec4 out_color;


uniform sampler2D colorTex_;
uniform float textureZoom_;

float lengthSquared(vec2 v) {
    return dot(v,v);
}
float lengthSquared(vec3 v) {
    return dot(v,v);
}

void main() {
    float normX = abs(frag_normCoords.x);
    float normY = abs(frag_normCoords.y);
    vec2 xshift = vec2(1.0f,0.0f);
    float r = (0.5*(normX*normX)+0.5f);
    if(r<normY
//Cut the halos from the link area. This is necessary for textures with alpha
#ifdef USE_TEXTURES
        || lengthSquared((frag_normCoords+xshift)*frag_halfLineDimensions)<frag_radiusLeftSquared
        || lengthSquared((frag_normCoords-xshift)*frag_halfLineDimensions)<frag_radiusRightSquared
#endif
    ) {
        discard;
    }
#ifdef USE_TEXTURES
    vec2 textureCoords = ((frag_normCoords+frag_texPosSeed)*frag_halfLineDimensions)*textureZoom_;
    vec4 textureColor = texture(colorTex_, textureCoords);
    out_color = textureColor*frag_color;
#else
    float yNormOnRing = normY/r;
    float zNorm = sqrt(1.0f - yNormOnRing*yNormOnRing);
    float zreal = zNorm*r*(frag_halfLineDimensions.y);
    out_color = vec4(frag_color.xyz*zNorm,1);
    gl_FragDepth = 1-zreal*0.1f;
#endif
}
