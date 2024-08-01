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

#version 330

#include "modules/mod_transfunc.frag"

in vec2 frag_texcoord;

uniform sampler2D controlTex_;
uniform sampler2D intensityTex_;
uniform sampler2D controlInteractiveTex_;
uniform sampler2D intensityInteractiveTex_;
uniform float rwmOffset_;
uniform float rwmScale_;

// transfer functions
uniform sampler1D transFuncTex_;
uniform TransFuncParameters transFuncParams_;
uniform vec3 channelShift_;
#if NUM_CHANNELS > 1
uniform TF_SAMPLER_TYPE transFuncTex2_;
uniform TransFuncParameters transFuncParams2_;
uniform vec3 channelShift2_;
#endif
#if NUM_CHANNELS > 2
uniform TF_SAMPLER_TYPE transFuncTex3_;
uniform TransFuncParameters transFuncParams3_;
uniform vec3 channelShift3_;
#endif
#if NUM_CHANNELS > 3
uniform TF_SAMPLER_TYPE transFuncTex4_;
uniform TransFuncParameters transFuncParams4_;
uniform vec3 channelShift4_;
#endif

void main() {
    vec4 textureVal;
    if(texture(controlTex_, frag_texcoord).r == 0.0) {
        if(texture(controlInteractiveTex_, frag_texcoord).r == 0.0) {
            discard;
        } else {
            textureVal = texture(intensityInteractiveTex_, frag_texcoord);
        }
    } else {
        textureVal = texture(intensityTex_, frag_texcoord);
    }
    vec4 intensity = rwmOffset_ + rwmScale_ * textureVal;

    vec4 result;
#if NUM_CHANNELS == 1 // assuming red-texture
    result = applyTF(transFuncParams_, transFuncTex_, intensity.r);
#elif NUM_CHANNELS == 2 // assuming red-green texture
    vec4 channel1 = applyTF(transFuncParams_, transFuncTex_, intensity.r);
    vec4 channel2 = applyTF(transFuncParams2_, transFuncTex2_, intensity.g);
    //result = clamp(channel1 + channel2, vec4(0.0), vec4(1.0));
    result.rgb = clamp(channel1.rgb*channel1.a + channel2.rgb*channel2.a, vec3(0.0), vec3(1.0));
    result.a = max(channel1.a, channel2.a);
#elif NUM_CHANNELS == 3 // assuming RGB texture
    vec4 channel1 = applyTF(transFuncParams_, transFuncTex_, intensity.r);
    vec4 channel2 = applyTF(transFuncParams2_, transFuncTex2_, intensity.g);
    vec4 channel3 = applyTF(transFuncParams3_, transFuncTex3_, intensity.b);
    //result = clamp(channel1 + channel2 + channel3, vec4(0.0), vec4(1.0));
    result.rgb = clamp(channel1.rgb*channel1.a + channel2.rgb*channel2.a + channel3.rgb*channel3.a, vec3(0.0), vec3(1.0));
    result.a = max(channel1.a, max(channel2.a, channel3.a));
#elif NUM_CHANNELS == 4 // assuming RGBA texture
    vec4 channel1 = applyTF(transFuncParams_, transFuncTex_, intensity.r);
    vec4 channel2 = applyTF(transFuncParams2_, transFuncTex2_, intensity.g);
    vec4 channel3 = applyTF(transFuncParams3_, transFuncTex3_, intensity.b);
    vec4 channel4 = applyTF(transFuncParams4_, transFuncTex4_, intensity.a);
    //result = channel1 + channel2 + channel3 + channel4;
    result.rgb = clamp(channel1.rgb*channel1.a + channel2.rgb*channel2.a + channel3.rgb*channel3.a + channel4.rgb*channel4.a, vec3(0.0), vec3(1.0));
    result.a = max(channel1.a, max(channel2.a, max(channel3.a, channel4.a)));
#else // more than four channels => unknown texture type => only use alpha-value
    result = applyTF(transFuncParams_, transFuncTex_, intensity.a);
#endif

    vec4 fragColor = result;
    FragData0 = fragColor;

    if (result.a > 0.0)
        gl_FragDepth = gl_FragCoord.z;
    else
        gl_FragDepth = 1.0;
}
