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

#include "modules/mod_sampler2d.frag"
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"
#include "modules/mod_sliceviewer.frag"

uniform sampler2D entryPoints_;            // ray entry points
uniform sampler2D entryPointsDepth_;       // ray entry points depth
uniform TextureParameters entryParameters_;

// declare volume
#ifdef VOLUME_1_ACTIVE
uniform VolumeParameters volumeStruct1_;
uniform sampler3D volume1_;    // volume data set 1
uniform TransFuncParameters transferFunc1_;
uniform TF_SAMPLER_TYPE_1 transferFuncTex1_;
uniform float blendingFactor1_;
#endif

#ifdef VOLUME_2_ACTIVE
uniform VolumeParameters volumeStruct2_;
uniform sampler3D volume2_;    // volume data set 2
uniform TransFuncParameters transferFunc2_;
uniform TF_SAMPLER_TYPE_2 transferFuncTex2_;
uniform float blendingFactor2_;
#endif

#ifdef VOLUME_3_ACTIVE
uniform VolumeParameters volumeStruct3_;
uniform sampler3D volume3_;    // volume data set 3
uniform TransFuncParameters transferFunc3_;
uniform TF_SAMPLER_TYPE_3 transferFuncTex3_;
uniform float blendingFactor3_;
#endif

#ifdef VOLUME_4_ACTIVE
uniform VolumeParameters volumeStruct4_;
uniform sampler3D volume4_;    // volume data set 4
uniform TransFuncParameters transferFunc4_;
uniform TF_SAMPLER_TYPE_4 transferFuncTex4_;
uniform float blendingFactor4_;
#endif

vec4 lookup(in vec3 first) {
    vec4 result = vec4(0.0);

#ifdef VOLUME_1_ACTIVE
    vec3 samplePos1 = (volumeStruct1_.worldToTextureMatrix_ * vec4(first, 1.0)).xyz;

    if(inUnitCube(samplePos1)) {
        vec4 voxel = getVoxel(volume1_, volumeStruct1_, samplePos1);

        vec4 color;
        if (volumeStruct1_.numChannels_ == 1) {
            // Assuming Intensity only volume - get value via the transfer function.
            color = applyTF(transferFunc1_, transferFuncTex1_, voxel);
        }
        else if (volumeStruct1_.numChannels_ == 3) {
            // Assuming RGB volume - convert RGB to HSV and apply transfer function to hue.
            vec4 tfColor = applyTF(transferFunc1_, transferFuncTex1_, rgb2hsv(voxel.rgb).r);
            //if (applyColorModulation_)
                //voxel.rgb *= tfColor.rgb;
            color = vec4(voxel.rgb, tfColor.a);
        }
        else if (volumeStruct1_.numChannels_ == 4) {
            // Assuming RGBA volume - no transfer function lookup necessary
            color = voxel;
        }
        else {
            // something bad happened
            color = vec4(1.0);
        }

        result = BLENDINGMODE1(result, color, blendingFactor1_);
    }
#endif

#ifdef VOLUME_2_ACTIVE
    vec3 samplePos2 = (volumeStruct2_.worldToTextureMatrix_ * vec4(first, 1.0)).xyz;

    if(inUnitCube(samplePos2)) {
        vec4 voxel = getVoxel(volume2_, volumeStruct2_, samplePos2);

        vec4 color;
        if (volumeStruct2_.numChannels_ == 1) {
            // Assuming Intensity only volume - get value via the transfer function.
            color = applyTF(transferFunc2_, transferFuncTex2_, voxel);
        }
        else if (volumeStruct2_.numChannels_ == 3) {
            // Assuming RGB volume - convert RGB to HSV and apply transfer function to hue.
            vec4 tfColor = applyTF(transferFunc2_, transferFuncTex2_, rgb2hsv(voxel.rgb).r);
            //if (applyColorModulation_)
                //voxel.rgb *= tfColor.rgb;
            color = vec4(voxel.rgb, tfColor.a);
        }
        else if (volumeStruct2_.numChannels_ == 4) {
            // Assuming RGBA volume - no transfer function lookup necessary
            color = voxel;
        }
        else {
            // something bad happened
            color = vec4(1.0);
        }

        result = BLENDINGMODE2(result, color, blendingFactor2_);
    }
#endif

#ifdef VOLUME_3_ACTIVE
    vec3 samplePos3 = (volumeStruct3_.worldToTextureMatrix_ * vec4(first, 1.0)).xyz;

    if(inUnitCube(samplePos3)) {
        vec4 voxel = getVoxel(volume3_, volumeStruct3_, samplePos3);

        vec4 color;
        if (volumeStruct3_.numChannels_ == 1) {
            // Assuming Intensity only volume - get value via the transfer function.
            color = applyTF(transferFunc3_, transferFuncTex3_, voxel);
        }
        else if (volumeStruct3_.numChannels_ == 3) {
            // Assuming RGB volume - convert RGB to HSV and apply transfer function to hue.
            vec4 tfColor = applyTF(transferFunc3_, transferFuncTex3_, rgb2hsv(voxel.rgb).r);
            //if (applyColorModulation_)
                //voxel.rgb *= tfColor.rgb;
            color = vec4(voxel.rgb, tfColor.a);
        }
        else if (volumeStruct3_.numChannels_ == 4) {
            // Assuming RGBA volume - no transfer function lookup necessary
            color = voxel;
        }
        else {
            // something bad happened
            color = vec4(1.0);
        }

        result = BLENDINGMODE3(result, color, blendingFactor3_);
    }
#endif

#ifdef VOLUME_4_ACTIVE
    vec3 samplePos4 = (volumeStruct4_.worldToTextureMatrix_ * vec4(first, 1.0)).xyz;

    if(inUnitCube(samplePos4)) {
        vec4 voxel = getVoxel(volume4_, volumeStruct4_, samplePos4);

        vec4 color;
        if (volumeStruct4_.numChannels_ == 1) {
            // Assuming Intensity only volume - get value via the transfer function.
            color = applyTF(transferFunc4_, transferFuncTex4_, voxel);
        }
        else if (volumeStruct4_.numChannels_ == 3) {
            // Assuming RGB volume - convert RGB to HSV and apply transfer function to hue.
            vec4 tfColor = applyTF(transferFunc4_, transferFuncTex4_, rgb2hsv(voxel.rgb).r);
            //if (applyColorModulation_)
                //voxel.rgb *= tfColor.rgb;
            color = vec4(voxel.rgb, tfColor.a);
        }
        else if (volumeStruct4_.numChannels_ == 4) {
            // Assuming RGBA volume - no transfer function lookup necessary
            color = voxel;
        }
        else {
            // something bad happened
            color = vec4(1.0);
        }

        result = BLENDINGMODE4(result, color, blendingFactor4_);
    }
#endif

    return result;
}

void main() {
    vec4 frontPos = textureLookup2D(entryPoints_, entryParameters_, gl_FragCoord.xy);

    if(frontPos.a > 0.9) {
        vec4 result = lookup(frontPos.rgb);
        if(result.a == 0.0) discard;
        FragData0 = result;
    }else
        FragData0 = vec4(0.0);

    gl_FragDepth = textureLookup2D(entryPointsDepth_, entryParameters_, gl_FragCoord.xy).r;
}
