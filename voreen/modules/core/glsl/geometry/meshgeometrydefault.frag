/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "modules/mod_sampler2d.frag"
#include "modules/mod_shading.frag"
#include "mod_geometryoutput.frag"

in VertexData {
  vec4 posEye;
  vec4 color;
#if defined(USE_NORMAL)
  vec3 normal;
#endif
#if defined(USE_TEXCOORD)
  vec2  texCoord;
  flat ivec2 texIndex;
#endif
} fragData;

uniform vec3 lightPositionEye_;
uniform bool enableLighting_;


#if defined(USE_TEXCOORD)
uniform bool useTextureArray_;

uniform sampler2DArray textures_; //Used if useTextureArray_ == true.
uniform sampler2D texture_; //Used if useTextureArray_ == true.
#endif

void main() {
#if defined(USE_TEXCOORD)
    // Array-texturing was enabled, but no texture coordinates were defined -> discard fragment.
    if(useTextureArray_ && fragData.texIndex.y == 0) {
        discard;
    }
#endif

#if defined(USE_NORMAL)
    vec3 norm = normalize(fragData.normal);
#endif

#if defined(USE_TEXCOORD)
    // If this check turns out to be a performance problem,
    // one might consider changing it back to a if/def block
    // (see trianglemesh.frag).
    vec4 ambCol = vec4(0.0);
    vec4 difCol = vec4(0.0);
    vec4 spcCol = vec4(0.0);
    if(useTextureArray_) {
        bool ambColSet = false;
        bool difColSet = false;
        bool spcColSet = false;

        int usedTextures = fragData.texIndex.x;
        int textureApplications = fragData.texIndex.y;
        int counter = 0;

        while(usedTextures != 0 && counter < 4) {
            // calculate int where the lowest bit that was set in the original int is deleted
            int tmp = int(usedTextures & int(usedTextures - 1));
            // calculate int containing only the previously deleted bit, take log2 to calculate layer in texture stack
            float texNum = log2(float(usedTextures - tmp));
            //float texNum = float(hashPos[(usedTextures - tmp) % 37]);
            // look up corresponding color in the stack
            vec4 texCol = texture(textures_, vec3(fragData.texCoord, texNum));
            // make sure that tex color is visible if set
            if(texCol.rgb != vec3(0.0) && texCol.a == 0.0)
                texCol.a = 1.0;

            // find out what this texture does by looking up the 4 corresponding texture function bits
            int textureFunction = int(textureApplications >> (counter * 4)) % 16;
            // ambient color
            if(int(textureFunction & 1) > 0) {
                ambCol = combineColorAndTexture(fragData.color, texCol);
                ambColSet = true;
            }
            // diffuse color
            if(int(textureFunction & 2) > 0) {
                difCol = combineColorAndTexture(fragData.color, texCol);
                difColSet = true;
            }
            // specular color
            if(int(textureFunction & 4) > 0) {
                spcCol = texCol;
                spcColSet = true;
            }

            // normal map
            // TODO differ from / support simpler one-channel bump-mapping
            if(int(textureFunction & 8) > 0) {
                // rotate normal read from geometry map to local coordinate system defined by geometry normal
                vec3 localNorm = normalize(2.0 * texCol.xyz - 1.0);
                vec3 tang = normalize(cross(norm, vec3(0.0, 1.0, 0.0)));
                vec3 btang = normalize(cross(norm, tang));
                //mat3 rot = transpose(mat3(tang, btang, fragData.normal));
                mat3 rot = mat3(tang, btang, norm);
                norm = rot * localNorm;
            }

            usedTextures = tmp;
            counter++;
        }

        // use same color for diffuse and ambient term if only one was set
        if(!ambColSet && difColSet)
            ambCol = difCol;
        else if(!difColSet && ambColSet)
            difCol = ambCol;

        // use default specular color white if not set
        if(!spcColSet)
            spcCol = vec4(1.0);
    } else {
        vec4 texCol = vec4(texture(texture_, fragData.texCoord).xyz, 1.0);
        ambCol = combineColorAndTexture(fragData.color, texCol);
        difCol = combineColorAndTexture(fragData.color, texCol);
        spcCol = vec4(1.0);
    }
#else
    vec4 ambCol = fragData.color;
    vec4 difCol = fragData.color;
    vec4 spcCol = vec4(1.0);
#endif

    vec4 result = ambCol;

#if defined(USE_NORMAL)
    if(enableLighting_)
        result.xyz = phongShading(norm, fragData.posEye.xyz, lightPositionEye_, vec3(0.0), ambCol.rgb, difCol.rgb, spcCol.rgb);
#endif

    outputFragment(result);
}
