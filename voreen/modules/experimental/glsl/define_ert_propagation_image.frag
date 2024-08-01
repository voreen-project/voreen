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

#include "modules/mod_sampler2d.frag"
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"

uniform sampler2D entryPoints_;                      // ray entry points
uniform TextureParameters entryParameters_;
uniform sampler2D exitPoints_;                        // ray exit points
uniform TextureParameters exitParameters_;
uniform sampler2D ertExitPoints_;                       // the ert exit point texture
uniform TextureParameters ertExitParameters_;
uniform sampler2D ertTexture_;                       // our texture
uniform TextureParameters ertTextureParams_;

uniform vec2 projectedLightDir_;

//Return true/false if valid/invalid coordinates
bool checkValidCoordinates(in vec2 coords){
    return (coords.x >= 0.0 && coords.y >= 0.0 && coords.x <= 1.0 && coords.y <= 1.0);
}

/***
 * Defines correct propagation zone, and represented zone color
 ***/
vec4 propagateERT(vec4 first, vec4 last, vec2 p) {
    vec4 result = vec4(0.0, 0.0, 0.0, 1.0);
    vec2 p_light = projectedLightDir_ * screenDimRCP_;
    vec2 prevCoord = p + p_light;

    //Define last hit point based on distance to entry point
    vec3 direction = normalize(last.rgb - first.rgb);
    vec3 sb = first.rgb + last.a * direction;

    //Discard if not valid coordinate
     if(checkValidCoordinates(prevCoord)){
         vec4 prevLast = textureLookup2Dnormalized(ertTexture_, ertTextureParams_, prevCoord);
         result.rgb = min(sb, max(last.rgb, prevLast.rgb));
     }
     else{
        result.rgb = min(sb, last.rgb);
    }

    return result;
}

/***
 * The main method.
 ***/
void main() {
    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    vec4 frontPos = textureLookup2Dnormalized(entryPoints_, entryParameters_, p);
    vec4 backPos = textureLookup2Dnormalized(exitPoints_, exitParameters_, p);

    //determine whether the ray has to be casted
    if (frontPos.rgb == backPos.rgb) {
        //background need no raycasting
        FragData0 = frontPos;
    } else {
        //fragCoords are lying inside the boundingbox
        vec4 ertPos = textureLookup2Dnormalized(ertExitPoints_, ertExitParameters_, p);
        FragData0 = propagateERT(frontPos, ertPos, p);
    }
}
