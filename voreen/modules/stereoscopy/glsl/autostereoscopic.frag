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

uniform vec2 screenDimRCP_;
uniform sampler2D colorTex_;
uniform sampler2D depthTex_;
uniform int interleaveCode_;

in vec2 frag_texcoord;

void main() {
    vec2 fragCoord = gl_FragCoord.xy * screenDimRCP_;
    vec2 fragCoordLeft = vec2(fragCoord.x/2.f,fragCoord.y);
    vec2 fragCoordRight= vec2(fragCoordLeft.x+0.5f,fragCoord.y);
    if(interleaveCode_ == 0) { //vertical
        if(mod(int(gl_FragCoord.x), 2) == 0){ //left
            FragData0 = texture(colorTex_, fragCoordLeft);
            gl_FragDepth = texture(depthTex_, fragCoordLeft).z;
        } else { //right
            FragData0 = texture(colorTex_, fragCoordRight);
            gl_FragDepth = texture(depthTex_, fragCoordRight).z;
        }
    } else if(interleaveCode_ == 1) { //horizontal
        if(mod(int(gl_FragCoord.y), 2) == 0){ //left
            FragData0 = texture(colorTex_, fragCoordLeft);
            gl_FragDepth = texture(depthTex_, fragCoordLeft).z;
        } else { //right
            FragData0 = texture(colorTex_, fragCoordRight);
            gl_FragDepth = texture(depthTex_, fragCoordRight).z;
        }
    } else { //checker
        if(mod(int(gl_FragCoord.x), 2) + mod(int(gl_FragCoord.y), 2) != 1) {
           FragData0 = texture(colorTex_, fragCoordLeft);
           gl_FragDepth = texture(depthTex_, fragCoordLeft).z;
        } else { //right
            FragData0 = texture(colorTex_, fragCoordRight);
            gl_FragDepth = texture(depthTex_, fragCoordRight).z;
        }
    }
}
