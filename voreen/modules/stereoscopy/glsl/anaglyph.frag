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

uniform vec2 screenDimRCP_;

uniform sampler2D colorTex_;
uniform sampler2D depthTex_;

uniform int colorCode_;

in vec2 frag_texcoord;

void main() {
        vec2 fragCoord = gl_FragCoord.xy * screenDimRCP_;
        vec2 fragCoordLeft = vec2(fragCoord.x/2.f,fragCoord.y);
        vec2 fragCoordRight= vec2(fragCoordLeft.x+0.5f,fragCoord.y);
        vec4 fragColor = texture(colorTex_, fragCoordLeft);
        FragData0.r = fragColor.r;
        FragData0.a = fragColor.a;

        fragColor = texture(colorTex_, fragCoordRight);

        if(colorCode_ == 2){ //RED_CYAN
            FragData0.gb = fragColor.gb;
            FragData0.a = max(FragData0.a,fragColor.a);
        } else if(colorCode_ == 1){ //RED_BLUE
            FragData0.b = fragColor.b;
            FragData0.g = 0.0;
            FragData0.a = max(FragData0.a,fragColor.a);
        } else { //RED_GREEN
            FragData0.g = fragColor.g;
            FragData0.b = 0.0;
            FragData0.a = max(FragData0.a,fragColor.a);
        }

        float FragDepthLeft = texture(depthTex_, fragCoordLeft).z;
        float FragDepthRight = texture(depthTex_, fragCoordRight).z;
        gl_FragDepth = max(FragDepthLeft,FragDepthRight);
}
