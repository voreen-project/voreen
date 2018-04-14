/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

uniform sampler2D colorTexLeft_;
uniform sampler2D colorTexRight_;

uniform sampler2D depthTexLeft_;
uniform sampler2D depthTexRight_;

uniform bool useDepthTexLeft_;
uniform bool useDepthTexRight_;

in vec2 frag_texcoord;

void main() {
    /*gl_FragColor = vec4(frag_texcoord,1.0,1.0);
    return;*/ //DEBUG
    if(frag_texcoord.x < 0.5){//left part of screen
        vec2 fragCoord = frag_texcoord;
        fragCoord.x *= 2;
        vec4 fragColor = texture(colorTexLeft_, fragCoord);
        gl_FragColor = fragColor;
        if(useDepthTexLeft_)
            gl_FragDepth = texture(depthTexLeft_, fragCoord).z;

    } else {//right part of screen
        vec2 fragCoord = frag_texcoord;
        fragCoord.x = (fragCoord.x-0.5)*2;
        vec4 fragColor = texture(colorTexRight_, fragCoord);
        gl_FragColor = fragColor;
        if(useDepthTexRight_)
            gl_FragDepth = texture(depthTexRight_, fragCoord).z;
    }
}
