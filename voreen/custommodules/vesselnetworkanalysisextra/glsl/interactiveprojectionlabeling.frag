/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

uniform sampler2D tex_;

in vec2 frag_texcoord;

out vec4 color;

void main() {
    vec2 fragCoord = frag_texcoord.xy;

    float texValue = texture(tex_, fragCoord).x;
    float alpha = texture(tex_, fragCoord).y;

    vec3 col;
    if(alpha > 0.0) {
        col = vec3(texValue);
    } else {
        //Checkerboard pattern
        if(((int(gl_FragCoord.x)/40) & 1) - ((int(gl_FragCoord.y)/40) & 1) == 0) {
            col = vec3(0.3, 0.3, 0.3);
        } else {
            col = vec3(0.7, 0.7, 0.7);
        }
    }
    color = vec4(col, 1.0);
}
