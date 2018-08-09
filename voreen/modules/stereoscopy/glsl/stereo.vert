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

layout(location=0)in vec2 vert_position;
layout(location=1)in vec2 tex_coord;

out vec2 frag_texcoord;

uniform bool useQuadBuffer_;
uniform bool useLeftQuadBuffer_;

void main() {
    gl_Position = vec4(vert_position,0,1);
    if(!useQuadBuffer_) {
        frag_texcoord = tex_coord;
    } else {
        vec2 tmp_tex = vec2(tex_coord.x/2.0, tex_coord.y);
        if(useLeftQuadBuffer_)
            frag_texcoord = tmp_tex;
        else
            frag_texcoord = vec2(tmp_tex.x+0.5,tmp_tex.y);
    }
}
