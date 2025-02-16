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

#define VERTEX_FLAG_NONE        0u
#define VERTEX_FLAG_SELECTED   (1u<<0)
#define VERTEX_FLAG_MOUSE_OVER (1u<<1)
#define VERTEX_FLAG_ON_PATH    (1u<<2)
#define VERTEX_FLAG_ORPHAN     (1u<<3)

layout(location = 0) in vec3 in_vertex;
layout(location = 1) in uint in_flags;

uniform mat4 viewMatrix;
uniform float radius_;
uniform vec4 colorOrdinary_;
uniform vec4 colorSelected_;
uniform vec4 colorMouseOver_;
uniform vec4 colorOnPath_;
uniform vec4 colorOrphan_;

out vec4 geom_position;
out vec4 geom_color;
out float geom_radius;

void main() {
    vec4 pos = viewMatrix * vec4(in_vertex, 1.0f);
    gl_Position = pos;
    geom_position = pos;
    float r = radius_;
    if((in_flags&VERTEX_FLAG_SELECTED)!=0u) {
        geom_color = colorSelected_;
        r *= 1.5f;
    } else if((in_flags&VERTEX_FLAG_MOUSE_OVER)!=0u) {
        geom_color = colorMouseOver_;
        r *= 1.5f;
    } else if((in_flags&VERTEX_FLAG_ORPHAN)!=0u) {
        geom_color = colorOrphan_;
    } else if((in_flags&VERTEX_FLAG_ON_PATH)!=0u) {
        geom_color = colorOnPath_;
    } else {
        geom_color = colorOrdinary_;
    }
    geom_radius = r;
}
