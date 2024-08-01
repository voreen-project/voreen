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

layout(location = 0) out vec4 color_;

uniform mat4 P_;

in VertexData{
    vec2 coord;
    float z;
    vec3 zdir;
    vec3 extdir;
    float radius;
    float depth;
    vec3 dir;
    vec3 color;
} VertexOut;

void main(){
    float t = sqrt(1-VertexOut.z*VertexOut.z)/VertexOut.z;
    float dx = VertexOut.radius*VertexOut.coord.x/(t-VertexOut.radius);
    float x = dx+VertexOut.coord.x;
    x = VertexOut.coord.x;

    if (x < 0 || x > 1)
        discard;
    float radius = VertexOut.radius*x;

    float zz = 1.0/((VertexOut.radius/radius)*(VertexOut.radius/radius)) - VertexOut.coord.y * VertexOut.coord.y;
    if (zz < 0) discard;
    float z = sqrt(zz);


    vec2  Pdepth = (P_*vec4(0, 0, VertexOut.depth+z*radius*inversesqrt(1-VertexOut.z*VertexOut.z), 1.0)).zw;
    float depth  = Pdepth.x/Pdepth.y;
    gl_FragDepth    = depth*0.5+0.5;

    vec3 normal = -VertexOut.coord.y*normalize(VertexOut.extdir)+z*normalize(VertexOut.zdir);//+dx*VertexOut.dir;
    normal = normalize(normal);
    color_ = vec4((0.5+0.5*clamp(dot(normal, vec3(0, 0, 1)), 0, 1))*VertexOut.color, 1);
    //color_ = vec4(gl_FragDepth);
}
