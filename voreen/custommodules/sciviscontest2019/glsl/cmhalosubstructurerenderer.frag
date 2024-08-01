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
#define GLSL_VERSION_130
#include "modules/mod_transfunc.frag"

out vec4 out_color;
in vec2 frag_coord;
in vec3 frag_vel;
in vec3 frag_view_pos;
in float frag_radius;
in float frag_mass;
flat in int  frag_vertexID;


uniform TransFuncParameters transferFunc_;
uniform sampler1D transferFuncTex_;
uniform float alphaFactor_;

uniform float scale_;
uniform float offset_;
uniform float hostHaloRingThickness_;

uniform mat4 projectionMatrix;

void main(){
    float oneMinusThickness = 1-hostHaloRingThickness_;
    float fdot = dot(frag_coord, frag_coord);
    if (fdot > 1 || (frag_vertexID==0 && fdot < (oneMinusThickness*oneMinusThickness)))
        discard;
    float dz= sqrt(1 - frag_coord.x*frag_coord.x - frag_coord.y*frag_coord.y);
    vec3 normal = vec3(frag_coord.x,frag_coord.y,dz);
    vec3 surfacePos = frag_view_pos + normal*frag_radius;
    vec3 halfwayVec = normalize(-surfacePos);

    vec4 view_z = vec4(0,0,frag_view_pos.z+dz*frag_radius,1);
    vec4 predevide_z = (projectionMatrix * view_z);

    gl_FragDepth = (0.5*predevide_z.z/predevide_z.w) + 0.5;


    float dcont=max(0.0,dot(normal,halfwayVec));
    //float fallOff = 1/(1+dot(surfacePos, surfacePos));
    float fallOff = 0.5;
    vec3 diffuse=pow(dcont,3)*vec3(1,1,1)*fallOff;
    out_color = vec4(0.5*(applyTF(transferFunc_, transferFuncTex_, frag_mass)).xyz+diffuse,alphaFactor_);
    out_color.a = alphaFactor_;
}
