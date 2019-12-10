/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

in vec3  geom_color;
in vec2  geom_coord;
in float geom_depth;
in float geom_active;
in float geom_id;

out vec4 frag_color;

uniform float radius_;
uniform float selectionid_;

void main(){
    float coord_magnitude = dot(geom_coord, geom_coord); 
    if ( coord_magnitude > 1)
        discard;
    float z = sqrt(1- coord_magnitude);
    vec3 normal = vec3(geom_coord, z);
    
    vec3 lightdir = normalize(vec3(0.5, 0.5, 1));
    
    float light_intensity = 0.5*clamp(dot(normal, lightdir), 0, 1)+0.5;
    
    vec3 color = geom_color;
    if (geom_active == 1.0 && coord_magnitude > 0.5){
        color = vec3(1)-geom_color;
    }
    if (geom_id == selectionid_ && coord_magnitude < 0.5){
        color = vec3(1)-geom_color;
    }
    frag_color = vec4(color*light_intensity, 1);
}