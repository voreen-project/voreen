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

#version 330

in vec3  geom_color;
in vec2  geom_coord;
in float geom_depth;
in float geom_selected;
out vec4 frag_color;

uniform mat4 P_;
uniform float radius_;
uniform float selectionId_;


bool eq(float a, float b){
    return abs(a-b) < 0.1;
}

void main(){
    vec3 color;
    float coord_magnitude = dot(geom_coord, geom_coord); 
    if ( coord_magnitude > 1)
        discard;
    float z = sqrt(1- coord_magnitude);
    vec3 normal = vec3(geom_coord, z);
    
    vec3 lightdir = normalize(vec3(0.5, 0.5, 1));
    
    float light_intensity = 0.5*clamp(dot(normal, lightdir), 0, 1)+0.5;
    
    
    color = geom_color*light_intensity;

    if (coord_magnitude > 0.5 && (eq(geom_selected, 1) || eq(geom_selected, 3)) ) {
        color = (vec3(1)-geom_color)*light_intensity;
    }
    if (coord_magnitude < 0.5 && (eq(geom_selected, 2) || eq(geom_selected, 3)) ){
        color = (vec3(1)-geom_color)*light_intensity;
    }

    frag_color = vec4(color, 1.0);
    
    float depth = radius_*z;
    vec2 depth_transformed_zw = (P_*vec4(0, 0, geom_depth+depth, 1)).zw;
    float depth_transformed = depth_transformed_zw.x/depth_transformed_zw.y; 
    
    gl_FragDepth = depth_transformed*0.5+0.5;
    //gl_FragDepth = gl_FragCoord.z;
}