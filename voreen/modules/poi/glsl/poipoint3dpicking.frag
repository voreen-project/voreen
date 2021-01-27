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

#version 330

in vec2 geom_coord;
in float geom_depth;
in float geom_id;


out uint out_color;

uniform mat4 P_;
uniform float radius_;

void main(){
    
    float coord_magnitude = dot(geom_coord, geom_coord); 
    if ( coord_magnitude > 1)
        discard;
    float z = sqrt(1- coord_magnitude);
    
    float depth = radius_*z;
    vec2 depth_transformed_zw = (P_*vec4(0, 0, geom_depth+z, 1)).zw;
    float depth_transformed = depth_transformed_zw.x/depth_transformed_zw.y; 
    
    gl_FragDepth = depth_transformed*0.5+0.5;
    out_color = uint(geom_id);
}