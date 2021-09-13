/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2012 University of Muenster, Germany.                        *
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

  //#extension GL_EXT_geometry_shader4 : enable
    
  layout(points) in;
  layout(triangle_strip, max_vertices=7) out; 
  
  float size =  0.1;
  
  uniform sampler1D ColorTexture;
  uniform float maxVelocity;
  uniform float hight_;
  uniform bool colorFromDir_;
  uniform mat3 MVPtinv_;

uniform mat4 MVP_;
uniform mat4 MV_;
uniform mat4 P_;


  in VertexData {
    vec3 normal;
	vec4 color;
  } VertexIn[];

  out VertexData {
    vec4 color;
  } VertexOut;
  
vec4 dirToColor(vec3 dir) {
    return vec4(dir*vec3(0.5)+vec3(0.5), 1.0);
}

void main() {
	vec3 pos = gl_in[0].gl_Position.xyz;
	vec3 dir = normalize(VertexIn[0].normal)*hight_;
	vec3 ortho = normalize(cross(pos, vec3(0, 0, 1)))*0.25*hight_;
	
	//vec3 normal = cross(dir, ortho);
	
	
	vec4 color = VertexIn[0].color;//*dot(vec3(0, 0, 1),s normal);
	
	
    VertexOut.color = color;
    gl_Position = P_*vec4(pos+dir, 1.0); 
	EmitVertex();
	
	VertexOut.color = color;
    gl_Position = P_*vec4(pos+ortho+0.7*dir, 1.0); 
	EmitVertex();
	
	VertexOut.color = color;
    gl_Position = P_*vec4(pos-ortho+0.7*dir, 1.0); 
	EmitVertex();
    EndPrimitive();
	
	gl_Position = P_*vec4(pos+ortho*0.3+0.7*dir, 1.0); 
	EmitVertex();
	
	gl_Position = P_*vec4(pos-ortho*0.3+0.7*dir, 1.0); 
	EmitVertex();
	
	gl_Position = P_*vec4(pos+ortho*0.3, 1.0); 
	EmitVertex();
	
	gl_Position = P_*vec4(pos-ortho*0.3, 1.0); 
	EmitVertex();
	EndPrimitive();
}


    
