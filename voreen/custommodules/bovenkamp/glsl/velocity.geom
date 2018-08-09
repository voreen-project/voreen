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

  //#extension GL_EXT_geometry_shader4 : enable
    
  layout(points) in;
  layout(triangle_strip, max_vertices=30) out; 
    
  uniform sampler1D ColorTexture;
  uniform float maxVelocity;
  uniform float hight_;
  uniform bool colorFromDir_;
  uniform vec3 sliceAlignement_;
  
  in VertexData {
    vec3 normal;
  } VertexIn[];

  out VertexData {
    //vec3 normal;
    //vec3 vertex_to_light;
    vec4 color;
  } VertexOut;
  
vec3 dirToColor(in vec3 dir) { //normalized dir
    vec3 tmp = dir/vec3(2.0)+vec3(0.5);
    return tmp;// tmp / max(tmp.x, max(tmp.y, tmp.z)); is not correct, since (0,0,0) is not hit
}

void main() {

    //get color
    VertexOut.color = texture(ColorTexture,length(VertexIn[0].normal)/maxVelocity);
    if(colorFromDir_)
        VertexOut.color = vec4(dirToColor(normalize(VertexIn[0].normal)), (VertexOut.color.a > 0.0 ? 1.0 : 0.0));

    vec3 dir = normalize(VertexIn[0].normal*(vec3(1.0,1.0,1.0)-sliceAlignement_));
    
        
    if(dir.x == 0 && dir.y == 0 && dir.z == 0){
            //no displacement
    }
    else{
        //change direction order
        if(sliceAlignement_ == vec3(0.0,1.0,0.0)){
            dir.y = dir.z;
            dir.z = 0.0;
        } else if(sliceAlignement_ == vec3(1.0,0.0,0.0)) {
            dir.x = dir.y;
            dir.y = dir.z;
            dir.z = 0.0;
        }
        
    
        //getAxis
        vec3 x_dir = cross(dir,vec3(0.0,0.0,1.0))*0.5*hight_;
        vec3 y_dir = dir*0.5*hight_;
        vec3 pos = gl_in[0].gl_Position.xyz;
        vec3 pos01, pos02, pos03, pos04, pos05, pos06, pos07;

        pos01 = pos - 0.5*x_dir - y_dir;
        pos02 = pos + 0.5*x_dir - y_dir;
        pos03 = pos - 0.5*x_dir;
        pos04 = pos + 0.5*x_dir;
        pos05 = pos - x_dir;
        pos06 = pos + x_dir;
        pos07 = pos + y_dir;
        
        /*ver01 = gl_ModelViewProjectionMatrix*vec4(pos01,1);
        ver02 = gl_ModelViewProjectionMatrix*vec4(pos02,1);
        ver03 = gl_ModelViewProjectionMatrix*vec4(pos03,1);
        ver04 = gl_ModelViewProjectionMatrix*vec4(pos04,1);
        ver05 = gl_ModelViewProjectionMatrix*vec4(pos05,1);
        ver06 = gl_ModelViewProjectionMatrix*vec4(pos06,1);
        ver07 = gl_ModelViewProjectionMatrix*vec4(pos07,1);
        ver08 = gl_ModelViewProjectionMatrix*vec4(pos08,1);
        ver09 = gl_ModelViewProjectionMatrix*vec4(pos09,1);
        ver10 = gl_ModelViewProjectionMatrix*vec4(pos10,1);
        ver11 = gl_ModelViewProjectionMatrix*vec4(pos11,1);
        ver12 = gl_ModelViewProjectionMatrix*vec4(pos12,1);
        ver13 = gl_ModelViewProjectionMatrix*vec4(pos13,1);
        
        nor01 = normalize(gl_NormalMatrix*(-dir));
        nor02 = normalize(gl_NormalMatrix*(    - x_dir + y_dir));
        nor03 = normalize(gl_NormalMatrix*(    + x_dir + y_dir));
        nor04 = normalize(gl_NormalMatrix*(dir - x_dir + y_dir));
        nor05 = normalize(gl_NormalMatrix*(dir + x_dir + y_dir));
        nor06 = normalize(gl_NormalMatrix*(dir - x_dir - y_dir));
        nor07 = normalize(gl_NormalMatrix*(dir + x_dir - y_dir));
        nor08 = normalize(gl_NormalMatrix*dir);
        nor09 = normalize(gl_NormalMatrix*(dir         + y_dir));
        nor10 = normalize(gl_NormalMatrix*(dir + x_dir        ));
        nor11 = normalize(gl_NormalMatrix*(dir         - y_dir));
        nor12 = normalize(gl_NormalMatrix*(dir - x_dir        ));*/
        
        //bottem
        gl_Position = gl_ModelViewProjectionMatrix*vec4(pos01,1.0); EmitVertex();
        gl_Position = gl_ModelViewProjectionMatrix*vec4(pos02,1.0); EmitVertex();
        gl_Position = gl_ModelViewProjectionMatrix*vec4(pos03,1.0); EmitVertex();
        gl_Position = gl_ModelViewProjectionMatrix*vec4(pos04,1.0); EmitVertex();
        EndPrimitive();
        //front
        gl_Position = gl_ModelViewProjectionMatrix*vec4(pos05,1.0); EmitVertex();
        gl_Position = gl_ModelViewProjectionMatrix*vec4(pos06,1.0); EmitVertex();
        gl_Position = gl_ModelViewProjectionMatrix*vec4(pos07,1.0); EmitVertex();
        EndPrimitive();
    }
}


    
