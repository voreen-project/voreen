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

  //#extension GL_EXT_geometry_shader4 : enable

  layout(points) in;
  layout(triangle_strip, max_vertices=30) out;

  float size =  0.1;

  uniform sampler1D ColorTexture;
  uniform float maxVelocity;
  uniform float hight_;
  uniform bool colorFromDir_;

  in VertexData {
    vec3 normal;
  } VertexIn[];

  out VertexData {
    vec3 normal;
    vec3 vertex_to_light;
    vec4 color;
  } VertexOut;

vec3 dirToColor(in vec3 dir) {
    vec3 tmp = dir/vec3(2.0)+vec3(0.5);
    return tmp;// tmp / max(tmp.x, max(tmp.y, tmp.z)); not correct since 0,0,0 can not be hit
}

void main() {

    if(VertexIn[0].normal.x == 0 && VertexIn[0].normal.y == 0 && VertexIn[0].normal.z == 0){
            //no displacement
    } else {

    vec3 dir = normalize(VertexIn[0].normal);

    VertexOut.color = texture(ColorTexture,length(VertexIn[0].normal)/maxVelocity);

    if(colorFromDir_)
        VertexOut.color = vec4(dirToColor(dir), (VertexOut.color.a > 0.0 ? 1.0 : 0.0));

        //getAxis
        vec3 x_dir;
        vec3 y_dir;
        vec3 pos = (gl_in[0].gl_Position -(0.5*vec4(dir*hight_,0))).xyz;
        vec3 pos01, pos02, pos03, pos04, pos05, pos06, pos07, pos08,
             pos09, pos10, pos11, pos12, pos13;
        vec4 ver01, ver02, ver03, ver04, ver05, ver06, ver07, ver08,
             ver09, ver10, ver11, ver12, ver13;
        vec3 nor01, nor02, nor03, nor04, nor05, nor06, nor07, nor08,
             nor09, nor10, nor11, nor12;

        if(dir.x == 0){
            x_dir = vec3(1.0,0.0,0.0);
        } else {
            if(dir.y == 0){
                x_dir = vec3(0.0,1.0,0.0);
            } else {
                x_dir = normalize(vec3(-dir.y,dir.x,0));
            }
        }

        y_dir = cross(dir,x_dir);
        x_dir = x_dir*size*hight_;
        y_dir = y_dir*size*hight_;
        dir *= hight_;

        pos01 = pos - x_dir + y_dir;
        pos02 = pos + x_dir + y_dir;
        pos03 = pos - x_dir - y_dir;
        pos04 = pos + x_dir - y_dir;
        pos05 = pos01 + dir*0.75;
        pos06 = pos02 + dir*0.75;
        pos07 = pos03 + dir*0.75;
        pos08 = pos04 + dir*0.75;
        pos09 = pos05 - x_dir + y_dir;
        pos10 = pos06 + x_dir + y_dir;
        pos11 = pos07 - x_dir - y_dir;
        pos12 = pos08 + x_dir - y_dir;
        pos13 = pos + dir;

        ver01 = gl_ModelViewProjectionMatrix*vec4(pos01,1);
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
        nor12 = normalize(gl_NormalMatrix*(dir - x_dir        ));

        //bottem
        VertexOut.normal = nor01;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos04,1)).xyz;
        gl_Position = ver04; EmitVertex();
        VertexOut.normal = nor01;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos02,1)).xyz;
        gl_Position = ver02; EmitVertex();
        VertexOut.normal = nor01;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos03,1)).xyz;
        gl_Position = ver03; EmitVertex();
        VertexOut.normal = nor01;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos01,1)).xyz;
        gl_Position = ver01; EmitVertex();
        EndPrimitive();
        //front
        VertexOut.normal = nor02;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos01,1)).xyz;
        gl_Position = ver01; EmitVertex();
        VertexOut.normal = nor02;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos05,1)).xyz;
        gl_Position = ver05; EmitVertex();
        VertexOut.normal = nor03;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos02,1)).xyz;
        gl_Position = ver02; EmitVertex();
        VertexOut.normal = nor03;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos06,1)).xyz;
        gl_Position = ver06; EmitVertex();
        //right
        VertexOut.normal = -nor02;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos04,1)).xyz;
        gl_Position = ver04; EmitVertex();
        VertexOut.normal = -nor02;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos08,1)).xyz;
        gl_Position = ver08; EmitVertex();
        //back
        VertexOut.normal = -nor03;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos03,1)).xyz;
        gl_Position = ver03; EmitVertex();
        VertexOut.normal = -nor03;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos07,1)).xyz;
        gl_Position = ver07; EmitVertex();
        //left
        VertexOut.normal = nor02;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos01,1)).xyz;
        gl_Position = ver01; EmitVertex();
        VertexOut.normal = nor02;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos05,1)).xyz;
        gl_Position = ver05; EmitVertex();
        EndPrimitive();
        //arrow bottem
        VertexOut.normal = nor01;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos09,1)).xyz;
        gl_Position = ver09; EmitVertex();
        VertexOut.normal = nor01;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos10,1)).xyz;
        gl_Position = ver10; EmitVertex();
        VertexOut.normal = nor01;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos11,1)).xyz;
        gl_Position = ver11; EmitVertex();
        VertexOut.normal = nor01;
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos12,1)).xyz;
        gl_Position = ver12; EmitVertex();
        EndPrimitive();

        //pyramide
        VertexOut.normal = nor09; //04
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos09,1)).xyz;
        gl_Position = ver09; EmitVertex();
        VertexOut.normal = nor09; //05
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos10,1)).xyz;
        gl_Position = ver10; EmitVertex();
        VertexOut.normal = nor09; //08
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos13,1)).xyz;
        gl_Position = ver13; EmitVertex();
        EndPrimitive();
        VertexOut.normal = nor10; //05
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos10,1)).xyz;
        gl_Position = ver10; EmitVertex();
        VertexOut.normal = nor10; //07
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos12,1)).xyz;
        gl_Position = ver12; EmitVertex();
        VertexOut.normal = nor10; //08
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos13,1)).xyz;
        gl_Position = ver13; EmitVertex();
        EndPrimitive();
        VertexOut.normal = nor11; //07
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos12,1)).xyz;
        gl_Position = ver12; EmitVertex();
        VertexOut.normal = nor11; //06
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos11,1)).xyz;
        gl_Position = ver11; EmitVertex();
        VertexOut.normal = nor11; //08
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos13,1)).xyz;
        gl_Position = ver13; EmitVertex();
        EndPrimitive();
        VertexOut.normal = nor12; //06
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos11,1)).xyz;
        gl_Position = ver11; EmitVertex();
        VertexOut.normal = nor12; //12
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos13,1)).xyz;
        gl_Position = ver13; EmitVertex();
        VertexOut.normal = nor12; //04
        VertexOut.vertex_to_light = (gl_LightSource[0].position - gl_ModelViewMatrix * vec4(pos09,1)).xyz;
        gl_Position = ver09; EmitVertex();
        EndPrimitive();
    }
}



