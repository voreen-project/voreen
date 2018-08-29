/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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
  //#define ARRAY_SIZE_ xxx;

  layout(lines) in;
  layout(line_strip, max_vertices=17) out;

  uniform int iteration_;
  uniform float maximum_;
  uniform float factor_;

  in VertexData {
    vec3 center;
    float timestep;
  } VertexIn[];

  out VertexData {
    float color;
  } VertexOut;


void main() {
    // load start setup
    int ind_a = 0;
    int ind_b = ARRAY_SIZE_-1;

    vec3 arr[ARRAY_SIZE_];

    float timestep = VertexIn[0].timestep;
    float timestep_step = 1.f/(ARRAY_SIZE_-1);
    float color = length(gl_in[0].gl_Position.xyz)/maximum_;
    float color_step = (length(gl_in[1].gl_Position.xyz)/maximum_ - color)/(ARRAY_SIZE_-1);

    arr[ind_a] = normalize(gl_in[0].gl_Position.xyz);
    arr[ind_b] = normalize(gl_in[1].gl_Position.xyz);

    // fill array
    int internalIteration = 1;
    int tmp_a, tmp_b;
    int tmp_start_c, tmp_c = 0;
    int tmp_old_c = 0;
    for(int i = 0; i < iteration_; i++) {
        tmp_start_c = ind_b/(2*(i+1));
        tmp_c = tmp_start_c;
        tmp_a = 0;
        tmp_b = 2*tmp_start_c;
        for(int j = 1; j <= internalIteration; j++) {
            arr[tmp_c] = normalize(arr[tmp_a] + arr[tmp_b]);
            tmp_c += tmp_old_c;
            tmp_a = tmp_b;
            tmp_b += tmp_old_c;
        }
        tmp_old_c = tmp_start_c;
        internalIteration *= 2;
    }

    for(int i = 0; i < ARRAY_SIZE_; i++){
        VertexOut.color = color;
        gl_Position = gl_ModelViewProjectionMatrix*vec4(arr[i]*timestep*factor_+VertexIn[0].center,1.0);
        EmitVertex();
        timestep += timestep_step;
        color += color_step;
    }
    EndPrimitive();
}


