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

uniform float radius;

varying vec3 vertex_light_position;
varying vec4 eye_position;

void main()
{
    // r^2 = (x - x0)^2 + (y - y0)^2 + (z - z0)^2
    float x = gl_TexCoord[0].x;
    float y = gl_TexCoord[0].y;
    float zz = 1.0 - x*x - y*y;

    if (zz <= 0.0) {
        discard;
    } else {
        float z = sqrt(zz);

        vec3 normal = vec3(x, y, z);
        vec3 light = vec3(0, 0, 1);
        // Lighting
        // use extern lightsource
        //float diffuse_value = max(dot(normal, vertex_light_position), 0.0);
        // use eye-position as lightsource
        float diffuse_value = max(dot(normal, light), 0.0f);

        vec4 pos = eye_position;
        //pos.z += z * (radius / 2.0f);
        pos.z += z * radius;
        pos = gl_ProjectionMatrix * pos;
        gl_FragDepth = (pos.z / pos.w + 1.0) / 2.0;

        vec4 col = gl_Color;
        col = col * diffuse_value;
        col[3] = 1.0f;

        if(z < 0.3)
            FragData0 = vec4(0.0, 0.0, 0.0, 0.0);
        else
            FragData0 = vec4(col.xyz, gl_Color.a);
    }
}
