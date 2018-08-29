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

varying vec3 normalVec;
varying vec3 lightVec;
varying vec3 viewVec;
varying float depth;

uniform vec3 camPos_;
uniform vec3 lightDir_;

void main()
{
    // The uniforms and gl_xxx were given in world coordinates, transform them
    // into eye space
    //
    normalVec = normalize(gl_NormalMatrix * gl_Normal);
    lightVec = normalize((gl_ModelViewMatrix * vec4(-lightDir_, 0.0)).xyz);
    viewVec = normalize((gl_ModelViewMatrix * (vec4(camPos_, 1.0) - gl_Vertex)).xyz);

    gl_FrontColor = gl_Color;
    gl_Position = ftransform();
    depth = gl_Position.z;
    gl_TexCoord[0] = gl_MultiTexCoord0;
}
