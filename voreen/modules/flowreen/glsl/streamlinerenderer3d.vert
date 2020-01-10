/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

layout(location=0) in vec3 position_;
layout(location=1) in vec3 velocity_;
layout(location=2) in float radius_;

//set from camera in setGloalShaderParameters
uniform mat4 projectionMatrix_;
uniform mat4 viewMatrix_;
//set manuell in process
uniform mat4 voxelToWorldMatrix_;
uniform mat4 velocityTransformMatrix_;

out vData
{
    vec3 position;
    vec3 velocity;
    float radius;
} vertex;

void main()
{
    gl_Position = projectionMatrix_ * viewMatrix_ * voxelToWorldMatrix_ * vec4(position_,1);
    vertex.position = gl_Position.xyz;
    vertex.velocity = (velocityTransformMatrix_* vec4(velocity_,1.0)).xyz;
    vertex.radius = radius_;
}
