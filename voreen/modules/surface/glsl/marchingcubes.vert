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

out uvec3 voxelPos;

uniform ivec3 volumeDimensions_;

void main() {
    // We extract surfaces in the spaces between the voxels.
    // Thus we actually actually only consider volumeDimensions_ - (1,1,1) cuboid spaces.
    ivec3 vd_minus_one = volumeDimensions_ - ivec3(1,1,1);
    int x = gl_VertexID%vd_minus_one.x;
    int y = ((gl_VertexID-x)/vd_minus_one.x)%vd_minus_one.y;
    int z = ((gl_VertexID-x-y*vd_minus_one.x)/vd_minus_one.x/vd_minus_one.y)%vd_minus_one.z;
    voxelPos = uvec3(x, y, z);
}
