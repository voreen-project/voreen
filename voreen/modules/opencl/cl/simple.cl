/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

__constant sampler_t smpNorm = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP | CLK_FILTER_LINEAR;

// Simple OpenCL kernel.  Extend it's arguments to dynamically generate ports and properties.
__kernel void standardKernel(read_only image2d_t inCol, //$ color inport
                             write_only image2d_t outCol, //$ color outport
                             float test //$ float property to test @min=0.0 @max=1.0 @value=0.5
                            )
{
    //output image pixel coordinates
    int2 target = (int2)(get_global_id(0), get_global_id(1));
    float2 targetNorm = convert_float2(target) / convert_float2((int2)(get_global_size(0), get_global_size(1)));

    float4 color = read_imagef(inCol, smpNorm, targetNorm);
    write_imagef(outCol, target, color * test);
}


