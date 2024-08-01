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

#include "mod_volumewrite.cl"
#include "mod_gradients.cl"

const sampler_t smp = CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

bool border(uint4 pos, uint4 dim) {
    return (pos.x == 0 || pos.x >= dim.x-1 || pos.y == 0 || pos.y >= dim.y-1 || pos.z == 0 || pos.z >= dim.z-1);
}

// The functions calcFilteredGradient and calcGradient are defined in "mod_gradients.cl".
// The function write_volumef and the type volume_t are defined in "mod_volumewrite.cl".
// The argument "volume_t volume" initially contains only information (like dimension and datatype) about the volume which we want to write into.
// The pointer to the actual data buffer has to be passed separately at the moment, but can (and should) be assigned to the volume_t by hand.
__kernel void gradient( read_only image3d_t inputVolume, __global volume_t* volume, __global void* volData, int copyIntensityChannel) {

    // This is (apparently) currently the only way to have a pointer to GPU memory within a struct.
    volume->data_ = (intptr_t)volData;

    uint4 pos = (uint4)(get_global_id(0), get_global_id(1), get_global_id(2), 1);
    uint4 dim = (uint4)(get_global_size(0), get_global_size(1), get_global_size(2), 1);
    float4 normPos = convert_float4(pos) / convert_float4(dim);

    int numChannels = (copyIntensityChannel == 1) ? 4 : 3;

    if (border(pos, dim))
        write_volumef(volume, pos, (float4)(0.5, 0.5, 0.5, 0.0));
    else {
        // The function write_volumef handles the datatype and indexing of the volume data internally.
#if defined(GRADIENT_FILTERED)
        write_volumef(volume, pos, (float4)(0.5f * calcFilteredGradient(normPos, inputVolume).xyz + 0.5f, read_imagef(inputVolume, smp, normPos).w));
#elif defined(GRADIENT_CENTRAL)
        write_volumef(volume, pos, (float4)(0.5f * calcGradient(normPos, inputVolume).xyz - 0.5f, read_imagef(inputVolume, smp, normPos).w));
#endif
    }
}

