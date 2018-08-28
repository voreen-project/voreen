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

const sampler_t smp = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

float4 readFloat(image2d_t img, sampler_t smp, int2 p) {
    if(get_image_channel_data_type(img) == CLK_UNSIGNED_INT16) {
        uint4 tmp = read_imageui(img, smp, p);
        float4 ret;
        ret.x = tmp.x / 65535.0f;
        ret.y = tmp.y / 65535.0f;
        ret.z = tmp.z / 65535.0f;
        ret.w = tmp.w / 65535.0f;
        return ret;
    }
    else {
        return read_imagef(img, smp, p);
    }
}

void writeFloat(write_only image2d_t img, int2 p, float4 value) {
    if(get_image_channel_data_type(img) == CLK_UNSIGNED_INT16) {
        uint4 write;
        write.x = value.x * 65535.0f;
        write.y = value.y * 65535.0f;
        write.z = value.z * 65535.0f;
        write.w = value.w * 65535.0f;
        write_imageui(img, p, write);
    }
    else {
        write_imagef(img, p, value);
    }
}

float4 toGrayScale(float4 color, float saturation) {
    float brightness = ((0.30 * color.x) + (0.59 * color.y) + (0.11 * color.z));
    float4 grayscale;
    grayscale.xyz = brightness;
    grayscale.w = color.w;

    return mix(grayscale, color, saturation);
}

__kernel void gr(read_only image2d_t srcimg, write_only image2d_t output, float saturation)
{
    int tidX = get_global_id(0), tidY = get_global_id(1);
    float4 color = readFloat(srcimg, smp, (int2)(tidX, tidY));
    color = toGrayScale(color, saturation);

    writeFloat(output, (int2)( tidX, tidY ), color);
}

