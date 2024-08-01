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

typedef struct __attribute__((packed)) {
    uint4 dimensions_;
    uchar numChannels_;
    uchar numBitsPerChannel_;
    uchar isSigned_;
    uchar isFloat_;
    intptr_t data_;
} volume_t;

void write_volumef(__global volume_t* volume, uint4 coord, float4 value) {

    if(volume->data_ == 0 || coord.x > volume->dimensions_.x || coord.y > volume->dimensions_.y || coord.z > volume->dimensions_.z)
        return;

    __global void* data = (__global void*)(volume->data_);

    uint scale = (1 << volume->numBitsPerChannel_) - 1;

    uint coordLinear = coord.x + coord.y * volume->dimensions_.x + coord.z * volume->dimensions_.x * volume->dimensions_.y;

    if(volume->isFloat_ > 0) {
        if(volume->numChannels_ == 1)
            ((__global float*)(data))[coordLinear] = value.x;
        else if(volume->numChannels_ == 3) {
            ((__global float*)(data))[coordLinear] = value.x;
            ((__global float*)(data))[coordLinear + 1] = value.y;
            ((__global float*)(data))[coordLinear + 2] = value.z;
        } else if(volume->numChannels_ == 4)
            ((__global float4*)(data))[coordLinear] = value;

        return;
    }

    if(volume->numBitsPerChannel_ == 8) {
        if(volume->numChannels_ == 1)
            ((__global uchar*)(data))[coordLinear] = convert_uchar(value.x * scale);
        else if(volume->numChannels_ == 3) {
            ((__global uchar*)(data))[coordLinear    ] = convert_uchar(value.x * scale);
            ((__global uchar*)(data))[coordLinear + 1] = convert_uchar(value.y * scale);
            ((__global uchar*)(data))[coordLinear + 2] = convert_uchar(value.z * scale);
        } else if(volume->numChannels_ == 4)
            ((__global uchar4*)(data))[coordLinear] = convert_uchar4(value * scale);

        return;
    }

    if(volume->numBitsPerChannel_ == 16) {
        if(volume->numChannels_ == 1)
            ((__global ushort*)(data))[coordLinear] = convert_ushort(value.x * scale);
        else if(volume->numChannels_ == 3) {
            ((__global ushort*)(data))[coordLinear    ] = convert_ushort(value.x * scale);
            ((__global ushort*)(data))[coordLinear + 1] = convert_ushort(value.y * scale);
            ((__global ushort*)(data))[coordLinear + 2] = convert_ushort(value.z * scale);
        } else if(volume->numChannels_ == 4)
            ((__global ushort4*)(data))[coordLinear] = convert_ushort4(value * scale);

        return;
    }
}

void write_volumef_norm(__global volume_t* volume, float4 coord, float4 value) {
    write_volumef(volume, convert_uint4(coord * convert_float4(volume->dimensions_)), value);
}

void write_volumeui(__global volume_t* volume, uint4 coord, uint4 value) {

    if(volume->data_ == 0 || coord.x > volume->dimensions_.x || coord.y > volume->dimensions_.y || coord.z > volume->dimensions_.z)
        return;

    __global void* data = (__global void*)(volume->data_);

    uint coordLinear = coord.x + coord.y * volume->dimensions_.x + coord.z * volume->dimensions_.x * volume->dimensions_.y;

    if(volume->isFloat_ > 0) {

        if(volume->numChannels_ == 1)
            ((__global float*)(data))[coordLinear] = convert_float(value.x);
        else if(volume->numChannels_ == 3) {
            ((__global float*)(data))[coordLinear] = convert_float(value.x);
            ((__global float*)(data))[coordLinear + 1] = convert_float(value.y);
            ((__global float*)(data))[coordLinear + 2] = convert_float(value.z);
        } else if(volume->numChannels_ == 4)
            ((__global float4*)(data))[coordLinear] = convert_float4(value);

        return;
    }

    if(volume->numBitsPerChannel_ == 8) {
        if(volume->numChannels_ == 1)
            ((__global uchar*)(data))[coordLinear] = convert_uchar(value.x);
        else if(volume->numChannels_ == 3) {
            ((__global uchar*)(data))[coordLinear    ] = convert_uchar(value.x);
            ((__global uchar*)(data))[coordLinear + 1] = convert_uchar(value.y);
            ((__global uchar*)(data))[coordLinear + 2] = convert_uchar(value.z);
        } else if(volume->numChannels_ == 4)
            ((__global uchar4*)(data))[coordLinear] = convert_uchar4(value);

        return;
    }

    if(volume->numBitsPerChannel_ == 16) {
        if(volume->numChannels_ == 1)
            ((__global ushort*)(data))[coordLinear] = convert_ushort(value.x);
        else if(volume->numChannels_ == 3) {
            ((__global ushort*)(data))[coordLinear    ] = convert_ushort(value.x);
            ((__global ushort*)(data))[coordLinear + 1] = convert_ushort(value.y);
            ((__global ushort*)(data))[coordLinear + 2] = convert_ushort(value.z);
        } else if(volume->numChannels_ == 4)
            ((__global ushort4*)(data))[coordLinear] = convert_ushort4(value);

        return;
    }
}

void write_volumeui_norm(__global volume_t* volume, float4 coord, uint4 value) {
    write_volumeui(volume, convert_uint4(coord * convert_float4(volume->dimensions_)), value);
}

