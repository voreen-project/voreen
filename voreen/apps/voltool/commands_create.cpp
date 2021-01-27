/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "commands_create.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "tgt/spline.h"

using tgt::ivec3;
using tgt::vec3;

namespace voreen {

CommandGenerateMask::CommandGenerateMask() /*:
    Command("--generatemask", "", "", "<DX DY DZ RADIUS CX CY CZ OUT>", 8) */
{
//    loggerCat_ += "." + name_;
}

bool CommandGenerateMask::checkParameters(const std::vector<std::string>& parameters) {
    /*bool correctType =
        is<int>(parameters[0]) && is<int>(parameters[1]) && is<int>(parameters[2]) && is<float>(parameters[3]) &&
        is<float>(parameters[4]) && is<float>(parameters[5]) && is<float>(parameters[6]);
    return (parameters.size() == 8) && correctType; */
    return false;
}

bool CommandGenerateMask::execute(const std::vector<std::string>& parameters) {
    /*vec3 center;
    ivec3 dimensions;
    dimensions.x = cast<int>(parameters[0]);
    dimensions.y = cast<int>(parameters[1]);
    dimensions.z = cast<int>(parameters[2]);
    float radius = cast<float>(parameters[3]);
    center.x = cast<float>(parameters[4]);
    center.y = cast<float>(parameters[5]);
    center.z = cast<float>(parameters[6]);

    LINFO("Generating mask dataset with dimensions: " << dimensions);

    VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dimensions);
    Volume* targetHandle = new Volume(target, vec3(1.0f), vec3(0.0f));
    oldVolumePosition(targetHandle);

    for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
        for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
            for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {
                vec3 diff = center - static_cast<vec3>(ivec3(voxel_x, voxel_y, voxel_z));
                if (length(diff) <= radius)
                    target->voxel(voxel_x, voxel_y, voxel_z) = (uint8_t)255;
                else
                    target->voxel(voxel_x, voxel_y, voxel_z)= (uint8_t)0;
            }
        }
    }

    if (target) {
        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        serializer->write(parameters.back(), targetHandle);
        delete targetHandle;
        return true;
    }
    else */
        return false;
}

//-----------------------------------------------------------------------------

CommandCreate::CommandCreate() /*:
    Command("--create", "", "Create datasets. \n\
\t\tWrites 8-bit dataset with dimensions SIZExSIZExSIZE\n\
\t\tPossible datasets:\n\
\t\tcornell: creates a cornell box\n\
\t\tcube:\n\
\t\tsynth: \n\
\t\tcloud: \n\
\t\taotestbox: \n\
\t\tshadowtestvol: \n\
\t\taorticarch: generates a bezier-shaped aorta with three arches",
"<[cornell|cube|synth|blobs|blobs2|blobs3|sphere|doublesphere|spherecoord|cloud|aotestbox|shadowtestvol|aorticarch] SIZE OUT8>", 3) */

{
//    loggerCat_ += "." + name_;
}

bool CommandCreate::checkParameters(const std::vector<std::string>& parameters) {
    /*std::set<std::string> set;
    set.insert("cornell");
    set.insert("cube");
    set.insert("blobs");
    set.insert("blobs2");
    set.insert("blobs3");
    set.insert("sphere");
    set.insert("doublesphere");
    set.insert("spherecoord");
    set.insert("synth");
    set.insert("cloud");
    set.insert("aotestbox");
    set.insert("shadowtestvol");
    set.insert("aorticarch");
    return (parameters.size() == 3) && is<int>(parameters[1]) && isValueInSet(parameters[0], set); */
    return false;
}

bool CommandCreate::execute(const std::vector<std::string>& parameters) {
/*    std::string operation = parameters[0];

    ivec3 dimensions;
    dimensions.x = cast<int>(parameters[1]);
    dimensions.y = dimensions.x;
    dimensions.z = dimensions.x;

    VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dimensions);

    if (operation == "cornell") {
        vec3 center = dimensions / 2;

        int s = dimensions.x;
        int thickness = s / 15;
        int border = 3;
        int is = s - ((thickness+border)*2); //the inner size of the box

        uint8_t white = 200;
        uint8_t red = 70;
        uint8_t green = 150;

        uint8_t box0 = 100;
        uint8_t box1 = 50;
        uint8_t sphere = 120;

//         uint8_t light = 240;
//         uint8_t fog = 20;

//         float radius = dimensions.x / 4;

        LINFO("Generating cornell dataset with dimensions: " << dimensions);

        //clear dataset:
        fillBox(target, ivec3(0,0,0), dimensions, 0);



        //white floor
        fillBox(target, ivec3(border, border, border), ivec3(s-border, border+thickness, s-border), white);
        //white back wall:
        fillBox(target, ivec3(border, border, border), ivec3(s-border, s-border, border+thickness), white);
        //white ceiling:
        fillBox(target, ivec3(border, s-border-thickness, border), ivec3(s-border, s-border, s-border), white);

        //THE NEXT 2 ARE OPTIONAL:
        //white front wall:
//         fillBox(target, ivec3(border, border, border+thickness+is), ivec3(s-border, s-border, s-border), white);
        //cut hole in ceiling (lightsource)
//         fillBox(target, ivec3(center.x-thickness*2, s-border-thickness, center.x-thickness*2), ivec3(center.x+thickness*2, s-border, center.x+thickness*2), light);
        fillBox(target,
                ivec3(static_cast<int>(center.x-thickness*2), s-border-thickness, static_cast<int>(center.x-thickness*2)),
                ivec3(static_cast<int>(center.x+thickness*2), s-border, static_cast<int>(center.x+thickness*2)),
                0);

        //green right wall:
        fillBox(target, ivec3(s-thickness-border,thickness+border,thickness+border), ivec3(s-border,s-thickness-border,s-border), green);

        //red right wall:
        fillBox(target, ivec3(border,thickness+border,thickness+border), ivec3(border+thickness,s-thickness-border,s-border), red);

        //box0: oriented box, left, back
        fillOrientedBox(target, vec3(static_cast<float>(border+thickness+(is*0.3f)), 0.0f, static_cast<float>(border+thickness+(is*0.3f))),
            vec3(1.0f, 0.0f, 2.0f), static_cast<float>(is*0.15), static_cast<float>(is*0.2), static_cast<float>(border+thickness),
            static_cast<float>(border+thickness+(is*0.6f)), box0);

        //box1: oriented box, right, front
        fillOrientedBox(target, vec3(static_cast<float>(border+thickness+(is*0.7f)), 0.0f, static_cast<float>(border+thickness+(is*0.8f))),
            vec3(2.0f, 0.0f, 1.0f), static_cast<float>(is*0.15), static_cast<float>(is*0.15), static_cast<float>(border+thickness),
            static_cast<float>(border+thickness+(is*0.3f)), box1);

        //sphere, left, front
        fillSphere(target, vec3(static_cast<float>(border+thickness+(is*0.3f)), static_cast<float>(border+thickness+is/6),
            static_cast<float>(border+thickness+(is*0.8f))), static_cast<float>(is/6), sphere);
//         fillBox(targetDataset_, ivec3(), ivec3(), );

//         void fillSphere(VolumeDataset* vds, vec3 center, float radius, uint8_t value);
//         void fillBox(VolumeDataset* vds, ivec3 start, ivec3 end, uint8_t value);

    }
    else if (operation == "cube") {
        vec3 center;
        center = dimensions / 2;

        int s = dimensions.x;
        int thickness = s / 15;
        int border = 3;
//         int is = s - ((thickness+border)*2); //the inner size of the box

        uint8_t white = 200;
//         uint8_t red = 70;
//         uint8_t green = 150;

//         uint8_t box0 = 100;
//         uint8_t box1 = 50;
//         uint8_t sphere = 120;

//         uint8_t light = 240;
//         uint8_t fog = 20;

//         float radius = (float)dimensions.x / 4.0f;

        LINFO("Generating cornell dataset with dimensions: " << dimensions);

        //clear dataset:
        fillBox(target, ivec3(0,0,0), dimensions, 0);

        //white floor
        fillBox(target, ivec3(border, border, border), ivec3(s-border, border+thickness, s-border), white);
        //white back wall:
        /*
        fillBox(targetDataset_, ivec3(border, border, border), ivec3(s-border, s-border, border+thickness), white);
        //white ceiling:
        fillBox(targetDataset_, ivec3(border, s-border-thickness, border), ivec3(s-border, s-border, s-border), white);

        //THE NEXT 2 ARE OPTIONAL:
        //white front wall:
//         fillBox(targetDataset_, ivec3(border, border, border+thickness+is), ivec3(s-border, s-border, s-border), white);
        //cut hole in ceiling (lightsource)
        fillBox(targetDataset_, ivec3(center.x-thickness*2, s-border-thickness, center.x-thickness*2), ivec3(center.x+thickness*2, s-border, center.x+thickness*2), light);

        //green right wall:
        fillBox(targetDataset_, ivec3(s-thickness-border,thickness+border,thickness+border), ivec3(s-border,s-thickness-border,s-border), green);

        //red right wall:
        fillBox(targetDataset_, ivec3(border,thickness+border,thickness+border), ivec3(border+thickness,s-thickness-border,s-border), red);

        //box0: oriented box, left, back
        fillOrientedBox(targetDataset_, vec3(border+thickness+(int)(is*0.3f), 0.0f, border+thickness+(int)(is*0.3f)), vec3(1.0f, 0.0f, 2.0f), (int)(is*0.15), (int)(is*0.2), border+thickness, border+thickness+(int)is*0.6f, box0);

        //box1: oriented box, right, front
        fillOrientedBox(targetDataset_, vec3(border+thickness+(int)(is*0.7f), 0.0f, border+thickness+(int)(is*0.8f)), vec3(2.0f, 0.0f, 1.0f), (int)(is*0.15), (int)(is*0.15), border+thickness, border+thickness+(int)is*0.3f, box1);

        //sphere, left, front
        fillSphere(targetDataset_, vec3(border+thickness+(int)(is*0.3f), border+thickness+is/6, border+thickness+(int)(is*0.8f)), is/6, sphere);
//         fillBox(targetDataset_, ivec3(), ivec3(), );

//         void fillSphere(VolumeDataset* vds, vec3 center, float radius, uint8_t value);
//         void fillBox(VolumeDataset* vds, ivec3 start, ivec3 end, uint8_t value);
    }
    else if (operation == "blobs") {
        delete target;
        target = 0;
        VolumeRAM_UInt8* target2 = new VolumeRAM_UInt8(dimensions);

        tgt::vec3 center = dimensions / 2;

        LINFO("Generating blobs dataset with dimensions: " << dimensions);

        //clear dataset:
        fillBox(target2, tgt::ivec3(0,0,0), dimensions, 0);
        tgt::vec3 blob1 = tgt::vec3(center.x, center.y, dimensions.z * 0.25f);
        tgt::vec3 blob2 = tgt::vec3(center.x, center.y, dimensions.z * 0.75f);

        float radius = 0.22f / (float)dimensions.x * 32.0f;
        tgt::vec3 m  = blob1 + tgt::vec3(1.0f,1.0f,1.0f);
        float max = (1.0f / distance(blob1, m)) + (1.0f / distance(blob1, m));

        max -= radius;

        for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
            for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
                for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {
                    float b1 = 1.0f / distance(blob1, tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z));
                    float b2 = 1.0f / distance(blob2, tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z));

                    float v = b1+b2;
                    v -= radius;

                    if(v >= max)
                        target2->voxel(voxel_x, voxel_y, voxel_z) = 255;
                    else if(v > 0.0)
                        target2->voxel(voxel_x, voxel_y, voxel_z) = static_cast<uint8_t>(255.0f * (v / max));
                    else
                        target2->voxel(voxel_x, voxel_y, voxel_z) = 0;
                }
            }
        }

        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        Volume h(target2, vec3(1.0f), vec3(0.0f));
        oldVolumePosition(&h);
        serializer->write(parameters.back(), &h);
        //delete target2;
    }
    else if (operation == "blobs2") {
        delete target;
        target = 0;
        VolumeRAM_UInt8* target2 = new VolumeRAM_UInt8(dimensions);

        tgt::vec3 center = dimensions / 2;

        LINFO("Generating blobs2 dataset with dimensions: " << dimensions);

        //clear dataset:
        fillBox(target2, tgt::ivec3(0,0,0), dimensions, 0);
        tgt::vec3 blob1 = tgt::vec3(center.x * 0.9f, center.y * 0.8f, dimensions.z * 0.3f);
        tgt::vec3 blob2 = tgt::vec3(center.x * 1.1f, center.y * 1.2f, dimensions.z * 0.7f);

        float radius = 0.22f / (float)dimensions.x * 32.0f;
        tgt::vec3 m  = blob1 + tgt::vec3(1.0f,1.0f,1.0f);
        float max = (1.0f / distance(blob1, m)) + (1.0f / distance(blob1, m));

        max -= radius;

        for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
            for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
                for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {
                    float b1 = 1.0f / distance(blob1, tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z));
                    float b2 = 1.0f / distance(blob2, tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z));

                    float v = b1+b2;
                    v -= radius;

                    if(v >= max)
                        target2->voxel(voxel_x, voxel_y, voxel_z) = 255;
                    else if(v > 0.0)
                        target2->voxel(voxel_x, voxel_y, voxel_z) = static_cast<uint8_t>(255.0f * (v / max));
                    else
                        target2->voxel(voxel_x, voxel_y, voxel_z) = 0;
                }
            }
        }

        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        Volume h(target2, vec3(1.0f), vec3(0.0f));
        oldVolumePosition(&h);
        serializer->write(parameters.back(), &h);
        //delete target2;
    }
    else if (operation == "blobs3") {
        delete target;
        target = 0;
        VolumeRAM_UInt8* target2 = new VolumeRAM_UInt8(dimensions);

        tgt::vec3 center = dimensions / 2;

        LINFO("Generating blobs3 dataset with dimensions: " << dimensions);

        //clear dataset:
        fillBox(target2, tgt::ivec3(0,0,0), dimensions, 0);
        tgt::vec3 blob1 = tgt::vec3(center.x * 0.9f, center.y * 0.8f, dimensions.z * 0.3f);
        tgt::vec3 blob2 = tgt::vec3(center.x * 1.1f, center.y * 1.2f, dimensions.z * 0.7f);
        tgt::vec3 blob3 = tgt::vec3(dimensions.x * 0.7f, center.y * 1.0f, dimensions.z * 0.5f);

        float radius = 0.32f / (float)dimensions.x * 32.0f;
        tgt::vec3 m  = blob1 + tgt::vec3(1.0f,1.0f,1.0f);
        float max = (1.0f / distance(blob1, m)) + (1.0f / distance(blob1, m)) + (1.0f / distance(blob3, m));

        max -= radius;

        for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
            for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
                for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {
                    float b1 = 1.0f / distance(blob1, tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z));
                    float b2 = 1.0f / distance(blob2, tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z));
                    float b3 = 1.0f / distance(blob3, tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z));

                    float v = b1+b2+b3;
                    v -= radius;

                    if(v >= max)
                        target2->voxel(voxel_x, voxel_y, voxel_z) = 255;
                    else if(v > 0.0)
                        target2->voxel(voxel_x, voxel_y, voxel_z) = static_cast<uint8_t>(255.0f * (v / max));
                    else
                        target2->voxel(voxel_x, voxel_y, voxel_z) = 0;
                }
            }
        }

        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        Volume h(target2, vec3(1.0f), vec3(0.0f));
        serializer->write(parameters.back(), &h);
        //delete target2;
    }
    else if (operation == "sphere") {
        delete target;
        target = 0;
        VolumeRAM_UInt8* target2 = new VolumeRAM_UInt8(dimensions);

        tgt::vec3 center = dimensions / 2;

        LINFO("Generating sphere dataset with dimensions: " << dimensions);

        //clear dataset:
//         fillBox(target2, tgt::ivec3(0,0,0), dimensions, 0);

        float radius = center.x * 0.8f;

        for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
            for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
                for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {
                    tgt::vec3 diff = center - tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z);
                    float r = length(diff);

                    if(r <= (float)radius) {
                        target2->voxel(voxel_x, voxel_y, voxel_z) = 255;
                    }
                    else {
                        target2->voxel(voxel_x, voxel_y, voxel_z) = 0;
                    }
                }
            }
        }

        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        Volume h(target2, vec3(1.0f), vec3(0.0f));
        serializer->write(parameters.back(), &h);
        //delete target2;
    }
    else if (operation == "doublesphere") {
        //two spheres, NOT concentric
        delete target;
        target = 0;
        VolumeRAM_UInt8* target2 = new VolumeRAM_UInt8(dimensions);

        tgt::vec3 center = dimensions / 2;
        tgt::vec3 center2 = center;
        center2.x *= 1.2;

        LINFO("Generating sphere dataset with dimensions: " << dimensions);

        //clear dataset:
//         fillBox(target2, tgt::ivec3(0,0,0), dimensions, 0);

        float radius = center.x * 0.8f;
        float radius2 = center.x * 0.4f;

        for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
            for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
                for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {
                    tgt::vec3 diff = center2 - tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z);
                    float r = length(diff);

                    if(r <= radius2) {
                        target2->voxel(voxel_x, voxel_y, voxel_z) = 255;
                    }
                    else {
                        diff = center - tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z);
                        float r = length(diff);

                        if(r <= radius) {
                            target2->voxel(voxel_x, voxel_y, voxel_z) = 128;
                        }
                        else {
                            target2->voxel(voxel_x, voxel_y, voxel_z) = 0;
                        }
                    }
                }
            }
        }

        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        Volume h(target2, vec3(1.0f), vec3(0.0f));
        serializer->write(parameters.back(), &h);
        //delete target2;
    }
    else if (operation == "spherecoord") {
        //the same sphere but as 32bit DS with natural param.
        delete target;
        target = 0;
        VolumeRAM_4xUInt8* target2 = new VolumeRAM_4xUInt8(dimensions);

        tgt::vec3 center = dimensions / 2;

        LINFO("Generating sphere dataset with dimensions: " << dimensions);

        //clear dataset:
//         fillBox(target2, tgt::ivec3(0,0,0), dimensions, 0);

        float radius = center.x * 0.8f;

        for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
            for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
                for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {
                    tgt::vec3 diff = center - tgt::vec3((float)voxel_x, (float)voxel_y, (float)voxel_z);
                    float r = length(diff);

                    if(r <= (float)radius) {
                        //y,x
                        double theta = atan2(voxel_y-center.y, voxel_x-center.x);
                        theta += tgt::PI;

                        float phi = acos((voxel_z-center.z) / r);

                        target2->voxel(voxel_x, voxel_y, voxel_z).r = static_cast<uint8_t>(theta / (tgt::PI * 2.0) * 255.0);
                        target2->voxel(voxel_x, voxel_y, voxel_z).g = static_cast<uint8_t>(phi / tgt::PI * 255.0);
                        target2->voxel(voxel_x, voxel_y, voxel_z).b = static_cast<uint8_t>((float) r / (float) radius * 255.0f);
                        if(target2->voxel(voxel_x, voxel_y, voxel_z).b == 0)
                            target2->voxel(voxel_x, voxel_y, voxel_z).b = 1;

                        target2->voxel(voxel_x, voxel_y, voxel_z).a = static_cast<uint8_t>(255.0f - (((float) r / (float) radius) * 255.0f));
//                         target2->voxel(voxel_x, voxel_y, voxel_z).a = value;
                    }
                    else {
                        target2->voxel(voxel_x, voxel_y, voxel_z).r = 0;
                        target2->voxel(voxel_x, voxel_y, voxel_z).g = 0;
                        target2->voxel(voxel_x, voxel_y, voxel_z).b = 0;
                        target2->voxel(voxel_x, voxel_y, voxel_z).a = 0;
                    }
                }
            }
        }

        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        Volume h(target2, vec3(1.0f), vec3(0.0f));
        serializer->write(parameters.back(), &h);
        //delete target2;
    }
     else if (operation == "synth") {
        vec3 center;
        center = dimensions / 2;

        float radius = dimensions.x / 4.0f;

        LINFO("Generating synth dataset with dimensions: " << dimensions);

        fillBox(target, ivec3(0,0,0), dimensions, 0);

        fillSphere(target, center, radius, 96);

        /*
        for (int voxel_z=5; voxel_z<dimensions.z-5; voxel_z++) {
            for (int voxel_y=5; voxel_y<dimensions.y-5; voxel_y++) {
                for (int voxel_x=5; voxel_x<(5+(dimensions.x/10)); voxel_x++) {
                   targetDataset_->setVoxel(voxel_x, voxel_y, voxel_z, (uint8_t)192);
                }
            }
        }
        for (int voxel_x=5; voxel_x<dimensions.x-5; voxel_x++) {
            for (int voxel_y=5; voxel_y<dimensions.y-5; voxel_y++) {
                for (int voxel_z=5; voxel_z<(5+(dimensions.z/10)); voxel_z++) {
                   target->voxel(voxel_x, voxel_y, voxel_z) = 128;
                }
            }
        }
    }
    else if (operation == "cloud") {
        dimensions.x *= 2;
        vec3 center;
        center = dimensions / 2;

        LINFO("Generating cloud dataset with dimensions: " << dimensions);

        fillBox(target, ivec3(0,0,0), dimensions, 0);

        // generate cloud pattern
        fillEllipsoid(target, tgt::dvec3(dimensions.x*0.5, dimensions.y*0.5, dimensions.z*0.5), tgt::dvec3(dimensions.x*0.3, dimensions.y*0.35, dimensions.z*0.3), 128);
        fillEllipsoid(target, tgt::dvec3(dimensions.x*0.2, dimensions.y*0.3, dimensions.z*0.45), tgt::dvec3(dimensions.x*0.2, dimensions.y*0.2, dimensions.z*0.15), 128);
        fillEllipsoid(target, tgt::dvec3(dimensions.x*0.8, dimensions.y*0.4, dimensions.z*0.6), tgt::dvec3(dimensions.x*0.15, dimensions.y*0.15, dimensions.z*0.15), 128);
        fillEllipsoid(target, tgt::dvec3(dimensions.x*0.6, dimensions.y*0.3, dimensions.z*0.8), tgt::dvec3(dimensions.x*0.1, dimensions.y*0.15, dimensions.z*0.12), 128);
        fillEllipsoid(target, tgt::dvec3(dimensions.x*0.3, dimensions.y*0.6, dimensions.z*0.55), tgt::dvec3(dimensions.x*0.15, dimensions.y*0.15, dimensions.z*0.2), 128);

        // perform noise perturbation
        vec3 lowFrequency = vec3(10, 10, 10);
        vec3 highAmplitude = vec3(static_cast<float>(dimensions.x)/10.0f, static_cast<float>(dimensions.y)/10.0f, static_cast<float>(dimensions.z)/10.0f);
        LINFO("applying low frequency high amplitude perturbation...");
        applyPerturbation(target, dimensions, lowFrequency, highAmplitude);

        vec3 highFrequency = vec3(1, 1, 1);
        vec3 lowAmplitude = vec3(static_cast<float>(dimensions.x)/20.0f, static_cast<float>(dimensions.y)/20.0f, static_cast<float>(dimensions.z)/20.0f);
        LINFO("applying high frequency low amplitude perturbation...");
        applyPerturbation(target, dimensions, highFrequency, lowAmplitude);
    }
    else if (operation == "aotestbox") {
        LINFO("Generating aotest");

        fillBox(target, ivec3(0,0,0), dimensions, 0);

        for (int i=0 ; i<dimensions.z; ++i) {
            LINFO("Generating aotest slice: "<< i);
            fillBox(target, ivec3(0,0,i), ivec3(dimensions.x,dimensions.y,i+1), i);
        }
    }
    else if (operation == "shadowtestvol") {
        LINFO("Generating shadowTestVolume");
        delete target;
        target = new VolumeRAM_UInt8(ivec3(64,64,64));

        fillBox(target, ivec3(0,0,0), ivec3(64,64,64), 0);

        fillBox(target, ivec3(0,0,0), ivec3(64,64,1), 25);
        fillBox(target, ivec3(0,0,1), ivec3(64,64,2), 25);

        fillBox(target, ivec3(25,25,45), ivec3(39,39,46), 25);
        fillBox(target, ivec3(25,25,46), ivec3(39,39,47), 25);
    }
    else if (operation == "aorticarch") {
        // --------------------
        // ------Settings------
        // --------------------
        // The granularity of the bezier-sweeping
        const int stepping = 50;

        const float arch1Position = 0.5f;
        const float arch2Position = 0.6f;
        const float arch3Position = 0.7f;

        const float aortaWidth = 0.117f * dimensions.x;
        const float aortaThickness = 0.0234f * dimensions.x;
        const float archWidth = 0.039f * dimensions.x;
        const float archThickness = 0.0117f * dimensions.x;

        // A vessel is built from an "interior"-intensity surrounded by a "base"-intensity
        const uint8_t arch1BaseIntensity = 100;
        const uint8_t arch2BaseIntensity = arch1BaseIntensity;
        const uint8_t arch3BaseIntensity = arch1BaseIntensity;
        const uint8_t arch1InteriorIntensity = 200;
        const uint8_t arch2InteriorIntensity = arch1InteriorIntensity;
        const uint8_t arch3InteriorIntensity = arch1InteriorIntensity;
        const uint8_t aortaBaseIntensity = 100;
        const uint8_t aortaInteriorIntensity = 200;

        // If we want a plaque somewhere in the aorta, specify the details here
        const bool plaqueEnable = true;
        const float plaqueBegin = 0.55f;
        const float plaqueEnd = 0.85f;
        const float plaqueThickness = 0.03125f * dimensions.x;
        const uint8_t plaqueIntensity = 1;
        const vec3 plaqueDeviation = vec3(0.0f, 0.0f, 0.0f);//vec3(0.f, 0.f, 0.03125f * dimensions.x);


        // --------------------
        // -------Drawing------
        // --------------------

        // Empty the volume first
        fillBox(target, ivec3(0,0,0), dimensions, 0);

        // Create the main aorta
        std::vector<vec3> aortaPoints;
        aortaPoints.push_back(vec3::zero);
        aortaPoints.push_back(vec3(dimensions.x*0.25f,dimensions.y*0.5f,dimensions.z*0.65f));
        aortaPoints.push_back(vec3(dimensions.x*0.75f,dimensions.y*0.5f,dimensions.z*0.65f));
        aortaPoints.push_back(vec3(static_cast<float>(dimensions.x), static_cast<float>(dimensions.y),0));
        tgt::BSpline aortaSpline(aortaPoints);

        // The first arch
        std::vector<vec3> arch1Points;
        arch1Points.push_back(aortaSpline.getPoint(arch1Position));
        const float arch1pos1 = (arch1Position-0.1f < 0) ? 0 : arch1Position-0.1f;
        const float arch1pos2 = (arch1Position-0.2f < 0) ? 0 : arch1Position-0.2f;
        const float arch1pos3 = (arch1Position-0.3f < 0) ? 0 : arch1Position-0.3f;
        arch1Points.push_back(vec3(aortaSpline.getPoint(arch1pos1).x , aortaSpline.getPoint(arch1pos1).y , 0.4f*dimensions.z));
        arch1Points.push_back(vec3(aortaSpline.getPoint(arch1pos2).x , aortaSpline.getPoint(arch1pos2).y , 0.7f*dimensions.z));
        arch1Points.push_back(vec3(aortaSpline.getPoint(arch1pos3).x , aortaSpline.getPoint(arch1pos3).y , 1.0f*dimensions.z));
        tgt::BSpline arch1Spline(arch1Points);

        // The second arch
        std::vector<vec3> arch2Points;
        arch2Points.push_back(aortaSpline.getPoint(arch2Position));
        arch2Points.push_back(vec3(aortaSpline.getPoint(arch2Position).x , aortaSpline.getPoint(arch2Position).y , 0.4f*dimensions.z));
        arch2Points.push_back(vec3(aortaSpline.getPoint(arch2Position).x , aortaSpline.getPoint(arch2Position).y , 0.6f*dimensions.z));
        arch2Points.push_back(vec3(aortaSpline.getPoint(arch2Position).x , aortaSpline.getPoint(arch2Position).y , 1.f*dimensions.z));
        tgt::BSpline arch2Spline(arch2Points);

        // The third arch
        std::vector<vec3> arch3Points;
        arch3Points.push_back(aortaSpline.getPoint(arch3Position));
        const float arch3pos1 = (arch3Position+0.1f > 1) ? 1 : arch3Position+0.1f;
        const float arch3pos2 = (arch3Position+0.2f > 1) ? 1 : arch3Position+0.2f;
        const float arch3pos3 = (arch3Position+0.3f > 1) ? 1 : arch3Position+0.3f;
        arch3Points.push_back(vec3(aortaSpline.getPoint(arch3pos1).x , aortaSpline.getPoint(arch3pos1).y , 0.4f*dimensions.z));
        arch3Points.push_back(vec3(aortaSpline.getPoint(arch3pos2).x , aortaSpline.getPoint(arch3pos2).y , 0.7f*dimensions.z));
        arch3Points.push_back(vec3(aortaSpline.getPoint(arch3pos3).x , aortaSpline.getPoint(arch3pos3).y , 1.0f*dimensions.z));
        tgt::BSpline arch3Spline(arch3Points);

        // Draw the base interior
        for (int i = 0; i <= stepping; ++i) {
            float progress = static_cast<float>(i)/static_cast<float>(stepping);
            std::cout << progress*50 << "%" << std::endl;
            fillSphere(target, arch1Spline.getPoint(progress), archWidth, arch1BaseIntensity);
            fillSphere(target, arch2Spline.getPoint(progress), archWidth, arch2BaseIntensity);
            fillSphere(target, arch3Spline.getPoint(progress), archWidth, arch3BaseIntensity);
            fillSphere(target, aortaSpline.getPoint(progress), aortaWidth, aortaBaseIntensity);
        }

        std::cout << "Baseline created..." << std::endl << "Starting drilling..." << std::endl;

        // Draw the interior above the base
        for (int i = 0; i <= stepping; ++i) {
            float progress = static_cast<float>(i)/static_cast<float>(stepping);
            std::cout << progress*50+50 << "%" << std::endl;
            fillSphere(target, arch1Spline.getPoint(progress), archWidth - archThickness, arch1InteriorIntensity);
            fillSphere(target, arch2Spline.getPoint(progress), archWidth - archThickness, arch2InteriorIntensity);
            fillSphere(target, arch3Spline.getPoint(progress), archWidth - archThickness, arch3InteriorIntensity);
            if ( plaqueEnable && (progress >= plaqueBegin) && (progress <= plaqueEnd) ) {
                fillSphere(target, aortaSpline.getPoint(progress) + plaqueDeviation, aortaWidth - (aortaThickness + plaqueThickness), aortaInteriorIntensity);
            }
            else
                fillSphere(target, aortaSpline.getPoint(progress), aortaWidth - aortaThickness, aortaInteriorIntensity);
        }

        if (plaqueEnable) {
            VolumeRAM_UInt8* petTarget = new VolumeRAM_UInt8(dimensions);
            fillBox(petTarget, ivec3(0,0,0), dimensions, 0);

            // Draw pet spot
            fillSphere(petTarget, aortaSpline.getPoint((plaqueBegin+plaqueEnd)/2.f) - (plaqueDeviation + plaqueThickness), plaqueThickness - 2.f, plaqueIntensity);

            VolumeSerializerPopulator volLoadPop;
            const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

            Volume h(petTarget, vec3(1.0f), vec3(0.0f));
            serializer->write(parameters.back() + "-pet.dat", &h);
        }

    }

    if (target) {
        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        Volume h(target, vec3(1.0f), vec3(0.0f));
        serializer->write(parameters.back(), &h);
    } */
    return true;
}



void CommandCreate::fillPlane(VolumeRAM_UInt8* vds, vec3 center, vec3 normal, uint8_t value) {
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
        for (size_t voxel_y=0; voxel_y<vds->getDimensions().y; voxel_y++) {
            for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
                vec3 pos = static_cast<vec3>(ivec3(voxel_x, voxel_y, voxel_z));
                vec3 diff = pos - center;

                if (fabsf(dot(normal, diff)) <= 1) {
                    vds->voxel(voxel_x, voxel_y, voxel_z) = value;
                }
            }
        }
    }
}

void CommandCreate::fillCircle(VolumeRAM_UInt8* vds, vec3 center, float radius, uint8_t value) {
    for (size_t voxel_y = 0; voxel_y < vds->getDimensions().y; ++voxel_y) {
        for (size_t voxel_x = 0; voxel_x < vds->getDimensions().x; ++voxel_x) {
            vec3 diff = center - vec3(static_cast<float>(voxel_x), static_cast<float>(voxel_y), center.z);

            if (length(diff) <= radius)
                vds->voxel(voxel_x, voxel_y, static_cast<size_t>(center.z)) = value;
        }
    }
}

void CommandCreate::fillOrientedCircle(VolumeRAM_UInt8* vds, vec3 center, vec3 normal, float radius, uint8_t value) {
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
        for (size_t voxel_y=0; voxel_y<vds->getDimensions().y; voxel_y++) {
            for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
                vec3 pos = static_cast<vec3>(ivec3(voxel_x, voxel_y, voxel_z));
                vec3 diff = pos - center;

                if ( (fabsf(dot(normal, diff)) <= 1) && (length(diff) <= radius) ) {
                    vds->voxel(voxel_x, voxel_y, voxel_z) = value;
                }
            }
        }
    }

}

void CommandCreate::fillSphere(VolumeRAM_UInt8* vds, vec3 center, float radius, uint8_t value) {
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
        for (size_t voxel_y=0; voxel_y<vds->getDimensions().y; voxel_y++) {
            for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
                vec3 diff = center - static_cast<vec3>(ivec3(voxel_x, voxel_y, voxel_z));

                if (length(diff) < radius)
                    vds->voxel(voxel_x, voxel_y, voxel_z) = value;
            }
        }
    }
}

void CommandCreate::fillEllipsoid(VolumeRAM_UInt8* vds, vec3 center, vec3 radius, uint8_t value) {
    // ugly code - don't do this at home
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
        for (size_t voxel_y=0; voxel_y<vds->getDimensions().y; voxel_y++) {
            for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
                if ((((voxel_x-center.x)*(voxel_x-center.x))/((radius.x*radius.x))) +
                    (((voxel_y-center.y)*(voxel_y-center.y))/((radius.y*radius.y))) +
                    (((voxel_z-center.z)*(voxel_z-center.z))/((radius.z*radius.z))) <= 1)
                {
                    vds->voxel(voxel_x, voxel_y, voxel_z) = value;
                }
            }
        }
    }
}

void CommandCreate::fillBox(VolumeRAM_UInt8* vds, ivec3 start, ivec3 end, uint8_t value) {
    ivec3 i;
    for (i.x = start.x; i.x < end.x; i.x++) {
        for (i.y = start.y; i.y < end.y; i.y++) {
            for (i.z = start.z; i.z < end.z; i.z++) {
                vds->voxel(i.x, i.y, i.z) = value;
            }
        }
    }
}

void CommandCreate::fillOrientedBox(VolumeRAM_UInt8* vds, vec3 center, vec3 dir, float lengthA, float lengthB,
                                    float yStart, float yEnd, uint8_t value)
{
    dir = normalize(dir);
    center.y = 0.0f;
    vec3 dir2 = cross(dir, vec3(0.0f, 1.0f, 0.0f));
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
       for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
           vec3 diff = static_cast<vec3>(ivec3(voxel_x, 0, voxel_z)) - center;
           float l = dot(dir, diff);
           float l2 = dot(dir2, diff);
           if ((fabsf(l) < lengthA) && (fabsf(l2) < lengthB))
               for (int voxel_y = static_cast<int>(yStart); static_cast<float>(voxel_y) < yEnd; voxel_y++) {
                   vds->voxel(voxel_x, voxel_y, voxel_z) = value;
               }
       }
    }
}

void CommandCreate::applyPerturbation(VolumeRAM* /*vds*/, ivec3 /*dimensions*/, vec3 /*frequency*/, vec3 /*amplitude*/) {
//     srand( static_cast<unsigned int>(42) );
//     LWARNING("ATTENTION: set random seed to fixed value.");
//
//     vec3 perturbVec;
//     for (int voxel_z=5; voxel_z<dimensions.z-5; voxel_z+=frequency.z) {
//         for (int voxel_y=5; voxel_y<dimensions.y-5; voxel_y+=frequency.y) {
//             for (int voxel_x=5; voxel_x<dimensions.x-5; voxel_x+=frequency.x) {
//                 perturbVec = vec3((rand()/(float)RAND_MAX)*amplitude.x,
//                                        (rand()/(float)RAND_MAX)*amplitude.y,
//                                        (rand()/(float)RAND_MAX)*amplitude.z);
//                 for (int offset_z=0; offset_z<frequency.z; offset_z++) {
//                     for (int offset_y=0; offset_y<frequency.y; offset_y++) {
//                         for (int offset_x=0; offset_x<frequency.x; offset_x++) {
//                             if ((voxel_x+offset_x + perturbVec.x < dimensions.x) &&
//                                 (voxel_y+offset_y + perturbVec.y < dimensions.y) &&
//                                 (voxel_z+offset_z + perturbVec.z < dimensions.z)) {
//                                     uint8_t* oldIntensity = targetDataset_->voxel(voxel_x+offset_x, voxel_y+offset_y, voxel_z+offset_z);
//                                     targetDataset_->setVoxel(voxel_x+offset_x, voxel_y+offset_y, voxel_z+offset_z, targetDataset_->voxel(voxel_x+offset_x+perturbVec.x, voxel_y+offset_y+perturbVec.y, voxel_z+offset_z+perturbVec.z));
//                                     //targetDataset_->setVoxel(voxel_x+offset_x+perturbVec.x, voxel_y+offset_y+perturbVec.y, voxel_z+offset_z+perturbVec.z, oldIntensity);
//                                 }
//                         }
//                     }
//                 }
//             }
//         }
//     }
}

}   //namespace voreen
