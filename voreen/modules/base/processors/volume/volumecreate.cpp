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

#include "volumecreate.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures//volume/volumeoperator.h"
#include "tgt/bspline.h"

using std::vector;
using tgt::BSpline;
using tgt::svec3;
using tgt::vec2;
using tgt::vec3;

namespace voreen {

const std::string VolumeCreate::loggerCat_("voreen.VolumeCreate");

VolumeCreate::VolumeCreate()
    : VolumeCreateBase()
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , currentSeed_(0)
    , operation_("operation", "Operation")
    , dimension_("dimension", "Dimension", tgt::ivec3(64), tgt::ivec3(2), tgt::ivec3(1024), Processor::VALID)
    , regenerate_("regenerate", "Regenerate Volume", Processor::INVALID_RESULT)
    , numShapes_("numShapes", "Number of shapes", 4, 1, 10)
    , keepCurrentShapes_("keepShapes", "Keep current shapes", false, Processor::INVALID_RESULT)
    , numSubdivisions_("numSubdivisions", "Number of Subdivisions", 1, 1, 1000)
{
    addPort(outport_);

    operation_.addOption("cornell",         "Cornell");
    operation_.addOption("cube",            "Cube");
    operation_.addOption("blobs",           "Blobs");
    operation_.addOption("blobs2",          "Blobs2");
    operation_.addOption("blobs3",          "Blobs3");
    operation_.addOption("sphere",          "Sphere");
    operation_.addOption("torus",           "Torus");
    operation_.addOption("doubletorus",     "Doubletorus");
    operation_.addOption("doublepartialtorus","Doubletorus (partial)");
    operation_.addOption("bumpysphere",     "Bumpy Sphere");
    operation_.addOption("doublesphere",    "Doublesphere");
    operation_.addOption("spherecoord",     "Spherecoord");
    operation_.addOption("synth",           "Synth");
    operation_.addOption("cloud",           "Cloud");
    operation_.addOption("aotestbox",       "Aotestbox");
    // TODO: seems to be bugged, deactivated
    //operation_.addOption("shadowtestvol",   "Shadowtestvol");
    operation_.addOption("aorticarch",      "Aorticarch");
    operation_.addOption("randomshapes",    "Random Shapes");
    operation_.addOption("dect", "Dual-Energy CT");
    operation_.addOption("stack", "Stack of different intensity values");
    addProperty(operation_);

    addProperty(dimension_);
    addProperty(regenerate_);
    addProperty(numShapes_);
    addProperty(keepCurrentShapes_);

    addProperty(numSubdivisions_);
}

void VolumeCreate::process() {
    VolumeRAM* outputVolume = 0;

    std::string operation = operation_.get();

    tgt::ivec3 dimensions(dimension_.get());

    VolumeRAM_UInt8* target = new VolumeRAM_UInt8(dimensions);

    outputVolume = target;

    setProgress(0.f);

    if (operation == "cornell") {
        createCornell(dimensions, target);
    }
    else if (operation == "cube") {
        createCube(dimensions, target);
    }
    else if (operation == "blobs") {
        createBlobs(dimensions, target);
    }
    else if (operation == "blobs2") {
        createBlobs2(dimensions, target);
    }
    else if (operation == "blobs3") {
        createBlobs3(dimensions, target);
    }
    else if (operation == "sphere") {
        createSphere(dimensions, target);
    }
    else if (operation == "doublesphere") {
        createDoubleSphere(dimensions, target);
    }
    else if (operation == "torus") {
        createTorus(dimensions, target);
    }
    else if (operation == "doubletorus") {
        createDoubleTorus(dimensions, target);
    }
    else if (operation == "doublepartialtorus") {
        createDoublePartialTorus(dimensions, target);
    }
    else if (operation == "bumpysphere") {
        createBumpySphere(dimensions, target);
    }
    else if (operation == "spherecoord") {
        delete target;
        VolumeRAM_4xUInt8* target2 = new VolumeRAM_4xUInt8(dimensions);
        createSphereCoord(dimensions, target2);
        outputVolume = target2;
    }
     else if (operation == "synth") {
         createSynth(dimensions, target);
    }
    else if (operation == "cloud") {
        createCloud(dimensions, target);
    }
    else if (operation == "aotestbox") {
        createAOTestBox(dimensions, target);
    }
    else if (operation == "shadowtestvol") {
        createShadowTestVolume(dimensions, target);
    }
    else if (operation == "aorticarch") {
        createAorticArch(dimensions, target);
    }
    else if (operation == "randomshapes") {
        createRandomShapes(dimensions, target);
    }
    else if (operation == "dect") {
        createDECT(dimensions, target);
    }
    else if (operation == "stack") {
        createStack(dimensions, target);
    }
    else {
        tgtAssert(false, "Should not get here");
    }

    setProgress(1.f);
    Volume* vh = new Volume(outputVolume, vec3(1.0f), vec3(0.0f));
    oldVolumePosition(vh);

    outport_.setData(vh);
}

void VolumeCreate::createCornell(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    vec3 center = dimensions / 2;

    int size = dimensions.x;
    int thickness = size / 15;
    int border = 3;
    int innerSize = size - ((thickness+border)*2);

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
    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);


    //fill with fog
    //fillBox(target, ivec3(border+thickness, border+thickness, border+thickness),
    //ivec3(border+thickness+is, border+thickness+is, border+thickness+is), fog);



    //white floor
    fillBox(target, tgt::ivec3(border, border, border), tgt::ivec3(size-border, border+thickness, size-border), white);
    //white back wall:
    fillBox(target, tgt::ivec3(border, border, border), tgt::ivec3(size-border, size-border, border+thickness), white);
    //white ceiling:
    fillBox(target, tgt::ivec3(border, size-border-thickness, border), tgt::ivec3(size-border, size-border, size-border), white);

    //THE NEXT 2 ARE OPTIONAL:
    //white front wall:
    //         fillBox(target, ivec3(border, border, border+thickness+is), ivec3(s-border, s-border, s-border), white);
    //cut hole in ceiling (lightsource)
    //         fillBox(target, ivec3(center.x-thickness*2, s-border-thickness, center.x-thickness*2), ivec3(center.x+thickness*2, s-border, center.x+thickness*2), light);
    fillBox(target,
        tgt::ivec3(static_cast<int>(center.x-thickness*2), size-border-thickness, static_cast<int>(center.x-thickness*2)),
        tgt::ivec3(static_cast<int>(center.x+thickness*2), size-border, static_cast<int>(center.x+thickness*2)),
        0);

    //green right wall:
    fillBox(target, tgt::ivec3(size-thickness-border,thickness+border,thickness+border),
        tgt::ivec3(size-border, size-thickness-border, size-border), green);

    //red right wall:
    fillBox(target, tgt::ivec3(border, thickness+border, thickness+border),
        tgt::ivec3(border+thickness, size-thickness-border, size-border), red);

    //box0: oriented box, left, back
    fillOrientedBox(target, vec3(static_cast<float>(border+thickness+(innerSize*0.3f)), 0.0f, static_cast<float>(border+thickness+(innerSize*0.3f))),
        vec3(1.0f, 0.0f, 2.0f), static_cast<float>(innerSize*0.15), static_cast<float>(innerSize*0.2), static_cast<float>(border+thickness),
        static_cast<float>(border+thickness+(innerSize*0.6f)), box0);

    //box1: oriented box, right, front
    fillOrientedBox(target, vec3(static_cast<float>(border+thickness+(innerSize*0.7f)), 0.0f, static_cast<float>(border+thickness+(innerSize*0.8f))),
        vec3(2.0f, 0.0f, 1.0f), static_cast<float>(innerSize*0.15), static_cast<float>(innerSize*0.15), static_cast<float>(border+thickness),
        static_cast<float>(border+thickness+(innerSize*0.3f)), box1);

    //sphere, left, front
    fillSphere(target, vec3(static_cast<float>(border+thickness+(innerSize*0.3f)), static_cast<float>(border+thickness+innerSize/6),
        static_cast<float>(border+thickness+(innerSize*0.8f))), static_cast<float>(innerSize/6), sphere);
}

void VolumeCreate::createCube(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    LINFO("Generating cube dataset with dimensions: " << dimensions);

    int border = 1;

    uint8_t white = 255;

    //clear dataset:
    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);

    fillBox(target, tgt::ivec3(border), dimensions - tgt::ivec3(border), white);
}

void VolumeCreate::createBlobs(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target ) {
    tgt::vec3 center = dimensions / 2;

    LINFO("Generating blobs dataset with dimensions: " << dimensions);

    //clear dataset:
    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);
    tgt::vec3 blob1 = tgt::vec3(center.x, center.y, dimensions.z * 0.25f);
    tgt::vec3 blob2 = tgt::vec3(center.x, center.y, dimensions.z * 0.75f);

    float radius = 0.22f / (float)dimensions.x * 32.0f;
    tgt::vec3 m  = blob1 + tgt::vec3(1.0f,1.0f,1.0f);
    float max = (1.0f / distance(blob1, m)) + (1.0f / distance(blob1, m));

    max -= radius;

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        float b1 = 1.0f / distance(blob1, vox);
        float b2 = 1.0f / distance(blob2, vox);

        float v = b1 + b2;
        v -= radius;

        if (v >= max)
            target->voxel(pos) = 255;
        else if (v > 0.f)
            target->voxel(pos) = static_cast<uint8_t>(255.f * (v / max));
        else
            target->voxel(pos) = 0;

    }
}

void VolumeCreate::createBlobs2(const tgt::ivec3 &dimensions, VolumeRAM_UInt8* target) {
    tgt::vec3 center = dimensions / 2;

    LINFO("Generating blobs2 dataset with dimensions: " << dimensions);

    //clear dataset:
    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);
    tgt::vec3 blob1 = tgt::vec3(center.x * 0.9f, center.y * 0.8f, dimensions.z * 0.3f);
    tgt::vec3 blob2 = tgt::vec3(center.x * 1.1f, center.y * 1.2f, dimensions.z * 0.7f);

    float radius = 0.22f / (float)dimensions.x * 32.0f;
    tgt::vec3 m  = blob1 + tgt::vec3(1.0f,1.0f,1.0f);
    float max = (1.0f / distance(blob1, m)) + (1.0f / distance(blob1, m));

    max -= radius;

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        float b1 = 1.0f / distance(blob1, vox);
        float b2 = 1.0f / distance(blob2, vox);

        float v = b1 + b2;
        v -= radius;

        if (v >= max)
            target->voxel(pos) = 255;
        else if (v > 0.f)
            target->voxel(pos) = static_cast<uint8_t>(255.f * (v / max));
        else
            target->voxel(pos) = 0;
    }
}

void VolumeCreate::createBlobs3(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    tgt::vec3 center = dimensions / 2;

    LINFO("Generating blobs3 dataset with dimensions: " << dimensions);

    //clear dataset:
    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);
    vec3 blob1 = vec3(center.x * 0.9f, center.y * 0.8f, dimensions.z * 0.3f);
    vec3 blob2 = vec3(center.x * 1.1f, center.y * 1.2f, dimensions.z * 0.7f);
    vec3 blob3 = vec3(dimensions.x * 0.7f, center.y * 1.0f, dimensions.z * 0.5f);

    float radius = 0.32f / (float)dimensions.x * 32.0f;
    tgt::vec3 m  = blob1 + tgt::vec3(1.0f,1.0f,1.0f);
    float max = (1.0f / distance(blob1, m)) + (1.0f / distance(blob1, m)) + (1.0f / distance(blob3, m));

    max -= radius;
    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        float b1 = 1.0f / distance(blob1, vox);
        float b2 = 1.0f / distance(blob2, vox);
        float b3 = 1.0f / distance(blob3, vox);

        float v = b1 + b2 + b3;
        v -= radius;

        if (v >= max)
            target->voxel(pos) = 255;
        else if (v > 0.f)
            target->voxel(pos) = static_cast<uint8_t>(255.f * (v / max));
        else
            target->voxel(pos) = 0;
    }
}

void VolumeCreate::createSphere(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    tgt::vec3 center = dimensions / 2;

    LINFO("Generating sphere dataset with dimensions: " << dimensions);

    float radius = center.x * 0.8f;

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        float r = length(center - vox);

        if (r <= radius)
            target->voxel(pos) = 255;
        else
            target->voxel(pos) = 0;
    }
}

void VolumeCreate::createDoubleSphere(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    tgt::vec3 center = dimensions / 2;
    tgt::vec3 center2 = center;
    center2.x *= 1.2;

    LINFO("Generating sphere dataset with dimensions: " << dimensions);

    float radius = center.x * 0.8f;
    float radius2 = center.x * 0.4f;

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        float r = length(center2 - vox);

        if (r <= radius2)
            target->voxel(pos) = 255;
        else {
            r = length(center2 - vox);

            if(r <= radius)
                target->voxel(pos) = 128;
            else
                target->voxel(pos) = 0;
        }
    }
}

void VolumeCreate::createTorus(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    tgt::vec3 center = dimensions / 2;
    tgt::vec3 center2 = center;

    LINFO("Generating torus dataset with dimensions: " << dimensions);

    float radius = center.x * 0.55f;
    float radius2 = center.x * 0.4f;

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        vec3 diff = center2 - vox;
        float r = pow(radius - sqrt((diff.x * diff.x) + (diff.y * diff.y)), 2.0f) + (diff.z*diff.z);

        if (r <= (radius2*radius2))
            target->voxel(pos) = 255;
        else
            target->voxel(pos) = 0;
    }
}

void VolumeCreate::createDoubleTorus(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    tgt::vec3 center = dimensions / 2;
    tgt::vec3 center2 = center;

    LINFO("Generating double torus dataset with dimensions: " << dimensions);

    float radius = center.x * 0.55f;
    float radius2 = center.x * 0.4f;
    float radius3 = center.x * 0.2f;

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        vec3 diff = center2 - vox;
        float r = pow(radius - sqrt((diff.x * diff.x) + (diff.y * diff.y)), 2.0f) + (diff.z*diff.z);

        if (r <= (radius2*radius2)) {
            if (r <= (radius3*radius3))
                target->voxel(pos) = 255;
            else
                target->voxel(pos) = 128;
        }
        else
            target->voxel(pos) = 0;
    }
}

void VolumeCreate::createDoublePartialTorus(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    tgt::vec3 center = dimensions / 2;
    tgt::vec3 center2 = center;

    LINFO("Generating double partial torus dataset with dimensions: " << dimensions);

    float radius = center.x * 0.55f;
    float radius2 = center.x * 0.4f;
    float radius3 = center.x * 0.2f;

    float angleTh = numShapes_.get() / (float)numShapes_.getMaxValue(); //using numShapes to avoid separate property

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        vec3 diff = center2 - vox;

        if ((diff.x < 0.f) || (fabs(diff.y) >= 1.f)) {
            float r = pow(radius - sqrt((diff.x * diff.x) + (diff.y * diff.y)), 2.0f) + (diff.z*diff.z);

            if (r <= (radius2*radius2)) {
                if (r <= (radius3*radius3)) {
                    float a = 1.0f - (atan2(fabs(diff.y), diff.x) / tgt::PIf);

                    if ((a < angleTh))
                        target->voxel(pos) = 255;
                    else
                        target->voxel(pos) = 128;
                }
                else
                    target->voxel(pos) = 128;
            }
            else
                target->voxel(pos) = 0;
        }
        else
            target->voxel(pos) = 0;
    }
}

void VolumeCreate::createBumpySphere(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    tgt::vec3 center = dimensions / 2;

    LINFO("Generating bumpysphere dataset with dimensions: " << dimensions);

    //clear dataset:
    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);

    float offset = 0.5f * numShapes_.get() / (float)numShapes_.getMaxValue(); //using numShapes to avoid separate property

    vec3 blob1 = center; //main blob (sphere)
    vec3 blob2 = vec3(center.x, dimensions.y * 0.9f, center.z);
    vec3 blob3 = vec3(dimensions.y * (0.5f - offset), center.y, center.z);
    vec3 blob4 = vec3(dimensions.y * (0.5f + offset), center.y, center.z);

    float radius = 0.32f / static_cast<float>(dimensions.x) * 32.0f;
    vec3 m = blob1 + vec3(1.0f,1.0f,1.0f);
    float max = (1.f / distance(blob1, m)) + (1.f / distance(blob1, m)) + (1.f / distance(blob3, m));

    max -= radius;

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        float b1 = 1.f / distance(blob1, vox);
        float b2 = 1.f / distance(blob2, vox);
        float b3 = 1.f / distance(blob3, vox);
        float b4 = 1.f / distance(blob4, vox);

        float v = b1 + b2 + b3 + b4;
        v -= radius;

        if (v >= max)
            target->voxel(pos) = 255;
        else if (v > 0.0)
            target->voxel(pos) = static_cast<uint8_t>(255.0f * (v / max));
        else
            target->voxel(pos) = 0;
    }
}

void VolumeCreate::createSphereCoord(const tgt::ivec3& dimensions, VolumeRAM_4xUInt8* target) {
    //the same sphere but as 32bit DS with natural param.
    tgt::vec3 center = dimensions / 2;

    LINFO("Generating sphere dataset with dimensions: " << dimensions);

    float radius = center.x * 0.8f;

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(0,0,0), svec3(dimensions), this) {
        vec3 vox(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));
        float r = length(center - vox);

        if (r <= radius) {
            //y,x
            float theta = atan2(vox.y-center.y, vox.x-center.x);
            theta += tgt::PIf;

            float phi = acos((vox.z-center.z) / r);

            target->voxel(pos).r = static_cast<uint8_t>(theta / (tgt::PIf * 2.f) * 255.f);
            target->voxel(pos).g = static_cast<uint8_t>(phi / tgt::PIf * 255.f);
            target->voxel(pos).b = static_cast<uint8_t>(r / radius * 255.f);
            if (target->voxel(pos).b == 0)
                target->voxel(pos).b = 1;

            target->voxel(pos).a = static_cast<uint8_t>(255.f - (r / radius * 255.f));
        }
        else {
            target->voxel(pos).r = 0;
            target->voxel(pos).g = 0;
            target->voxel(pos).b = 0;
            target->voxel(pos).a = 0;
        }
    }
}

void VolumeCreate::createSynth(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    vec3 center;
    center = dimensions / 2;

    float radius = dimensions.x / 4.0f;

    LINFO("Generating synth dataset with dimensions: " << dimensions);

    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);

    fillSphere(target, center, radius, 96);

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, svec3(5,5,5), svec3(dimensions.x - 5, dimensions.y - 5, 5+dimensions.z/10), this) {
        target->voxel(pos) = 128;
    }
}

void VolumeCreate::createCloud(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    // TODO: this might be a dirty hack to just remove this line
    // there was an out of bounds error
    //dimensions.x *= 2;
    vec3 center;
    center = dimensions / 2;

    LINFO("Generating cloud dataset with dimensions: " << dimensions);

    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);

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

void VolumeCreate::createAOTestBox(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    LINFO("Generating aotest");

    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);

    for (int i=0 ; i<dimensions.z; ++i) {
        LINFO("Generating aotest slice: "<< i);
        fillBox(target, tgt::ivec3(0,0,i), tgt::ivec3(dimensions.x,dimensions.y,i+1), i);
    }
}

void VolumeCreate::createShadowTestVolume(const tgt::ivec3& /*dimensions*/, VolumeRAM_UInt8* target) {
    LINFO("Generating shadowTestVolume");
    delete target;
    target = new VolumeRAM_UInt8(tgt::ivec3(64,64,64));

    fillBox(target, tgt::ivec3(0,0,0), tgt::ivec3(64,64,64), 0);

    fillBox(target, tgt::ivec3(0,0,0), tgt::ivec3(64,64,1), 25);
    fillBox(target, tgt::ivec3(0,0,1), tgt::ivec3(64,64,2), 25);

    fillBox(target, tgt::ivec3(25,25,45), tgt::ivec3(39,39,46), 25);
    fillBox(target, tgt::ivec3(25,25,46), tgt::ivec3(39,39,47), 25);
}

void VolumeCreate::createAorticArch(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
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
    //const uint8_t plaqueIntensity = 1;
    const vec3 plaqueDeviation = vec3(0.0f, 0.0f, 0.0f);//vec3(0.f, 0.f, 0.03125f * dimensions.x);


    // --------------------
    // -------Drawing------
    // --------------------

    // Empty the volume first
    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);

    // Create the main aorta
    vector<vec3> aortaPoints;
    aortaPoints.push_back(vec3::zero);
    aortaPoints.push_back(vec3(dimensions.x*0.25f,dimensions.y*0.5f,dimensions.z*0.65f));
    aortaPoints.push_back(vec3(dimensions.x*0.75f,dimensions.y*0.5f,dimensions.z*0.65f));
    aortaPoints.push_back(vec3(static_cast<float>(dimensions.x), static_cast<float>(dimensions.y),0));
    BSpline aortaSpline(aortaPoints);

    // The first arch
    vector<vec3> arch1Points;
    arch1Points.push_back(aortaSpline.getPoint(arch1Position));
    const float arch1pos1 = (arch1Position-0.1f < 0) ? 0 : arch1Position-0.1f;
    const float arch1pos2 = (arch1Position-0.2f < 0) ? 0 : arch1Position-0.2f;
    const float arch1pos3 = (arch1Position-0.3f < 0) ? 0 : arch1Position-0.3f;
    arch1Points.push_back(vec3(aortaSpline.getPoint(arch1pos1).x , aortaSpline.getPoint(arch1pos1).y , 0.4f*dimensions.z));
    arch1Points.push_back(vec3(aortaSpline.getPoint(arch1pos2).x , aortaSpline.getPoint(arch1pos2).y , 0.7f*dimensions.z));
    arch1Points.push_back(vec3(aortaSpline.getPoint(arch1pos3).x , aortaSpline.getPoint(arch1pos3).y , 1.0f*dimensions.z));
    BSpline arch1Spline(arch1Points);

    // The second arch
    vector<vec3> arch2Points;
    arch2Points.push_back(aortaSpline.getPoint(arch2Position));
    arch2Points.push_back(vec3(aortaSpline.getPoint(arch2Position).x , aortaSpline.getPoint(arch2Position).y , 0.4f*dimensions.z));
    arch2Points.push_back(vec3(aortaSpline.getPoint(arch2Position).x , aortaSpline.getPoint(arch2Position).y , 0.6f*dimensions.z));
    arch2Points.push_back(vec3(aortaSpline.getPoint(arch2Position).x , aortaSpline.getPoint(arch2Position).y , 1.f*dimensions.z));
    BSpline arch2Spline(arch2Points);

    // The third arch
    vector<vec3> arch3Points;
    arch3Points.push_back(aortaSpline.getPoint(arch3Position));
    const float arch3pos1 = (arch3Position+0.1f > 1) ? 1 : arch3Position+0.1f;
    const float arch3pos2 = (arch3Position+0.2f > 1) ? 1 : arch3Position+0.2f;
    const float arch3pos3 = (arch3Position+0.3f > 1) ? 1 : arch3Position+0.3f;
    arch3Points.push_back(vec3(aortaSpline.getPoint(arch3pos1).x , aortaSpline.getPoint(arch3pos1).y , 0.4f*dimensions.z));
    arch3Points.push_back(vec3(aortaSpline.getPoint(arch3pos2).x , aortaSpline.getPoint(arch3pos2).y , 0.7f*dimensions.z));
    arch3Points.push_back(vec3(aortaSpline.getPoint(arch3pos3).x , aortaSpline.getPoint(arch3pos3).y , 1.0f*dimensions.z));
    BSpline arch3Spline(arch3Points);

    // Draw the base interior
    for (int i = 0; i <= stepping; ++i) {
        float progress = static_cast<float>(i)/static_cast<float>(stepping);
        fillSphere(target, arch1Spline.getPoint(progress), archWidth, arch1BaseIntensity);
        fillSphere(target, arch2Spline.getPoint(progress), archWidth, arch2BaseIntensity);
        fillSphere(target, arch3Spline.getPoint(progress), archWidth, arch3BaseIntensity);
        fillSphere(target, aortaSpline.getPoint(progress), aortaWidth, aortaBaseIntensity);
        setProgress(progress /2.f);
    }

    std::cout << "Baseline created..." << std::endl << "Starting drilling..." << std::endl;

    // Draw the interior above the base
    for (int i = 0; i <= stepping; ++i) {
        float progress = static_cast<float>(i)/static_cast<float>(stepping);
        fillSphere(target, arch1Spline.getPoint(progress), archWidth - archThickness, arch1InteriorIntensity);
        fillSphere(target, arch2Spline.getPoint(progress), archWidth - archThickness, arch2InteriorIntensity);
        fillSphere(target, arch3Spline.getPoint(progress), archWidth - archThickness, arch3InteriorIntensity);
        if (plaqueEnable && (progress >= plaqueBegin) && (progress <= plaqueEnd))
            fillSphere(target, aortaSpline.getPoint(progress) + plaqueDeviation, aortaWidth - (aortaThickness + plaqueThickness), aortaInteriorIntensity);
        else
            fillSphere(target, aortaSpline.getPoint(progress), aortaWidth - aortaThickness, aortaInteriorIntensity);
        setProgress(0.5f + progress / 2.f);
    }

    if (plaqueEnable) {
        // TODO what to do with the pet volume? Removed until solution is found
        /*
        VolumeRAM_UInt8* petTarget = new VolumeRAM_UInt8(dimensions);
        fillBox(petTarget, ivec3(0,0,0), dimensions, 0);

        // Draw pet spot
        fillSphere(petTarget, aortaSpline.getPoint((plaqueBegin+plaqueEnd)/2.f) - (plaqueDeviation + plaqueThickness), plaqueThickness - 2.f, plaqueIntensity);

        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

        serializer->save(parameters.back() + "-pet.dat", petTarget);
        delete petTarget;
        */
    }
}

void VolumeCreate::createRandomShapes(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    int s = dimensions.x;
    //int thickness = s / 15;
    int thickness = s / 30;
    int border = 3;
    int is = s - ((thickness+border)*2); //the inner size of the box

    //uint8_t white = 200;
    //uint8_t red = 70;
    //uint8_t green = 150;

    //uint8_t box0 = 100;
    //uint8_t box1 = 50;
    //uint8_t sphere = 120;

    //         uint8_t light = 240;
    //         uint8_t fog = 20;

    //         float radius = dimensions.x / 4;

    LINFO("Generating dataset containing random shapes with dimensions: " << dimensions);

    //clear dataset:
    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);

    //fill with fog
    /*
    fillBox(target, ivec3(border+thickness, border+thickness, border+thickness),
    ivec3(border+thickness+is, border+thickness+is, border+thickness+is), fog);
    */

    //white floor
    fillBox(target, tgt::ivec3(border, border, border), tgt::ivec3(s-border, border+thickness, s-border), 255);
    //white back wall:
    fillBox(target, tgt::ivec3(border, border, border), tgt::ivec3(s-border, s-border, border+thickness), 255);

    //cut hole in ceiling (lightsource)
    //         fillBox(target, ivec3(center.x-thickness*2, s-border-thickness, center.x-thickness*2), ivec3(center.x+thickness*2, s-border, center.x+thickness*2), light);

    int numshapes = numShapes_.get();
    const float max = float(RAND_MAX);

    if (!keepCurrentShapes_.get() || currentSeed_ == 0)
        currentSeed_ = static_cast<unsigned int>(time(0));

    srand(currentSeed_);

    for (int i = 0; i < numshapes; ++i) {
        int color = (i + 1) * (255 / numshapes);

        int r = int(rand()) % 3;

        if (r == 0) {
            float r0 = static_cast<float>(rand()) / (max * 1.2f);
            float r1 = static_cast<float>(rand()) / (max * 1.2f);
            float r2 = static_cast<float>(rand()) / (max * 1.2f);
            float r3 = static_cast<float>(rand()) / max;
            float r4 = static_cast<float>(rand()) / max;
            float r5 = std::max(static_cast<float>(rand()) / (max * 4.0f), 0.1f);
            float r6 = std::max(static_cast<float>(rand()) / (max * 4.0f), 0.1f);
            float r7 = static_cast<float>(rand()) / (max * 1.5f);
            float r8 = static_cast<float>(rand()) / (max * 1.5f);

            float sign1 = (static_cast<float>(rand()) / max) > 0.5f ? 1.f : -1.f;
            float sign2 = (static_cast<float>(rand()) / max) > 0.5f ? 1.f : -1.f;

            bool end = false;
            do {
                vec2 corner = vec2(r0, r2) + normalize(vec2(sign1 * r3, sign2 * r4)) * vec2(r5, r6);
                if (corner.x > 0.95f || corner.y > 0.95f || corner.x < 0.05f || corner.y < 0.05f) {
                    r0 = static_cast<float>(rand()) / (max * 1.2f);
                    r2 = static_cast<float>(rand()) / (max * 1.2f);
                    r3 = static_cast<float>(rand()) / max;
                    r4 = static_cast<float>(rand()) / max;
                    r5 = static_cast<float>(rand()) / (max * 3.0f);
                    r6 = static_cast<float>(rand()) / (max * 3.0f);
                }
                else
                    end = true;
            } while(!end);

            fillOrientedBox(target,
                vec3(border+thickness+is*r0, border+thickness+r1*is, border+thickness+is*r2),
                vec3(sign1*r3, 0.f, sign2*r4),
                is*r5,
                is*r6,
                border+thickness+is*r7,
                border+thickness+is*(std::min(r7+r8, 0.9f)),
                color);
        } else if (r == 1) {
            float r0 = static_cast<float>(rand()) / (max * 1.2f) + 0.15f;
            float r1 = static_cast<float>(rand()) / (max * 1.2f) + 0.15f;
            float r2 = static_cast<float>(rand()) / (max * 1.2f) + 0.15f;
            float r3 = std::max(static_cast<float>(rand()) / (max * 4.f), 0.1f);
            float r4 = std::max(static_cast<float>(rand()) / (max * 4.f), 0.1f);
            float r5 = std::max(static_cast<float>(rand()) / (max * 4.f), 0.1f);
            fillEllipsoid(target,
                vec3(border+thickness+is*r0, border+thickness+r1*is, border+thickness+is*r2),
                vec3(dimensions.x*r3, dimensions.y*r4, dimensions.z*r5),
                color);
        } else if (r == 2) {
            float r0 = static_cast<float>(rand()) / (max * 1.2f);
            float r1 = static_cast<float>(rand()) / (max * 1.2f);
            float r2 = static_cast<float>(rand()) / (max * 1.2f);
            float r3 = static_cast<float>(rand()) / max;
            float r4 = static_cast<float>(rand()) / max;
            float r5 = static_cast<float>(rand()) / (max * 4.0f);
            float r6 = static_cast<float>(rand()) / max;
            float r7 = static_cast<float>(rand()) / max;
            float r8 = static_cast<float>(rand()) / (max * 4.0f);
            float r9 = static_cast<float>(rand()) / (max * 4.0f);

            if (fabs(r6 - r7) < 0.1f) {
                if (r6 > 0.1f)
                    r6 -= 0.1f;
                if (r7 < 0.9f)
                    r7 += 0.1f;
            }

            float sign1 = (static_cast<float>(rand()) / max) > 0.5f ? 1.f : -1.f;
            float sign2 = (static_cast<float>(rand()) / max) > 0.5f ? 1.f : -1.f;

            fillCone(target,
                vec3(border+thickness+is*r0, border+thickness+r1*is, border+thickness+is*r2),
                vec3(sign1*r3, 0.f, sign2*r4),
                is*r5,
                border+thickness+is*std::min(r6, r7),
                border+thickness+is*std::max(r6, r7),
                is*r8,
                is*r9,
                color);
        }
    }
}

void VolumeCreate::createDECT(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    fillBox(target, tgt::ivec3(0,0,0), dimensions, 0);
}

void VolumeCreate::createStack(const tgt::ivec3& dimensions, VolumeRAM_UInt8* target) {
    int previousHeight = 0;
    int height;
    for (int i = 1; i <= numSubdivisions_.get(); ++i) {
        height = static_cast<int>(static_cast<float>(i) * dimensions.z / numSubdivisions_.get());
        uint8_t color = static_cast<uint8_t>(static_cast<float>(i) / static_cast<float>(numSubdivisions_.get()) * 255.f);
        fillBox(target, tgt::ivec3(0,0,previousHeight), tgt::ivec3(dimensions.x, dimensions.y, height), color);
        previousHeight = height;
    }
}

}   // namespace
