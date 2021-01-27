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

#include "volumecreatebase.h"

#include "voreen/core/datastructures/volume/volumeram.h"

using tgt::ivec3;
using tgt::length;
using tgt::vec2;
using tgt::vec3;

namespace voreen {

void VolumeCreateBase::fillPlane(VolumeRAM_UInt8* vds, vec3 center, vec3 normal, uint8_t value) {
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
        for (size_t voxel_y=0; voxel_y<vds->getDimensions().y; voxel_y++) {
            for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
                vec3 pos = static_cast<vec3>(tgt::svec3(voxel_x, voxel_y, voxel_z));
                vec3 diff = pos - center;

                if (fabsf(dot(normal, diff)) <= 1) {
                    vds->voxel(voxel_x, voxel_y, voxel_z) = value;
                }
            }
        }
    }
}

void VolumeCreateBase::fillCircle(VolumeRAM_UInt8* vds, vec3 center, float radius, uint8_t value) {
    for (size_t voxel_y = 0; voxel_y < vds->getDimensions().y; ++voxel_y) {
        for (size_t voxel_x = 0; voxel_x < vds->getDimensions().x; ++voxel_x) {
            vec3 diff = center - vec3(static_cast<float>(voxel_x), static_cast<float>(voxel_y), center.z);

            if (length(diff) <= radius)
                vds->voxel(voxel_x, voxel_y, static_cast<size_t>(center.z)) = value;
        }
    }
}

void VolumeCreateBase::fillOrientedCircle(VolumeRAM_UInt8* vds, vec3 center, vec3 normal, float radius, uint8_t value) {
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
        for (size_t voxel_y=0; voxel_y<vds->getDimensions().y; voxel_y++) {
            for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
                vec3 pos = static_cast<vec3>(tgt::svec3(voxel_x, voxel_y, voxel_z));
                vec3 diff = pos - center;

                if ( (fabsf(dot(normal, diff)) <= 1) && (length(diff) <= radius) ) {
                    vds->voxel(voxel_x, voxel_y, voxel_z) = value;
                }
            }
        }
    }

}

void VolumeCreateBase::fillSphere(VolumeRAM_UInt8* vds, vec3 center, float radius, uint8_t value) {
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
        for (size_t voxel_y=0; voxel_y<vds->getDimensions().y; voxel_y++) {
            for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
                vec3 diff = center - static_cast<vec3>(tgt::svec3(voxel_x, voxel_y, voxel_z));

                if (length(diff) < radius)
                    vds->voxel(voxel_x, voxel_y, voxel_z) = value;
            }
        }
    }
}

void VolumeCreateBase::fillEllipsoid(VolumeRAM_UInt8* vds, vec3 center, vec3 radius, uint8_t value) {
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

void VolumeCreateBase::fillBox(VolumeRAM_UInt8* vds, ivec3 start, ivec3 end, uint8_t value) {
    ivec3 i;
    for (i.x = start.x; i.x < end.x; i.x++) {
        for (i.y = start.y; i.y < end.y; i.y++) {
            for (i.z = start.z; i.z < end.z; i.z++) {
                vds->voxel(i.x, i.y, i.z) = value;
            }
        }
    }
}

void VolumeCreateBase::fillOrientedBox(VolumeRAM_UInt8* vds, vec3 center, vec3 dir, float lengthA, float lengthB,
    float yStart, float yEnd, uint8_t value)
{
    dir = normalize(dir);
    center.y = 0.0f;
    vec3 dir2 = cross(dir, vec3(0.0f, 1.0f, 0.0f));
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
        for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
            vec3 diff = static_cast<vec3>(tgt::svec3(voxel_x, 0, voxel_z)) - center;
            float l = dot(dir, diff);
            float l2 = dot(dir2, diff);
            if ((fabsf(l) < lengthA) && (fabsf(l2) < lengthB))
                for (int voxel_y = static_cast<int>(yStart); static_cast<float>(voxel_y) < yEnd; voxel_y++) {
                    vds->voxel(voxel_x, voxel_y, voxel_z) = value;
                }
        }
    }
}

void VolumeCreateBase::fillCone(VolumeRAM_UInt8* vds, vec3 center, vec3 /*dir*/, float /*lengthC*/, float yStart, float yEnd, float topRadius, float bottomRadius, uint8_t value) {
    //dir = normalize(dir);
    //center.y = 0.0f;
    //vec3 dir2 = cross(dir, vec3(0.0f, 1.0f, 0.0f));
    for (size_t voxel_z=0; voxel_z<vds->getDimensions().z; voxel_z++) {
        for (size_t voxel_y = static_cast<size_t>(yStart); static_cast<float>(voxel_y) < yEnd; voxel_y++) {
            for (size_t voxel_x=0; voxel_x<vds->getDimensions().x; voxel_x++) {
                float distXZ = length(vec2(center.x, center.z) - vec2((float)voxel_x, (float)voxel_z));
                float interpol = topRadius > bottomRadius ? ((float(voxel_y) - yStart) / (yEnd - yStart)) : 1.f - ((float(voxel_y) - yStart) / (yEnd - yStart));
                if(distXZ < std::min(topRadius, bottomRadius) + fabs(topRadius - bottomRadius) * interpol)
                    vds->voxel(voxel_x, voxel_y, voxel_z) = value;
            }
        }
    }
}

void VolumeCreateBase::applyPerturbation(VolumeRAM* /*vds*/, ivec3 /*dimensions*/, vec3 /*frequency*/, vec3 /*amplitude*/) {
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

tgt::vec3 VolumeCreateBase::getRandVec() const {
    return tgt::vec3(float(rand()), float(rand()), float(rand())) / (float)RAND_MAX;
}

float VolumeCreateBase::getRandFloat() const {
    return float(rand()) / RAND_MAX;
}

} // namespace
