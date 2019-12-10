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

#include "volumecenter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <sstream>

namespace voreen {

const std::string VolumeCenter::loggerCat_("voreen.experimental.VolumeCenter");

VolumeCenter::VolumeCenter()
    : VolumeProcessor(),
    inport_(Port::INPORT, "volume"),
    applyDatasetTransformationMatrix_("useDatasetTrafoMatrix", "Apply data set trafo matrix", true),
    centerWorld_("centerWorld", "Center (World Coordinates)", tgt::vec3(0.0f), tgt::vec3(-999.0f), tgt::vec3(999.0f), VALID),
    centerVoxel_("centerVoxel", "Center (Voxel Coordinates)", tgt::vec3(0.0f), tgt::vec3(-999.0f), tgt::vec3(999.0f), VALID),
    centerVolume_("centerVolume", "Center (Volume/Texture Coordinates)", tgt::vec3(0.0f), tgt::vec3(-999.0f), tgt::vec3(999.0f), VALID)
{
    addPort(inport_);

    addProperty(applyDatasetTransformationMatrix_);
    addProperty(centerWorld_);
    addProperty(centerVoxel_);
    addProperty(centerVolume_);

    centerWorld_.setReadOnlyFlag(true);
    centerVoxel_.setReadOnlyFlag(true);
    centerVolume_.setReadOnlyFlag(true);
}

VolumeCenter::~VolumeCenter() {}

void VolumeCenter::process() {
    const VolumeRAM_UInt8* inputVolume = dynamic_cast<const VolumeRAM_UInt8*>(inport_.getData()->getRepresentation<VolumeRAM>());
    if(!inputVolume)
        return;

    tgt::ivec3 dimensions = inputVolume->getDimensions();

    tgt::dvec3 center(0.0);
    long numVoxels = 0;

    for (int voxel_z=0; voxel_z<dimensions.z; voxel_z++) {
        for (int voxel_y=0; voxel_y<dimensions.y; voxel_y++) {
            for (int voxel_x=0; voxel_x<dimensions.x; voxel_x++) {
                uint8_t v = inputVolume->voxel(voxel_x, voxel_y, voxel_z);
                if(v != 0) {
                    center += tgt::dvec3(static_cast<double>(voxel_x), static_cast<double>(voxel_y), static_cast<double>(voxel_z));
                    numVoxels++;
                }
            }
        }
    }

    center /= (double)numVoxels;

    tgt::mat4 voxelToWorld = applyDatasetTransformationMatrix_.get() ? inport_.getData()->getVoxelToWorldMatrix() : inport_.getData()->getVoxelToPhysicalMatrix();

    centerWorld_.set(voxelToWorld * tgt::vec3(static_cast<float>(center.x), static_cast<float>(center.y), static_cast<float>(center.z)));
    centerVoxel_.set(center);
    centerVolume_.set(tgt::vec3((float)center.x / (float)dimensions.x, (float)center.y / (float)dimensions.y, (float)center.z /(float) dimensions.z));
}

}   // namespace
