/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "volumecubify.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

const std::string VolumeCubify::loggerCat_("voreen.base.VolumeCubify");

VolumeCubify::VolumeCubify()
    : CachingVolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input", "Volume Input"),
    outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false),
    enableProcessing_("enableProcessing", "Enable")

{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
}

void VolumeCubify::process() {

    if (!enableProcessing_.get()) {
        outport_.setData(inport_.getData(), false);
        return;
    }

    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    VolumeRAM* outputVolume = 0;

    tgt::ivec3 oldims = inputVolume->getDimensions();
    tgt::ivec3 newdims = tgt::ivec3(tgt::max(oldims));
    tgt::ivec3 llf = (newdims - oldims) / 2;
    tgt::ivec3 urb = (newdims + oldims) / 2;

    if(dynamic_cast<const VolumeRAM_UInt8*>(inputVolume))
        outputVolume = new VolumeRAM_UInt8(newdims);
    else if(dynamic_cast<const VolumeRAM_UInt16*>(inputVolume))
        outputVolume = new VolumeRAM_UInt16(newdims);
    else
        LERROR("Unsupported value for getBitsStored()");

    for (int voxel_z=0; voxel_z < newdims.z; voxel_z++) {
        for (int voxel_y=0; voxel_y < newdims.y; voxel_y++) {
            for (int voxel_x=0; voxel_x < newdims.x; voxel_x++) {

                tgt::ivec3 pos = tgt::ivec3(voxel_x, voxel_y, voxel_z);

                if(tgt::hor(tgt::lessThan(pos, llf)) || tgt::hor(tgt::greaterThanEqual(pos, urb)))
                    outputVolume->setVoxelNormalized(0.f, pos);
                else {
                    tgt::ivec3 oldPos = pos - llf;
                    outputVolume->setVoxelNormalized(inputVolume->getVoxelNormalized(oldPos), pos);
                }
            }
        }
    }

    // assign computed volume to outport
    if (outputVolume) {
        Volume* h = new Volume(outputVolume, inport_.getData());
        h->setOffset(h->getOffset() - (tgt::vec3(llf) * h->getSpacing()));
        outport_.setData(h);
    }
    else
        outport_.setData(0);
}

}   // namespace
