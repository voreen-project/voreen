/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

    // Clear old data.
    outport_.setData(nullptr);

    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    if(!inputVolume) {
        LERROR("No RAM representation available");
        return;
    }

    const tgt::svec3 oldDims = inputVolume->getDimensions();
    const tgt::svec3 newDims(tgt::max(oldDims));
    if(oldDims == newDims) {
        LINFO("Volume already is a cube");
        outport_.setData(inport_.getData(), false);
        return;
    }

    const tgt::svec3 llf = (newDims - oldDims) / static_cast<size_t>(2);
    const tgt::svec3 urb = (newDims + oldDims) / static_cast<size_t>(2);
    const size_t numChannels = inputVolume->getNumChannels();

    VolumeRAM* outputVolume = inputVolume->createNew(newDims, true);
    if(!outputVolume) {
        LERROR("Output volume could not be created");
        return;
    }

#ifdef VRN_MODULE_OPENMP
    long newDimsZ = static_cast<long>(newDims.z);
    #pragma omp parallel for
    for (long voxel_zP=0; voxel_zP < newDimsZ; voxel_zP++) {
        size_t voxel_z = static_cast<size_t>(voxel_zP);
#else
    for (size_t voxel_z=0; voxel_z < newDims.z; voxel_z++) {
#endif
        for (size_t voxel_y=0; voxel_y < newDims.y; voxel_y++) {
            for (size_t voxel_x=0; voxel_x < newDims.x; voxel_x++) {
                tgt::svec3 pos(voxel_x, voxel_y, voxel_z);
                if (tgt::hor(tgt::lessThan(pos, llf)) || tgt::hor(tgt::greaterThanEqual(pos, urb))) {
                    for(size_t channel=0; channel < numChannels; channel++) {
                        outputVolume->setVoxelNormalized(0.f, pos, channel);
                    }
                }
                else {
                    tgt::ivec3 oldPos = pos - llf;
                    for(size_t channel=0; channel < numChannels; channel++) {
                        outputVolume->setVoxelNormalized(inputVolume->getVoxelNormalized(oldPos), pos);
                    }
                }
            }
        }
    }

    // assign computed volume to outport
    Volume* volume = new Volume(outputVolume, inport_.getData());
    volume->setOffset(volume->getOffset() - (tgt::vec3(llf) * volume->getSpacing()));
    outport_.setData(volume);
}

}   // namespace
