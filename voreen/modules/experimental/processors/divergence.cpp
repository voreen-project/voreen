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

#include "divergence.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

const std::string Divergence::loggerCat_("voreen.Divergence");

Divergence::Divergence()
    : VolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input"),
    outport_(Port::OUTPORT, "volumehandle.output", "volumehandle.output", false)
{
    addPort(inport_);
    addPort(outport_);
}

Divergence::~Divergence() {
}

Processor* Divergence::create() const {
    return new Divergence();
}

template<class T, class U>
VolumeAtomic<U>* calcDivergence(const VolumeAtomic<T>* vol) {
    tgt::ivec3 dim = vol->getDimensions();
    tgt::ivec3 pos;
    VolumeAtomic<U>* out = new VolumeAtomic<U>(dim);

    for (pos.z = 0; pos.z < dim.z; pos.z++) {
        for (pos.y = 0; pos.y < dim.y; pos.y++) {
            for (pos.x = 0; pos.x < dim.x; pos.x++) {
                U div = 0.0;

                if((pos.x>0) && (pos.x<(dim.x-1)))
                    div += vol->voxel(pos.x+1, pos.y, pos.z).x - vol->voxel(pos.x-1, pos.y, pos.z).x;

                if((pos.y>0) && (pos.y<(dim.y-1)))
                    div += vol->voxel(pos.x, pos.y+1, pos.z).y -vol->voxel(pos.x, pos.y-1, pos.z).y;

                if((pos.z>0) && (pos.z<(dim.z-1)))
                    div += vol->voxel(pos.x, pos.y, pos.z+1).z - vol->voxel(pos.x, pos.y, pos.z-1).z;

                //out->voxel(pos) = abs(div) * 10000.0f;
                out->voxel(pos) = -div * 10000.0f;
            }
        }
    }
    return out;
}

void Divergence::process() {
    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    tgtAssert(inputVolume, "no input volume");
    VolumeRAM* outputVolume = 0;

    if (dynamic_cast<const VolumeRAM_3xFloat*>(inputVolume)) {
        outputVolume = calcDivergence<tgt::vec3, float>(static_cast<const VolumeRAM_3xFloat*>(inputVolume));
    }
    else {
        LWARNING("Unsupported volume type!");
    }

    // assign computed volume to outport
    if (outputVolume)
        outport_.setData(new Volume(outputVolume, inport_.getData()));
    else
        outport_.setData(0);
}

}   // namespace
