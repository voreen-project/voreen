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

#include "jacobian.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

const std::string Jacobian::loggerCat_("voreen.Jacobian");

Jacobian::Jacobian()
    : VolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input"),
    outport_(Port::OUTPORT, "volumehandle.output", "volumehandle.output", false)
{
    addPort(inport_);
    addPort(outport_);
}

Jacobian::~Jacobian() {
}

Processor* Jacobian::create() const {
    return new Jacobian();
}

template<class T, class U>
VolumeAtomic<U>* calcJacobian(const VolumeAtomic<T>* vol) {
    tgt::ivec3 dim = vol->getDimensions();
    tgt::ivec3 pos;
    VolumeAtomic<U>* out = new VolumeAtomic<U>(dim);

    T cv[6];
    for (pos.z = 0; pos.z < dim.z; pos.z++) {
        for (pos.y = 0; pos.y < dim.y; pos.y++) {
            for (pos.x = 0; pos.x < dim.x; pos.x++) {
                U jac;

                if((pos.x == 0) || (pos.x == dim.x-1)
                   || (pos.y == 0) || (pos.y == dim.y-1)
                   || (pos.z == 0) || (pos.z == dim.z-1))
                    jac = U::identity;
                else {
                    cv[0] = vol->voxel(pos.x + 1, pos.y, pos.z);
                    cv[1] = vol->voxel(pos.x - 1, pos.y, pos.z);
                    cv[2] = vol->voxel(pos.x, pos.y + 1, pos.z);
                    cv[3] = vol->voxel(pos.x, pos.y - 1, pos.z);
                    cv[4] = vol->voxel(pos.x, pos.y, pos.z + 1);
                    cv[5] = vol->voxel(pos.x, pos.y, pos.z - 1);

                    jac[0][0] = (cv[0].x - cv[1].x);
                    jac[0][1] = (cv[2].x - cv[3].x);
                    jac[0][2] = (cv[4].x - cv[5].x);

                    jac[1][0] = (cv[0].y - cv[1].y);
                    jac[1][1] = (cv[2].y - cv[3].y);
                    jac[1][2] = (cv[4].y - cv[5].y);

                    jac[2][0] = (cv[0].z - cv[1].z);
                    jac[2][1] = (cv[2].z - cv[3].z);
                    jac[2][2] = (cv[4].z - cv[5].z);
                }

                out->voxel(pos) = jac;
            }
        }
    }
    return out;
}

void Jacobian::process() {
    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    tgtAssert(inputVolume, "no input volume");
    VolumeRAM* outputVolume = 0;

    if (dynamic_cast<const VolumeRAM_3xFloat*>(inputVolume)) {
        outputVolume = calcJacobian<tgt::vec3, tgt::mat3>(static_cast<const VolumeRAM_3xFloat*>(inputVolume));
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
