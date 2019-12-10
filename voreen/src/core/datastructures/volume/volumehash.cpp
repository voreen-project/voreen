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

#include "voreen/core/datastructures/volume/volumehash.h"

#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/io/serialization/serialization.h"

#include "voreen/core/utils/hashing.h"

namespace voreen {

VolumeHash::VolumeHash() :
    VolumeDerivedData(),
    hash_("")
{}

VolumeHash::VolumeHash(const std::string& hash) :
    VolumeDerivedData()
{
    setHash(hash);
}

VolumeDerivedData* VolumeHash::create() const {
    return new VolumeHash();
}

VolumeDerivedData* VolumeHash::createFrom(const VolumeBase* handle) const {
    tgtAssert(handle, "no volume");

    // prefer disk representation over ram representation
    if (handle->hasRepresentation<VolumeDisk>()) {
        return new VolumeHash(handle->getRepresentation<VolumeDisk>()->getHash());
    }
    else if (handle->hasRepresentation<VolumeRAM>()) { //< do not load volume into ram for hashing
        const VolumeRAM* v = handle->getRepresentation<VolumeRAM>();
        tgtAssert(v, "no volume");

        size_t s = v->getNumVoxels() * v->getBytesPerVoxel();

        std::string h = VoreenHash::getHash(v->getData(), s);
        return new VolumeHash(h);
    }
    else {
        LWARNING("Unable to compute volume hash: neither disk nor ram representation available");
        return new VolumeHash("");
    }
}

void VolumeHash::serialize(Serializer& s) const  {
    s.serialize("hash", hash_);
}

void VolumeHash::deserialize(Deserializer& s) {
    s.deserialize("hash", hash_);
}

} // namespace voreen
