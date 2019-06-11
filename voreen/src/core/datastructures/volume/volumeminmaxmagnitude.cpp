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

#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"

namespace voreen {

VolumeMinMaxMagnitude::VolumeMinMaxMagnitude()
    : VolumeDerivedData()
{}

VolumeMinMaxMagnitude::VolumeMinMaxMagnitude(const VolumeMinMaxMagnitude& other)
    : VolumeDerivedData()
    , minMagnitude_(other.minMagnitude_)
    , maxMagnitude_(other.maxMagnitude_)
{}

VolumeMinMaxMagnitude::VolumeMinMaxMagnitude(float minMag, float maxMag)
    : VolumeDerivedData()
    , minMagnitude_(minMag)
    , maxMagnitude_(maxMag)
{}

VolumeDerivedData* VolumeMinMaxMagnitude::create() const {
    return new VolumeMinMaxMagnitude();
}

VolumeDerivedData* VolumeMinMaxMagnitude::createFrom(const VolumeBase* handle) const {
    tgtAssert(handle, "no volume");

    float minMagnitude = 0.f, maxMagnitude = 1.f;
    if (handle->hasRepresentation<VolumeRAM>()) {
        const VolumeRAM* v = handle->getRepresentation<VolumeRAM>();
        tgtAssert(v, "no volume");

        minMagnitude = v->minMagnitude();
        maxMagnitude = v->maxMagnitude();
    }
    else if (handle->hasRepresentation<VolumeDisk>()) {
        const VolumeDisk* volumeDisk = handle->getRepresentation<VolumeDisk>();
        tgtAssert(volumeDisk, "no disk volume");

        minMagnitude = FLT_MAX;
        maxMagnitude = -FLT_MAX;

        // compute min/max values slice-wise
        size_t numSlices = handle->getDimensions().z;
        tgtAssert(numSlices > 0, "empty volume");
        for (size_t i=0; i<numSlices; i++) {

            // interruption point after each slice!
            boost::this_thread::interruption_point();

            try {
                VolumeRAM* sliceVolume = volumeDisk->loadSlices(i, i);
                minMagnitude = std::min(minMagnitude, sliceVolume->minMagnitude());
                maxMagnitude = std::max(maxMagnitude, sliceVolume->maxMagnitude());
                delete sliceVolume;
            }
            catch (tgt::Exception& e) {
                LWARNING("Unable to compute min/max values: failed to load slice from disk volume: " << e.what());
                break;
            }
        }
    }
    else {
        LWARNING("Unable to compute min/max values: neither disk nor ram representation available");
    }
    tgtAssert(minMagnitude <= maxMagnitude, "invalid min/max values");

    return new VolumeMinMaxMagnitude(minMagnitude,maxMagnitude);
}

float VolumeMinMaxMagnitude::getMinMagnitude() const {
    return minMagnitude_;
}

float VolumeMinMaxMagnitude::getMaxMagnitude() const {
    return maxMagnitude_;
}

void VolumeMinMaxMagnitude::serialize(Serializer& s) const  {
    s.serialize("minMagnitude", minMagnitude_);
    s.serialize("maxMagnitude", maxMagnitude_);
}

void VolumeMinMaxMagnitude::deserialize(Deserializer& s) {
    s.deserialize("minMagnitude", minMagnitude_);
    s.deserialize("maxMagnitude", maxMagnitude_);
}

} // namespace voreen
