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

VolumeMinMaxMagnitude* createFromVolumeRAM(const VolumeRAM* volume, RealWorldMapping rwm) {

    float minMagnitude = 0.f, maxMagnitude = 0.f;
    float minNormMagnitude = 0.f, maxNormMagnitude = 0.f;

    for (size_t i=0; i<volume->getNumVoxels(); i++) {
        float magnitude = 0.0f, normalizedMagnitude = 0.0f;
        for(size_t c=0; c<volume->getNumChannels(); c++) {

            // Update normalized magnitude.
            float value = volume->getVoxelNormalized(i, c);
            normalizedMagnitude += value * value;

            // Update real world magnitude.
            value = rwm.normalizedToRealWorld(value);
            magnitude += value * value;
        }

        minMagnitude = std::min(minMagnitude, magnitude);
        maxMagnitude = std::max(maxMagnitude, magnitude);

        minNormMagnitude = std::min(minNormMagnitude, normalizedMagnitude);
        maxNormMagnitude = std::max(maxNormMagnitude, normalizedMagnitude);
    }

    minMagnitude = std::sqrt(minMagnitude);
    maxMagnitude = std::sqrt(maxMagnitude);
    tgtAssert(minMagnitude <= maxMagnitude, "invalid min/max magnitude values");

    minNormMagnitude = std::sqrt(minNormMagnitude);
    maxNormMagnitude = std::sqrt(maxNormMagnitude);
    tgtAssert(minNormMagnitude <= maxNormMagnitude, "invalid min/max normalized magnitude values");

    return new VolumeMinMaxMagnitude(minMagnitude, maxMagnitude, minNormMagnitude, maxNormMagnitude);
}

VolumeMinMaxMagnitude::VolumeMinMaxMagnitude()
    : VolumeDerivedData()
{}

VolumeMinMaxMagnitude::VolumeMinMaxMagnitude(const VolumeMinMaxMagnitude& other)
    : VolumeDerivedData()
    , minMagnitude_(other.minMagnitude_)
    , maxMagnitude_(other.maxMagnitude_)
    , minNormalizedMagnitude_(other.minNormalizedMagnitude_)
    , maxNormalizedMagnitude_(other.maxNormalizedMagnitude_)
{}

VolumeMinMaxMagnitude::VolumeMinMaxMagnitude(float minMag, float maxMag, float minNormMag, float maxNormMag)
    : VolumeDerivedData()
    , minMagnitude_(minMag)
    , maxMagnitude_(maxMag)
    , minNormalizedMagnitude_(minNormMag)
    , maxNormalizedMagnitude_(maxNormMag)
{}

VolumeDerivedData* VolumeMinMaxMagnitude::create() const {
    return new VolumeMinMaxMagnitude();
}

VolumeDerivedData* VolumeMinMaxMagnitude::createFrom(const VolumeBase* handle) const {
    tgtAssert(handle, "no volume");

    RealWorldMapping rwm = handle->getRealWorldMapping();

    if (handle->hasRepresentation<VolumeRAM>()) {
        const VolumeRAM* v = handle->getRepresentation<VolumeRAM>();
        tgtAssert(v, "no volume");
        return createFromVolumeRAM(v, rwm);
    }
    else if (handle->hasRepresentation<VolumeDisk>()) {
        const VolumeDisk* volumeDisk = handle->getRepresentation<VolumeDisk>();
        tgtAssert(volumeDisk, "no disk volume");

        float minMagnitude = std::numeric_limits<float>::max();
        float maxMagnitude = 0.0f;
        float minNormalizedMagnitude = std::numeric_limits<float>::max();
        float maxNormalizedMagnitude = 0.0f;

        // compute min/max values slice-wise
        size_t numSlices = handle->getDimensions().z;
        tgtAssert(numSlices > 0, "empty volume");
        for (size_t i=0; i<numSlices; i++) {

            // interruption point after each slice!
            boost::this_thread::interruption_point();

            try {
                std::unique_ptr<VolumeRAM> sliceVolume(volumeDisk->loadSlices(i, i));
                std::unique_ptr<VolumeMinMaxMagnitude> sliceMinMaxMagnitude(createFromVolumeRAM(sliceVolume.get(), rwm));

                minMagnitude = std::min(minMagnitude, sliceMinMaxMagnitude->getMinMagnitude());
                maxMagnitude = std::max(maxMagnitude, sliceMinMaxMagnitude->getMaxMagnitude());

                minNormalizedMagnitude = std::min(minNormalizedMagnitude, sliceMinMaxMagnitude->getMinNormalizedMagnitude());
                maxNormalizedMagnitude = std::max(maxNormalizedMagnitude, sliceMinMaxMagnitude->getMaxNormalizedMagnitude());
            }
            catch (tgt::Exception& e) {
                LWARNING("Unable to compute min/max values: failed to load slice from disk volume: " << e.what());
                return nullptr;
            }
        }

        return new VolumeMinMaxMagnitude(minMagnitude, maxMagnitude, minNormalizedMagnitude, maxNormalizedMagnitude);;
    }

    LWARNING("Unable to compute min/max magnitude values: neither disk nor ram representation available");
    return nullptr;
}

float VolumeMinMaxMagnitude::getMinMagnitude() const {
    return minMagnitude_;
}

float VolumeMinMaxMagnitude::getMaxMagnitude() const {
    return maxMagnitude_;
}

float VolumeMinMaxMagnitude::getMinNormalizedMagnitude() const {
    return minNormalizedMagnitude_;
}

float VolumeMinMaxMagnitude::getMaxNormalizedMagnitude() const {
    return maxNormalizedMagnitude_;
}

void VolumeMinMaxMagnitude::serialize(Serializer& s) const  {
    s.serialize("minMagnitude", minMagnitude_);
    s.serialize("maxMagnitude", maxMagnitude_);
    s.serialize("minNormalizedMagnitude", minNormalizedMagnitude_);
    s.serialize("maxNormalizedMagnitude", maxNormalizedMagnitude_);
}

void VolumeMinMaxMagnitude::deserialize(Deserializer& s) {
    s.deserialize("minMagnitude", minMagnitude_);
    s.deserialize("maxMagnitude", maxMagnitude_);
    try {
        // Normalized values have been added later on.
        // Note, that former versions of min and max magnitude do not exist any longer.
        // This breaks compatibility with magnitudes calculated on other than
        // (vec) float or (vec) double but those haven't been used so far in voreen.
        s.deserialize("minNormalizedMagnitude", minNormalizedMagnitude_);
        s.deserialize("maxNormalizedMagnitude", maxNormalizedMagnitude_);
    }
    catch (SerializationException&) {
        LWARNING("Couldn't find normalized magnitudes, assuming default real world mapping");
        minNormalizedMagnitude_ = minMagnitude_;
        maxNormalizedMagnitude_ = maxMagnitude_;
    }
}

} // namespace voreen
