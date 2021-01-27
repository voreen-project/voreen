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

#include "voreen/core/datastructures/volume/volumeminmax.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"

namespace voreen {

VolumeMinMax::VolumeMinMax()
    : VolumeDerivedData()
{}

VolumeMinMax::VolumeMinMax(float min, float max, float minNorm, float maxNorm)
    : VolumeDerivedData()
{
    minValues_.push_back(min);
    maxValues_.push_back(max);
    minNormValues_.push_back(minNorm);
    maxNormValues_.push_back(maxNorm);
}

VolumeMinMax::VolumeMinMax(const std::vector<float>& minValues, const std::vector<float>& maxValues,
    const std::vector<float>& minNormValues, const std::vector<float>& maxNormValues )
{
    tgtAssert(!minValues.empty() && !maxValues.empty() && !minNormValues.empty() && !minNormValues.empty(), "empty vector passed");
    tgtAssert(minValues.size() == maxValues.size() && minValues.size() == minNormValues.size() && minValues.size() == maxNormValues.size(),
        "passed vectors differ in size");

    minValues_ = minValues;
    maxValues_ = maxValues;
    minNormValues_ = minNormValues;
    maxNormValues_ = maxNormValues;
}

VolumeMinMax::VolumeMinMax(const VolumeMinMax& other)
    : minValues_(other.minValues_)
    , maxValues_(other.maxValues_)
    , minNormValues_(other.minNormValues_)
    , maxNormValues_(other.maxNormValues_)
{}

VolumeDerivedData* VolumeMinMax::create() const {
    return new VolumeMinMax();
}

VolumeDerivedData* VolumeMinMax::createFrom(const VolumeBase* handle) const {
    tgtAssert(handle, "no volume");

    std::vector<float> tMinValues;
    std::vector<float> tMaxValues;
    std::vector<float> tMinNormValues;
    std::vector<float> tMaxNormValues;

    for (size_t i=0; i<handle->getNumChannels(); i++) {
        // use RAM representation only if already present,
        // otherwise use disk representation, if available
        float minNorm = 0.f;
        float maxNorm = 0.f;
        if (handle->hasRepresentation<VolumeRAM>()) {
            VolumeRAMRepresentationLock v(handle);
            tgtAssert(*v, "no volume");

            minNorm = v->minNormalizedValue(i);
            maxNorm = v->maxNormalizedValue(i);
        }
        else if (handle->hasRepresentation<VolumeDisk>()) {
            const VolumeDisk* volumeDisk = handle->getRepresentation<VolumeDisk>();
            tgtAssert(volumeDisk, "no disk volume");

            minNorm = FLT_MAX;
            maxNorm = -FLT_MAX;

            // compute min/max values slice-wise
            size_t numSlices = handle->getDimensions().z;
            tgtAssert(numSlices > 0, "empty volume");
            for (size_t slice=0; slice<numSlices; slice++) {

                // interruption point after each slice!
                boost::this_thread::interruption_point();

                // computation for the slice
                try {
                    VolumeRAM* sliceVolume = volumeDisk->loadSlices(slice, slice);
                    minNorm = std::min(minNorm, sliceVolume->minNormalizedValue(i));
                    maxNorm = std::max(maxNorm, sliceVolume->maxNormalizedValue(i));
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
        tgtAssert(minNorm <= maxNorm, "invalid min/max values");

        float min = handle->getRealWorldMapping().normalizedToRealWorld(minNorm);
        float max = handle->getRealWorldMapping().normalizedToRealWorld(maxNorm);

        tMinValues.push_back(min);
        tMaxValues.push_back(max);
        tMinNormValues.push_back(minNorm);
        tMaxNormValues.push_back(maxNorm);
    }

    return new VolumeMinMax(tMinValues, tMaxValues, tMinNormValues, tMaxNormValues);
}

size_t VolumeMinMax::getNumChannels() const {
    return minValues_.size();
}

float VolumeMinMax::getMin(size_t channel /*= 0*/) const {
    tgtAssert(channel < minValues_.size(), "invalid channel");
    return minValues_.at(channel);
}

float VolumeMinMax::getMax(size_t channel /*= 0*/) const {
    tgtAssert(channel < maxValues_.size(), "invalid channel");
    return maxValues_.at(channel);
}

float VolumeMinMax::getMinNormalized(size_t channel /*= 0*/) const {
    tgtAssert(channel < minNormValues_.size(), "invalid channel");
    return minNormValues_.at(channel);
}

float VolumeMinMax::getMaxNormalized(size_t channel /*= 0*/) const {
    tgtAssert(channel < maxNormValues_.size(), "invalid channel");
    return maxNormValues_.at(channel);
}

void VolumeMinMax::serialize(Serializer& s) const  {
    s.serialize("minValues", minValues_, "channel");
    s.serialize("maxValues", maxValues_, "channel");
    s.serialize("minNormValues", minNormValues_, "channel");
    s.serialize("maxNormValues", maxNormValues_, "channel");
}

void VolumeMinMax::deserialize(Deserializer& s) {
    try {
        s.deserialize("minValues", minValues_, "channel");
        s.deserialize("maxValues", maxValues_, "channel");
        s.deserialize("minNormValues", minNormValues_, "channel");
        s.deserialize("maxNormValues", maxNormValues_, "channel");
    }
    catch (SerializationException& /*e*/) {
        // try to deserialize old format (single channel)
        float min, max, minNorm, maxNorm;
        s.deserialize("min", min);
        s.deserialize("max", max);
        s.deserialize("minNorm", minNorm);
        s.deserialize("maxNorm", maxNorm);

        minValues_.clear();
        maxValues_.clear();
        minNormValues_.clear();
        maxNormValues_.clear();

        minValues_.push_back(min);
        maxValues_.push_back(max);
        minNormValues_.push_back(minNorm);
        maxNormValues_.push_back(maxNorm);
    }
}

} // namespace voreen
