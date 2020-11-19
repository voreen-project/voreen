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

#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgradient.h"

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

using tgt::vec3;
using tgt::ivec3;

Histogram1D createHistogram1DFromVolume(const VolumeBase* handle, size_t bucketCount, size_t channel /*= 0*/) {
    RealWorldMapping rwm = handle->getRealWorldMapping();

    VolumeMinMax* volumeMinMax = 0;
    while(!volumeMinMax) {
        volumeMinMax = handle->getDerivedData<VolumeMinMax>();
        boost::this_thread::interruption_point();
    }

    tgtAssert(volumeMinMax, "null pointer in volume min max derived data");
    tgtAssert(channel < volumeMinMax->getNumChannels(), "invalid channel");


    float min = volumeMinMax->getMinNormalized(channel);
    float max = volumeMinMax->getMaxNormalized(channel);
    min = rwm.normalizedToRealWorld(min);
    max = rwm.normalizedToRealWorld(max);

    return createHistogram1DFromVolume(handle, bucketCount, min, max, channel);
}

VRN_CORE_API Histogram1D createHistogram1DFromVolume(const VolumeBase* handle, size_t bucketCount,
    float realWorldMin, float realWorldMax, size_t channel /*= 0*/)
{
    tgtAssert(realWorldMin <= realWorldMax, "invalid real world range");

    RealWorldMapping rwm = handle->getRealWorldMapping();
    Histogram1D h(realWorldMin, realWorldMax, (int)bucketCount);

    // prefer RAM over disk representation, but only if RAM volume is already present
    const VolumeRAM* volumeRam = 0;
    const VolumeDisk* volumeDisk = 0;
    if (handle->hasRepresentation<VolumeRAM>())
        volumeRam = handle->getRepresentation<VolumeRAM>();
    else if (handle->hasRepresentation<VolumeDisk>())
        volumeDisk = handle->getRepresentation<VolumeDisk>();
    else {
        LWARNINGC("voreen.Histogram", "Unable to compute 1D histogram: neither disk nor RAM representation available");
        return h;
    }
    tgtAssert(volumeRam || volumeDisk, "no representation");

    // iterate over slices
    tgt::svec3 dims = handle->getDimensions();
    tgt::svec3 pos;
    for (pos.z = 0; pos.z < dims.z; ++pos.z) {
        boost::this_thread::interruption_point();

        if (volumeRam) {
            // access volume data in RAM directly
            for (pos.y = 0; pos.y < dims.y; ++pos.y) {
                for (pos.x = 0; pos.x < dims.x; ++pos.x) {
                    float val = volumeRam->getVoxelNormalized(pos, channel);
                    val = rwm.normalizedToRealWorld(val);
                    h.addSample(val);
                }
            }
        }
        else if (volumeDisk) {
            try {
                // temporarily load current slice into RAM
                VolumeRAM* sliceVolume = volumeDisk->loadSlices(pos.z, pos.z);
                tgtAssert(sliceVolume, "null pointer returned (exception expected)");
                for (pos.y = 0; pos.y < dims.y; ++pos.y) {
                    for (pos.x = 0; pos.x < dims.x; ++pos.x) {
                        float val = sliceVolume->getVoxelNormalized(tgt::svec3(pos.x, pos.y, 0), channel);
                        val = rwm.normalizedToRealWorld(val);
                        h.addSample(val);
                    }
                }
                delete sliceVolume;
            }
            catch (tgt::Exception& e) {
                LWARNINGC("voreen.Histogram", "Unable to compute 1D histogram: failed to load slice from disk volume: " + std::string(e.what()));
                return h;
            }
        }
        else {
            tgtAssert(false, "should never get here");
        }
    }

    return h;
}

//-----------------------------------------------------------------------------

Histogram2D createHistogram2DFromVolume(const VolumeBase* handle, int bucketCountIntensity, int bucketCountGradient, size_t channel) {

    RealWorldMapping rwm = handle->getRealWorldMapping();
    ivec3 dims = handle->getDimensions();
    vec3 sp = handle->getSpacing();

    VolumeMinMax* vmm = 0;
    while(!vmm) {
        vmm = handle->getDerivedData<VolumeMinMax>();
        boost::this_thread::interruption_point();
    }

    tgtAssert(vmm, "null pointer in volume min max derived data");
    tgtAssert(channel < vmm->getNumChannels(), "invalid channel");

    float min = vmm->getMinNormalized(channel);
    float max = vmm->getMaxNormalized(channel);
    min = rwm.normalizedToRealWorld(min);
    max = rwm.normalizedToRealWorld(max);

    float minGradLength = 0.0f; // always 0
    float maxGradLength = 0.0f;

    VolumeRAMRepresentationLock vol(handle);
    if(!*vol) {
        LWARNINGC("voreen.Histogram", "Unable to compute 2D histogram: no RAM representation available");
        Histogram2D h(min, max, bucketCountIntensity, minGradLength, maxGradLength, bucketCountGradient);
        return h;
    }

    //TODO: improve performance
    ivec3 pos;
    for (pos.z = 0; pos.z < dims.z; ++pos.z) {
        boost::this_thread::interruption_point();

        for (pos.y = 0; pos.y < dims.y; ++pos.y) {
            for (pos.x = 0; pos.x < dims.x; ++pos.x) {
                //vec3 grad = VolumeOperatorGradient::calcGradientCentralDifferences(vol, sp, pos);
                vec3 grad = VolumeOperatorGradient::calcGradientSobel(*vol, sp, pos, channel);

                float nlength = tgt::length(grad) * rwm.getScale();

                if (nlength > maxGradLength)
                    maxGradLength = nlength;
            }
        }
    }

    Histogram2D h(min, max, bucketCountIntensity, minGradLength, maxGradLength, bucketCountGradient);
    for (pos.z = 0; pos.z < dims.z; ++pos.z) {
        boost::this_thread::interruption_point();

        for (pos.y = 0; pos.y < dims.y; ++pos.y) {
            for (pos.x = 0; pos.x < dims.x; ++pos.x) {
                //vec3 grad = VolumeOperatorGradient::calcGradientCentralDifferences(vol, sp, pos);
                vec3 grad = VolumeOperatorGradient::calcGradientSobel(*vol, sp, pos, channel);

                float nlength = tgt::length(grad) * rwm.getScale();

                float v = vol->getVoxelNormalized(pos, channel);
                v = rwm.normalizedToRealWorld(v);

                h.addSample(v, nlength);
            }
        }
    }

    return h;
}

//-----------------------------------------------------------------------------

VolumeHistogramIntensity::VolumeHistogramIntensity() :
    VolumeDerivedData(),
    histograms_()
{}

VolumeHistogramIntensity::VolumeHistogramIntensity(const VolumeHistogramIntensity& h)
    : VolumeDerivedData()
    , histograms_(h.histograms_)
{}

VolumeHistogramIntensity::VolumeHistogramIntensity(const Histogram1D& h) {
    histograms_.push_back(h);
}

VolumeHistogramIntensity::VolumeHistogramIntensity(const std::vector<Histogram1D>& histograms)
    : VolumeDerivedData()
    , histograms_(histograms)
{}

VolumeDerivedData* VolumeHistogramIntensity::createFrom(const VolumeBase* handle) const {
    tgtAssert(handle, "no volume");
    VolumeHistogramIntensity* h = new VolumeHistogramIntensity();
    for (size_t channel = 0; channel < handle->getNumChannels(); channel++) {
        Histogram1D hist = createHistogram1DFromVolume(handle, 256, channel);
        h->histograms_.push_back(hist);
    }
    return h;
}

size_t VolumeHistogramIntensity::getNumChannels() const {
    return histograms_.size();
}

size_t VolumeHistogramIntensity::getBucketCount(size_t channel) const {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    return histograms_.at(channel).getNumBuckets();
}

uint64_t VolumeHistogramIntensity::getValue(int bucket, size_t channel) const {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    return histograms_.at(channel).getBucket(bucket);
}

uint64_t VolumeHistogramIntensity::getValue(size_t bucket, size_t channel) const {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    return getValue(static_cast<int>(bucket), channel);
}

uint64_t VolumeHistogramIntensity::getValue(float i, size_t channel) const {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    size_t bucketCount = histograms_.at(channel).getNumBuckets();
    float m = (bucketCount - 1.f);
    int bucket = static_cast<int>(floor(i * m));
    return getValue(bucket, channel);
}

float VolumeHistogramIntensity::getNormalized(int i, size_t channel) const {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    return histograms_.at(channel).getBucketNormalized(i);
}

float VolumeHistogramIntensity::getNormalized(float i, size_t channel) const {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    size_t bucketCount = histograms_.at(channel).getNumBuckets();
    float m = (bucketCount - 1.f);
    int bucket = static_cast<int>(floor(i * m));
    return getNormalized(bucket, channel);
}

float VolumeHistogramIntensity::getLogNormalized(int i, size_t channel) const {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    const Histogram1D& hist = histograms_.at(channel);
    return hist.getBucketLogNormalized(i);
}

float VolumeHistogramIntensity::getLogNormalized(float i, size_t channel) const {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    size_t bucketCount = histograms_.at(channel).getNumBuckets();
    float m = (bucketCount - 1.f);
    int bucket = static_cast<int>(floor(i * m));
    return getLogNormalized(bucket, channel);
}

void VolumeHistogramIntensity::serialize(Serializer& s) const  {
    s.serialize("histograms", histograms_, "channel");
}

void VolumeHistogramIntensity::deserialize(Deserializer& s) {
    try {
        s.deserialize("histograms", histograms_, "channel");
    }
    catch (SerializationException& /*e*/) {
        // try to deserialize old format (single channel)
        Histogram1D hist;
        s.deserialize("histogram", hist);

        histograms_.clear();
        histograms_.push_back(hist);
    }
}

const Histogram1D& VolumeHistogramIntensity::getHistogram(size_t channel) const {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    return histograms_.at(channel);
}

Histogram1D& VolumeHistogramIntensity::getHistogram(size_t channel) {
    tgtAssert(channel < histograms_.size(), "invalid channel");
    return histograms_.at(channel);
}

VolumeDerivedData* VolumeHistogramIntensity::create() const {
    return new VolumeHistogramIntensity();
}

//-----------------------------------------------------------------------------

VolumeHistogramIntensityGradient::VolumeHistogramIntensityGradient()
    : VolumeDerivedData()
    //, hist_(0.0f, 1.0f, 1, 0.0f, 1.0f, 1)
{}

VolumeHistogramIntensityGradient::VolumeHistogramIntensityGradient(const Histogram2D& h, uint64_t maxBucket)
    : VolumeDerivedData()
{
    hist_.push_back(h);
    maxBucket_.push_back(maxBucket);
}

VolumeHistogramIntensityGradient::VolumeHistogramIntensityGradient(const std::vector<Histogram2D>& histograms, const std::vector<uint64_t>& maxBuckets)
    : VolumeDerivedData()
    , hist_(histograms)
    , maxBucket_(maxBuckets)
{
    tgtAssert(histograms.size() == maxBuckets.size(), "Number of histograms and number of maxBuckets differ.");
}

VolumeDerivedData* VolumeHistogramIntensityGradient::createFrom(const VolumeBase* handle) const {
    tgtAssert(handle, "no volume");
    VolumeHistogramIntensityGradient* h = new VolumeHistogramIntensityGradient();
    for(size_t ch = 0; ch < handle->getNumChannels(); ch++) {
        Histogram2D hist = createHistogram2DFromVolume(handle, 256, 256, ch);
        h->hist_.push_back(hist);
        h->maxBucket_.push_back(hist.getMaxBucket());
    }
    return h;
}

size_t VolumeHistogramIntensityGradient::getNumChannels() const {
    return hist_.size();
}

size_t VolumeHistogramIntensityGradient::getBucketCountIntensity(size_t channel) const {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return hist_[channel].getNumBuckets(0);
}

size_t VolumeHistogramIntensityGradient::getBucketCountGradient(size_t channel) const {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return hist_[channel].getNumBuckets(1);
}

int VolumeHistogramIntensityGradient::getValue(int i, int g, size_t channel) const {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return static_cast<int>(hist_[channel].getBucket(i, g));
}

float VolumeHistogramIntensityGradient::getNormalized(int i, int g, size_t channel) const {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return (static_cast<float>(getValue(i, g, channel)) / static_cast<float>(getMaxBucket(channel)));
}

float VolumeHistogramIntensityGradient::getLogNormalized(int i, int g, size_t channel) const {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return (logf(static_cast<float>(1+getValue(i, g, channel)) ) / log(static_cast<float>(1+getMaxBucket(channel))));
}

int VolumeHistogramIntensityGradient::getMaxBucket(size_t channel) const {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return static_cast<int>(maxBucket_[channel]);
}

float VolumeHistogramIntensityGradient::getMinValue(int dim, size_t channel) const {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return hist_[channel].getMinValue(dim);
}

float VolumeHistogramIntensityGradient::getMaxValue(int dim, size_t channel) const {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return hist_[channel].getMaxValue(dim);
}

void VolumeHistogramIntensityGradient::serialize(Serializer& s) const {
    s.serialize("histograms", hist_, "channel");
    s.serialize("maxBuckets", maxBucket_, "channel");
}

void VolumeHistogramIntensityGradient::deserialize(Deserializer& s) {
   try {
        s.deserialize("histograms", hist_, "channel");
        s.deserialize("maxBuckets", maxBucket_, "channel");
    }
    catch (SerializationException& /*e*/) {
        // try to deserialize old format (single channel)
        Histogram2D hist;
        s.deserialize("histogram", hist);

        hist_.clear();
        hist_.push_back(hist);
    }
}

VolumeDerivedData* VolumeHistogramIntensityGradient::create() const {
    return new VolumeHistogramIntensityGradient();
}

const Histogram2D& VolumeHistogramIntensityGradient::getHistogram(size_t channel) const {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return hist_.at(channel);
}

Histogram2D& VolumeHistogramIntensityGradient::getHistogram(size_t channel) {
    tgtAssert(channel < hist_.size(), "invalid channel");
    return hist_.at(channel);
}

} // namespace voreen
