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

#ifndef VRN_HISTOGRAM_H
#define VRN_HISTOGRAM_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumederiveddata.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

#include <string>
#include <iostream>
#include <fstream>

namespace voreen {

    /** Base class */
    class Histogram {
    };

template<typename T, int ND>
    class HistogramGeneric : public Histogram, public Serializable {
    public:
        int getNumBuckets(int dim) const {
            if((dim >= 0) && (dim < ND))
                return bucketCounts_[dim];
            else {
                tgtAssert(false, "Dimension-index out of range!");
                return 0;
            }
        }

        size_t getNumBuckets() const {
            return buckets_.size();
        }

        /// Returns the number of samples inserted into this histogram.
        uint64_t getNumSamples() const {
            return numSamples_;
        }

        /// Returns the number of samples in this bucket.
        uint64_t getBucket(size_t b) const {
            if(b < getNumBuckets())
                return buckets_[b];
            else {
                tgtAssert(false, "Index out of range!");
                return 0;
            }
        }

        uint64_t getBucket(size_t b1, size_t b2, ... /* size_t for each further dimension */) const {
            if(ND >= 2) {
                int c[ND];

                c[0] = static_cast<int>(b1);
                c[std::min(1,ND-1)] = static_cast<int>(b2); //Workaround to fix WARNING C4789 UNDER VC11

                va_list args;
                va_start(args, b2);
                for(int i=2; i<ND; i++)
                    c[i] = static_cast<int>(va_arg(args, size_t));
                va_end(args);

                int b = getBucketNumber(c);
                return getBucket(b);
            } else {
                tgtAssert(false, "histogram needs at least 2 dimensions for this function!");
                return 0;
            }
        }

        uint64_t getMaxBucket() const {
            //TODO: optimize?
            uint64_t max = 0;
            for(size_t i=0; i<getNumBuckets(); i++)
                if(buckets_[i] > max)
                    max = buckets_[i];

            return max;
        }

        /// Returns the normalized count in the bucket (getBucket(b) / getMaxBucket())
        float getBucketNormalized(int b) const {
            return (float) getBucket(b) / (float) getMaxBucket();
        }

        float getBucketLogNormalized(int b) const {
            return (logf(static_cast<float>(1+getBucket(b)) ) / logf( static_cast<float>(1+getMaxBucket())));
        }


        void increaseBucket(size_t b) {
            if(b < getNumBuckets()) {
                buckets_[b]++;
                numSamples_++;
            }
            else
                tgtAssert(false, "Index out of range!");
        }

        void increaseBucket(size_t bucket, uint64_t value) {
            if(bucket < getNumBuckets()) {
                buckets_[bucket] += value;
                numSamples_ += value;
            }
            else
                tgtAssert(false, "Index out of range!");
        }

        void decreaseBucket(size_t b) {
            if(b < getNumBuckets() && buckets_[b] > 0) {
                buckets_[b]--;
                numSamples_--;
            }
            else
                tgtAssert(false, "Index out of range!");
        }

        T getMinValue(int dim) const {
            return minValues_[dim];
        }

        T getMaxValue(int dim) const {
            return maxValues_[dim];
        }

        virtual void serialize(Serializer& s) const {
            std::vector<T> temp;

            for(int i=0; i<ND; i++)
                temp.push_back(minValues_[i]);
            s.serialize("minValues", temp);
            temp.clear();

            for(int i=0; i<ND; i++)
                temp.push_back(maxValues_[i]);
            s.serialize("maxValues", temp);
            temp.clear();

            std::vector<int> temp2;
            for(int i=0; i<ND; i++)
                temp2.push_back(bucketCounts_[i]);
            s.serialize("bucketCounts", temp2);

            s.serializeBinaryBlob("binaryBuckets", buckets_);
            //s.serialize("buckets", buckets_);
        }

        virtual void deserialize(Deserializer& s) {
            std::vector<T> temp;

            s.deserialize("minValues", temp);
            if(temp.size() != ND)
                throw tgt::CorruptedFileException("Dimension mismatch!");
            for(int i=0; i<ND; i++)
                minValues_[i] = temp[i];
            temp.clear();

            s.deserialize("maxValues", temp);
            if(temp.size() != ND)
                throw tgt::CorruptedFileException("Dimension mismatch!");
            for(int i=0; i<ND; i++)
                maxValues_[i] = temp[i];
            temp.clear();

            std::vector<int> temp2;
            s.deserialize("bucketCounts", temp2);
            if(temp2.size() != ND)
                throw tgt::CorruptedFileException("Dimension mismatch!");
            for(int i=0; i<ND; i++)
                bucketCounts_[i] = temp2[i];

            try {
                s.deserializeBinaryBlob("binaryBuckets", buckets_);
            } catch(...) {
                //try old
                s.deserialize("buckets", buckets_);
            }

            size_t numBuckets = 1;
            for(int i=0; i<ND; i++)
                numBuckets *= (size_t)bucketCounts_[i];

            if(numBuckets != buckets_.size())
                throw tgt::CorruptedFileException("Bucket number mismatch!");

            numSamples_ = 0;
            for(size_t i=0; i<buckets_.size(); i++)
                numSamples_ += buckets_[i];
        }

    protected:
        void addSample(double sample, ... /* double sample for each further dimension */) {
            T values[ND];
            values[0] = static_cast<T>(sample);

            va_list args;
            va_start(args, sample);
            for(int i=1; i<ND; i++) {
                //values[i] = va_arg(args, T);
                values[i] = static_cast<T>(va_arg(args, double));
            }

            va_end(args);

            int c[ND];
            for(int i=0; i<ND; i++)
                c[i] = mapValueToBucket(values[i], i);

            int b = getBucketNumber(c);
            increaseBucket(b);
        }

        int mapValueToBucket(T v, int dim) const {
            int bucket =  static_cast<int>(bucketCounts_[dim] * ((v - minValues_[dim])/ (maxValues_[dim] - minValues_[dim])));
            if(bucket < 0) {
                return 0;
            } else if(bucket >= bucketCounts_[dim]) {
                return (bucketCounts_[dim] - 1);
            } else {
                return bucket;
            }
        }

        int getBucketNumber(const int* c) const {
            int n = 0;
            int helper = 1;
            for(int i=0; i<ND; i++) {
                if((c[i] >= 0) && (c[i] < bucketCounts_[i])) {
                    n += helper * c[i];
                    helper *= bucketCounts_[i];
                }
                else {
                    //TODO
                }
            }
            return n;
        }

        HistogramGeneric(double minValue, double maxValue, int bucketCount, ... /*, double minValue, double maxValue, int bucketCount for each further dimension */) : numSamples_(0) {
            va_list args;
            va_start(args, bucketCount);

            minValues_[0] = static_cast<T>(minValue);
            maxValues_[0] = static_cast<T>(maxValue);
            bucketCounts_[0] = bucketCount;

            for(int i=1; i<ND; i++) {
                //minValues_[i] = va_arg(args, T);
                //maxValues_[i] = va_arg(args, T);
                minValues_[i] = static_cast<T>(va_arg(args, double));
                maxValues_[i] = static_cast<T>(va_arg(args, double));
                bucketCounts_[i] = va_arg(args, int);
            }

            va_end(args);

            int numBuckets = 1;
            for(int i=0; i<ND; i++)
                numBuckets *= bucketCounts_[i];

            buckets_.assign(numBuckets, 0);
        }

    private:
        T minValues_[ND];
        T maxValues_[ND];
        int bucketCounts_[ND];

        std::vector<uint64_t> buckets_;
        uint64_t numSamples_;
    };

template <typename T>
class Histogram1DGeneric : public HistogramGeneric<T, 1> {
    public:
        Histogram1DGeneric(T minValue, T maxValue, int bucketCount) : HistogramGeneric<T, 1>(static_cast<double>(minValue), static_cast<double>(maxValue), bucketCount) {}

        void addSample(T value) {
            HistogramGeneric<T, 1>::addSample(static_cast<double>(value));
        }
        T getMinValue() const {
            return HistogramGeneric<T, 1>::getMinValue(0);
        }

        T getMaxValue() const {
            return HistogramGeneric<T, 1>::getMaxValue(0);
        }
    private:
};

class VRN_CORE_API Histogram1D : public Histogram1DGeneric<float> {
    public:
        Histogram1D(float minValue, float maxValue, int bucketCount) : Histogram1DGeneric<float>(minValue, maxValue, bucketCount) {}
        Histogram1D() : Histogram1DGeneric<float>(0.f, 1.f, 256) {}
};

VRN_CORE_API Histogram1D createHistogram1DFromVolume(const VolumeBase* handle, size_t bucketCount, size_t channel = 0);
VRN_CORE_API Histogram1D createHistogram1DFromVolume(const VolumeBase* handle, size_t bucketCount, float realWorldMin, float realWorldMax, size_t channel = 0);


//--------------------------------------------------------------------------

template <typename T>
class Histogram2DGeneric : public HistogramGeneric<T, 2> {
    public:
        Histogram2DGeneric(T minValue1, T maxValue1, int bucketCount1, T minValue2, T maxValue2, int bucketCount2) : HistogramGeneric<T, 2>(static_cast<double>(minValue1), static_cast<double>(maxValue1), bucketCount1, static_cast<double>(minValue2), static_cast<double>(maxValue2), bucketCount2) {}

        void addSample(T value1, T value2) {
            HistogramGeneric<T, 2>::addSample(static_cast<double>(value1), static_cast<double>(value2));
        }

    private:
};

class VRN_CORE_API Histogram2D : public Histogram2DGeneric<float> {
    public:
        Histogram2D(float minValue1, float maxValue1, int bucketCount1, float minValue2, float maxValue2, int bucketCount2) : Histogram2DGeneric<float>(minValue1, maxValue1, bucketCount1, minValue2, maxValue2, bucketCount2) {}
        Histogram2D() : Histogram2DGeneric<float>(0.f, 1.f, 256, 0.f, 1.f, 256) {}
};

VRN_CORE_API Histogram2D createHistogram2DFromVolume(const VolumeBase* handle, int bucketCountIntensity, int bucketCountGradient, size_t channel = 0);

//--------------------------------------------------------------------------

/// 1D Intensity Histogram.
class VRN_CORE_API VolumeHistogramIntensity : public VolumeDerivedData {
public:
    /// Copy constructor.
    VolumeHistogramIntensity(const VolumeHistogramIntensity& h);
    VolumeHistogramIntensity(const Histogram1D& h);
    VolumeHistogramIntensity(const std::vector<Histogram1D>& histograms);

    /// Empty default constructor required by VolumeDerivedData interface.
    VolumeHistogramIntensity();
    virtual std::string getClassName() const { return "VolumeHistogramIntensity"; }
    virtual VolumeDerivedData* create() const;

    /**
     * Creates a histogram with a bucket count of 256.
     *
     * @see VolumeDerivedData
     */
    virtual VolumeDerivedData* createFrom(const VolumeBase* handle) const;

    /// Returns the number of channel histograms stored in this derived data.
    size_t getNumChannels() const;

    /// Returns the number buckets of the specified channel histogram.
    size_t getBucketCount(size_t channel = 0) const;

    /// get value in bucket i
    uint64_t getValue(int i, size_t channel = 0) const;

    /// get value in bucket i
    uint64_t getValue(size_t i, size_t channel = 0) const;

    /// get value in bucket nearest to i
    uint64_t getValue(float i, size_t channel = 0) const;

    /// Returns normalized (with max.) histogram value at bucket i
    float getNormalized(int i, size_t channel = 0) const;

    /// Returns normalized (with max.) histogram value at bucket nearest to i
    float getNormalized(float i, size_t channel = 0) const;

    /// Returns normalized logarithmic histogram value at bucket i
    float getLogNormalized(int i, size_t channel = 0) const;

    /// Returns normalized logarithmic histogram value at bucket nearest to i
    float getLogNormalized(float i, size_t channel = 0) const;

    /// @see VolumeDerivedData
    virtual void serialize(Serializer& s) const;

    /// @see VolumeDerivedData
    virtual void deserialize(Deserializer& s);

    const Histogram1D& getHistogram(size_t channel = 0) const;
    Histogram1D& getHistogram(size_t channel = 0);

protected:
    std::vector<Histogram1D> histograms_;
};

// ----------------------------------------------------------------------------

/// 2D histogram using intensity and gradient length.
class VRN_CORE_API VolumeHistogramIntensityGradient : public VolumeDerivedData {
public:
    /// Empty default constructor required by VolumeDerivedData interface.
    VolumeHistogramIntensityGradient();

    /**
     * Create a single channeled 2D histogram.
     * @param h 2D histogram
     * @param maxBucket max bucket value.
     */
    VolumeHistogramIntensityGradient(const Histogram2D& h, uint64_t maxBucket);

    /**
     * Create a histogram from the given histograms and maxBuckets.
     * @param h 2D histogram
     * @param maxBucket max bucket value.
     *
     * @note histograms and maxBuckets must have the same size.
     */
    VolumeHistogramIntensityGradient(const std::vector<Histogram2D>& histograms, const std::vector<uint64_t>& maxBuckets);

    virtual std::string getClassName() const { return "VolumeHistogramIntensityGradient"; }
    virtual VolumeDerivedData* create() const;

    virtual VolumeDerivedData* createFrom(const VolumeBase* handle) const;

    /// Returns the number of channel histograms stored in this derived data.
    size_t getNumChannels() const;

    /// Returns voxels in bucket.
    int getValue(int i, int g, size_t channel = 0) const;

    /// Returns normalized (with max.) histogram value
    float getNormalized(int i, int g, size_t channel = 0) const;

    /// Returns normalized logarithmic histogram value
    float getLogNormalized(int i, int g, size_t channel = 0) const;

    /// Returns the maximal bucket value in the histogram
    int getMaxBucket(size_t channel = 0) const;

    float getMinValue(int dim, size_t channel = 0) const;
    float getMaxValue(int dim, size_t channel = 0) const;

    size_t getBucketCountIntensity(size_t channel = 0) const;
    size_t getBucketCountGradient(size_t channel = 0) const;

    /// @see VolumeDerivedData (currently unimplemented)
    virtual void serialize(Serializer& s) const;

    /// @see VolumeDerivedData (currently unimplemented)
    virtual void deserialize(Deserializer& s);

    const Histogram2D& getHistogram(size_t channel = 0) const;
    Histogram2D& getHistogram(size_t channel = 0);
protected:
    std::vector<Histogram2D> hist_;
    std::vector<uint64_t> maxBucket_;
};

} // namespace voreen

#endif // VRN_HISTOGRAM_H
