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

#include "voreen/core/datastructures/volume/histogramutils.h"

namespace voreen {

    Histogram1D* subHistogram(Histogram1D* histogram, float minVal, float maxVal) {

        tgtAssert(minVal < maxVal, "minVal must be smaller than maxVal");

        float oldMin = histogram->getMinValue();
        float oldMax = histogram->getMaxValue();
        float oldWidth = oldMax - oldMin;
        int oldCount = static_cast<int>(histogram->getNumBuckets());

        int firstBucket = static_cast<int>((minVal - oldMin) / oldWidth * oldCount);
        int lastBucket = static_cast<int>((maxVal - oldMin) / oldWidth * oldCount);
        int bucketCount = lastBucket - firstBucket + 1;

        Histogram1D* h = new Histogram1D(minVal, maxVal, bucketCount);

        // check the bucket limits of the histograms against each other
        if (static_cast<size_t>(firstBucket) >= histogram->getNumBuckets() || lastBucket < 0)
            return h;

        int bucket = 0;
        if (firstBucket < 0) {
            bucket = -firstBucket;
            firstBucket = 0;
        }
        if (static_cast<size_t>(lastBucket) >= histogram->getNumBuckets()) {
            lastBucket = static_cast<int>(histogram->getNumBuckets()) - 1;
        }

        // copy all relevant buckets to the new histogram
        for (int i = firstBucket; i < lastBucket; i++) {
            h->increaseBucket(bucket, histogram->getBucket(i));
            bucket++;
        }

        return h;
    }

    GaussianCurve* otsuTreshold(Histogram1D* histogram) {

        // this threshold method is proposed in Otsu 1979
        // We're classifying the underlying values of the histogram
        // in two optimal gaussian distributions.

        uint64_t sampleCount = histogram->getNumSamples();
        tgtAssert(sampleCount >= 2, "The histogram contains too few samples for an Otsu classification.")
        size_t bucketCount = histogram->getNumBuckets();

        // normalize the histogram to obtain the probability distribution
        // also check if there are more than one buckets with samples in it
        size_t nonEmptyBuckets = 0;
        float* p = new float[bucketCount];
        for (size_t i = 0; i < bucketCount; i++) {
            p[i] = histogram->getBucket(i) / (float)sampleCount;

            if (p[i] > 0)
                nonEmptyBuckets++;
        }
        tgtAssert(nonEmptyBuckets >= 2, "The histogram contains too few non empty buckets for an Otsu classification.")

        GaussianCurve* c = new GaussianCurve[2];

        float mK = 0; // used for efficiently summing i*p(i) from 0 to k
        float w0 = 0, w1 = 0;
        float m0, m1, mTotal = 0;
        bool success = false;
        // the total mean value mTotal does not change over the course of the algorithm
        // and can therefore be calculated beforehand
        for (size_t k = 1; k < bucketCount; k++)
            mTotal += k * p[k];

        float max = 0.0f;   // previous best (maximal) varianceB variance value
        float varianceB = 0.0f; // current varianceB variance
        size_t t = 0;     // previous best threshold
        float bestW0 = 0, bestW1 = 0;   // w0 and w1 for threshold t

        // try all possible thresholds k and save the values for the optimum k
        // We're searching for the k that maximizes the between class variance
        for (size_t k = 0; k < bucketCount; k++) {

            mK += static_cast<float>(k) * p[k];

            w0 += p[k];
            // skip empty buckets at the margin of the histogram
            if (w0 <= 0 || w0 >=1)
                continue;
            w1 = 1.f - w0;

            m0 = mK / w0;
            m1 = (mTotal - mK) / w1;

            // the between-class variance shall be maximized
            varianceB = w0 * w1 * (m1 - m0) * (m1 - m0);
            if (varianceB > max) {
                // new best threshold t is found
                t = k;
                max = varianceB;

                // save the current means
                c[0].mean = m0;
                c[1].mean = m1;
                // .. and w0, w1
                bestW0 = w0;
                bestW1 = w1;

                // a classifiaction was possible
                success = true;
            }
        }

        tgtAssert(success, "Otsu classification failed!");

        // calculate the resulting variances
        m0 = c[0].mean;
        m1 = c[1].mean;
        float v0 = 0, v1 = 0;
        for (size_t k = 0; k <= t; k++)
            v0 += (static_cast<float>(k) - m0)*(static_cast<float>(k) - m0)*p[k] / bestW0;
        for (size_t k = (t+1); k < bucketCount; k++)
            v1 += (static_cast<float>(k) - m1)*(static_cast<float>(k) - m1)*p[k] / bestW1;
        c[0].variance = v0;
        c[1].variance = v1;

        // map the curves from 'histogram bucket-space' to the real space represented by the histogram
        c[0] = mapCurve(c[0], tgt::vec2(0, static_cast<float>(histogram->getNumBuckets())), tgt::vec2(histogram->getMinValue(), histogram->getMaxValue()));
        c[1] = mapCurve(c[1], tgt::vec2(0, static_cast<float>(histogram->getNumBuckets())), tgt::vec2(histogram->getMinValue(), histogram->getMaxValue()));

        // cleanup
        delete[] p;

        return c;
    }


    //================================
    //       Helper Functions
    //================================

    GaussianCurve mapCurve(GaussianCurve& c, tgt::vec2 from, tgt::vec2 to) {
        float widthF = from.y - from.x;
        float widthT = to.y - to.x;

        GaussianCurve r;
        r.mean = (c.mean - from.x) / widthF * widthT + to.x;    // new mean
        r.variance = c.variance * widthT * widthT / widthF / widthF; // new variance

        return r;
    }


    /*
     * Method which computes otsu but only returns the threshold
     * TODO: remove code duplication with the method above
     */
    float computeOptimalOtsuThreshold(const Histogram1D* histogram) {
        uint64_t sampleCount = histogram->getNumSamples();
        tgtAssert(sampleCount >= 2, "The histogram contains too few samples for an Otsu classification.")
        size_t bucketCount = histogram->getNumBuckets();

        // normalize the histogram to obtain the probability distribution
        // also check if there are more than one buckets with samples in it
        size_t nonEmptyBuckets = 0;
        float* p = new float[bucketCount];
        for (size_t i = 0; i < bucketCount; i++) {
            p[i] = histogram->getBucket(i) / (float)sampleCount;

            if (p[i] > 0)
                nonEmptyBuckets++;
        }
        tgtAssert(nonEmptyBuckets >= 2, "The histogram contains too few non empty buckets for an Otsu classification.")

        float mK = 0; // used for efficiently summing i*p(i) from 0 to k
        float w0 = 0, w1 = 0;
        float m0, m1, mTotal = 0;
        bool success = false;
        // the total mean value mTotal does not change over the course of the algorithm
        // and can therefore be calculated beforehand
        for (size_t k = 1; k < bucketCount; k++)
            mTotal += static_cast<float>(k) * p[k];

        float max = 0.0f;   // previous best (maximal) varianceB variance value
        float varianceB = 0.0f; // current varianceB variance
        size_t t = 0;     // previous best threshold

        // try all possible thresholds k and save the values for the optimum k
        // We're searching for the k that maximizes the between class variance
        for (size_t k = 0; k < bucketCount; k++) {

            mK += static_cast<float>(k) * p[k];

            w0 += p[k];
            // skip empty buckets at the margin of the histogram
            if (w0 <= 0 || w0 >=1)
                continue;
            w1 = 1.f - w0;

            m0 = mK / w0;
            m1 = (mTotal - mK) / w1;

            // the between-class variance shall be maximized
            varianceB = w0 * w1 * (m1 - m0) * (m1 - m0);
            if (varianceB > max) {
                // new best threshold t is found
                t = k;
                max = varianceB;

                // a classifiaction was possible
                success = true;
            }
        }

        tgtAssert(success, "Otsu classification failed!");

        // map the threshold from 'histogram bucket-space' to the real space represented by the histogram
        float widthF = static_cast<float>(histogram->getNumBuckets());
        float widthT = histogram->getMaxValue() - histogram->getMinValue(); 
        float tTransformed = static_cast<float>(t) / widthF * widthT + histogram->getMinValue();

        return tTransformed;
    }


} // namespace voreen
