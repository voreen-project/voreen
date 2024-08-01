/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "voreen/core/utils/statistics.h"
#include "tgt/assert.h"
#include "float.h"

#include <algorithm>

namespace voreen {

Statistics::Statistics(bool collectSamples)
   : min_(FLT_MAX)
   , max_(FLT_MIN)
   , runningMean_(0.0f)
   , m2_(0.0f)
   , sum_(0.0f)
   , collectSamples_(collectSamples)
   , numSamples_(0)
{
}

void Statistics::reset() {
   min_ = FLT_MAX;
   max_ = FLT_MIN;
   runningMean_ = 0.0f;
   m2_ = 0.0f;
   sum_ = 0.0f;
   samples_.clear();
   numSamples_ = 0;
}

void Statistics::addSample(float v) {
    numSamples_++;

    sum_ += v;

    if(collectSamples_)
        samples_.push_back(v);

    if(v > max_)
        max_ = v;
    if(v < min_)
        min_ = v;

    double delta = v - runningMean_;
    runningMean_ = runningMean_ + delta/(double)numSamples_;
    if (numSamples_ > 1)
        m2_ = m2_ + delta*(v - runningMean_);
}


float Statistics::getMin() const { return min_; }
float Statistics::getMax() const { return max_; }
float Statistics::getSum() const { return static_cast<float>(sum_); }
size_t Statistics::getNumSamples() const { return numSamples_; }

float Statistics::getMean() const { return static_cast<float>(runningMean_); }
float Statistics::getVariance() const { return static_cast<float>(m2_/(double)numSamples_); }// std dev for population
float Statistics::getStdDev() const { return sqrt(getVariance()); }
float Statistics::getRelStdDev() const { return getStdDev() / getMean(); }

float Statistics::getMedian() const {
    float median = 0.0f;
    if(!collectSamples_) {
        tgtAssert(false, "collectSamples needs to be true in order to calculate the median");
        return median;
    }

    if((numSamples_ % 2) == 1) {
        size_t n = numSamples_/2;
        nth_element(samples_.begin(), samples_.begin()+n, samples_.end());
        median = *(samples_.begin()+n);
    }
    else {
        //number is even:
        size_t n2 = numSamples_/2;
        size_t n1 = n2-1;
        nth_element(samples_.begin(), samples_.begin()+n1, samples_.end());
        nth_element(samples_.begin()+n2, samples_.begin()+n2, samples_.end());
        float m1 = *(samples_.begin()+n1);
        float m2 = *(samples_.begin()+n2);
        median = (m1 + m2) * 0.5f;
    }

    return median;
}

float Statistics::getQ1() const {
    if(!collectSamples_) {
        tgtAssert(false, "collectSamples needs to be true in order to calculate Q3");
        return 0.0f;
    }

    if((numSamples_ % 4) != 0) {
        size_t n = numSamples_/4;
        nth_element(samples_.begin(), samples_.begin()+n, samples_.end());
        return *(samples_.begin()+n);
    }
    else {
        //number is even:
        size_t n2 = numSamples_/4;
        size_t n1 = n2-1;

        nth_element(samples_.begin(), samples_.begin()+n1, samples_.end());
        nth_element(samples_.begin()+n2, samples_.begin()+n2, samples_.end());
        float q1_1 = *(samples_.begin()+n1);
        float q1_2 = *(samples_.begin()+n2);
        return (q1_1 + q1_2) * 0.5f;
    }
}

float Statistics::getQ3() const {
    if(!collectSamples_) {
        tgtAssert(false, "collectSamples needs to be true in order to calculate Q3");
        return 0.0f;
    }

    if(((numSamples_*3) % 4) != 0) {
        size_t n = (numSamples_*3)/4;
        nth_element(samples_.begin(),
                samples_.begin()+n,
                samples_.end());
        return *(samples_.begin()+n);
    }
    else {
        //number is even:
        size_t n2 = (numSamples_*3)/4;
        size_t n1 = n2-1;

        nth_element(samples_.begin(),
                samples_.begin()+n1,
                samples_.end());
        nth_element(samples_.begin()+n2,
                samples_.begin()+n2,
                samples_.end());
        float q3_1 = *(samples_.begin()+n1);
        float q3_2 = *(samples_.begin()+n2);
        return (q3_1 + q3_2) * 0.5f;
    }
}

} // namespace
