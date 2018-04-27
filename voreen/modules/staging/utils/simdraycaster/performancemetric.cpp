/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "performancemetric.h"

#include "tgt/assert.h"

#include <algorithm>
#include <cmath>
#include <sstream>


#ifdef WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <Windows.h>
#endif 
#ifdef UNIX
#include <sys/time.h>
#endif

namespace voreen{
double PerformanceMetric::getHighPrecisionTimer(){
#ifdef WIN32
    LARGE_INTEGER counter;
    LARGE_INTEGER freq;
    QueryPerformanceCounter(&counter);
    QueryPerformanceFrequency(&freq);
    return (double)counter.QuadPart/(double)freq.QuadPart;
#elif UNIX
    // start using CLOCK_MONOTONIC
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double f =  1.0*tv.tv_sec+tv.tv_usec/1000000.0;
    return f;
#else
#error "getHighPrecissionTimer not implemented for this platform"
#endif
}

PerformanceMetric::PerformanceMetric(int samples /*= 100*/){
    runs_.resize(samples);
    std::fill(runs_.begin(), runs_.end(), 0.0);
    currentRunBegin_ = 0.0;
    clearRuns();

}

void PerformanceMetric::beginRun(){
    tgtAssert(currentRunBegin_ == 0.0f, "voreen.PerformanceMetric: Start new run while last run not finished");
    currentRunBegin_ = getHighPrecisionTimer();
}

double PerformanceMetric::endRun(){
    tgtAssert(currentRunBegin_ != 0.0, "voreen.PerformanceMetric:Ends run before run is started");
    double end = getHighPrecisionTimer();
    double diff = end-currentRunBegin_;
    tgtAssert(diff, "voreen.PerformanceMetric: Time of run is less than 0");
    runs_[end_] = diff;
    end_ = (end_+1)%runs_.size();
    if(begin_ == end_){
        begin_ = (begin_+1)%runs_.size();
    }
    currentRunBegin_ = 0.0f;
    return diff;

}

double PerformanceMetric::getMeanTime() const{
    tgtAssert(begin_ != end_, "voreen.PerformanceMetric: No runs yet");
    double sum = 0.0;
    for(int i = begin_; i != end_; i = (i+1)%runs_.size()){
        sum += runs_[i];
    }
    return sum/getSampleCount();
}

double PerformanceMetric::getMedianTime() const{
    tgtAssert(begin_ != end_, "voreen.PerformanceMetric: No runs yet");
    std::vector<double> runs = getAllRuns();
    std::sort(runs.begin(), runs.end(), std::less<double>());
    if (runs.size() % 2 == 0){
        double mid1 = runs[runs.size()/2];
        double mid2 = runs[(runs.size()+1)/2];
        return (mid1+mid2)/2;
    }else{
        return runs[runs.size()/2];
    }
}

double PerformanceMetric::getStandardDeviation() const{
    tgtAssert(begin_ != end_, "voreen.PerformanceMetric: No runs yet");
    int N = getSampleCount();
    if (N == 1) return 0;

    double sum = 0;
    double mean = getMeanTime();
    for(int i = begin_; i != end_; i = (i+1)%runs_.size()){
        double d = runs_[i]-mean;
        sum += d*d;
    }
    return std::sqrt(sum/(N-1));
}

double PerformanceMetric::getLastRun() const{
    tgtAssert(begin_ != end_, "voreen.PerformanceMetric: No runs yet");
    size_t idx = end_ > 0? end_-1: runs_.size()-1;
    return runs_[idx];
}

std::vector<double> PerformanceMetric::getAllRuns() const{
    std::vector<double> runs;
    // build linear buffer from the internal ringbuffer
    for(int i = begin_; i != end_; i = (i+1)%runs_.size()){
        runs.push_back(runs_[i]);
    }
    return runs;
}

std::string PerformanceMetric::getTextInfo() const{
    std::stringstream ss;
    ss << "Statistics form last " << getSampleCount() << " runs" << std::endl;
    ss << "Time of last run:    " << getLastRun()*1000 << "ms" << std::endl;
    ss << "Mean time of runs:   " << getMeanTime()*1000 << "ms" << std::endl;
    ss << "Median time of runs: " << getMedianTime()*1000 << "ms" << std::endl;
    ss << "Standard derivation: " << getStandardDeviation()*1000 << "ms" << std::endl;
    return ss.str();
}

int PerformanceMetric::getSampleCount() const{
    if (begin_ == 0) 
        return end_;
    else
        return static_cast<int>(runs_.size());
}

void PerformanceMetric::clearRuns(){
    std::fill(runs_.begin(), runs_.end(), 0.0);
    
    begin_ = 0;
    end_ = 0;
}

}
