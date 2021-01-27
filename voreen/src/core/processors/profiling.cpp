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

#include "voreen/core/processors/profiling.h"
#include <iomanip>
#include "tgt/tgt_gl.h"
#include "tgt/logmanager.h"

#include <algorithm>

#ifdef _MSC_VER  // high-performance counter
#include <windows.h>
#include <winbase.h>
#endif

namespace voreen {

//const std::string ProfilingBlock::loggerCat_ = "voreen.ProfilingBlock";
const std::string PerformanceSample::loggerCat_ = "voreen.PerformanceSample";

PerformanceSample::PerformanceSample() : parent_(NULL), name_(""), time_(-1.0f) {
}

PerformanceSample::PerformanceSample(PerformanceSample* sample) : parent_(sample->getParent()), name_(sample->getName()), time_(sample->getTime()), measurements_(sample->getMeasurementCount()) {
}

PerformanceSample::PerformanceSample(PerformanceSample* parent, std::string name) : parent_(parent), name_(name), time_(-1.0f), measurements_(1) {
}

PerformanceSample* PerformanceSample::addChild(std::string name) {
    children_.push_back(PerformanceSample(this, name));
    return &(children_.back());
}

PerformanceSample* PerformanceSample::addChild(PerformanceSample sample) {
    children_.push_back(sample);
    return &(children_.back());
}

void PerformanceSample::print(int level, std::string recordName) const {
    std::string spaces = "";
    for(int j=0; j<level; j++)
        spaces += ' ';

    float percent = 100.0f;
    if(parent_)
        percent = 100.0f * time_ / parent_->getTime();

    float percentChildren = 100.0f * (getChildrenTime()) / time_;

    if (children_.size() > 0 && parent_)
        LINFO(spaces << recordName << name_ << " (" << children_.size() << " children): " << std::setprecision(10) << time_ << " secs (" << percent << "% of parent) (" << percentChildren << "% in children)");
    else if (children_.size() > 0 && !parent_)
        LINFO(spaces << recordName << name_ << " (" << children_.size() << " children): " << std::setprecision(10) << time_ << " secs (" << percentChildren << "% in children)");
    else if (children_.size() == 0 && parent_)
        LINFO(spaces << recordName << name_ << ": " << std::setprecision(10) << time_ << " secs (" << percent << "% of parent)");
    else
        LINFO(spaces << recordName << name_ << ": " << std::setprecision(10) << time_ << " secs");

    for(size_t i=0; i<children_.size(); i++)
        children_[i].print(level+1);
}

float PerformanceSample::getChildTime(std::string name) const {
    float time = 0.0f;
    for(size_t i=0; i<children_.size(); i++){
        if(children_[i].getName() == name)
            return children_[i].getTime();
    }
    return time;
}

float PerformanceSample::getChildrenTime() const {
    float time = 0.0f;
    for(size_t i=0; i<children_.size(); i++)
        time += children_[i].getTime();
    return time;
}

PerformanceSample* PerformanceSample::getChild(std::string name) {
    for (size_t i = 0; i < children_.size(); ++i) {
        if(children_[i].getName() == name)
            return &(children_[i]);
    }
    return NULL;
}

std::vector<PerformanceSample*> PerformanceSample::getChildren() {
    std::vector<PerformanceSample*> sampleChildren;
    for (size_t i = 0; i < children_.size(); ++i) {
        sampleChildren.push_back(&(children_[i]));
    }
    return sampleChildren;
}

void PerformanceSample::setTime(float time) {
    time_ = time;
}

float PerformanceSample::getTime() const {
    return time_;
}

std::string PerformanceSample::getName() const {
    return name_;
}

PerformanceSample PerformanceSample::operator+ (PerformanceSample other_sample) {
    PerformanceSample sum(this);
    std::vector<PerformanceSample*> this_children = this->getChildren();
    std::vector<PerformanceSample*> other_children = other_sample.getChildren();
    size_t min_children = std::min(this_children.size(), other_children.size());
    for(size_t i=0; i<min_children; i++){
        sum.addChild((*this_children[i]) + (*other_children[i]));
    }
    sum.setTime(this->getTime() + other_sample.getTime());
    sum.increaseMeasurementCount();
    return sum;
}

int PerformanceSample::getMeasurementCount(){
    return measurements_;
}

void PerformanceSample::setMeasurementCount(int val){
    measurements_ = val;
}

void PerformanceSample::increaseMeasurementCount(){
    measurements_++;
}

void PerformanceSample::normalize(){
    std::vector<PerformanceSample*> this_children = this->getChildren();
    for(size_t i=0; i<this_children.size(); i++){
        this_children[i]->normalize();
    }
    this->setTime(this->getTime() / (float)this->getMeasurementCount());
    this->setMeasurementCount(1);
}

//----------------------------------------------------------------

PerformanceRecord::PerformanceRecord() : current_(0) {
}

PerformanceRecord::~PerformanceRecord() {
    deleteSamples();
}

void PerformanceRecord::startBlock(const ProfilingBlock* const pb) {
    if(current_ == 0)
        current_ = new PerformanceSample(current_, pb->getName());
    else {
        current_ = current_->addChild(pb->getName());
    }
}

void PerformanceRecord::endBlock(const ProfilingBlock* const pb) {
    current_->setTime(pb->getTime());
    if(current_->getParent() == 0)
        samples_.push_back(current_);

    current_ = current_->getParent();
}

const std::vector<PerformanceSample*> PerformanceRecord::getSamples() const {
    return samples_;
}

PerformanceSample* PerformanceRecord::getLastSample() const {
    if(samples_.empty())
        return 0;
    else
        return samples_.back();
}

void PerformanceRecord::deleteSamples() {
    while(!samples_.empty()) {
        delete samples_.back();
        samples_.pop_back();
    }
}

void PerformanceRecord::setName(std::string str){
    name_ = str;
}

std::string PerformanceRecord::getName() const {
    return name_;
}

//----------------------------------------------------------------

ProfilingBlock::ProfilingBlock(std::string name, PerformanceRecord& pr) : name_(name), pr_(pr) {
    //LINFO("Starting Block " << name);
    pr_.startBlock((const ProfilingBlock* const) this);
    //glFinish();

#ifdef _MSC_VER  // use high-performance counter on windows, if available
    if (QueryPerformanceFrequency(&hpfTicksPerSecond_) && QueryPerformanceCounter(&hpfStart_)) {
        useHpf_ = true;
    }
    else { // hpf not available
        useHpf_ = false;
        start_ = clock();
    }
#else
    start_ = clock();
#endif
}

ProfilingBlock::~ProfilingBlock() {
    //glFinish();
#ifdef _MSC_VER
    if (useHpf_)
        QueryPerformanceCounter(&hpfEnd_);
    else
        end_ = clock();
#else
    end_ = clock();
#endif
    //LINFO("Finishing Block " << name_);
    pr_.endBlock((const ProfilingBlock* const)this);
}

float ProfilingBlock::getTime() const {
#ifdef _MSC_VER
    if (useHpf_)
        return ((hpfEnd_.QuadPart-hpfStart_.QuadPart) / static_cast<float>(hpfTicksPerSecond_.QuadPart));
    else
        return ((end_-start_) / static_cast<float>(CLOCKS_PER_SEC));
#else
    return ((end_-start_) / static_cast<float>(CLOCKS_PER_SEC));
#endif
}

std::string ProfilingBlock::getName() const {
    return name_;
}

} // namespace
