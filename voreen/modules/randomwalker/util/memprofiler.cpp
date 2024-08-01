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

#include "memprofiler.h"
#include <algorithm>
#include "tgt/assert.h"

namespace voreen {

#ifdef VRN_MODULE_RANDOMWALKER_PROFILING_ENABLED
ProfileDataCollector::ProfileDataCollector()
    : current_(0)
    , maximum_(0)
{
}


void ProfileDataCollector::allocate(size_t size) const {
    size_t oldval = current_.fetch_add(size);
    size_t newval = oldval + size;

    size_t prev = maximum_;
    while(prev < newval && !maximum_.compare_exchange_weak(prev, newval)) {}
}
void ProfileDataCollector::deallocate(size_t size) const {
    tgtAssert(current_ >= size, "Invalid deallocation");
    current_.fetch_sub(size);
}
size_t ProfileDataCollector::peak() const {
    return maximum_;
}
void ProfileDataCollector::reset() const {
    current_ = 0;
    maximum_ = 0;
}

ProfileAllocation::ProfileAllocation(const ProfileDataCollector& collector, size_t size)
    : collector_(collector)
    , size_(size)
{
    collector_.allocate(size);
}

ProfileAllocation::ProfileAllocation(ProfileAllocation&& other)
    : collector_(other.collector_)
    , size_(other.size_)
{
    other.size_ = 0;
}
ProfileAllocation& ProfileAllocation::operator=(ProfileAllocation&& other) {
    if(&other != this) {
        this->~ProfileAllocation();
        new(this) ProfileAllocation(std::move(other));
    }
    return *this;

}

ProfileAllocation::~ProfileAllocation() {
    release();
}

void ProfileAllocation::release() {
    collector_.deallocate(size_);
    size_ = 0;
}

void ProfileAllocation::increaseBy(size_t size) {
    collector_.deallocate(size_);
    size_ += size;
    collector_.allocate(size_);
}

ProfileAllocation cgSystemFloatEllpack(const ProfileDataCollector& collector, size_t numEquations) {
    size_t numEntries = 7 * numEquations;

    size_t size = 0;

    size += numEntries * sizeof(float); //matbuf
    size += numEntries * sizeof(size_t); //matindicesbuf

    size += numEquations * sizeof(float); //tmpbuf
    size += numEquations * sizeof(float); //xbuf
    size += numEquations * sizeof(float); //rbuf
    size += numEquations * sizeof(float); //pbuf

    // Preconditioning
    size += numEquations * sizeof(float); //procond matbuf
    size += numEquations * sizeof(size_t); //precond matindicesbuf
    size += numEquations * sizeof(size_t); //zbuf (preconditioning)

    return ProfileAllocation(collector, size);
}

ProfileAllocation cgSystemFloatCSR(const ProfileDataCollector& collector, size_t numEquations, size_t numEntries) {
    size_t size = 0;

    size += numEntries * sizeof(float); //matbuf
    size += numEntries * sizeof(size_t); //colindex
    size += (numEquations+1) * sizeof(size_t); //rowindex

    size += numEquations * sizeof(float); //tmpbuf
    size += numEquations * sizeof(float); //xbuf
    size += numEquations * sizeof(float); //rbuf
    size += numEquations * sizeof(float); //pbuf

    // Preconditioning
    size += numEquations * sizeof(float); //matbuf
    size += numEquations * sizeof(size_t); //colindex
    size += (numEquations+1) * sizeof(size_t); //rowindex
    size += numEquations * sizeof(size_t); //zbuf (preconditioning)

    return ProfileAllocation(collector, size);
}

#endif

}
