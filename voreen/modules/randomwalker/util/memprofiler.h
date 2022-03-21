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

#ifndef VRN_RANDOM_WALKER_MEM_PROFILER_H
#define VRN_RANDOM_WALKER_MEM_PROFILER_H


#include <cstddef>
#include <atomic>

#define VRN_MODULE_RANDOMWALKER_PROFILING_ENABLED 1


namespace voreen {

#ifdef VRN_MODULE_RANDOMWALKER_PROFILING_ENABLED
struct ProfileDataCollector {
    ProfileDataCollector();

    void allocate(size_t size) const;
    void deallocate(size_t size) const;
    size_t peak() const;
    void reset() const;

private:
    mutable std::atomic<size_t> current_;
    mutable std::atomic<size_t> maximum_;
};

struct ProfileAllocation {
    ProfileAllocation(const ProfileDataCollector& collector, size_t size);
    ProfileAllocation(ProfileAllocation&& other);
    ProfileAllocation& operator=(ProfileAllocation&& other);
    ProfileAllocation(ProfileAllocation& other) = delete;
    ProfileAllocation& operator=(ProfileAllocation& other) = delete;
    ~ProfileAllocation();
    void release();
    void increaseBy(size_t size);

private:
    const ProfileDataCollector& collector_;
    size_t size_;
};

ProfileAllocation cgSystemFloatEllpack(const ProfileDataCollector& collector, size_t numEquations);
ProfileAllocation cgSystemFloatCSR(const ProfileDataCollector& collector, size_t numEquations, size_t numEntries);

#else

struct ProfileDataCollector {
    ProfileDataCollector() {}

    void allocate(size_t size) const {}
    void deallocate(size_t size) const {}
    void reset() const {}
};

struct ProfileAllocation {
    ProfileAllocation(const ProfileDataCollector& collector, size_t size) {}
    ProfileAllocation(ProfileAllocation&& other) {}
    ProfileAllocation& operator=(ProfileAllocation&& other) {}
    ProfileAllocation(ProfileAllocation& other) = delete;
    ProfileAllocation& operator=(ProfileAllocation& other) = delete;
    ~ProfileAllocation() {}
    void increaseBy(size_t size) {}

    void release() {}
};

ProfileAllocation cgSystemFloatEllpack(const ProfileDataCollector& collector, size_t numEquations) {}
ProfileAllocation cgSystemFloatCSR(const ProfileDataCollector& collector, size_t numEquations, size_t numEntries) {}

#endif

}

#endif
