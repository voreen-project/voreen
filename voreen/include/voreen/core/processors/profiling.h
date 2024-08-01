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

#ifndef VRN_PROFILING_H
#define VRN_PROFILING_H

#ifdef WIN32
#include <windows.h>
#include <winbase.h>
#endif

#include <string>
#include <vector>
#include <stack>
#include <time.h>

#include "voreen/core/voreencoreapi.h"

namespace voreen {

class VRN_CORE_API PerformanceSample  {
public:
    PerformanceSample();
    PerformanceSample(PerformanceSample* sample);
    PerformanceSample(PerformanceSample* parent, std::string name);

    PerformanceSample* addChild(std::string name);
    PerformanceSample* addChild(PerformanceSample sample);
    void print(int level = 0, std::string recordName = "") const;

    void setTime(float time);
    float getTime() const;

    std::string getName() const;

    float getChildTime(std::string) const;
    float getChildrenTime() const;

    PerformanceSample* getChild(std::string);
    std::vector<PerformanceSample*> getChildren();
    PerformanceSample* getParent() const { return parent_; }

    PerformanceSample operator+ (PerformanceSample);

    int getMeasurementCount();
    void setMeasurementCount(int);
    void increaseMeasurementCount();

    void normalize();

protected:
    PerformanceSample* parent_;
    std::vector<PerformanceSample> children_;
    std::string name_;
    float time_;
    int measurements_;

    static const std::string loggerCat_;
};

class ProfilingBlock;

/**
 * @brief Holds profiling info for an object.
 */
class VRN_CORE_API PerformanceRecord {
    friend class ProfilingBlock;
public:
    PerformanceRecord();
    ~PerformanceRecord();

    const std::vector<PerformanceSample*> getSamples() const;
    PerformanceSample* getLastSample() const;
    void deleteSamples();
    void setName(std::string);
    std::string getName() const;

protected:
    void startBlock(const ProfilingBlock* const pb);
    void endBlock(const ProfilingBlock* const pb);

    PerformanceSample* current_;
    std::string name_;

    //history:
    std::vector<PerformanceSample*> samples_;
};

/**
 * @brief The constructor/destructor of this class is used to time blocks.
 *
 * Use macro PROFILING_BLOCK("NAME") to measure runtime of a block.
 */
class VRN_CORE_API ProfilingBlock {
public:
    ProfilingBlock(std::string name, PerformanceRecord& pr);
    ~ProfilingBlock();

    float getTime() const;
    std::string getName() const;
protected:
    std::string name_;
    PerformanceRecord& pr_;

    clock_t start_;
    clock_t end_;

#ifdef _MSC_VER  // use high-performance counter on windows systems
    bool useHpf_;
    LARGE_INTEGER hpfTicksPerSecond_;
    LARGE_INTEGER hpfStart_;
    LARGE_INTEGER hpfEnd_;
#endif

    //static const std::string loggerCat_;
};

#ifdef VRN_PRINT_PROFILING
#define PROFILING_BLOCK(name) \
    ProfilingBlock block(name, performanceRecord_);
#else
#define PROFILING_BLOCK(name)
#endif

} // namespace

#endif //VRN_PROFILING_H
