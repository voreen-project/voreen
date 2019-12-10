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

#ifndef VRN_PERFORMANCEMETRIC_H
#define VRN_PERFORMANCEMETRIC_H

#include <vector>
#include <string>

#include "tgt/vector.h"

namespace voreen{
/**
 * Messures performance for some events.
 * Multiple events are saved and statistics over these are used.
 */
class PerformanceMetric{
public:
    /**
     * Get high precision timing. This is not an absolute time, but only for
     * each program iteration. Timing should be in usec range
     * \return Time in secounds
     */
    static double getHighPrecisionTimer();


    explicit PerformanceMetric(int samples = 100);

    /**
     * Begins a measurement
     */
    void beginRun();
    /**
     * Ends a mesurement
     */
    double endRun();

    /**
     * Clear out all saved runs (for example if conditions have changed)
     */
    void clearRuns();


    // Standard staticial measures
    double getMeanTime() const;
    double getMedianTime() const;
    double getStandardDeviation() const;
    double getLastRun() const;

    /**
     * Number of current samples used for statistics
     */ 
    int getSampleCount() const;

    /**
     * All saved datapoints
     */
    std::vector<double> getAllRuns() const;
    
    /**
     * Summary as a human readable text
     */
    std::string getTextInfo() const;

private:
    std::vector<double> runs_;
    int begin_;
    int end_;
    double currentRunBegin_;
    
};
};
#endif 
