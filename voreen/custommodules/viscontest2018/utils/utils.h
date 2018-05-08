/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_ENSEMBLE_UTILS_H
#define VRN_ENSEMBLE_UTILS_H

#include "voreen/core/voreencoreapi.h"

#include "../datastructures/ensembledataset.h"

#include <vector>
#include <map>

namespace voreen {

/**
 * Utility function mapping a value within range A to the equivalent value in range B.
 */
template<typename T, typename S>
S mapRange(const T& valA, const T& minA, const T& maxA, const S& minB, const S& maxB) {
    //tgtAssert(valA >= minA && valA <= maxA, "value out of range"); // value may lay outside intentionally!
    // Cast in receiver type.
    return S(minB + (maxB - minB) * (valA - minA) / (maxA - minA));
}

class VRN_CORE_API TimeStepMapper {
public:
    TimeStepMapper(const EnsembleDataset* dataset);
    ~TimeStepMapper();
    std::map<float, std::vector<std::pair<std::string, EnsembleDataset::TimeStep>>> getTimeSteps();

    float getMinDuration() const { return minDuration_; }
    float getMaxTime() const { return maxTime_; }
    int getMaxNumTimeSteps()  const { return maxNumTimeSteps_; }
    void iterate(const std::function <void(std::string, EnsembleDataset::TimeStep, float)>& callback);

private:
    const EnsembleDataset* dataset_;
    float minDuration_;
    float maxTime_;
    int maxNumTimeSteps_;
    std::map<float, std::vector<std::pair<std::string, EnsembleDataset::TimeStep>>> timeStepMap_;

    void init();
};

}

#endif
