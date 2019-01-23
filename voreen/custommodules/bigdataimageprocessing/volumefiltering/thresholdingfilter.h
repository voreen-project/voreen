/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_THRESHOLDINGFILTER_H
#define VRN_THRESHOLDINGFILTER_H

#include "volumefilter.h"

#include <functional>

namespace voreen {

enum ThresholdingStrategyType {
    LOWER_T,
    UPPER_T,
};

class ThresholdingFilter : public VolumeFilter {
public:

    ThresholdingFilter(float threshold, float replacement, ThresholdingStrategyType thresholdingStrategyType, const std::string& sliceBaseType);
    ThresholdingFilter(float threshold, ThresholdingStrategyType thresholdingStrategyType);
    virtual ~ThresholdingFilter();

    int zExtent() const;
    const std::string& getSliceBaseType() const;

    // For now we only support single channel volumes
    size_t getNumInputChannels() const { return 1; };
    size_t getNumOutputChannels() const { return 1; };

    std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const;

    ThresholdingStrategyType getThresholdingStrategyType() const;

private:

    std::function<bool(float, float)> strategy_;

    bool binarize_;
    const float threshold_;
    const float replacement_;
    const ThresholdingStrategyType thresholdingStrategyType_;
    const std::string sliceBaseType_;
};

} // namespace voreen

#endif // VRN_THRESHOLDINGFILTER_H
