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

#ifndef VRN_BINARYMEDIANFILTER_H
#define VRN_BINARYMEDIANFILTER_H

#include "volumefilter.h"

namespace voreen {

class BinaryMedianFilter : public VolumeFilter {
public:
    // Careful! SamplingStrategy::Set(x) might not yield the expected result für x > binarizationThreshold (due to separation)
    BinaryMedianFilter(const tgt::ivec3& extent, float binarizationThreshold, uint32_t objectVoxelThreshold, const SamplingStrategy<float>& samplingStrategy);
    BinaryMedianFilter(const tgt::ivec3& extent, float binarizationThreshold, const SamplingStrategy<float>& samplingStrategy);

    virtual ~BinaryMedianFilter();

    std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const;
    int zExtent() const;
    const std::string& getSliceBaseType() const;

    // For now we only support single channel volumes
    size_t getNumInputChannels() const { return 1; };
    size_t getNumOutputChannels() const { return 1; };

protected:
    bool isObjectValue(float input) const;

private:
    const tgt::ivec3 extent_;
    const tgt::ivec3 kernelDimensions_;

    const float binarizationThreshold_;
    const uint32_t objectVoxelThreshold_;

    const SamplingStrategy<float> samplingStrategy_;
    const static std::string SLICE_BASE_TYPE;
};


} // namespace voreen

#endif // VRN_BINARYMEDIANFILTER_H
