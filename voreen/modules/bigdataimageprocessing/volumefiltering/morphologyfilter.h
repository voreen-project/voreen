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

#ifndef VRN_MORPHOLOGYFILTER_H
#define VRN_MORPHOLOGYFILTER_H

#include "volumefilter.h"

namespace voreen {

enum MorphologyOperatorType {
    DILATION_T,
    EROSION_T,
};

enum MorphologyOperatorShape {
    CUBE_T,
    SPHERE_T,
};

class MorphologyFilter : public VolumeFilter {
public:
    MorphologyFilter(const tgt::ivec3& extent, MorphologyOperatorType type, MorphologyOperatorShape shape, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType);
    virtual ~MorphologyFilter() {}

    int zExtent() const;
    const std::string& getSliceBaseType() const;

    std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const;

    // For now we only support single channel volumes
    size_t getNumInputChannels() const { return 1; };
    size_t getNumOutputChannels() const { return 1; };

    MorphologyOperatorType getMorphologyOperatorType() const;
    MorphologyOperatorShape getMorphologyOperatorShape() const;

private:

    std::unique_ptr<VolumeRAM> getFilteredSliceCubeMorphology(const CachingSliceReader* src, int z) const;
    std::unique_ptr<VolumeRAM> getFilteredSliceSphereMorphology (const CachingSliceReader* src, int z) const;

    std::function<float(float, float)> morphFunc_;

    const tgt::ivec3 extent_;
    const MorphologyOperatorType type_;
    const MorphologyOperatorShape shape_;
    const SamplingStrategy<float> samplingStrategy_;
    const std::string sliceBaseType_;
};

} // namespace voreen

#endif // VRN_MORPHOLOGYFILTER_H
