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

#ifndef VRN_GAUSSIANFILTER_H
#define VRN_GAUSSIANFILTER_H

#include "volumefilter.h"

namespace voreen {

class GaussianFilter : public VolumeFilter {
public:
    static int suitableExtent(float standardDeviation);
    static tgt::ivec3 suitableExtent(const tgt::vec3& standardDeviation);

    static float suitableStandardDeviation(int extent);
    static tgt::vec3 suitableStandardDeviation(const tgt::ivec3& extent);


    GaussianFilter(const tgt::vec3& standardDeviation, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels = 1);
    GaussianFilter(float standardDeviation, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels = 1);

    GaussianFilter(const tgt::ivec3& extent, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels = 1);
    GaussianFilter(int extent, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels = 1);

    GaussianFilter(const tgt::vec3& standardDeviation, const tgt::ivec3& zExtent, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels = 1);
    GaussianFilter(float standardDeviation, int extent, const SamplingStrategy<float>& samplingStrategy, const std::string& sliceBaseType, size_t numChannels = 1);

    virtual ~GaussianFilter();


    std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const;
    int zExtent() const;
    const std::string& getSliceBaseType() const;

    size_t getNumInputChannels() const;
    size_t getNumOutputChannels() const;

private:
    float getKernelValX(int centeredPos) const;
    float getKernelValY(int centeredPos) const;
    float getKernelValZ(int centeredPos) const;

    const tgt::ivec3 neighborhoodDimensions_;
    const tgt::ivec3 kernelDimensions_;
    float* halfKernelX_;
    float* halfKernelY_;
    float* halfKernelZ_;

    const SamplingStrategy<float> samplingStrategy_;
    const std::string sliceBaseType_;
    size_t numChannels_;
};


} // namespace voreen

#endif // VRN_GAUSSIANFILTER_H
