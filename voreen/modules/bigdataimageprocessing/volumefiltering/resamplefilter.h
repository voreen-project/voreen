/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_RESAMPLEFILTER_H
#define VRN_RESAMPLEFILTER_H

#include "volumefilter.h"

namespace voreen {

class ResampleFilter : public VolumeFilter {
public:
    ResampleFilter(const tgt::svec3& targetDimensions, size_t numChannels);

    virtual ~ResampleFilter();


    std::unique_ptr<VolumeRAM> getFilteredSlice(const CachingSliceReader* src, int z) const;
    int zExtent() const;
    boost::optional<tgt::svec3> getOverwrittenDimensions() const;

    size_t getNumInputChannels() const;
    size_t getNumOutputChannels() const;

private:
    const tgt::svec3 dimensions_;
    //const SamplingStrategy<float> samplingStrategy_;
    const size_t numChannels_;
};


} // namespace voreen

#endif // VRN_RESAMPLEFILTER_H
