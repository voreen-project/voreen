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

#include "resamplefilter.h"
#include "slicereader.h"

namespace voreen {

ResampleFilter::ResampleFilter(const tgt::svec3& targetDimensions, size_t numChannels)
    : dimensions_(targetDimensions)
    , numChannels_(numChannels)
{
}

ResampleFilter::~ResampleFilter() {
}
boost::optional<tgt::svec3> ResampleFilter::getOverwrittenDimensions() const {
    return dimensions_;
}

std::unique_ptr<VolumeRAM> ResampleFilter::getFilteredSlice(const CachingSliceReader* src, int z) const {
    tgt::ivec3 dim = dimensions_;
    tgtAssert(z >= 0 && z<dim.z, "Invalid z pos in slice request");

    VolumeFactory volumeFactory;
    std::string format = volumeFactory.getFormat(src->getMetaData().getBaseType(), numChannels_);
    std::unique_ptr<VolumeRAM> outputSlice(volumeFactory.create(format, tgt::svec3(dim.xy(), 1)));

    tgt::vec3 thisToBaseScale(tgt::vec3(src->getDimensions()) / tgt::vec3(dimensions_));
    tgt::vec3 thisToBaseOffset(thisToBaseScale * tgt::vec3(0.5) - tgt::vec3(0.5));

    for(int y = 0; y < dim.y; ++y) {
        for(int x = 0; x < dim.x; ++x) {
            tgt::ivec3 srcPos = tgt::round(thisToBaseScale * tgt::vec3(x, y, z) + thisToBaseOffset);
            for(size_t channel = 0; channel < numChannels_; channel++) {
                float val = src->getVoxelNormalized(srcPos, channel);
                outputSlice->setVoxelNormalized(val, x, y, 0, channel);
            }
        }
    }

    return outputSlice;
}

int ResampleFilter::zExtent() const {
    return 0;
}

size_t ResampleFilter::getNumInputChannels() const {
    return numChannels_;
}

size_t ResampleFilter::getNumOutputChannels() const {
    return numChannels_;
}


} // namespace voreen
