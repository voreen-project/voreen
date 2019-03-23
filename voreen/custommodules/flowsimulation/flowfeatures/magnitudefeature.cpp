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

#include "magnitudefeature.h"

#include "custommodules/bigdataimageprocessing/volumefiltering/slicereader.h"

namespace voreen {

MagnitudeFeature::MagnitudeFeature(const std::string& sliceBaseType)
    : sliceBaseType_(sliceBaseType)
{
}

MagnitudeFeature::~MagnitudeFeature() {
}

int MagnitudeFeature::zExtent() const {
    return 1;
}

const std::string& MagnitudeFeature::getSliceBaseType() const {
    return sliceBaseType_;
}

std::unique_ptr<VolumeRAM> MagnitudeFeature::getFilteredSlice(const CachingSliceReader* src, int z) const {
    tgtAssert(z >= 0 && z < src->getSignedDimensions().z, "Invalid z pos in slice request");

    const tgt::ivec3& dim = src->getSignedDimensions();
    std::unique_ptr<VolumeRAM> outputSlice(VolumeFactory().create(sliceBaseType_, tgt::svec3(dim.xy(), 1)));

    #pragma omp parallel for
    for (int y = 0; y < dim.y; ++y) {
        for (int x = 0; x < dim.x; ++x) {
            tgt::vec3 velocity;
            for(size_t channel = 0; channel < getNumInputChannels(); channel++) {
                velocity.elem[channel] = src->getVoxelNormalized(tgt::ivec3(x, y, z), channel);
            }
            outputSlice->setVoxelNormalized(tgt::length(velocity), x, y, 0);
        }
    }
    return outputSlice;
}

} // namespace voreen
