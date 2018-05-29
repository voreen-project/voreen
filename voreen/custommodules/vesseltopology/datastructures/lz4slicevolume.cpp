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

#include "lz4slicevolume.h"
#include "voreen/core/voreenapplication.h"

#include <fstream>

namespace voreen {

const static std::string METADATA_ROOT_NODE_STRING = "metadata";

LZ4SliceVolumeMetadata::LZ4SliceVolumeMetadata(tgt::svec3 dimensions)
    : dimensions_(dimensions)
{
}

LZ4SliceVolumeMetadata::LZ4SliceVolumeMetadata(const std::string& xmlfile)
    : dimensions_(tgt::svec3::one)
{
    XmlDeserializer ds;
    std::ifstream filestream(xmlfile);
    ds.read(filestream);
    ds.deserialize(METADATA_ROOT_NODE_STRING, *this);
}

void LZ4SliceVolumeMetadata::serialize(Serializer& s) const {
    s.serialize("dimensions", dimensions_);
}
void LZ4SliceVolumeMetadata::deserialize(Deserializer& s) {
    s.deserialize("dimensions", dimensions_);
}

void LZ4SliceVolumeMetadata::save(const std::string& xmlfile) const {
    XmlSerializer ser;
    std::ofstream filestream(xmlfile);
    ser.serialize(METADATA_ROOT_NODE_STRING, *this);
}

LZ4SliceVolume<uint8_t> binarizeVolume(const VolumeBase& volume, float binarizationThresholdSegmentationNormalized) {
    const auto dimensions = volume.getDimensions();

    LZ4SliceVolumeBuilder<uint8_t> builder(VoreenApplication::app()->getUniqueTmpFilePath(".lz4vol"), LZ4SliceVolumeMetadata(dimensions));

    for(size_t z = 0; z<dimensions.z; ++z) {
        //progress.setProgress(static_cast<float>(z)/dimensions.z);

        std::unique_ptr<const VolumeRAM> inSlice(volume.getSlice(z));
        auto outSlice(builder.getNextWritableSlice());

        for(size_t y = 0; y<dimensions.y; ++y) {
            for(size_t x = 0; x<dimensions.x; ++x) {
                if(inSlice->getVoxelNormalized(x,y,0) > binarizationThresholdSegmentationNormalized) {
                    outSlice->voxel(x,y,0) = 1;
                } else {
                    outSlice->voxel(x,y,0) = 0;
                }
            }
        }
    }

    return LZ4SliceVolumeBuilder<uint8_t>::finalize(std::move(builder));
}

}
