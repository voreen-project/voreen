/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "rawvoxvolumereader.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <limits>

#include "tgt/exception.h"
#include "tgt/logmanager.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

using tgt::svec3;
using tgt::vec3;

namespace voreen {

const std::string RawVoxVolumeReader::loggerCat_ = "voreen.io.VolumeReader.rawvox";

struct RawVoxHeader {
    int32_t magic_;

    int32_t sizeX_;
    int32_t sizeY_;
    int32_t sizeZ_;

    int32_t bitsPerVoxel_;
};

VolumeList* RawVoxVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    LINFO("Reading file " << fileName);

    std::fstream fin(fileName.c_str(), std::ios::in | std::ios::binary);
    if (!fin.good())
        throw tgt::IOException();

    RawVoxHeader header;
    fin.read(reinterpret_cast<char*>(&header), sizeof(header));
    svec3 dimensions = svec3(header.sizeX_, header.sizeY_, header.sizeZ_);

    if(header.magic_ != 1381388120) {
        throw tgt::CorruptedFileException("Wrong magic number.");
    }

    VolumeRAM* dataset;
    switch(header.bitsPerVoxel_) {
        case 8:
            LINFO("Reading 8 bit dataset");
            dataset = new VolumeRAM_UInt8(dimensions);
            break;
        case 16:
            LINFO("Reading 16 bit dataset");
            dataset = new VolumeRAM_UInt16(dimensions);
            break;
        case 32:
            LINFO("Reading 32 bit (float) dataset");
            dataset = new VolumeRAM_Float(dimensions);
            break;
        default:
            LERROR("Unknown bpp!");
            throw tgt::UnsupportedFormatException("Unexpected bpp.");
    }

    fin.read(reinterpret_cast<char*>(dataset->getData()), dataset->getNumBytes());

    if ( fin.eof() ) {
        delete dataset;
        throw tgt::CorruptedFileException();
    }

    fin.close();

    VolumeList* volumeList = new VolumeList();
    Volume* volumeHandle = new Volume(dataset, vec3(1.0f), vec3(0.0f));
    oldVolumePosition(volumeHandle);
    volumeHandle->setOrigin(fileName);
    volumeHandle->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(volumeHandle->getBaseType()));
    volumeList->add(volumeHandle);

    return volumeList;
}

VolumeReader* RawVoxVolumeReader::create(ProgressBar* /*progress*/) const {
    return new RawVoxVolumeReader(/*progress*/);
}

} // namespace voreen
