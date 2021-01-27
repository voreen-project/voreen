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

#include "tuvvolumereader.h"

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

using tgt::ivec3;
using tgt::vec3;

namespace voreen {

const std::string TUVVolumeReader::loggerCat_ = "voreen.io.VolumeReader.tuv";

VolumeList* TUVVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    LINFO("Reading file " << fileName);

    std::fstream fin(fileName.c_str(), std::ios::in | std::ios::binary);
    if (!fin.good())
        throw tgt::IOException();

    unsigned short dim[3];
    fin.read(reinterpret_cast<char*>(dim),6);
    ivec3 dimensions = ivec3(dim[0], dim[1], dim[2]);

    LINFO("Read 16 bit dataset");
    VolumeRAM_UInt16* dataset;
    try {
        dataset = new VolumeRAM_UInt16(dimensions);
    } catch (std::bad_alloc&) {
        throw; // throw it to the caller
    }

    fin.read(reinterpret_cast<char*>(dataset->voxel()), dataset->getNumBytes());

    if ( fin.eof() )
        throw tgt::CorruptedFileException();

    fin.close();

    VolumeList* volumeList = new VolumeList();
    Volume* volumeHandle = new Volume(dataset, vec3(1.0f), vec3(0.0f));
    oldVolumePosition(volumeHandle);
    volumeHandle->setOrigin(fileName);
    volumeHandle->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(volumeHandle->getBaseType()));
    volumeList->add(volumeHandle);

    return volumeList;
}

VolumeReader* TUVVolumeReader::create(ProgressBar* /*progress*/) const {
    return new TUVVolumeReader(/*progress*/);
}

} // namespace voreen
