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

#include "synth2dreader.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/vector.h"

#include "voreen/core/io/textfilereader.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "modules/core/io/rawvolumereader.h"

using tgt::ivec3;

namespace voreen {

const std::string Synth2DReader::loggerCat_ = "voreen.io.VolumeReader.synth2d";

Synth2DReader::Synth2DReader(ProgressBar* progress)
    : VolumeReader(progress) {
    extensions_.push_back("vol");
}

VolumeList* Synth2DReader::read(const std::string &fileName) {
    RawVolumeReader rawReader(getProgressBar());

    FILE* fin = fopen(fileName.c_str(), "rb");
    char buf[4096];
    fread(buf, 4096, 1, fin);
    fclose(fin);

    Synth2DVolumeHeader* header = (Synth2DVolumeHeader *)buf;
    if ((header->magic[0] != 'V') || (header->magic[1] != 'O') || (header->magic[2] != 'L') || (header->magic[3] != 'U') || (header->version != 4)
                                  || (header->bytesPerChannel != 1) || (header->numChannels != 1 && header->numChannels != 3 && header->numChannels != 4))
        throw tgt::CorruptedFileException("error while reading data", fileName);

    std::string om("I");
    if(header->numChannels == 3)
        om = std::string("RGB");
    else if(header->numChannels == 4)
        om = std::string("RGBA");

    ivec3 dims(header->volSize);
    rawReader.setReadHints(dims,                     // dimensions of the volume
                           ivec3(1, 1, 1),           // thickness of one slice
                           om,                       // intensity, rgb or rgba image
                           "UCHAR",                  // one unsigned char per voxel
                           1,                        // number of time frames
                           4096);                    // header skip

    VolumeList* volumeList = rawReader.read(fileName);

    return volumeList;
}

VolumeReader* Synth2DReader::create(ProgressBar* progress) const {
    return new Synth2DReader(progress);
}

} // namespace voreen
