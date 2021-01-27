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

// include this before any windows headers
#include "ddsbase.h"

#include "pvmvolumereader.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/texturemanager.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

using std::string;
using tgt::Texture;

namespace voreen {

const std::string PVMVolumeReader::loggerCat_ = "voreen.pvm.PVMVolumeReader";

PVMVolumeReader::PVMVolumeReader(ProgressBar* progress)
  : VolumeReader(progress)
{
    extensions_.push_back("pvm");
}

namespace {

inline void endian_swap(uint16_t& x) {
    x = (x>>8) | (x<<8);
}

} // namespace

VolumeList* PVMVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    uint8_t* data;
    uint8_t* tmpData;

    unsigned int width, height, depth, components;
    float scalex, scaley, scalez;
    unsigned char *description;
    unsigned char *courtesy;
    unsigned char *parameter;
    unsigned char *comment;

    LINFO("Reading PVM volume " << fileName);

    /*
        TODO This subroutine returns an array created with malloc but it should
        be created with 'new[]' because this chunk of data will be deleted with 'delete[]'.
        This can cause hard to find errors.

        As a temporary workaround the data are copied over into a new array
        and the c-array is deleted with 'free'.
        Because of some c-pointer vodoo done in ddsbase.cpp free must be invoked
        after the use of all other returned pointers. (roland)
    */
    tmpData = readPVMvolume(const_cast<char*>(fileName.c_str()), getProgressBar(),
                            &width, &height, &depth, &components,
                            &scalex, &scaley, &scalez, &description, &courtesy,
                            &parameter, &comment);

    if (!tmpData) {
        LERROR("PVM Reading failed");
        return 0;
    }

    data = new uint8_t[width * height * depth * components];
    memcpy(data, tmpData, width * height * depth * components);

    VolumeRAM* dataset = 0;

    if (!data) {
        throw tgt::IOException();
    }
    else {
        LINFO("Size: " << width << " x " << height << " x " << depth);
        LINFO("Spacing: " << scalex << " x " << scaley << " x " << scalez);
        LINFO("Components: " << components);
        if (description)
            LINFO("Description: " << description);
        if (courtesy)
            LINFO("Courtesy: " << courtesy);
        if (parameter)
            LINFO("Parameter: " << parameter);
        if (comment)
            LINFO("Comment: " << comment);
        if (components == 1) {
            LINFO("Create 8 bit data set.");
            dataset = new VolumeRAM_UInt8(data, tgt::ivec3(width, height, depth));
        }
        else if (components == 2) {
            // the endianness conversion in ddsbase.cpp seem to be broken,
            // so we perform it here instead
            uint16_t* data16 = reinterpret_cast<uint16_t*>(data);
            int numElements = width * height * depth;
            for (int i=0; i < numElements; i++) {
                endian_swap(data16[i]);
            }

            dataset = new VolumeRAM_UInt16((uint16_t*)data,
                                       tgt::ivec3(width, height, depth));
        }
        else LERROR("Bit depth not supported.");
    }

    // TODO now it is safe to free
    free(tmpData);

    VolumeList* volumeList = new VolumeList();
    if (dataset) {
        Volume* volumeHandle = new Volume(dataset, tgt::vec3(scalex, scaley, scalez), tgt::vec3(0.0f));
        oldVolumePosition(volumeHandle);
        volumeHandle->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(volumeHandle->getBaseType()));
        volumeHandle->setOrigin(VolumeURL(fileName));
        volumeList->add(volumeHandle);
    }

    return volumeList;
}

VolumeReader* PVMVolumeReader::create(ProgressBar* progress) const {
    return new PVMVolumeReader(progress);
}

} // namespace voreen
