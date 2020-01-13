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

#include "inveonvolumereader.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

#include "voreen/core/datastructures/volume/operators/volumeoperatorswapendianness.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormirror.h"

using tgt::ivec3;
using std::string;

namespace voreen {

const std::string InveonVolumeReader::loggerCat_ = "voreen.io.VolumeReader.InveonVolumeReader";

//-----------------------------------------------------------------------------

InveonVolumeReader::InveonVolumeReader(ProgressBar* progress)
    : VolumeReader(progress) {
    extensions_.push_back("hdr");
    protocols_.push_back("inveon");
}

VolumeBase* InveonVolumeReader::read(const VolumeURL& origin) {
    VolumeBase* result = 0;

    int volumeId = -1;
    std::string tmp = origin.getSearchParameter("volumeId");
    if (! tmp.empty())
        volumeId = stoi(tmp);

    VolumeList* collection = read(origin.getPath(), volumeId);

    if (collection && collection->size() == 1) {
        result = collection->first();
    }
    else if (collection && collection->size() > 1) {
        while(!collection->empty()) {
           VolumeBase* vh = collection->first();
           collection->remove(vh);
           delete vh;
        }
        delete collection;
        throw tgt::FileException("Only one volume expected", origin.getPath());
    }

    delete collection;

    return result;
}

VolumeList* InveonVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);

    int volumeId = -1;
    std::string tmp = origin.getSearchParameter("volumeId");
    if (! tmp.empty())
        volumeId = stoi(tmp);

    return read(url, volumeId);
}

VolumeList* InveonVolumeReader::read(const std::string &url, int volumeId) {
    LINFO("Loading dataset " << url << " vid: " << volumeId);

    VolumeURL origin(url);
    std::string fileName = origin.getPath();
    InveonHeader hdr = readHeader(fileName);
    string rawFilename = FileSys.fullBaseName(fileName);

    VolumeList* vc = new VolumeList();

    VolumeRepresentation* volume = new VolumeDiskRaw(rawFilename, hdr.voreenVoxelType_, hdr.dims_, 0, hdr.bigEndian_);
    Volume* vh = new Volume(volume, hdr.spacing_, tgt::vec3(0.0f));

    {
        VolumeURL origin(fileName);
        origin.addSearchParameter("volumeId", itos(0));
        vh->setOrigin(origin);
        vh->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(vh->getBaseType()));

        vc->add(vh);
    }
    return vc;
}

InveonVolumeReader::InveonHeader InveonVolumeReader::readHeader(const std::string filename) {
    tgt::File* file = FileSys.open(filename);

    InveonHeader hdr;

    std::string line;
    while((line = file->getLine()) != "") {
        if(line[0] == '#')
            continue;

        size_t sep = line.find(' ');
        if(sep != string::npos) {
            string cmd = line.substr(0, sep);
            string data = line.substr(sep+1);

            if(cmd == "x_dimension")
                hdr.dims_.x = stoi(data);
            else if(cmd == "y_dimension")
                hdr.dims_.y = stoi(data);
            else if(cmd == "z_dimension")
                hdr.dims_.z = stoi(data);
            else if(cmd == "pixel_size_x")
                hdr.spacing_.x = stof(data);
            else if(cmd == "pixel_size_y")
                hdr.spacing_.y = stof(data);
            else if(cmd == "pixel_size_z")
                hdr.spacing_.z = stof(data);
            else if(cmd == "data_type") {
                int type = stoi(data);
                switch(type) {
                    case 0: throw tgt::UnsupportedFormatException("Unknown data_type", filename);
                        break;
                    case 1: hdr.voreenVoxelType_ = "int8";
                            hdr.bigEndian_ = false;
                            break;
                    case 2: hdr.voreenVoxelType_ = "int16";
                            hdr.bigEndian_ = false;
                            break;
                    case 3: hdr.voreenVoxelType_ = "int32";
                            hdr.bigEndian_ = false;
                            break;
                    case 4: hdr.voreenVoxelType_ = "float";
                            hdr.bigEndian_ = false;
                            break;
                    case 5: hdr.voreenVoxelType_ = "float";
                            hdr.bigEndian_ = true;
                            break;
                    case 6: hdr.voreenVoxelType_ = "int16";
                            hdr.bigEndian_ = true;
                            break;
                    case 7: hdr.voreenVoxelType_ = "int32";
                            hdr.bigEndian_ = true;
                            break;
                    default: throw tgt::CorruptedFileException("Unknown data_type", filename);
                            break;
                }
            }
            else if(cmd == "file_type") {
                if(stoi(data) != 5)
                    throw tgt::UnsupportedFormatException("Unknown file_type. Only image data (5) supported.", filename);
            }
            else if(cmd == "total_frames") {
                hdr.numFrames_ = stoi(data);
            }
            else if(cmd == "number_of_dimensions") {
                if(stoi(data) > 3)
                    throw tgt::UnsupportedFormatException("Unsupported number of dimensions.", filename);
            }
            else if(cmd == "end_of_header") {
                goto main_header_finished;
            }
        }
    }
main_header_finished:
    //TODO: frame headers
    delete file;
    return hdr;
}

std::vector<VolumeURL> InveonVolumeReader::listVolumes(const std::string& url) const {
    VolumeURL urlOrigin(url);
    std::string fileName = urlOrigin.getPath();

    InveonHeader hdr = readHeader(fileName);

    std::vector<VolumeURL> result;

    // temp hack for one volume:
    VolumeURL origin("inveon", fileName);
    origin.addSearchParameter("volumeId", itos(0));
    origin.getMetaDataContainer().addMetaData("volumeId", new IntMetaData(0));
    origin.getMetaDataContainer().addMetaData("Spacing", new Vec3MetaData(hdr.spacing_));
    origin.getMetaDataContainer().addMetaData("Dimensions", new IVec3MetaData(hdr.dims_));
    result.push_back(origin);

    //for(size_t i=0; i<s.subVolumes_.size(); i++) {
        //VolumeURL origin("ecat7", fileName);
        //origin.addSearchParameter("volumeId", itos(s.subVolumes_[i].de_.id_));
        //origin.getMetaDataContainer().addMetaData("volumeId", new IntMetaData(s.subVolumes_[i].de_.id_));
        //origin.getMetaDataContainer().addMetaData("Spacing", new Vec3MetaData(s.subVolumes_[i].getSpacing()));
        //origin.getMetaDataContainer().addMetaData("Dimensions", new IVec3MetaData(s.subVolumes_[i].getDimensions()));
        //s.transformMetaData(origin.getMetaDataContainer(), static_cast<int>(i));
        //result.push_back(origin);
    //}

    return result;
}

VolumeReader* InveonVolumeReader::create(ProgressBar* progress) const {
    return new InveonVolumeReader(progress);
}

} // namespace voreen
