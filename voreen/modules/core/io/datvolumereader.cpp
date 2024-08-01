/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "datvolumereader.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/vector.h"

#include "modules/core/io/rawvolumereader.h"
#include "voreen/core/io/textfilereader.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

using tgt::vec3;
using tgt::ivec3;
using tgt::hor;
using tgt::lessThanEqual;

namespace voreen {

const std::string DatVolumeReader::loggerCat_ = "voreen.DatVolumeReader";

DatVolumeReader::DatVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("dat");
}

std::string DatVolumeReader::getRelatedRawFileName(const std::string& fileName) {
    TextFileReader reader(fileName);
    if (! reader)
        return "";

    std::string type;
    std::string objectFilename = "";
    std::istringstream args;
    while (reader.getNextLine(type, args, false) == true) {
        if (type == "ObjectFileName:") {
            args >> objectFilename;
            break;
        }
    }
    return objectFilename;
}

VolumeBase* DatVolumeReader::read(const VolumeURL& origin) {
    VolumeBase* result = 0;

    int timeframe = -1;
    std::string tmp = origin.getSearchParameter("timeframe");
    if (! tmp.empty())
        timeframe = stoi(tmp);

    VolumeList* collection = read(origin.getPath(), timeframe);

    if (collection && collection->size() == 1) {
        result = collection->first();
    }
    else if (collection && collection->size() > 1) {
        delete collection;
        throw tgt::FileException("Only one volume expected", origin.getPath());
    }

    delete collection;

    return result;
}

VolumeList* DatVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    int timeframe = -1;
    std::string tmp = origin.getSearchParameter("timeframe");
    if (! tmp.empty())
        timeframe = stoi(tmp);

    return readSlices(url, 0, 0, timeframe);
}

VolumeList* DatVolumeReader::read(const std::string &url, int timeframe) {
    return readSlices(url, 0, 0, timeframe);
}

VolumeList* DatVolumeReader::readMetaFile(const std::string &fileName, size_t firstSlice, size_t lastSlice, int timeframe) {
    RawVolumeReader::ReadHints h;
    std::string objectFilename;
    vec3 sliceThickness = vec3(1.f, 1.f, 1.f);
    int numFrames = 1;
    std::string taggedFilename;
    int nbrTags;
    std::string objectType;
    std::string gridType;
    bool error = false;

    LINFO("Loading dat file " << fileName);
    TextFileReader reader(fileName);

    if (!reader)
        throw tgt::FileNotFoundException("reading dat file", fileName);

    std::string type;
    std::istringstream args;

    while (reader.getNextLine(type, args, false)) {
        if (type == "ObjectFileName:") {
            args >> objectFilename;
            LDEBUG(type << " " << objectFilename);
        } else if (type == "TaggedFileName:") {
            args >> taggedFilename;
            LDEBUG(type << " " << taggedFilename);
        } else if (type == "Resolution:") {
            args >> h.dimensions_.x;
            args >> h.dimensions_.y;
            args >> h.dimensions_.z;
            LDEBUG(type << " " << h.dimensions_);
        } else if (type == "NumFrames:") {
            args >> numFrames;
            LDEBUG(type << " " << numFrames);
        } else if (type == "SliceThickness:") {
            args >> sliceThickness.x >> sliceThickness.y >> sliceThickness.z;
            LDEBUG(type << " " << sliceThickness);
        } else if (type == "Format:") {
            args >> h.format_;
            LDEBUG(type << " " << h.format_);
        } else if (type == "NbrTags:") {
            args >> nbrTags;
            LDEBUG(type << " " << nbrTags);
        } else if (type == "ObjectType:") {
            args >> objectType;
            LDEBUG(type << " " << objectType);
        } else if (type == "ObjectModel:") {
            args >> h.objectModel_;
            LDEBUG(type << " " << h.objectModel_);
        } else if (type == "GridType:") {
            args >> gridType;
            LDEBUG(type << " " << gridType);
        } else if (type == "ByteOrder:") {
            std::string s;
            args >> s;
            if (s == "big-endian" || s == "bigendian" || s == "bigEndian")
                h.bigEndianByteOrder_ = true;
            else if (s == "little-endian" || s == "littleendian" || s == "littleEndian")
                h.bigEndianByteOrder_ = false;
            else
                LERROR("Unknown byte order: " << s);
        } else if (type == "TransformMatrix:") {
            // first argument contains number of row
            std::string row;
            args >> row;
            // last element is the number
            int index = atoi((row.substr(row.length()-1, row.length()-1)).c_str());
            args >> h.transformation_[index].x;
            args >> h.transformation_[index].y;
            args >> h.transformation_[index].z;
            args >> h.transformation_[index][3];
            if (index == 3) {
                LDEBUG("TransformMatrix:");
                LDEBUG("  " << h.transformation_[0]);
                LDEBUG("  " << h.transformation_[1]);
                LDEBUG("  " << h.transformation_[2]);
                LDEBUG("  " << h.transformation_[3]);
            }
        }
        else if (type == "Modality:") {
            std::string modalityStr;
            args >> modalityStr;
            LDEBUG(type << " " << modalityStr);
            h.modality_ = Modality(modalityStr);
        }
        else if (type == "Checksum:") {
            std::string checksumStr;
            args >> checksumStr;
            LDEBUG(type << " " << checksumStr);
            if(checksumStr.length() == 32)
                h.hash_ = checksumStr;
        }
        else if (type == "TimeStep:") {
            args >> h.timeStep_;
            LDEBUG(type << " " << h.timeStep_);
        } else if (type == "SliceOrder:") {
            args >> h.sliceOrder_;
            LDEBUG(type << " " << h.sliceOrder_);
        } else {
            LERROR("Unknown type: " << type);
        }

        if (args.fail()) {
            LERROR("Format error");
            error = true;
        }
    }

    // check whether necessary meta-data could be read
    if (objectFilename.empty()) {
        LERROR("No raw file specified");
        error = true;
    }

    if (hor(lessThanEqual(h.dimensions_, ivec3(0)))) {
        LERROR("Invalid resolution or resolution not specified: " << h.dimensions_);
        error = true;
    }

    h.spacing_ = sliceThickness;

    if (!error) {
        RawVolumeReader rawReader(getProgressBar());

        // do we have a relative path?
        if ((objectFilename.substr(0, 1) != "/")  && (objectFilename.substr(0, 1) != "\\") &&
            (objectFilename.substr(1, 2) != ":/") && (objectFilename.substr(1, 2) != ":\\"))
        {
            size_t p = fileName.find_last_of("\\/");
            // construct path relative to dat file
            objectFilename = fileName.substr(0, p + 1) + objectFilename;
        }

        int start = 0;
        int end = numFrames;
        if (timeframe != -1) {
            if (timeframe >= numFrames)
                throw tgt::FileException("Specified time frame not in volume", fileName);

            start = timeframe;
            end = timeframe+1;
        }

        VolumeList* toReturn = new VolumeList();
        for (int frame = start; frame < end; ++frame) {
            h.timeframe_ = frame;
            rawReader.setReadHints(h);

            VolumeList* volumeList = rawReader.readSlices(objectFilename, firstSlice, lastSlice);
            if (!volumeList->empty()) {
                VolumeURL origin(fileName);
                origin.addSearchParameter("timeframe", itos(frame));

                Volume* vh = static_cast<Volume*>(volumeList->first());
                vh->setOrigin(origin);
                vh->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(vh->getBaseType()));
                vh->setTimestep(static_cast<float>(frame));

                oldVolumePosition(vh);

                if(!h.hash_.empty())
                    vh->setHash(h.hash_);

                toReturn->add(volumeList->first());
            }
            delete volumeList;
        }
        return toReturn;
    }
    else {
        throw tgt::CorruptedFileException("error while reading data", fileName);
    }
}

VolumeList* DatVolumeReader::readSlices(const std::string &url, size_t firstSlice, size_t lastSlice, int timeframe) {
    VolumeURL origin(url);
    return readMetaFile(origin.getPath(), firstSlice, lastSlice, timeframe);
}

VolumeList* DatVolumeReader::readBrick(const std::string& url, tgt::ivec3 brickStartPos, int brickSize) {
    VolumeURL origin(url);
    return readMetaFileBrick(origin.getPath(), brickStartPos, brickSize);
}

VolumeList* DatVolumeReader::readVolumeFileBrick(const std::string &fileName, const tgt::ivec3& dims,
    tgt::ivec3 brickStartPos, int brickSize) {
    RawVolumeReader rawReader(getProgressBar());

    rawReader.setReadHints(dims,               // dimensions of the volume
                           ivec3(1, 1, 1),     // thickness of one slice
                           "I",                // intensity image
                           "USHORT",           // one unsigned short per voxel
                           1,                  // number of time frames
                           6);                 // header skip

    VolumeList* volumeList = rawReader.readBrick(fileName,brickStartPos, brickSize);

    return volumeList;
}

VolumeList* DatVolumeReader::readMetaFileBrick(const std::string &fileName,tgt::ivec3 brickStartPos,
    int brickSize) {
    RawVolumeReader::ReadHints h;
    std::string objectFilename;
    std::string taggedFilename;
    vec3 sliceThickness = vec3(1.f, 1.f, 1.f);

    int nbrTags;
    std::string objectType;
    std::string gridType;
    bool error = false;

    //LINFO("Loading file " << fileName);
    TextFileReader reader(fileName);

    if (!reader)
        throw tgt::FileNotFoundException("reading dat file", fileName);

    std::string type;
    std::istringstream args;
    while (reader.getNextLine(type, args, false)) {

        if (type == "ObjectFileName:") {
            args >> objectFilename;
        } else if (type == "TaggedFileName:") {
            args >> taggedFilename;
        } else if (type == "Resolution:") {
            args >> h.dimensions_.x >> h.dimensions_.y >> h.dimensions_.z;
        } else if (type == "SliceThickness:") {
            args >> sliceThickness.x >> sliceThickness.y >> sliceThickness.z;
        } else if (type == "Format:") {
            args >> h.format_;
        } else if (type == "NbrTags:") {
            args >> nbrTags;
        } else if (type == "ObjectType:") {
            args >> objectType;
        } else if (type == "ObjectModel:") {
            args >> h.objectModel_;
        } else if (type == "GridType:") {
            args >> gridType;
        } else if (type == "TransformMatrix:") {
            // first argument contains number of row
            std::string row;
            args >> row;
            // last element is the number
            int index = atoi((row.substr(row.length()-1, row.length()-1)).c_str());
            args >> h.transformation_[index].x;
            args >> h.transformation_[index].y;
            args >> h.transformation_[index].z;
            args >> h.transformation_[index][3];
        }
        else if (type == "Modality:") {
            std::string modalityStr;
            args >> modalityStr;
            h.modality_ = Modality(modalityStr);
        }
        else if (type == "TimeStep:") {
            args >> h.timeStep_;
        } else if (type == "ByteOrder:") {
            std::string s;
            args >> s;
            if (s == "big-endian" || s == "bigendian" || s == "bigEndian")
                h.bigEndianByteOrder_ = true;
            else if (s == "little-endian" || s == "littleendian" || s == "littleEndian")
                h.bigEndianByteOrder_ = false;
            else
                LERROR("Unknown byte order: " << s);
        } else {
            LERROR("Unknown type: " << type);
        }

        if (args.fail()) {
            LERROR("Format error");
            error = true;
        }
    }

    // check whether necessary meta-data could be read
    if (objectFilename.empty()) {
        LERROR("No raw file specified");
        error = true;
    }

    h.spacing_ = sliceThickness;

    if (hor(lessThanEqual(h.dimensions_, ivec3(0)))) {
        LERROR("Invalid resolution or resolution not specified: " << h.dimensions_);
        error = true;
    }

    if (!error) {
        RawVolumeReader rawReader(getProgressBar());
        rawReader.setReadHints(h);

        // do we have a relative path?
        if ((objectFilename.substr(0, 1) != "/")  && (objectFilename.substr(0, 1) != "\\") &&
            (objectFilename.substr(1, 2) != ":/") && (objectFilename.substr(1, 2) != ":\\"))
        {
            size_t p = fileName.find_last_of("\\/");
            // construct path relative to dat file
            objectFilename = fileName.substr(0, p + 1) + objectFilename;
        }

        VolumeList* volumeList = rawReader.readBrick(objectFilename, brickStartPos, brickSize);
        return volumeList;
    }
    else {
        throw tgt::CorruptedFileException("error while reading data", fileName);
    }
}

VolumeReader* DatVolumeReader::create(ProgressBar* progress) const {
    return new DatVolumeReader(progress);
}

} // namespace voreen
