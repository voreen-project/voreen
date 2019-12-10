/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "mhdvolumereader.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/vector.h"
#include "tgt/filesystem.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

using tgt::vec3;
using tgt::ivec3;
using tgt::hor;
using tgt::lessThanEqual;
using std::string;

namespace voreen {

const std::string MhdVolumeReader::loggerCat_ = "voreen.io.VolumeReader.dat";

class MhdLineParser {
    public:
        MhdLineParser(string line, string fileName) : name_(""), fileName_(fileName) {
            std::vector<string> exploded = strSplit(line, ' ');
            if(exploded.size() >= 3) {
                if(exploded[1] != "=")
                    return;

                name_ = exploded[0];

                for(size_t i=2; i<exploded.size(); i++) {
                    if((exploded[i] != "") && (exploded[i] != " "))
                        data_.push_back(exploded[i]);
                }
            }
        }

        string getName() {
            return name_;
        }

        void checkNumData(size_t numExpected) {
            if(data_.size() != numExpected)
                throw tgt::CorruptedFileException("Error parsing " + getName() + ": Expected " + itos(numExpected) + " arguments, got " + itos(data_.size()) + ".", fileName_);
        }

        int getDataInt() {
            checkNumData(1);
            return stoi(data_[0]);
        }

        string getDataString() {
            checkNumData(1);
            return data_[0];
        }

        bool getDataBool() {
            checkNumData(1);
            string st = data_[0];
            if(st == "False" || st == "0" || st == "false")
                return false;
            else if(st == "True" || st == "1" || st == "true")
                return true;
            else
                throw tgt::CorruptedFileException("Error parsing " + getName() + ": Expected bool (False/True), got " + st + ".", fileName_);
        }

        ivec3 getDataIVec3() {
            checkNumData(3);
            return ivec3(stoi(data_[0]), stoi(data_[1]), stoi(data_[2]));
        }

        vec3 getDataVec3() {
            checkNumData(3);
            return vec3(stof(data_[0]), stof(data_[1]), stof(data_[2]));
        }

    protected:
        string name_;
        string fileName_;
        std::vector<string> data_;
};

MhdVolumeReader::MhdVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("mhd");
}

VolumeList* MhdVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    LINFO("Opening " << fileName);

    string rawFilename = "";
    ivec3 dimensions(0,0,0);
    vec3 offset(0.0f);
    vec3 spacing(0.0f);
    string voxelType = "";
    int numChannels = 1;
    int64_t headerSkip = 0;

    tgt::File* file = FileSys.open(fileName);
    if ((!file) || (!file->isOpen())) {
        LERROR("Cannot open file '" << fileName);
    }
    else {
        while(!file->eof()) {
            std::string line = file->getLine();
            line = strReplaceAll(line, "\r", "");
            line = strReplaceAll(line, "\n", "");

            MhdLineParser parsedLine(line, fileName);

            //Raw data info:
            if(parsedLine.getName() == "ObjectType") {
                if(parsedLine.getDataString() != "Image")
                    throw tgt::FileException("Only images supported.");
            }
            else if(parsedLine.getName() == "NDims") {
                if(parsedLine.getDataInt() != 3)
                    throw tgt::FileException("Only 3D supported.");
            }
            else if(parsedLine.getName() == "ElementNumberOfChannels") {
                numChannels = parsedLine.getDataInt();
            }
            else if(parsedLine.getName() == "DimSize") {
                LINFO("Dimensions: " << parsedLine.getDataIVec3());
                dimensions = parsedLine.getDataIVec3();
            }
            else if(parsedLine.getName() == "ElementType") {
                LINFO("ElementType: '" << parsedLine.getDataString() << "'");
                voxelType = parsedLine.getDataString();
            }
            else if(parsedLine.getName() == "ElementDataFile") {
                rawFilename = parsedLine.getDataString();

                if(rawFilename == "LIST")
                    throw tgt::FileException("Slice-list mode not supported!");
                if(rawFilename == "LOCAL") {
                    rawFilename = fileName;
                    headerSkip = -1;
                }
            }
            else if(parsedLine.getName() == "ElementByteOrderMSB") {
                if(parsedLine.getDataBool())
                    throw tgt::FileException("ElementByteOrderMSB: Expected False.");
            }
            else if(parsedLine.getName() == "HeaderSize") {
                headerSkip = parsedLine.getDataInt();
            }
            //Metadata:
            else if(parsedLine.getName() == "Offset") {
                offset = parsedLine.getDataVec3();
            }
            else if(parsedLine.getName() == "ElementSpacing") {
                if(spacing != vec3(0.0f)) {
                    LWARNING("ElementSize and ElementSpacing specified, using ElementSpacing as spacing.");
                }

                spacing = parsedLine.getDataVec3();
            }
            else if(parsedLine.getName() == "ElementSize") {
                //use only if ElementSpacing has not been specified:
                if(spacing == vec3(0.0f))
                    spacing = parsedLine.getDataVec3();
                else {
                    LWARNING("ElementSize and ElementSpacing specified, using ElementSpacing as spacing.");
                }

            }
            //TODO:
            //Modality
            //MET_STRING
            //One of enum type: MET_MOD_CT, MET_MOD_MR, MET_MOD_US… See MetaImageTypes.h

            //SequenceID
            //MET_INT_ARRAY[4]
            //Four values comprising a DICOM sequence: Study, Series, Image numbers

            //ElementMin
            //MET_FLOAT
            //Minimum value in the data

            //ElementMax
            //MET_FLOAT
            //Maximum value in the data

            //Error handling:
            else if(parsedLine.getName() == "") {
                LWARNING("Failure to parse line: " << line);
            }
            else {
                LWARNING("Unkown attribute: " << parsedLine.getName());
            }
        }

        file->close();
    }
    delete file;

    if(spacing == vec3(0.0f))
        spacing = vec3(1.0f);

    string voreenVoxelType = convertVoxelTypeString(voxelType, numChannels);
    if(voreenVoxelType == "")
        throw tgt::FileException("Unsupported voxel type: " + voxelType + ", " + itos(numChannels) + " channels");

    VolumeRepresentation* volume;
    string directory = tgt::FileSystem::dirName(fileName);
    string fullRawFilename = directory+"/"+rawFilename;
    if(!FileSys.fileExists(fullRawFilename))
        throw tgt::FileException("Raw file '" + fullRawFilename + "' does not exist!");

    volume = (VolumeRAM*) new VolumeDiskRaw(fullRawFilename , voreenVoxelType, dimensions, headerSkip);

    Volume* vh = new Volume(volume, spacing, offset);
    vh->setOrigin(origin);
    vh->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(vh->getBaseType()));

    VolumeList* vc = new VolumeList();
    vc->add(vh);

    return vc;
}

string MhdVolumeReader::convertVoxelTypeString(string voxelType, int numChannels) {
    if(numChannels == 1) {
        if(voxelType == "MET_USHORT")
            return "uint16";
        else if(voxelType == "MET_SHORT")
            return "int16";

        else if(voxelType == "MET_UCHAR")
            return "uint8";
        else if(voxelType == "MET_CHAR")
            return "int8";

        else if(voxelType == "MET_UINT")
            return "uint32";
        else if(voxelType == "MET_INT")
            return "int32";

        else if(voxelType == "MET_ULONG")
            return "uint64";
        else if(voxelType == "MET_LONG")
            return "int64";

        else if(voxelType == "MET_FLOAT")
            return "float";
        else if(voxelType == "MET_DOUBLE")
            return "double";
        else {
            LERROR("Unsupported voxel type: " << voxelType << ", " << numChannels << " channels");
            return "";
        }
    }
    else {
        LERROR("Unsupported voxel type: " << voxelType << ", " << numChannels << " channels");
        return "";
    }
}

VolumeReader* MhdVolumeReader::create(ProgressBar* progress) const {
    return new MhdVolumeReader(progress);
}

} // namespace voreen
