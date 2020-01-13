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

#include "nrrdvolumereader.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/vector.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/io/textfilereader.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "modules/core/io/rawvolumereader.h"

using tgt::vec3;
using tgt::ivec3;

namespace voreen {

const std::string NrrdVolumeReader::loggerCat_ = "voreen.base.NrrdVolumeReader";


VolumeBase* NrrdVolumeReader::read(const VolumeURL& origin) {
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

VolumeList* NrrdVolumeReader::read(const std::string &url) {
    VolumeURL origin(url);
    int timeframe = -1;
    std::string tmp = origin.getSearchParameter("timeframe");
    if (! tmp.empty())
        timeframe = stoi(tmp);

    return read(url, timeframe);
}

VolumeList* NrrdVolumeReader::read(const std::string &url, int timeframe) {
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    std::string objectFilename = "";
    ivec3 resolution;
    vec3 sliceThickness(1.f);
    std::string format;
    std::string model;
    bool error = false;
    int numFrames = 1;
    int dimension = 3;

    // NNRD fields for voxel to world matrix
    vec3 spaceOrigin(0.f);
    vec3 dirX(1.f, 0.f, 0.f);
    vec3 dirY(0.f, 1.f, 0.f);
    vec3 dirZ(0.f, 0.f, 1.f);

    LINFO("NrrdVolumeReader: " << fileName);
    TextFileReader reader(fileName);
    reader.setSeparators(":\t\n\r");
    if (!reader) {
        throw tgt::FileNotFoundException("Unable to open", fileName);
    }

    std::string type;
    std::istringstream args;

    std::string magic = reader.getMagicNumber();

    if (magic != "NRRD0001" && magic != "NRRD0002" && magic != "NRRD0003" && magic != "NRRD0004" && magic != "NRRD0005") {
        LERROR("Format error");
        return 0;
    }

    while (reader.getNextLine(type, args, false)) {
        LINFO("Type : " << type);

        if (type == "datafile" || type == "data file") {
            args >> objectFilename;
            LINFO("Value: " << objectFilename);
        } else if (type == "dimension") {
            args >> dimension;
            LINFO("Value: " << dimension);
        } else if (type == "sizes") {
            args >> resolution[0];
            args >> resolution[1];
            args >> resolution[2];
            if (dimension == 4)
                args >> numFrames;
            LINFO("Value: " << resolution[0] << " x " <<
                  resolution[1] << " x " << resolution[2]);
        } else if (type == "spacings" || type == "thicknesses") {
            args >> sliceThickness[0] >> sliceThickness[1] >> sliceThickness[2];
            LINFO("Value: " << sliceThickness[0] << " " <<
                  sliceThickness[1] << " " << sliceThickness[2]);
        } else if (type == "type") {
            char whitespace;
            args.get(whitespace);
            getline(args, format);
            if (format == "uchar" ||
                format == "unsigned char" ||
                format == "uint8" ||
                format == "uint8_t") {
                format = "UCHAR";
                model = "I";
            } else if (format == "ushort" ||
                         format == "unsigned short" ||
                       format == "unsigned short int" ||
                       format == "uint16" ||
                       format == "uint16_t") {
                format = "USHORT";
                model = "I";
            } else if (format == "short" ||
                       format == "signed short" ||
                       format == "signed short int" ||
                       format == "short int" ||
                       format == "int16" ||
                       format == "int16_t") {
                    format = "SHORT";
                    model = "I";
            } else if (format == "float") {
                format = "FLOAT";
                model = "I";
            } else if (format == "uint" ||
                        format == "unsigned int" ||
                       format == "uint32" ||
                       format == "uint32_t") {
                format = "UCHAR";
                model = "RGBA";
            }
            LINFO("Value: " << format);
        } else if (type == "space origin") {
            bool success = reader.getVec3(args, spaceOrigin);
            if (success)
                LINFO("Value: " << spaceOrigin[0] << " x " << spaceOrigin[1] << " x " << spaceOrigin[2]);
            else {
                LERROR("Could not parse Space Origin field from NRRD header.");
                error = true;
            }
        } else if (type == "space directions") {
            bool success = reader.getVec3(args, dirX) && reader.getVec3(args, dirY) && reader.getVec3(args, dirZ);
            if (success)
                LINFO("Value: " << dirX << ", " << dirY << ", " << dirZ);
            else {
                LERROR("Could not parse Space Origin field from NRRD header.");
                error = true;
            }
        } else {
            LWARNING("Unknown type: " << type);
        }
        //... "encoding:" and other value-types should be added
        if (args.fail()) {
            LERROR("Format error");
            error = true;
        }
    }

    if (!error) {
        RawVolumeReader rawReader(getProgressBar());
        RawVolumeReader::ReadHints h(resolution, sliceThickness, model, format);

        // TODO:  does not work as expected - fix
        //        handle space orientation (RAS/LPS/etc.) as well
        // calculate voxel to world matrix as described in
        // http://www.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
        /*h.transformation_ = tgt::mat4(dirX.x, dirY.x, dirZ.x, spaceOrigin.x,
                                      dirX.y, dirY.y, dirZ.y, spaceOrigin.y,
                                      dirX.z, dirY.z, dirZ.z, spaceOrigin.z,
                                      0.f   , 0.f   , 0.f   , 1.f          );*/

        // do we have a relative path?
        if ((objectFilename.substr(0,1) != "/")  && (objectFilename.substr(0,1) != "\\") &&
            (objectFilename.substr(1,2) != ":/") && (objectFilename.substr(1,2) != ":\\"))
        {
            size_t p = fileName.find_last_of("\\");
            if (p == std::string::npos)
                p = fileName.find_last_of("/");

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

            VolumeList* collection = rawReader.read(objectFilename);
            if (!collection->empty()) {
                VolumeURL origin(fileName);
                origin.addSearchParameter("timeframe", itos(frame));
                Volume* vh = static_cast<Volume*>(collection->first());
                vh->setOrigin(origin);
                vh->setTimestep(static_cast<float>(frame));
                vh->setRealWorldMapping(RealWorldMapping::createDenormalizingMapping(vh->getBaseType()));
                oldVolumePosition(vh);
                toReturn->add(vh);
            }
            delete collection;
        }
        return toReturn;

    } else
        return 0;
}

VolumeReader* NrrdVolumeReader::create(ProgressBar* progress) const {
    return new NrrdVolumeReader(progress);
}

} // namespace voreen
