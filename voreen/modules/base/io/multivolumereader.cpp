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

#include "multivolumereader.h"

#include <fstream>
#include <iostream>

#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"

using std::string;

namespace {
    std::string extractPurePathFromFileName(const std::string& fileName) {
        if (fileName.rfind("/") != string::npos)
            return fileName.substr(0, fileName.rfind("/") + 1);
        else if (fileName.rfind("\\") != string::npos)
            return fileName.substr(0, fileName.rfind("\\") + 1);
        else
            return "/";
    }
}
namespace voreen {

const string MultiVolumeReader::loggerCat_ = "voreen.base.MultiVolumeReader";

MultiVolumeReader::MultiVolumeReader(VolumeSerializerPopulator* populator, ProgressBar* progress)
    : VolumeReader(progress),
      populator_(populator)
{
    extensions_.push_back("mv");
    protocols_.push_back("mv");
}

VolumeReader* MultiVolumeReader::create(ProgressBar* progress) const {
    return new MultiVolumeReader(0, progress);
}

VolumeList* MultiVolumeReader::read(const std::string& url) {
    LINFO("Loading multi volume file " << url);
    VolumeURL urlOrigin(url);

    std::vector<VolumeURL> origins = listVolumes(url);
    if (origins.empty())
        throw tgt::FileException("No volumes listed in multi-volume file", url);

    VolumeList* volumeList = new VolumeList();

    std::string refFile = urlOrigin.getSearchParameter("file");
    if (refFile == "") {
        // no particular file specified in URL => load all listed ones
        for (size_t i=0; i<origins.size(); i++) {
            VolumeBase* handle = read(origins.at(i));
            if (handle)
                volumeList->add(handle);
        }
    }
    else {
        // load specified file
        for (size_t i=0; i<origins.size(); i++) {
            if (origins.at(i).getSearchParameter("file") == refFile) {
                VolumeBase* handle = read(origins.at(i));
                if (handle) {
                    volumeList->add(handle);
                    break;
                }
            }
        }

        if (volumeList->empty()) {
            delete volumeList;
            throw tgt::FileException("File '" + refFile + "' not listed in multi-volume file", urlOrigin.getPath());
        }
    }

    return volumeList;
}

VolumeBase* MultiVolumeReader::read(const VolumeURL& origin) {
    std::string path = extractPurePathFromFileName(origin.getPath());
    std::string refFile = origin.getSearchParameter("file");
    if (refFile == "")
        throw tgt::FileException("'file' parameter missing in URL", origin.getURL());

    VolumeURL refOrigin(path + refFile);
    LINFO("Loading file " + refOrigin.getPath());
    VolumeBase* handle = 0;
    if (populator_)
        handle = populator_->getVolumeSerializer()->read(refOrigin);
    else {
        VolumeSerializerPopulator populator;
        handle = populator.getVolumeSerializer()->read(refOrigin);
    }
    if (handle)
        handle->setOrigin(origin);

    return handle;
}

std::vector<VolumeURL> MultiVolumeReader::listVolumes(const std::string& url) const {
    std::vector<VolumeURL> result;

    VolumeURL urlOrigin(url);
    std::string fileName = urlOrigin.getPath();

    std::fstream fin(fileName.c_str(), std::ios::in | std::ios::binary);
    if (!fin.good() || fin.eof() || !fin.is_open())
        throw tgt::FileNotFoundException("Unable to open multi volume file for reading", fileName);

    do {
        string line;
        getline(fin, line);

        if (line.empty())
            continue;

        // If the line was delimited by a '\r\n' the '\r' will still be the last character
        // so remove it
        if (line[line.length()-1] == char(13))
            line = line.substr(0, line.length()-1);

        VolumeURL origin("mv", fileName);
        origin.addSearchParameter("file", line);
        origin.getMetaDataContainer().addMetaData("File Name", new StringMetaData(line));
        result.push_back(origin);
    } while (!fin.eof());

    fin.close();

    return result;
}

} // namespace voreen
