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

#include "lz4slicevolumefilereader.h"
#include "../datastructures/lz4slicevolume.h"


#include "tgt/exception.h"
#include "tgt/assert.h"
#include "tgt/vector.h"

#include <boost/thread.hpp>


// ---------------------------------------------------------------------------------------------------------------------

namespace voreen {

const std::string LZ4SliceVolumeFileReader::loggerCat_ = "voreen.lz4.LZ4SliceVolumeFileReader";

LZ4SliceVolumeFileReader::LZ4SliceVolumeFileReader() : VolumeReader()
{
    extensions_.push_back(LZ4SliceVolumeBase::FILE_EXTENSION);

    protocols_.push_back(LZ4SliceVolumeBase::FILE_EXTENSION);
}

VolumeList* LZ4SliceVolumeFileReader::read(const std::string &url) {
    std::vector<VolumeURL> urls = listVolumes(url);
    VolumeList* volumeList = new VolumeList();

    try {
        for(const auto& url : urls) {
            volumeList->add(read(url));
        }
    } catch(tgt::IOException e) {
        delete volumeList;
        throw e;
    }
    return volumeList;
}

VolumeBase* LZ4SliceVolumeFileReader::read(const VolumeURL& origin) {
    std::string fileName = origin.getPath();
    LINFO("Loading " << fileName);

    std::unique_ptr<LZ4SliceVolumeBase> fileVolume = LZ4SliceVolumeBase::open(fileName);
    return std::move(*fileVolume).toVolume().release();
}

std::vector<VolumeURL> LZ4SliceVolumeFileReader::listVolumes(const std::string& urlStr) const {
    std::vector<VolumeURL> urls;
    urls.push_back(VolumeURL(urlStr));
    return urls;
}

VolumeReader* LZ4SliceVolumeFileReader::create(ProgressBar* progress) const {
    return new LZ4SliceVolumeFileReader();
}

} // namespace voreen
