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

#include "voreen/core/io/volumereader.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/io/progressbar.h"

#include "tgt/filesystem.h"

#include <fstream>

namespace voreen {

const std::string VolumeReader::loggerCat_("voreen.VolumeReader");

VolumeReader::VolumeReader(ProgressBar* progress /*= 0*/)
  : progress_(progress)
{}

VolumeList* VolumeReader::readSlices(const std::string&, size_t, size_t) {
    throw tgt::FileException("This VolumeReader can not read slice wise.");
}

VolumeBase* VolumeReader::read(const VolumeURL& origin) {

    VolumeBase* result = 0;

    VolumeList* collection = read(origin.getPath());

    if (collection && collection->size() == 1) {
        result = collection->first();
    }
    else if (collection && collection->size() > 1) {

        // Delete all contained volumes.
        for(size_t i = 0; i < collection->size(); i++)
            delete collection->at(i);
        delete collection;

        throw tgt::FileException("Only one volume expected", origin.getPath());
    }

    delete collection;

    return result;
}

VolumeList* VolumeReader::readBrick(const std::string& url, tgt::ivec3, int) {
    throw tgt::FileException("This file format does not support brick-wise reading of volume data.", url);
}

void VolumeReader::read(VolumeRAM* volume, FILE* fin) {
    if (progress_) {

        size_t max = tgt::max(volume->getDimensions());

        // validate dimensions
        if (max == 0 || max > 1e5) {
            LERROR("Invalid dimensions: " << volume->getDimensions());
            std::ostringstream stream;
            stream << volume->getDimensions();
            throw VoreenException("Invalid dimensions: " + stream.str());
        }

        // no remainder possible because getNumBytes is a multiple of max
        size_t sizeStep = volume->getNumBytes() / static_cast<size_t>(max);

        for (size_t i = 0; i < size_t(max); ++i) {
            if (fread(reinterpret_cast<char*>(volume->getData()) + sizeStep * i, 1, sizeStep, fin) == 0)
                LWARNING("fread() failed");
            progress_->setProgress(static_cast<float>(i) / static_cast<float>(max));
        }
    }
    else {
        if (fread(reinterpret_cast<char*>(volume->getData()), 1, volume->getNumBytes(), fin) == 0)
        LWARNING("fread() failed");
    }
}

std::vector<VolumeURL> VolumeReader::listVolumes(const std::string& url) const {
    std::vector<VolumeURL> result;

    VolumeURL urlOrigin(url);
    if (tgt::FileSystem::fileExists(urlOrigin.getPath()))
        result.push_back(urlOrigin);

    return result;
}

void VolumeReader::reverseXSliceOrder(VolumeRAM* const volume) const {
    if (volume == 0)
        return;

    size_t bytesPerVoxel = volume->getBytesPerVoxel();
    const tgt::ivec3 dim = volume->getDimensions();

    typedef unsigned char BYTE;
    BYTE* const data = reinterpret_cast<BYTE* const>(volume->getData());
    BYTE* aux = new BYTE[bytesPerVoxel];

    for (int x = 0; x < (dim.x / 2); ++x) {
        int _x = (dim.x - (x + 1));
        for (int z = 0; z < dim.z; ++z) {
            int a = (z * (dim.x * dim.y));
            for (int y = 0; y < dim.y; ++y) {
                int b = (y * dim.x);
                int l = a + b + x;
                int r = a + b + _x;
                memcpy(aux, (data + (l * bytesPerVoxel)), bytesPerVoxel);
                memcpy((data + (l * bytesPerVoxel)), (data + (r * bytesPerVoxel)), bytesPerVoxel);
                memcpy((data + (r * bytesPerVoxel)), aux, bytesPerVoxel);
            }
        }
    }
    delete [] aux;
}

void VolumeReader::reverseYSliceOrder(VolumeRAM* const volume) const {
    if (volume == 0)
        return;

    size_t bytesPerVoxel = volume->getBytesPerVoxel();
    const tgt::ivec3 dim = volume->getDimensions();

    typedef unsigned char BYTE;
    BYTE* const data = reinterpret_cast<BYTE* const>(volume->getData());
    BYTE* aux = new BYTE[bytesPerVoxel];

    for (int y = 0; y < (dim.y / 2); ++y) {
        int _y = (dim.y - (y + 1));
        int b = (y * dim.x);
        int _b = (_y * dim.x);
        for (int z = 0; z < dim.z; ++z) {
            int a = (z * (dim.x * dim.y));
            for (int x = 0; x < dim.x; ++x) {
                int l = a + b + x;
                int r = a + _b + x;
                memcpy(aux, (data + (l * bytesPerVoxel)), bytesPerVoxel);
                memcpy((data + (l * bytesPerVoxel)), (data + (r * bytesPerVoxel)), bytesPerVoxel);
                memcpy((data + (r * bytesPerVoxel)), aux, bytesPerVoxel);
            }
        }
    }
    delete [] aux;
}

void VolumeReader::reverseZSliceOrder(VolumeRAM* const volume) const {
    if (volume == 0)
        return;

    size_t bytesPerVoxel = volume->getBytesPerVoxel();
    const tgt::ivec3 dim = volume->getDimensions();
    size_t sliceSize = dim.x * dim.y * bytesPerVoxel;

    typedef unsigned char BYTE;
    BYTE* const data = reinterpret_cast<BYTE* const>(volume->getData());
    BYTE* aux = new BYTE[sliceSize];
    for (int z = 0; z < dim.z / 2; ++z) {
        // copy slice i to auxiliary buffer
        //
        memcpy(aux, (data + (z * sliceSize)), sliceSize);

        // copy slice r = (dim.z - 1 - i) to position i
        //
        int r = (dim.z - (z + 1));
        memcpy((data + (z * sliceSize)), (data + (r * sliceSize)), sliceSize);

        // copy auxiliary buffer back to position r
        //
        memcpy((data + (r * sliceSize)), aux, sliceSize);
    }
    delete [] aux;
}

VolumeURL VolumeReader::convertOriginToRelativePath(const VolumeURL& origin, const std::string& basePath) const {

    return VolumeURL(origin.getProtocol(),
                        tgt::FileSystem::relativePath(origin.getPath(), basePath),
                        origin.getSearchString());

}

VolumeURL VolumeReader::convertOriginToAbsolutePath(const VolumeURL& origin, const std::string& basePath) const {

    // build new path only if this is not an absolute path
    if (!tgt::FileSystem::isAbsolutePath(origin.getPath())) {
        return VolumeURL(origin.getProtocol(),
                            tgt::FileSystem::absolutePath(basePath + "/" + origin.getPath()),
                            origin.getSearchString());
    }
    else
        return origin;
}

const std::vector<std::string>& VolumeReader::getSupportedExtensions() const {
    return extensions_;
}

const std::vector<std::string>& VolumeReader::getSupportedFilenames() const {
    return filenames_;
}

const std::vector<std::string>& VolumeReader::getSupportedProtocols() const {
    return protocols_;
}

bool VolumeReader::canHandleMissingExtension() const {
    return false;
}

bool VolumeReader::canSupportFileWatching() const {
    return false;
}

void VolumeReader::setProgressBar(ProgressBar* progressBar) {
    progress_ = progressBar;
}

ProgressBar* VolumeReader::getProgressBar() const {
    return progress_;
}

} // namespace voreen
