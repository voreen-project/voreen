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

#include "utils.h"

#include "voreen/core/io/volumeserializer.h"

#ifdef VRN_MODULE_HDF5
#include "modules/hdf5/io/hdf5volumereader.h"
#endif

namespace voreen {

VolumeReader* EnsembleVolumeReaderPopulator::getVolumeReader(const std::string& path) const {

#ifdef VRN_MODULE_HDF5
    // For HDF5 files we first try to use the multi-channel reader.
    std::string ext = tgt::FileSystem::fileExtension(path);
    if (ext == "h5" || ext == "hdf5") {
        static HDF5VolumeReaderOriginal hdf5Reader;
        return &hdf5Reader;
    }
#endif

    try {
        return volumeSerializerPopulator_.getVolumeSerializer()->getReaders(path).front();
    }
    catch (tgt::UnsupportedFormatException&) {
    }

    return nullptr;
}

}


