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

#ifndef VRN_VOLUMESERIALIZERPOPULATOR_H
#define VRN_VOLUMESERIALIZERPOPULATOR_H

#include "voreen/core/voreencoreapi.h"
#include <vector>
#include <string>

namespace voreen {

// forward declarations
class ProgressBar;
class VolumeReader;
class VolumeWriter;
class VolumeSerializer;
class ProgressBar;

/**
 * This class creates a VolumeSerializer and registers all known VolumeReader
 * and VolumeWriter objects. You can get a pointer to the VolumeSerializer.
 * The VolumeSerializer will be destroyed when this class is deleted.
 */
class VRN_CORE_API VolumeSerializerPopulator {
public:
    /**
     * Constructor.
     *
     * @param progressBar Optional progress bar to assign to the volume readers and writers.
     *          Is <emph>not</emph> deleted by the populator's destructor.
     */
    VolumeSerializerPopulator(ProgressBar* progressBar = 0);
    ~VolumeSerializerPopulator();

    /// get a pointer of the generated VolumeSerializer
    VolumeSerializer* getVolumeSerializer() const;

    /**
     * Returns all file extensions that are supported by the
     * available volume readers.
     */
    std::vector<std::string> getSupportedReadExtensions() const;

    /**
     * Returns all filenames that are supported by the
     * available volume readers.
     */
    std::vector<std::string> getSupportedReadFilenames() const;

    /**
     * Returns all file extensions that are supported by the
     * available volume writers.
     */
    std::vector<std::string> getSupportedWriteExtensions() const;

    /**
     * Returns all filenames that are supported by the
     * available volume writers.
     */
    std::vector<std::string> getSupportedWriteFilenames() const;

private:
    std::vector<VolumeReader*> readers_;
    std::vector<VolumeWriter*> writers_;
    VolumeSerializer* const vs_;
    ProgressBar* progressBar_;
};

} // namespace voreen

#endif // VRN_VOLUMESERIALIZERPOPULATOR_H
