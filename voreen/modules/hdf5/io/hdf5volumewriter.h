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

#ifndef VRN_HDF5VOLUMEWRITER_H
#define VRN_HDF5VOLUMEWRITER_H

#include "voreen/core/io/volumewriter.h"
#include "tgt/vector.h"

namespace voreen {

/**
 * Writes the volume into a .hdf5 file (hierarchical data format)
 *
 * The volume will be a single DataSet named "volume" at the
 * root of the file.
 * Additional Attributes (e.g. spacing) will be part of the
 * DataSet.
 */
class VRN_CORE_API HDF5VolumeWriter : public VolumeWriter {
public:
    HDF5VolumeWriter(ProgressBar* progress = 0);
    virtual VolumeWriter* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "HDF5VolumeWriter"; }
    virtual std::string getFormatDescription() const { return "Hierarchical Data Format 5"; }

    /**
     * Write the given volume to the file specified by filename using dataSet at dataSetLocation.
     * @throws tgt::IOException if the file could not be opened or created, the dataset could not be created or there is no volumerepresentation to read from.
     */
    void write(const std::string& filename, const VolumeBase* volumeHandle, const std::string& dataSetLocation, bool truncateFile, int deflateLevel, tgt::svec3 chunkDim, bool shuffleEnabled);
    /**
     * See VolumeWriter.
     */
    virtual void write(const std::string& filename, const VolumeBase* volumeHandle);

    static const std::string VOLUME_DATASET_NAME;

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif //VRN_HDF5VOLUMEWRITER_H
