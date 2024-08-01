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

#ifndef VRN_VOLUMEDISKDICOM_H
#define VRN_VOLUMEDISKDICOM_H

#include "voreen/core/datastructures/volume/volumedisk.h"

#include "../dicominfo.h"

namespace voreen {

/**
 * Disk volume representing a single DICOM data set.
 */
class VRN_CORE_API VolumeDiskDicom : public VolumeDisk {
public:

    /**
     * @param format @see VolumeFactory
     * @param dimensions voxel dimensions of the volume
     * @param info meta-information describing the data set
     * @param the (ordered) list of files the data of the volume is stored in
     */
    VolumeDiskDicom(const std::string& format, tgt::svec3 dimensions, const DicomInfo& info, const std::vector<std::string>& files);

    virtual ~VolumeDiskDicom();

    /// Computes a hash string from the data properties (including filenames).
    virtual std::string getHash() const;

    /**
     * Loads the data set from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @throw tgt::Exception if the volume could not be loaded
     */
    virtual VolumeRAM* loadVolume() const;

    /**
     * Loads a set of consecutive slices of the DICOM data set from disk and returns them as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param firstSlice first slice of the slice range to load (inclusive)
     * @param lastSlice last slice of the slice range (inclusive)
     *
     * @throw tgt::Exception if the slices could not be loaded
     */
    virtual VolumeRAM* loadSlices(const size_t firstSlice, const size_t lastSlice) const;

    /**
     * Loads a brick of the DICOM volume from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param offset lower-left-front corner voxel of the brick to load
     * @param dimensions dimension of the brick to load
     *
     * @throw tgt::Exception if the brick could not be loaded
     */
    virtual VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const;

protected:

    DicomInfo info_;
    std::vector<std::string> sliceFiles_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
