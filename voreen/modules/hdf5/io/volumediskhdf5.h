/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEDISK_HDF5_H
#define VRN_VOLUMEDISK_HDF5_H

#include "voreen/core/datastructures/volume/volumedisk.h"
#include "hdf5filevolume.h"

namespace voreen {

/**
 * Disk volume representing a single channel (3D volume) of a HDF5 file.
 */
class VRN_CORE_API VolumeDiskHDF5 : public VolumeDisk {
public:
    /**
     * Constructor.
     * @param volume The handle for the volume contained in a hdf5 file.
     * @param channel The channel within the file volume this volumedisk should represent
     */
    VolumeDiskHDF5(std::unique_ptr<HDF5FileVolume> volume, size_t channel=0);

    /**
     * Destructor.
     */
    virtual ~VolumeDiskHDF5();

    /**
     * Computes a hash string from the datastack properties: file name, volume location and channel.
     */
    virtual std::string getHash() const;

    /**
     * Loads the channel/timestep from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @throw tgt::Exception if the volume could not be loaded
     */
    virtual VolumeRAM* loadVolume() const;

    /**
     * Loads a set of consecutive z slices of the HDF5 channel/timestep from disk
     * and returns them as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param firstZSlice first slice of the slice range to load (inclusive)
     * @param lastZSlice last slice of the slice range (inclusive)
     *
     * @throw tgt::Exception if the slices could not be loaded
     */
    virtual VolumeRAM* loadSlices(const size_t firstZSlice, const size_t lastZSlice) const;

    /**
     * Loads a brick of the HDF5 channel/timestep volume from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param offset lower-left-front corner voxel of the brick to load
     * @param dimensions dimension of the brick to load
     *
     * @throw tgt::Exception if the brick could not be loaded
     */
    virtual VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const;

protected:

    /// The volume inside a HDF5 file
    const std::unique_ptr<HDF5FileVolume> volume_;
    /// Channel of the volume
    const size_t channel_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif //VRN_VOLUMEDISK_HDF5_H
