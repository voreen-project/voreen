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

#ifndef VRN_VOLUMEDISK_OMETIFF_H
#define VRN_VOLUMEDISK_OMETIFF_H

#include "voreen/core/datastructures/volume/volumedisk.h"

namespace voreen {

/// Meta-information about one OmeTiff file on the disk.
struct OMETiffFile {
    OMETiffFile(const std::string& filename, size_t numDirectories, size_t firstZ, size_t firstT, size_t firstC);
    OMETiffFile();
    std::string toString() const;

    std::string filename_;   ///< absolute path to tiff file
    size_t numDirectories_;  ///< number of directories, i.e. slices, contained by the file
    size_t firstZ_;          ///< z coordinate of the first slice in the file
    size_t firstT_;          ///< t coordinate of the first slice in the file
    size_t firstC_;          ///< c coordinate of the first slice in the file
};

//-------------------------------------------------------------------------------------------------

/// Meta-information about an entire OMETiff disk stack.
struct OMETiffStack {
    OMETiffStack(const std::vector<OMETiffFile>& files, const std::string& dimensionOrder, const std::string& datatype,
        tgt::svec3 volumeDim, tgt::vec3 voxelSpacing, size_t sizeC, size_t sizeT, size_t numSlices);
    OMETiffStack();

    std::vector<OMETiffFile> files_;    ///< list of OMETiff files constituting the stack

    std::string dimensionOrder_;  ///< dimension order of the stack, e.g., XYZTC
    std::string datatype_;        ///< data type of the pixel data, @see VolumeFactory
    tgt::vec3 voxelSpacing_;      ///< stacks voxel spacing in mm, i.e., physical size * 1000
    tgt::svec3 volumeDim_;        ///< spatial dimensions of one 3D volume of the stack
    size_t sizeC_;                ///< number of channels the stack contains
    size_t sizeT_;                ///< number of timesteps the stack contains
    size_t numSlices_;            ///< total number of 2D slices contained by the stack
};

//-------------------------------------------------------------------------------------------------

/**
 * Disk volume representing a single channel/timestep (3D volume) of an OmeTiff datastack.
 */
class VRN_CORE_API VolumeDiskOmeTiff : public VolumeDisk {
public:
    /**
     * @param format @see VolumeFactory
     * @param dimensions voxel dimensions of the volume
     * @param datastack meta-information describing the OMETiff datastack. Is passed to the OmeTiffVolumeReader for read operations.
     * @param channel channel the disk volume represents. Is passed to the OmeTiffVolumeReader for read operations.
     * @param timestep timestep the disk volume represents. Is passed to the OmeTiffVolumeReader for read operations.
     */
    VolumeDiskOmeTiff(const std::string& format, tgt::svec3 dimensions,
        const OMETiffStack& datastack, size_t channel, size_t timestep);
    virtual ~VolumeDiskOmeTiff();

    /**
     * Returns the datastack information.
     */
    const OMETiffStack& getDatastack() const;

    /// Computes a hash string from the datastack properties (including filenames) and the selected channel/timestep.
    virtual std::string getHash() const;

    /**
     * Loads the channel/timestep from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @throw tgt::Exception if the volume could not be loaded
     */
    virtual VolumeRAM* loadVolume() const;

    /**
     * Loads a set of consecutive z slices of the OmeTiff channel/timestep from disk
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
     * Loads a brick of the OmeTiff channel/timestep volume from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param offset lower-left-front corner voxel of the brick to load
     * @param dimensions dimension of the brick to load
     *
     * @throw tgt::Exception if the brick could not be loaded
     */
    virtual VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const;

protected:
    OMETiffStack datastack_;

    size_t channel_;
    size_t timestep_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
