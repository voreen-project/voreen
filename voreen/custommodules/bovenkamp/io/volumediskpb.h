/***********************************************************************************
*                                                                                 *
* Voreen - The Volume Rendering Engine                                            *
*                                                                                 *
* Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_VOLUMEDISK_PB_H
#define VRN_VOLUMEDISK_PB_H

#include "voreen/core/datastructures/volume/volumedisk.h"

namespace voreen {

/**
* Disk volume representing a single channel (3D volume) of a PB file.
*/
class VRN_CORE_API VolumeDiskPB : public VolumeDisk {
public:
    /**
     * Constructor for magnitude volumes.
     * @param magnitudeFilename Single filename pointing to a magnitude file.
     * @param invertPosition Vector holding invert flag for each position component
     * @param dimensions Dimensions of the volume
     * @param timeStep Actual time step represented by this volume
     */
    VolumeDiskPB(const std::string& magnitudeFilename,
                 const tgt::bvec3& invertPosition,
                 const tgt::svec3& dimensions,
                 int timeStep);

    /**
    * Constructor for velocity volumes.
    * @param velocityXFilename Filename pointing to the x component velocity file.
    * @param velocityYFilename Filename pointing to the y component velocity file.
    * @param velocityZFilename Filename pointing to the z component velocity file.
    * @param invertPosition Vector holding invert flag for each position component
    * @param invertVelocity Vector holding invert flag for each velocity component
    * @param timeStep Actual time step represented by this volume
    */
    VolumeDiskPB(const std::string& velocityXFilename,
                 const std::string& velocityYFilename,
                 const std::string& velocityZFilename,
                 const tgt::bvec3& invertPosition,
                 const tgt::bvec3& invertVelocity,
                 const tgt::svec3& dimensions,
                 int timeStep);

    /**
     * Destructor.
     */
    virtual ~VolumeDiskPB();

    /**
     * Computes a hash string from the following properties: file names and time step
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
     * Loads a set of consecutive z slices of the PB files from disk
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
     * Loads a brick of the PB files from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param offset lower-left-front corner voxel of the brick to load
     * @param dimensions dimension of the brick to load
     *
     * @throw tgt::Exception if the brick could not be loaded
     */
    virtual VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const;

private:

    /// Helper function to retrieve volume data from file.
    void readFile(const std::string& filename, VolumeRAM* volume, size_t channel, const tgt::svec3& brickOffset, const tgt::svec3& brickDimensions) const;

    /// List of filenames together defining the volume.
    std::vector<std::string> filenames_;

    /// Determines if axis should be inverted when loaded.
    const tgt::bvec3 invertPosition_;
    const tgt::bvec3 invertVelocity_;

    /// Time step of the volume
    const int timeStep_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif //VRN_VOLUMEDISK_PB_H
