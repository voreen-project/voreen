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

#ifndef VRN_VOLUMEDISK_H
#define VRN_VOLUMEDISK_H

#include <typeinfo>
#include <vector>

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumegl.h"

#include <stdexcept>
namespace voreen {

/**
 * Base class for representations storing the information to do a lazy loading of the volume data from disk.
 */
class VRN_CORE_API VolumeDisk : public VolumeRepresentation {
public:
    /**
     * @param format @see VolumeFactory
     * @param dimensions voxel dimensions of the volume
     */
    VolumeDisk(const std::string& format, tgt::svec3 dimensions);

    /// @see VolumeFactory
    virtual std::string getFormat() const { return format_; }

    /// Returns the base type (e.g., "float" for a representation of format "Vector3(float)").
    virtual std::string getBaseType() const;

    virtual size_t getNumChannels() const;
    virtual size_t getBytesPerVoxel() const;

    /**
     * Returns a hash string of the disk volume.
     * For performance reasons, the hash is usually not computed from the volume itself,
     * but rather from identifying information such as the disk filename(s).
     */
    virtual std::string getHash() const = 0;

    /**
     * Loads the entire volume from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @throw tgt::Exception if the volume could not be loaded
     */
    virtual VolumeRAM* loadVolume() const = 0;

    /**
     * Loads a set of consecutive z slices from disk and returns them as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param firstZSlice first slice of the slice range to load (inclusive)
     * @param lastZSlice last slice of the slice range to load (inclusive)
     *
     * @throw tgt::Exception if the slices could not be loaded
     */
    virtual VolumeRAM* loadSlices(const size_t firstZSlice, const size_t lastZSlice) const = 0;

    /**
     * Loads a brick of the volume from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param offset lower-left-front corner voxel of the brick to load
     * @param dimensions dimension of the brick to load
     *
     * @throw tgt::Exception if the brick could not be loaded
     */
    virtual VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const = 0;

protected:
    const std::string format_;

    /// Cache often used volume information
    const std::string baseType_;
    const size_t numChannels_;
    const size_t bytesPerVoxel_;

    static const std::string loggerCat_;
};

//-------------------------------------------------------------------------------------------------

/**
 * Disk volume that references a raw (binary) file.
 */
class VRN_CORE_API VolumeDiskRaw : public VolumeDisk {
public:
    /**
     * @param filename Absolute file name of the raw file.
     * @param format @see VolumeFactory
     * @param dimensions voxel dimensions of the volume
     * @param offset If offset >= 0: Bytes offset from beginning of file. If offset < 0: assume data is aligned to end of file.
     * @param swapEndian Swaps endianness if true.
     */
    VolumeDiskRaw(const std::string& filename, const std::string& format, tgt::svec3 dimensions, int64_t offset = 0, bool swapEndian = false);
    VolumeDiskRaw(const VolumeDiskRaw* diskrep);
    virtual ~VolumeDiskRaw();

    std::string getFileName() const { return filename_; }

    /// offset in the file (in bytes).
    int64_t getOffset() const { return offset_; }
    bool getSwapEndian() const { return swapEndian_; }

    /// Computes a hash string from the filename, the format, the dimensions, the offset and the swapEndian parameter.
    virtual std::string getHash() const;

    /**
     * Loads the entire volume from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @throw tgt::Exception if the volume could not be loaded
     */
    virtual VolumeRAM* loadVolume() const;

    /**
     * Loads a set of consecutive z slices from disk and returns them as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param firstZSlice first slice of the slice range to load (inclusive)
     * @param lastZSlice last slice of the slice range to load (inclusive)
     *
     * @throw tgt::Exception if the slices could not be loaded
     */
    virtual VolumeRAM* loadSlices(const size_t firstZSlice, const size_t lastZSlice) const;

    /**
     * Loads a brick of the volume from disk and returns it as VolumeRAM.
     * The caller is responsible for deleting the returned object.
     *
     * @param offset lower-left-front corner voxel of the brick to load
     * @param dimensions dimension of the brick to load
     *
     * @throw tgt::Exception if the brick could not be loaded
     */
    virtual VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const;

protected:
    std::string filename_;
    int64_t offset_;
    bool swapEndian_;

    static const std::string loggerCat_;
};

//-------------------------------------------------------------------------------------------------

/// Creates a VolumeRam from a VolumeDisk.
class VRN_CORE_API  RepresentationConverterLoadFromDisk : public RepresentationConverter<VolumeRAM> {
public:
    virtual bool canConvert(const VolumeRepresentation* source) const;
    virtual VolumeRepresentation* convert(const VolumeRepresentation* source) const;
};

//-------------------------------------------------------------------------------------------------

/// Creates a VolumeGL from a VolumeDisk.
class VRN_CORE_API  RepresentationConverterLoadFromDiskToGL : public RepresentationConverter<VolumeGL> {
public:
    virtual bool canConvert(const VolumeRepresentation* source) const;
    virtual VolumeRepresentation* convert(const VolumeRepresentation* source) const;
};

} // namespace voreen

#endif
