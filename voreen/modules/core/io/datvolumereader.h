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

#ifndef VRN_DATVOLUMEREADER_H
#define VRN_DATVOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

/**
 * Reader for <tt>.dat</tt> files.
 *
 * The following tags are understood by the reader:
 * - <tt>ObjectFileName</tt>: Name of the corresponding <tt>.raw</tt> file. The name may not
 *    contain whitespaces.
 * - <tt>TaggedFileName</tt>: ignored
 * - <tt>Resolution</tt>: Number of voxels in x-, y- and z-direction
 * - <tt>SliceThickness</tt>: Thickness for x-, y- and z-direction
 * - <tt>Format</tt>: Data format like \c UCHAR, \c USHORT, \c FLOAT. Must be supported by the
 *   RawVolumeReader, see RawVolumeReader::readHints() for details.
 * - <tt>NbrTags</tt>: ignored
 * - <tt>ObjectType</tt>: ignored
 * - <tt>ObjectModel</tt>: \c I (intensity) or \c RGBA or \c RGB
 * - <tt>GridType</tt>: only \c EQUIDISTANT is supported
 * - <tt>BitsStored</tt>: specifies logical data type, e.g. 12 for 12 bit CT data where each
 *   voxel is stored in a 16 bit \c short.
 * - <tt>ZeroPoint</tt>: offset for the voxel data for supporting signed data
 * - <tt>MetaString</tt>: a string containing arbitrary meta-data
 * - <tt>NumFrames</tt>: number of time frames, each time frame will be loaded into a separate volume (default: 1)
 *
 * Example file <tt>nucleon.dat</tt>:
 \verbatim
ObjectFileName: nucleon.raw
TaggedFileName: ---
Resolution:     41 41 41
SliceThickness: 1.0 1.0 1.0
Format:         UCHAR
NbrTags:        0
ObjectType:     TEXTURE_VOLUME_OBJECT
ObjectModel:    I
GridType:       EQUIDISTANT
 \endverbatim
 */
class VRN_CORE_API DatVolumeReader : public VolumeReader {
public:
    DatVolumeReader(ProgressBar* progress = 0);
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "DatVolumeReader"; }
    virtual std::string getFormatDescription() const { return "Voreen dat/raw format"; }

    /**
     * Reads the file name called 'ObjectFileName' from the .dat file and
     * returns it. This is particular useful to determine the correct name
     * of the .raw file related to a .dat file, because it might be different
     * than the .dat file name!
     */
    static std::string getRelatedRawFileName(const std::string& fileName);

    /**
     * Loads a single volume from the passed origin.
     *
     * \see VolumeReader::read
     **/
    virtual VolumeBase* read(const VolumeURL& origin);

    /**
     * Loads one or multiple volumes from the specified URL.
     *
     * \see VolumeReader::read
     **/
    virtual VolumeList* read(const std::string& url);

    /**
     * Loads one or multiple volumes from the specified URL.
     *
     * \param   url         url to load volume from
     * \param   timeframe   time frame to select from volume, if -1 all time frames will be selected
     **/
    virtual VolumeList* read(const std::string& url, int timeframe);

    /**
     * Loads the given slices of one or multiple volumes from the specified URL.
     *
     * \param   url         url to load volume from
     * \param   firstSlice  first slice to load
     * \param   lastSlice   last slice to load, if 0 all slices from volume will be loaded
     * \param   timeframe   time frame to select from volume, if -1 all time frames will be selected
     **/
    using VolumeReader::readSlices;
    virtual VolumeList* readSlices(const std::string& url, size_t firstSlice = 0, size_t lastSlice = 0, int timeframe = -1);

    virtual VolumeList* readBrick(const std::string& url, tgt::ivec3 brickStartPos, int brickSize);

private:
    static const std::string loggerCat_;

    VolumeList* readVolumeFile(const std::string& fileName, const tgt::ivec3& dims,size_t firstSlice, size_t lastSlice);

    VolumeList* readVolumeFileBrick(const std::string& fileName,const tgt::ivec3& dims, tgt::ivec3 brickStartPos,
        int brickSize);

    VolumeList* readMetaFile(const std::string& fileName,size_t firstSlice, size_t lastSlice, int timeframe = -1);

    VolumeList* readMetaFileBrick(const std::string& fileName, tgt::ivec3 brickStartPos, int brickSize);
};

} // namespace voreen

#endif // VRN_DATVOLUMEREADER_H
