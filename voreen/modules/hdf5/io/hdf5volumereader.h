/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_HDF5VOLUMEREADER_H
#define VRN_HDF5VOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

/**
 * Reader for volumes in the hierarchical data format 5 (HDF5).
 * http://www.hdfgroup.org/
 *
 * HDF5 is a generic data format. DataSets of 3 oder 4 dimensions
 * in the tree of a file are considered volumes. If present the
 * fourth dimension indicates the channel. Additional attributes
 * (e.g. spacing) have to be specified as an Attribute of the
 * DataSet they describe.
 *
 */
class VRN_CORE_API HDF5VolumeReaderBase : public VolumeReader {
public:
    HDF5VolumeReaderBase(bool separatedChannels);
    virtual ~HDF5VolumeReaderBase() {}

    /**
     * See VolumeReader.
     */
    virtual VolumeList* read(const std::string& url);

    /*
     * See VolumeReader
     */
    virtual VolumeBase* read(const VolumeURL& origin);

    /**
     * See VolumeReader.
     *
     * All 3 oder 4 dimensional DataSets within the file will be considered a
     * volume. Location and channel number are distinguish different
     * volumes in the same file.
     */
    virtual std::vector<VolumeURL> listVolumes(const std::string& url) const;

    /**
     * HDF5 format supports filewatching, thus this functions returns true.
     */
    virtual bool canSupportFileWatching() const;

private:
    bool separatedChannels_;
    static const std::string loggerCat_;
};

class VRN_CORE_API HDF5VolumeReader : public HDF5VolumeReaderBase {
public:
    HDF5VolumeReader();
    virtual VolumeReader* create(ProgressBar* progress = nullptr) const;

    virtual std::string getClassName() const   { return "HDF5VolumeReader"; }
    virtual std::string getFormatDescription() const { return "3D HDF5 format"; /*multi-channel volumes are separated into multiple single-channel volumes.*/ }
};

class VRN_CORE_API HDF5VolumeReaderOriginal : public HDF5VolumeReaderBase {
public:
    HDF5VolumeReaderOriginal();
    virtual VolumeReader* create(ProgressBar* progress = nullptr) const;

    virtual std::string getClassName() const   { return "HDF5VolumeReaderOriginal"; }
    virtual std::string getFormatDescription() const { return "3D HDF5 format; multi-channel volumes remain untouched."; }
};

} // namespace voreen

#endif // VRN_HDF5VOLUMEREADER_H
