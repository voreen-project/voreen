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

#ifndef VRN_LZ4VOLUMEREADER_H
#define VRN_LZ4VOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

/**
 * Reader for volumes in the hierarchical data format 5 (LZ4).
 * http://www.hdfgroup.org/
 *
 * LZ4 is a generic data format. DataSets of 3 oder 4 dimensions
 * in the tree of a file are considered volumes. If present the
 * fourth dimension indicates the channel. Additional attributes
 * (e.g. spacing) have to be specified as an Attribute of the
 * DataSet they describe.
 *
 */
class VRN_CORE_API LZ4SliceVolumeFileReader : public VolumeReader {
public:
    LZ4SliceVolumeFileReader();
    ~LZ4SliceVolumeFileReader() {}
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "LZ4SliceVolumeReader"; }
    virtual std::string getFormatDescription() const { return "LZ4 Slice Volume Format"; }

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

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_LZ4VOLUMEREADER_H
