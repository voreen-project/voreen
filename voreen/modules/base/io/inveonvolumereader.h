/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#ifndef VRN_INVEONVOLUMEREADER_H
#define VRN_INVEONVOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

/**
 * Reader for Siemens Inveon (.hdr/.img) files
 */
class VRN_CORE_API InveonVolumeReader : public VolumeReader {
    struct InveonHeader {
        tgt::ivec3 dims_;
        tgt::vec3 spacing_;
        bool bigEndian_;
        std::string voreenVoxelType_;
        int numFrames_;

        InveonHeader() : dims_(1)
                       , spacing_(0.0f)
                       , bigEndian_(false)
                       , voreenVoxelType_("")
                       , numFrames_(1)
                       {}
    };
public:
    InveonVolumeReader(ProgressBar* progress = 0);
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "InveonVolumeReader"; }
    virtual std::string getFormatDescription() const { return "Siemens Inveon files"; }

    VolumeBase* read(const VolumeURL& origin);

    virtual VolumeList* read(const std::string& url);

    /**
     * Loads one or multiple volumes from the specified URL.
     *
     * \param   url         url to load volume from
     * \param   volumeId    id to select the volume, if -1 all volumes will be selected
     **/
    virtual VolumeList* read(const std::string& url, int volumeId);

    std::vector<VolumeURL> listVolumes(const std::string& url) const;

    static InveonHeader readHeader(const std::string filename);
private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_INVEONVOLUMEREADER_H
