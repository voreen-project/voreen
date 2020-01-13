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

#ifndef VRN_ANALYZEVOLUMEREADER_H
#define VRN_ANALYZEVOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

class VRN_CORE_API AnalyzeVolumeReader : public VolumeReader {
public:
    AnalyzeVolumeReader(ProgressBar* progress = 0);
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "AnalyzeVolumeReader"; }
    virtual std::string getFormatDescription() const { return "Analyze format"; }

    static std::string getRelatedImgFileName(const std::string& fileName);

    virtual VolumeBase* read(const VolumeURL& origin);

    virtual VolumeList* read(const std::string& url);

    std::vector<VolumeURL> listVolumes(const std::string& url) const;

private:
    VolumeList* read(const std::string& url, int volId, int channel);

    ///@param volId -1: read all, else: read volume volId
    virtual VolumeList* readAnalyze(const std::string& fileName, int volId, int channel);

    ///@param volId -1: read all, else: read volume volId
    virtual VolumeList* readNifti(const std::string& fileName, bool standalone, int volId, int channel);

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_ANALYZEVOLUMEREADER_H
