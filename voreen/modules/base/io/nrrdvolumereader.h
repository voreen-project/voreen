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

#ifndef VRN_NRRDVOLUMEREADER_H
#define VRN_NRRDVOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

/**
 * Reader for <tt>.nrrd</tt> volume files (nearly raw raster data).
 * TODO: This reader is still incomplete and largely untested.
 *
 * See http://teem.sourceforge.net/nrrd/ for details about the file format.
 */
class VRN_CORE_API NrrdVolumeReader : public VolumeReader {
public:
    NrrdVolumeReader(ProgressBar* progress = 0) :
        VolumeReader(progress)
    {
        extensions_.push_back("nrrd");
        extensions_.push_back("nhdr");
    }

    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "NrrdVolumeReader"; }
    virtual std::string getFormatDescription() const { return "Nearly Raw Raster Data"; }

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

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_NRRDVOLUMEREADER_H
