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

#ifndef VRN_PHILIPSUSVOLUMEREADER_H
#define VRN_PHILIPSUSVOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

/**
 * Reader for Philips Sonos 7500 <tt>.dcm</tt> files. This reader is implemented based
 * on the example reader specified in the document 'Reading Philips 3D DICOM Files from
 * Sonos 7500 D.0'.
 */
class PhilipsUSVolumeReader : public VolumeReader {
public:
    PhilipsUSVolumeReader(ProgressBar* progress = 0);
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "PhilipsUSVolumeReader"; }
    virtual std::string getFormatDescription() const { return "Philips Sonos 7500 format"; }

    virtual VolumeList* read(const std::string& url);

    virtual VolumeList* readSlices(const std::string& url, size_t firstSlice=0, size_t lastSlice=0);

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_PHILIPSUSVOLUMEREADER_H
