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

#ifndef VRN_TIFFVOLUMEREADER_H
#define VRN_TIFFVOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

class IOProgress;

/**
 * Reads a multi-image TIFF file into a volume dataset.
 */
class VRN_CORE_API TiffVolumeReader : public VolumeReader {
public:
    TiffVolumeReader(ProgressBar* progress = 0);
    ~TiffVolumeReader() {}
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "TiffVolumeReader"; }
    virtual std::string getFormatDescription() const { return "3D TIFF format"; }

    virtual VolumeList* read(const std::string& url);

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_TIFFVOLUMEREADER_H
