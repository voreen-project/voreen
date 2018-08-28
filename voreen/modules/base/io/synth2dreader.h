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

#ifndef VRN_SYNTH2DREADER_H
#define VRN_SYNTH2DREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

struct VRN_CORE_API Synth2DVolumeHeader {
    char magic[4];
    int version;
    char texName[256];
    bool wrap;
    int volSize;
    int numChannels;
    int bytesPerChannel;
};

/**
 * Reader for .vol files by Johann Kopf
 */
class VRN_CORE_API Synth2DReader : public VolumeReader {
public:

    Synth2DReader(ProgressBar* progress = 0);
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "Synth2DReader"; }
    virtual std::string getFormatDescription() const { return "Synth2D .vol format by Johann Kopf"; }

    virtual VolumeList* read(const std::string& url);

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_SYNTH2DREADER_H
