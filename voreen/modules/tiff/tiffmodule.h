/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_TIFFMODULE_H
#define VRN_TIFFMODULE_H

#include "voreen/core/voreenmodule.h"

#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class TiffModule : public VoreenModule {

public:
    TiffModule(const std::string& modulePath);

    virtual std::string getDescription() const {
        return "Provides a volume reader for multi-image TIFF files, using the libtiff library.";
    }

private:
    /// Determine directory count of a (OME-)Tiff sequence only for the first and last file,
    /// and estimate it for the other files. This speeds up the initial opening of a (OME-)Tiff sequence.
    BoolProperty estimateDirectoryCount_;

};

} // namespace

#endif // VRN_TIFFMODULE_H
