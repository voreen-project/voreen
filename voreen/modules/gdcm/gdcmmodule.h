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

#ifndef VRN_GDCMMODULE_H
#define VRN_GDCMMODULE_H

#include "voreen/core/voreenmodule.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class GdcmModule : public VoreenModule
{
public:
    GdcmModule(const std::string& modulePath);

    std::string getDescription() const {
        return "Module for using DICOM data with Grassroots DICOM (GDCM)";
    }

    std::string getAeTitle() const {
        return aetProperty_.get();
    }

    int getIncomingPortNumber() const {
        return incomingPortProperty_.get();
    }

    std::string getDefaultScpUrl() const {
        return scpUrlProperty_.get();
    }

    std::string getDefaultScpAet() const {
        return scpAetProperty_.get();
    }

    int getDefaultScpPortNumber() const {
        return scpPortProperty_.get();
    }

    bool ignoreSliceSpacing() const {
        return ignoreSliceSpacing_.get();
    }

    bool useGdcmRescaling() const {
        return useGdcmRescaling_.get();
    }

private:
    StringProperty aetProperty_; ///< user (SCU) AE Title
    IntProperty incomingPortProperty_; ///< port for incoming connection

    StringProperty scpUrlProperty_; ///< default DICOM SCP (standard when connecting)
    StringProperty scpAetProperty_; ///< AET of the default SCP
    IntProperty scpPortProperty_; ///< port to connect to the default SCP

    BoolProperty ignoreSliceSpacing_;   ///< if selected the dataset is read despite the fact that the slice spacing differs more than 10% within the dataset
    BoolProperty useGdcmRescaling_;     ///< if enabled the rescaling of datasets with intercept and slope varying between the slices is handled by GDCM (might lead to problems very few cases), otherwise a slower method is used
};

} //end namespace

#endif // VRN_GDCMMODULE_H
