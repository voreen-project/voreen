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

// Don't forget to change the header guard!!!
#ifndef VRN_SAMPLEMODULE_H
#define VRN_SAMPLEMODULE_H

//include module base class
#include "voreen/core/voreenmodule.h"

//use namespace voreen
namespace voreen {

/**
 * Each Module has to inherit from VoreenModule.
 */
class SampleModule : public VoreenModule {

public:

    /**
     * Constructor of the module.
     * @param modulePath the path to the module
     */
    SampleModule(const std::string& modulePath);

    /**
     * Sets the description to be shown in the VoreenVE GUI.
     */
    virtual std::string getDescription() const;
};

} // namespace

#endif // VRN_SAMPLEMODULE_H
