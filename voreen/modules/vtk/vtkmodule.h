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

#ifndef VRN_VTKMODULE_H
#define VRN_VTKMODULE_H

#include "voreen/core/voreenmodule.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class VTKModule: public VoreenModule {

public:
    VTKModule(const std::string& modulePath);

    virtual std::string getDescription() const {
        return "Enables use of VTK library and adds support for additional volume formats such as vti and NetCDF";
    }

    /**
     * Enables or disables forcing a disk representation of loaded volumes, even
     * if they have no native one. If enabled, a hdf5 disk representation will be used.
     * Default: enabled
     */
    void setForceDiskRepresentation(bool enabled);

    /**
     * Returns whether disk representations are forced.
     */
    bool getForceDiskRepresentation() const;

    /**
     * Returns the global instance of this class.
     *
     * @note Does not create the instance. If the module class has not been
     *       instantiated yet, the null pointer is returned.
     */
    static VTKModule* getInstance();

private:

    static VTKModule* instance_;

    BoolProperty forceDiskRepresentation_;
};

} // namespace

#endif
