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

#ifndef VRN_INTERPOLATIONFUNCTIONBASE_H
#define VRN_INTERPOLATIONFUNCTIONBASE_H

#include <string>
#include "tgt/tgt_gl.h"
#include "tgt/vector.h"
#include "voreen/core/voreencoreapi.h"
#include "voreen/core/voreenobject.h"

namespace voreen {

/// This class is an non-templated interface for interpolation functions.
class VRN_CORE_API InterpolationFunctionBase : public VoreenSerializableObject {
public:
    virtual ~InterpolationFunctionBase() {}

    /// Returns the category of the interpolation.
    virtual std::string getCategory() const = 0;

    /// Returns the name of the interpolation to be displayed in the GUI.
    virtual std::string getGuiName() const = 0;

    virtual std::string getClassName() const = 0;

    virtual InterpolationFunctionBase* create() const = 0;
    //virtual VoreenSerializableObject* create() const = 0;

protected:
};

} // namespace voreen

#endif
