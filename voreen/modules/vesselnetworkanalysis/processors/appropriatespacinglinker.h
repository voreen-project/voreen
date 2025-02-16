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

#pragma once

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "tgt/vector.h"

namespace voreen {

class AppropriateSpacingLinker : public VolumeProcessor {
public:
    AppropriateSpacingLinker();
    virtual ~AppropriateSpacingLinker();

    virtual std::string getClassName() const         { return "AppropriateSpacingLinker";      }
    virtual std::string getCategory() const       { return "Util"; }
    virtual std::string setDescriptions() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() { setDescription( "Find an appropriate spacing to be used, for example, in a resampling processor This processor simply reads the minimum of the spacing in each dimension from the specified volume, multiplies it by a specified factor and displays the result in a property."); }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual void process();

private:
    VolumePort inport_;

    FloatVec3Property spacingDisplay_;
    FloatProperty multiplier_;

    static const std::string loggerCat_;
};

} // namespace voreen
