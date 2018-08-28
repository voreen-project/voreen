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

#ifndef VRN_ARBITRARYVOLUMECLIPPING_H
#define VRN_ARBITRARYVOLUMECLIPPING_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

class ArbitraryVolumeClipping : public Processor {
public:
    ArbitraryVolumeClipping();
    virtual ~ArbitraryVolumeClipping();

    virtual std::string getCategory() const { return "Volume"; }
    virtual std::string getClassName() const { return "ArbitraryVolumeClipping"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual Processor* create() const { return new ArbitraryVolumeClipping(); }

    virtual void process();

protected:
    virtual void setDescriptions() {
        setDescription("Converts spherical coordinates to properties of a tangential plane.");
    }

    FloatProperty depth_;
    FloatProperty azimuth_;
    FloatProperty elevation_;

    FloatVec3Property plane_;
    FloatProperty planeDist_;

    VolumePort inport_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
