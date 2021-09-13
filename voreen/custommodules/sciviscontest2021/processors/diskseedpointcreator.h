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

#ifndef VRN_DISKSEEDPOINTCREATOR_H
#define VRN_DISKSEEDPOINTCREATOR_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "voreen/core/ports/geometryport.h"

namespace voreen {

class DiskSeedPointCreator : public Processor {
public:
    DiskSeedPointCreator();
    virtual ~DiskSeedPointCreator();

    virtual std::string getCategory() const override{ return "Volume"; }
    virtual std::string getClassName() const override{ return "DiskSeedPointCreator"; }
    virtual Processor::CodeState getCodeState() const override{ return CODE_STATE_EXPERIMENTAL; }
    virtual Processor* create() const override{ return new DiskSeedPointCreator(); }
    /** @see processor */
    virtual void process() override;
protected:

    FloatVec3Property normal_;      ///< Clipping plane normal
    FloatProperty distance_;    ///< Clipping plane position in world space
    FloatProperty radiusMin_;
    FloatProperty radiusMax_;
    FloatProperty shiftFactor_;
    IntProperty numSeedPoints_;                             ///< number of seed points
    IntProperty seedTime_;                                  ///< seed

    GeometryPort geomOutport_;      ///< Output of the generated slice geometry

    static const std::string loggerCat_;    ///< meow
};

} // namespace voreen

#endif
