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

#ifndef VRN_SPHERICALGEOMETRYTRANSFORMATION_H
#define VRN_SPHERICALGEOMETRYTRANSFORMATION_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class VRN_CORE_API SphericalGeometryTransformation : public Processor {
public:
    SphericalGeometryTransformation();
    virtual ~SphericalGeometryTransformation();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Geometry";                        }
    virtual std::string getClassName() const  { return "SphericalGeometryTransformation"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;           }

protected:
    virtual void setDescriptions() {
        setDescription("Transform the input geometry into a spherical coordinate system.");
    }

    virtual void process();

private:
    GeometryPort inport_;
    VolumePort volumeport_;
    GeometryPort outport_;

    BoolProperty enableProcessing_;

    FloatProperty radiusMin_;
    FloatProperty radiusMax_;
    FloatProperty shiftFactor_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif // VRN_SPHERICALGEOMETRYTRANSFORMATION_H
