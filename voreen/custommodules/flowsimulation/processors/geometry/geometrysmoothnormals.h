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

#ifndef VRN_GEOMETRYSMOOTHNORMALS_H
#define VRN_GEOMETRYSMOOTHNORMALS_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/geometryport.h"

#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class GlMeshGeometryUInt32Normal;

/**
 * Calculates and sets smooth normals for the input geometry.
 */
class VRN_CORE_API GeometrySmoothNormals : public Processor {
public:
    GeometrySmoothNormals();
    virtual ~GeometrySmoothNormals();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "GeometrySmoothNormals";  }
    virtual std::string getCategory() const  { return "Geometry";               }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;  }

protected:
    virtual void setDescriptions() {
        setDescription("Calculates and sets smooth normals for the input geometry.");
    }

    virtual void process();

    GeometryPort inport_;        ///< Inport for a list of mesh geometries to close.
    GeometryPort outport_;       ///< Outport for a list of mesh geometries that were closed.

    BoolProperty enabled_;       ///< Determines whether the closing operation is performed.
    FloatProperty epsilon_;      ///< Epsilon for equality check.

    /// category used in logging
    static const std::string loggerCat_;
};

} //namespace

#endif
