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

#ifndef VRN_CLIPREGIONGEOMETRYCREATOR_H
#define VRN_CLIPREGIONGEOMETRYCREATOR_H


#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/boundingboxproperty.h"

namespace voreen {

/**
 * Creates a mesh corresponding to a bounding box (e.g., clip region) in volume coordinates.
 */
class VRN_CORE_API ClipRegionGeometryCreator : public Processor {

public:

    ClipRegionGeometryCreator();

    virtual Processor* create() const;

    virtual std::string getClassName() const { return "ClipRegionGeometryCreator"; }
    virtual std::string getCategory() const { return "Geometry"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_TESTING; }

protected:
    virtual void setDescriptions() {
        setDescription("Creates a mesh corresponding to the clipping region / bounding box of an input volume.");
    }

    virtual void process();

    virtual void adjustPropertiesToInput();

    VolumePort inport_;
    GeometryPort outport_;

    IntBoundingBoxProperty boundingBox_;

};

} // namespace

#endif
