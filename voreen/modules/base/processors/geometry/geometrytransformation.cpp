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

#include "geometrytransformation.h"

#include "voreen/core/datastructures/geometry/geometry.h"
#include "tgt/logmanager.h"

namespace voreen {

const std::string GeometryTransformation::loggerCat_("voreen.base.GeometryTransformation");

GeometryTransformation::GeometryTransformation()
    : Processor()
    , inport_(Port::INPORT, "geometry.input", "Geometry Input")
    , outport_(Port::OUTPORT, "geometry.output", "Geometry Output", false)
    , enableProcessing_("enableProcessing", "Enable")
    , transformMatrix_("transformMatrix", "Transformation Matrix", tgt::mat4::identity, tgt::mat4(-1e10), tgt::mat4(1e10))
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
    addProperty(transformMatrix_);
}

GeometryTransformation::~GeometryTransformation() {}

Processor* GeometryTransformation::create() const {
    return new GeometryTransformation();
}

void GeometryTransformation::process() {
    const Geometry* inputGeometry = inport_.getData();
    tgtAssert(inputGeometry, "no input geometry");
    if (!enableProcessing_.get()) {
        outport_.setData(inputGeometry, false);
        return;
    }

    // clone and transform input geometry
    std::unique_ptr<Geometry> outputGeometry = inputGeometry->clone();
    outputGeometry->transform(transformMatrix_.get());
    outport_.setData(outputGeometry.release());
}

}   // namespace
