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

#include "geometryclipping.h"

namespace voreen {

const std::string GeometryClipping::loggerCat_("voreen.base.GeometryClipping");

GeometryClipping::GeometryClipping()
    : Processor()
    , inport_(Port::INPORT, "geometry.geometry", "Geometry Input")
    , outport_(Port::OUTPORT, "geometry.clippedgeometry", "Clipped Geometry Output")
    , enabled_("enabled", "Enable", true)
    , invert_("invert_", "Invert", false)
    , normal_("planeNormal", "Plane Normal", tgt::vec3(0, 1, 0), tgt::vec3(-1), tgt::vec3(1))
    , position_("planePosition", "Plane Position", 0.0f, -/*1e4*/FLT_MAX, /*1e4*/FLT_MAX, Processor::INVALID_RESULT, NumericProperty<float>::STATIC)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enabled_);
    addProperty(invert_);
    addProperty(normal_);
    addProperty(position_);
}

GeometryClipping::~GeometryClipping()
{}

Processor* GeometryClipping::create() const {
    return new GeometryClipping();
}

void GeometryClipping::process() {
    const Geometry* inputGeometry = inport_.getData();
    tgtAssert(inport_.getData(), "no geometry");

    if (!enabled_.get()) {
        outport_.setData(inputGeometry, false);
        return;
    }

    tgt::vec3 n = normalize(normal_.get());
    float d = position_.get();

    if (invert_.get()) {
        n *= -1.f;
        d *= -1.f;
    }

    tgt::plane plane = tgt::plane(n, d);

    std::unique_ptr<Geometry> outputGeometry = inputGeometry->clone();
    double epsilon = static_cast<double>(tgt::length(outputGeometry->getBoundingBox(false).diagonal())) * 1e-6;
    epsilon = std::max(0.0, epsilon);
    outputGeometry->clip(plane, epsilon);

    outport_.setData(outputGeometry.release());
}

}  //namespace
