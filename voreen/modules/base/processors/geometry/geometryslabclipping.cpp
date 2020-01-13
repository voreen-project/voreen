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

#include "geometryslabclipping.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

namespace voreen {

GeometrySlabClipping::GeometrySlabClipping()
    : Processor()
    , inport_(Port::INPORT, "geometry.geometry", "Geometry Input")
    , outport_(Port::OUTPORT, "geometry.clippedgeometry", "Clipped Geometry Output")
    , enabled_("enabled", "Enabled", true)
    , normal_("slabNormal", "Slab Normal", tgt::vec3(0, 1, 0), tgt::vec3(-1), tgt::vec3(1))
    , position_("slabPosition", "Slab Position", 0.0f, -100.0f, 100.0f, Processor::INVALID_RESULT, NumericProperty<float>::DYNAMIC)
    , thickness_("slabThickness", "Slab Thickness", 0.5f, 0.0f, 200.0f, Processor::INVALID_RESULT, NumericProperty<float>::DYNAMIC)
    , wheelInteractionHandler_("wheelInteractionHandler", "Slab Position", &position_,
        tgt::Event::MODIFIER_NONE, false, false)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enabled_);
    thickness_.onChange(MemberFunctionCallback<GeometrySlabClipping>(this, &GeometrySlabClipping::thicknessChanged));
    addProperty(normal_);
    addProperty(position_);
    addProperty(thickness_);

    addInteractionHandler(wheelInteractionHandler_);
}

GeometrySlabClipping::~GeometrySlabClipping() {
}

Processor* GeometrySlabClipping::create() const {
    return new GeometrySlabClipping();
}

void GeometrySlabClipping::process() {
    const Geometry* inputGeometry = inport_.getData();
    tgtAssert(inport_.getData(), "no geometry");

    if (!enabled_.get()) {
        outport_.setData(inputGeometry, false);
        return;
    }

    std::unique_ptr<Geometry> outputGeometry = inputGeometry->clone();
    double epsilon = static_cast<double>(tgt::length(outputGeometry->getBoundingBox(false).diagonal())) * 1e-6;
    outputGeometry->clip(tgt::plane(tgt::normalize(normal_.get()),  position_.get() + 0.5f * thickness_.get()), epsilon);
    outputGeometry->clip(tgt::plane(tgt::normalize(-normal_.get()), -position_.get() + 0.5f * thickness_.get()), epsilon);

    outport_.setData(outputGeometry.release());
}

void GeometrySlabClipping::thicknessChanged() {
    position_.setStepping(std::min(thickness_.get(), 0.1f));
}

}  //namespace
