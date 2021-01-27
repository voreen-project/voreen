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

#include "geometrydelay.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"

#include <vector>

using tgt::vec3;
using tgt::ivec3;
using tgt::ivec2;
using std::vector;

namespace voreen {

const std::string GeometryDelay::loggerCat_("voreen.GeometryDelay");

GeometryDelay::GeometryDelay()
  : Processor(),
    geometryInport_(Port::INPORT, "geometry.inport"),
    geometryOutport_(Port::OUTPORT, "geometry.outport"),
    delayUpdates_("delayUpdates", "Delay Updates", true),
    update_("update", "Update")
{
    addPort(geometryInport_);
    addPort(geometryOutport_);

    addProperty(delayUpdates_);
    update_.onChange(MemberFunctionCallback<GeometryDelay>(this, &GeometryDelay::forceUpdate));
    addProperty(update_);
}

GeometryDelay::~GeometryDelay() {
}

Processor* GeometryDelay::create() const {
    return new GeometryDelay();
}

void GeometryDelay::invalidate(int inv) {
    if(!delayUpdates_.get())
        Processor::invalidate(inv);
}

void GeometryDelay::process() {
    forceUpdate();
}

void GeometryDelay::forceUpdate() {
    const Geometry* g = geometryInport_.getData();

    // assign result to outport
    geometryOutport_.setData(const_cast<Geometry*>(g), false);
}


} // namespace voreen
