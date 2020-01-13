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

#include "geometryeventblocker.h"

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

const std::string GeometryEventBlocker::loggerCat_("voreen.GeometryEventBlocker");

GeometryEventBlocker::GeometryEventBlocker()
  : Processor(),
    geometryInport_(Port::INPORT, "inport"),
    geometryOutport_(Port::OUTPORT, "outport")
{
    addPort(geometryInport_);
    addPort(geometryOutport_);
}

GeometryEventBlocker::~GeometryEventBlocker() {
}

Processor* GeometryEventBlocker::create() const {
    return new GeometryEventBlocker();
}

void GeometryEventBlocker::process() {
    const Geometry* g = geometryInport_.getData();
    geometryOutport_.setData(const_cast<Geometry*>(g), false);
}

void GeometryEventBlocker::onEvent(tgt::Event* /*e*/) {
    return;
}

} // namespace voreen
