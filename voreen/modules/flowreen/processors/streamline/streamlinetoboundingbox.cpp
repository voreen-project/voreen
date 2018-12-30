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

#include "streamlinetoboundingbox.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

namespace voreen {

const std::string StreamlineToBoundingBox::loggerCat_("voreen.flowreen.StreamlineToBoundingBox");

StreamlineToBoundingBox::StreamlineToBoundingBox()
    : Processor()
    // ports
    , inport_(Port::INPORT, "input", "Streamline Input", false)
    , outport_(Port::OUTPORT, "outport", "Bounding Box Output")
{
    addPort(inport_);
    addPort(outport_);
}

StreamlineToBoundingBox::~StreamlineToBoundingBox() {
}

void StreamlineToBoundingBox::process() {
    const StreamlineListBase* streamlines = inport_.getData();

    if(streamlines) {
        PointListGeometryVec3* list = new PointListGeometryVec3();
        list->addPoint(streamlines->getOriginalWorldBounds().getLLF());
        list->addPoint(streamlines->getOriginalWorldBounds().getURB());
        outport_.setData(list);
    }
}

}   // namespace
