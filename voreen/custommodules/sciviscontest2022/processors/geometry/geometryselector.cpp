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

#include "geometryselector.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"

namespace voreen {

    const std::string GeometrySelector::loggerCat_("voreen.core.GeometrySelector");

    GeometrySelector::GeometrySelector()
        : Processor(),
        inport_(Port::INPORT, "geometrysequence.in", "GeometrySequence Input", false),
        outport_(Port::OUTPORT, "geometrysequence.out", "GeometrySequence Output", false),
        geometryID_("imageID", "Selected Geometry", -1, -1, std::numeric_limits<int>::max() - 1)
    {
        addPort(inport_);
        addPort(outport_);

        addProperty(geometryID_);
    }

    Processor* GeometrySelector::create() const {
        return new GeometrySelector();
    }

    void GeometrySelector::process() {
        const Geometry* geometry = inport_.getData();

        if (geometryID_.get() == -1 || geometry == nullptr) {
            outport_.clear();
            outport_.invalidatePort();
            return;
        }

        // see if cast to geometrysequence is possible
        if (auto sequence = dynamic_cast<const GeometrySequence*>(geometry)) {
            std::unique_ptr<Geometry> selected = sequence->getGeometry(geometryID_.get())->clone();
            outport_.setData(selected.release());
        }
        else if (auto segmentlist = dynamic_cast<const PointSegmentListGeometryVec3*>(geometry)) {
            auto id = geometryID_.get();
            // I managed to crash voreen without this test
            if (id < segmentlist->getNumSegments()) {
                auto selected = std::unique_ptr<PointListGeometryVec3>(new PointListGeometryVec3());
                const auto& selectedSegment = segmentlist->getSegment(geometryID_.get());
                selected->setData(selectedSegment);
                outport_.setData(selected.release());
            }
        }
        else {
            outport_.setData(geometry->clone().release());
        }
    }

    void GeometrySelector::initialize() {
        Processor::initialize();

        adjustPropertiesToInput();
    }

    void GeometrySelector::deinitialize() {
        Processor::deinitialize();
    }

    void GeometrySelector::adjustPropertiesToInput() {

        const Geometry* geometry = inport_.getData();

        if (!geometry)
            return;

        // see if cast to geometrysequence is possible
        if (auto sequence = dynamic_cast<const GeometrySequence*>(geometry)) {
            //if we have a list, adapt min and max values
            int size = static_cast<int>(sequence->getNumGeometries());
            geometryID_.setMinValue(std::min(0, size - 1));
            geometryID_.setMaxValue(size - 1);

            // set to first image if no image was present earlier
            if (size > 0 && geometryID_.get() == -1)
                geometryID_.set(0);
        }
        else if (auto segmentlist = dynamic_cast<const PointSegmentListGeometryVec3*>(geometry)) {
            int size = static_cast<int>(segmentlist->getNumSegments());
            geometryID_.setMinValue(std::min(0, size - 1));
            geometryID_.setMaxValue(size - 1);

            // set to first image if no image was present earlier
            if (size > 0 && geometryID_.get() == -1)
                geometryID_.set(0);
        }
        else {
            // geometry is no sequence, thus its just one element
            geometryID_.setMinValue(0);
            geometryID_.setMaxValue(0);
            geometryID_.set(0);
        }
    }

} // namespace
