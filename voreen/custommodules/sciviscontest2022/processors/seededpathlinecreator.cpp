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

#include "seededpathlinecreator.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include <tuple>
#include <set>

namespace voreen {

const std::string SeededPathlineCreator::loggerCat_{ "flowanalysis.SeededPathlineCreator" };

SeededPathlineCreator::SeededPathlineCreator() :
    PathlineCreator(),
    seedingCurveInport_(Port::PortDirection::INPORT, "inport.seedingCurve", "Seeding Curve")
{
    addPort(seedingCurveInport_);
}

bool SeededPathlineCreator::isReady() const {
    return PathlineCreator::isReady() && seedingCurveInport_.hasData();
}

PathlineCreatorInput SeededPathlineCreator::prepareComputeInput() {
    auto result = PathlineCreator::prepareComputeInput();
    // Clear seeding points:
    result.seedPoints.clear();
    // Reset seed points:
    if (seedingCurveInport_.getData()->getClassName() == "PointListGeometryVec3") {
        auto geometry = static_cast<const PointListGeometryVec3*>(seedingCurveInport_.getData());
        const auto& transformation = geometry->getTransformationMatrix();
        result.seedPoints.reserve(geometry->getNumPoints());
        for (auto i = 0u; i < geometry->getNumPoints(); ++i) {
            result.seedPoints.emplace_back(transformation * geometry->getPoint(i));
        }
    }
    if (seedingCurveInport_.getData()->getClassName() == "PointSegmentListGeometryVec3") {
        auto geometry = static_cast<const PointSegmentListGeometryVec3*>(seedingCurveInport_.getData());
        const auto& transformation = geometry->getTransformationMatrix();
        const auto& points = geometry->getPoints();
        for (const auto& point : points) {
            result.seedPoints.emplace_back(transformation * point);
        }
    }
    
    LINFO("Seeding with " << result.seedPoints.size() << " seeding points");
    return result;
}

} // namespace voreen
