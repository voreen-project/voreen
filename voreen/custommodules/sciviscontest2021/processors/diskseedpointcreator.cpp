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

#include "diskseedpointcreator.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "tgt/glmath.h"

#include "modules/flowanalysis/utils/flowutils.h"

#include <random>

namespace voreen {

const std::string DiskSeedPointCreator::loggerCat_("voreen.DiskSeedPointCreator");

DiskSeedPointCreator::DiskSeedPointCreator()
    : Processor()
    , normal_("planeNormal", "Plane Normal", tgt::vec3(0.f, 1.f, 0.f), tgt::vec3(-1.f), tgt::vec3(1.f))
    , distance_("planeDistance", "Plane Distance", 0.0f, -1000.f, 1000.f)
    , radiusMin_("radiusMin", "Set Radius Min", 3485.0f, 0.0f, 10000.0f)
    , radiusMax_("radiusMax", "Set Radius Max", 6371.0f, 0.0f, 10000.0f)
    , shiftFactor_("shiftFactor", "Set shift factor", 100.0f, 0.0f, 1000.0f)
    , numSeedPoints_("numSeedPoints", "Number of Seed Points", 5000, 1, 200000)
    , seedTime_("seedTime", "Current Random Seed", static_cast<int>(time(0)), std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , geomOutport_(Port::OUTPORT, "geometry", "Geometry Output")
{
    addPort(geomOutport_);

    addProperty(normal_);
    normal_.setGroupID("plane");
    addProperty(distance_);
    distance_.setGroupID("plane");
    distance_.setNumDecimals(4);
    setPropertyGroupGuiName("plane","Plane Settings");

    addProperty(radiusMin_);
    addProperty(radiusMax_);
    addProperty(shiftFactor_);

    addProperty(numSeedPoints_);
    addProperty(seedTime_);
}

DiskSeedPointCreator::~DiskSeedPointCreator() {
}

void DiskSeedPointCreator::process() {

    if(normal_.get() == tgt::vec3::zero) {
        LERROR ("Plane normal of (0,0,0) is not allowed! No output gererated!");
        geomOutport_.setData(nullptr);
        return;
    }

    std::function<float()> rnd(std::bind(std::uniform_real_distribution<float>(0.0f, 1.0f), std::mt19937(seedTime_.get())));

    tgt::plane plane(normal_.get(), distance_.get()); // TODO: make adjustable
    tgt::mat4 m = createTransformationMatrix(tgt::vec3::zero, normal_.get());

    float rMin = radiusMin_.get() / shiftFactor_.get();
    float rMax = radiusMax_.get() / shiftFactor_.get();

    std::vector<tgt::vec3> points;
    for(int i=0; i<numSeedPoints_.get(); i++) {

        float angle = rnd() * 2.0f * 3.1415926535897932384626f;
        float r = rMin + (rMax - rMin) * std::sqrt(rnd());
        float x = r * std::cos(angle);
        float y = r * std::sin(angle);

        tgt::vec3 p(x, y, 0.0f);
        points.push_back(m * p);
    }

    auto disk = new PointSegmentListGeometryVec3();
    disk->addSegment(points);
    geomOutport_.setData(disk);
}


} // namespace voreen
