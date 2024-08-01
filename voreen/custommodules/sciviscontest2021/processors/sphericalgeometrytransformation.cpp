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

#include "sphericalgeometrytransformation.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "tgt/logmanager.h"


namespace {

void convSCtoCC(tgt::vec3& pos, const tgt::vec3& dimensions, float minRadius, float maxRadius, float shiftFactor) {

    float phi   = (pos.x + 180.0f) / 180.0f * tgt::PIf; //-(pos.x - 180.0f) / 180.0f * PI; TODO: According to NetCDF info.
    float r     = (minRadius + (pos.y / dimensions.y) * (maxRadius - minRadius)) / shiftFactor;
    float theta = pos.z / 180.0f * tgt::PIf; // TODO: check if we need to apply a phase shift or to negate?

    pos.x = r * std::cos(phi) * std::sin(theta);
    pos.y = r * std::sin(phi) * std::sin(theta);
    pos.z = r * std::cos(theta);
}

}


namespace voreen {

const std::string SphericalGeometryTransformation::loggerCat_("voreen.scivis2021.SphericalGeometryTransformation");

SphericalGeometryTransformation::SphericalGeometryTransformation()
    : Processor()
    , inport_(Port::INPORT, "geometry.input", "Geometry Input")
    , volumeport_(Port::INPORT, "volume.input", "Volume Input")
    , outport_(Port::OUTPORT, "geometry.output", "Geometry Output", false)
    , enableProcessing_("enableProcessing", "Enable")
    , radiusMin_("radiusMin", "Set Radius Min", 3485.0f, 0.0f, 10000.0f)
    , radiusMax_("radiusMax", "Set Radius Max", 6371.0f, 0.0f, 10000.0f)
    , shiftFactor_("shiftFactor", "Set shift factor", 1000.0f, 0.0f, 1000.0f)
{
    addPort(inport_);
    addPort(volumeport_);
    addPort(outport_);

    addProperty(enableProcessing_);
    addProperty(radiusMin_);
    addProperty(radiusMax_);
    addProperty(shiftFactor_);
}

SphericalGeometryTransformation::~SphericalGeometryTransformation() {}

Processor* SphericalGeometryTransformation::create() const {
    return new SphericalGeometryTransformation();
}

void SphericalGeometryTransformation::process() {
    auto inputGeometry = dynamic_cast<const PointSegmentListGeometryVec3*>(inport_.getData());
    /*
    //DEBUG:
    auto inputGeometry = new PointSegmentListGeometryVec3();
    auto meridian = std::vector<tgt::vec3>();
    for(size_t i=0; i<180; i++) {
        meridian.push_back(tgt::vec3(0.0f, volumeport_.getData()->getDimensions().y, i));
    }
    inputGeometry->addSegment(meridian);
    auto equator = std::vector<tgt::vec3>();
    for(size_t i=0; i<180; i++) {
        equator.push_back(tgt::vec3(i, volumeport_.getData()->getDimensions().y, volumeport_.getData()->getDimensions().z/2));
    }
    inputGeometry->addSegment(equator);
    // NOTE: Also set transformation matrix to identity!
     */

    if(!inputGeometry) {
        LERROR("Currently, only PointSegmentListGeometryVec3 supported");
        return;
    }

    if (!enableProcessing_.get()) {
        outport_.setData(inport_.getData(), false);
        return;
    }

    tgt::vec3 dimensions = volumeport_.getData()->getDimensions();

    // We transform the points back into voxel space, since the voxel position is encoding the spherical coordinates.
    tgt::mat4 transformation = volumeport_.getData()->getWorldToVoxelMatrix() * inputGeometry->getTransformationMatrix();

    std::unique_ptr<PointSegmentListGeometryVec3> outputGeometry(new PointSegmentListGeometryVec3);

    for(int i=0; i<inputGeometry->getNumSegments(); i++) {

        const auto& segment = inputGeometry->getSegment(i);

        std::vector<tgt::vec3> transformedSegment;
        transformedSegment.reserve(segment.size());

        for(const auto& p : segment) {
            auto transformedPoint = transformation * p;
            convSCtoCC(transformedPoint, dimensions, radiusMin_.get(), radiusMax_.get(), shiftFactor_.get());
            transformedSegment.emplace_back(transformedPoint);
        }

        outputGeometry->addSegment(transformedSegment);
    }

    outport_.setData(outputGeometry.release());
}

}   // namespace
