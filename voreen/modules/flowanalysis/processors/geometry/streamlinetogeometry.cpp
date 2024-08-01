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

#include "streamlinetogeometry.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

const std::string StreamlineToGeometry::loggerCat_("flowanalysis.StreamlineToGeometry");

StreamlineToGeometry::StreamlineToGeometry()
    : Processor()
    // ports
    , inport_(Port::INPORT, "input", "Streamline Input", false)
    , outport_(Port::OUTPORT, "outport", "Geometry Output")
    , geometryType_("geometryType", "Geometry Type")
    , targetCoordinateSystem_("targetCoordinateSystem", "Target Coordinate System")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(geometryType_);
    geometryType_.addOption("glmesh", "GLMeshGeometry", GEOMETRY_GLMESH);
    geometryType_.addOption("pointsegmentlist", "PointSegmentListGeometry", GEOMETRY_POINTSEGMENTLIST);
    addProperty(targetCoordinateSystem_);
    targetCoordinateSystem_.addOption("world", "World Coordinates", WORLD_COORDINATES);
    targetCoordinateSystem_.addOption("voxel", "Voxel Coordinates", VOXEL_COORDINATES);
}

StreamlineToGeometry::~StreamlineToGeometry() {
}

void StreamlineToGeometry::process() {
    const StreamlineListBase* streamlines = inport_.getData();

    if(streamlines) {
        if(geometryType_.getValue() == GEOMETRY_GLMESH) {
            GlMeshGeometryUInt32Color* geometry = new GlMeshGeometryUInt32Color();
            if(targetCoordinateSystem_.getValue() == VOXEL_COORDINATES) {
                geometry->setTransformationMatrix(streamlines->getOriginalWorldToVoxelMatrix());
            } // else: streamline data is in world space by default.
            geometry->setPrimitiveType(GL_LINE_STRIP);
            geometry->enablePrimitiveRestart();

            for (const Streamline& streamline : streamlines->getStreamlines()) {
                for (size_t i = 0; i < streamline.getNumElements(); i++) {
                    const Streamline::StreamlineElement& element = streamline.getElementAt(i);
                    VertexColor vertex;
                    vertex.pos_ = element.position_;
                    //TODO: Use adjustable transfer function.
                    vertex.color_ = tgt::vec4(tgt::vec3(tgt::length(element.velocity_) / streamlines->getMaxMagnitude()), 1.0f);
                    geometry->addIndex(static_cast<uint32_t>(geometry->getNumVertices()));
                    geometry->addVertex(vertex);
                }
                geometry->addIndex(geometry->getPrimitiveRestartIndex());
            }

            outport_.setData(geometry);
        }
        else if(geometryType_.getValue() == GEOMETRY_POINTSEGMENTLIST) {
            PointSegmentListGeometryVec3* geometry = new PointSegmentListGeometryVec3();
            if(targetCoordinateSystem_.getValue() == VOXEL_COORDINATES) {
                geometry->setTransformationMatrix(streamlines->getOriginalWorldToVoxelMatrix());
            } // else: streamline data is in world space by default.

            for (const Streamline& streamline : streamlines->getStreamlines()) {
                std::vector<tgt::vec3> segment(streamline.getNumElements());
                for (size_t i = 0; i < streamline.getNumElements(); i++) {
                    segment[i] = streamline.getElementAt(i).position_;
                }
                geometry->addSegment(segment);
            }

            outport_.setData(geometry);
        }
        else {
            tgtAssert(false, "Unsupported geometry type");
        }
    }
    else {
        outport_.clear();
    }
}

}   // namespace
