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

#include "streamlinetogeometry.h"

#include "voreen/core/datastructures/geometry/glmeshgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

namespace voreen {

const std::string StreamlineToGeometry::loggerCat_("voreen.flowreen.StreamlineToGeometry");

StreamlineToGeometry::StreamlineToGeometry()
    : Processor()
    // ports
    , inport_(Port::INPORT, "input", "Streamline Input", false)
    , outport_(Port::OUTPORT, "outport", "Geometry Output")
    , geometryType_("geometryType", "Geometrcy Type")
{
    addPort(inport_);
    addPort(outport_);
    addProperty(geometryType_);
    geometryType_.addOption("glmesh", "GLMeshGeometry", GEOMETRY_GLMESH);
    geometryType_.addOption("pointsegmentlist", "PointSegmentListGeometry", GEOMETRY_POINTSEGMENTLIST);
}

StreamlineToGeometry::~StreamlineToGeometry() {
}

void StreamlineToGeometry::process() {
    const StreamlineListBase* streamlines = inport_.getData();

    if(streamlines) {
        if(geometryType_.getValue() == GEOMETRY_GLMESH) {
            GlMeshGeometryUInt32Color* geometry = new GlMeshGeometryUInt32Color();
            geometry->setPrimitiveType(GL_LINE_STRIP);
            for (const Streamline& streamline : streamlines->getStreamlines()) {
                for (size_t i = 0; i < streamline.getNumElements(); i++) {
                    VertexColor vertex;
                    vertex.pos_ = streamline.getElementAt(i).position_;
                    //TODO: convert to proper color here
                    vertex.color_ = tgt::vec4(streamline.getElementAt(i).velocity_, 0.0f);
                    geometry->addVertex(vertex);
                    geometry->addIndex(geometry->getNumIndices());
                }
                geometry->addIndex(geometry->getPrimitiveRestartIndex());
            }
            outport_.setData(geometry);
        }
        else if(geometryType_.getValue() == GEOMETRY_POINTSEGMENTLIST) {
            PointSegmentListGeometryVec3* geometry = new PointSegmentListGeometryVec3();
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
