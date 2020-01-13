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

#include "vesselgraphcenterlineconverter.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include <vector>
#include <memory>

namespace voreen {

const std::string VesselGraphCenterlineConverter::loggerCat_("voreen.vesselnetworkanalysis.vesselgraphcenterlineextractor");

VesselGraphCenterlineConverter::VesselGraphCenterlineConverter()
    : graphInport_(Port::INPORT, "vesselgraphcenterlineextractor_graph.inport", "Graph", false, Processor::INVALID_RESULT)
    , outport_(Port::OUTPORT, "vesselgraphcenterlineextractor.outport", "Edges", false, Processor::VALID)
{
    addPort(graphInport_);
    addPort(outport_);
}

VesselGraphCenterlineConverter::~VesselGraphCenterlineConverter() {
}
VoreenSerializableObject* VesselGraphCenterlineConverter::create() const {
    return new VesselGraphCenterlineConverter();
}

void VesselGraphCenterlineConverter::process() {
    const VesselGraph* input = graphInport_.getData();
    if(!input) {
        outport_.setData(nullptr);
        return;
    }
    const VesselGraph& graph = *input;

    PointSegmentListGeometryVec3* edgeGeom = new PointSegmentListGeometryVec3();
    for(auto& edge : graph.getEdges()) {
        std::vector<tgt::vec3> segment(edge.getVoxels().size()+2);
        segment.front() = edge.getNode1().pos_;
        std::transform(edge.getVoxels().begin(), edge.getVoxels().end(), segment.begin()+1, [] (const VesselSkeletonVoxel& voxel) {
                return voxel.pos_;
                });
        segment.back() = edge.getNode2().pos_;
        edgeGeom->addSegment(segment);
    }
    outport_.setData(edgeGeom);
}

} // namespace voreen
