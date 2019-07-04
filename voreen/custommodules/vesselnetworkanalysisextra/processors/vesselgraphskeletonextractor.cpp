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

#include "vesselgraphskeletonextractor.h"

#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include <vector>
#include <memory>

namespace voreen {

const std::string VesselGraphSkeletonExtractor::loggerCat_("voreen.vesselnetworkanalysisextra.vesselgraphskeletonextractor");

VesselGraphSkeletonExtractor::VesselGraphSkeletonExtractor()
    : graphInport_(Port::INPORT, "vesselgraphskeletonextractor_graph.inport", "Graph", false, Processor::INVALID_RESULT)
    , nodeOutport_(Port::OUTPORT, "vesselgraphskeletonextractor_nodes.outport", "Nodes", false, Processor::VALID)
    , edgeOutport_(Port::OUTPORT, "vesselgraphskeletonextractor_edges.outport", "Edges", false, Processor::VALID)
{
    addPort(graphInport_);
    addPort(nodeOutport_);
    addPort(edgeOutport_);
}

VesselGraphSkeletonExtractor::~VesselGraphSkeletonExtractor() {
}
VoreenSerializableObject* VesselGraphSkeletonExtractor::create() const {
    return new VesselGraphSkeletonExtractor();
}

bool VesselGraphSkeletonExtractor::isReady() const {
    return isInitialized() && graphInport_.isReady() && (nodeOutport_.isReady() || edgeOutport_.isReady());
}

void VesselGraphSkeletonExtractor::process() {
    const VesselGraph* input = graphInport_.getData();
    if(!input) {
        nodeOutport_.setData(nullptr);
        edgeOutport_.setData(nullptr);
        return;
    }
    const VesselGraph& graph = *input;
    PointListGeometryVec3* nodeGeom = new PointListGeometryVec3();
    for(const auto& node : graph.getNodes()) {
        for(const auto& voxel: node.voxels_) {
            nodeGeom->addPoint(voxel);
        }
    }
    nodeOutport_.setData(nodeGeom);

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
    edgeOutport_.setData(edgeGeom);
}

} // namespace voreen
