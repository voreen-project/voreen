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

#include "subgraphextractor.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "../algorithm/vesselgraphnormalization.h"

#include <numeric>
#include <queue>
#include <unordered_set>
#include <unordered_map>

namespace voreen {



const std::string SubGraphExtractor::loggerCat_("voreen.vesseltopology.subgraphextractor");


SubGraphExtractor::SubGraphExtractor()
    : Processor()
    , inport_(Port::INPORT, "graph.input", "Graph Input", false, Processor::INVALID_RESULT)
    , startingPoint_(Port::INPORT, "pointlist.seedpoints", "Starting Point", false, Processor::INVALID_RESULT)
    , outport_(Port::OUTPORT, "graph.output", "Normalized Graph Output", false, Processor::VALID)
    , enabled_("enabled", "Enabled", true)
    , maxEdgeDistance_("maxEdgeDistance_", "Maximum Number of Edges to Starting Point")
    , keepBounds_("keepBounds", "Keep Bounds of Input Graph", true)
{

    addPort(inport_);
    addPort(startingPoint_);
    addPort(outport_);

    addProperty(enabled_);
    addProperty(maxEdgeDistance_);
    addProperty(keepBounds_);
}

SubGraphExtractor::~SubGraphExtractor() {
}

static bool tryExtractPoints(const Geometry* geometry, std::vector<tgt::vec3>& output) {
    auto transformation = geometry->getTransformationMatrix();
    auto transform = [&transformation] (const tgt::vec3& v) {
        return (transformation*tgt::vec4(v, 1)).xyz();
    };
    if(const PointSegmentListGeometryVec3* seedList = dynamic_cast<const PointSegmentListGeometryVec3*>(geometry)) {
        for(int i = 0; i < seedList->getNumSegments(); i++) {
            auto segment = seedList->getSegment(i);
            for(const tgt::vec3& v: segment) {
                output.push_back(transform(v));
            }
        }
    } else if(const PointListGeometryVec3* seeds = dynamic_cast<const PointListGeometryVec3*>(geometry)) {
        for(const tgt::vec3& v: *seeds) {
            output.push_back(transform(v));
        }
    } else {
        return false;
    }
    return true;
}

static bool tryExtractCenterPoint(const GeometryPort& g, tgt::vec3& output) {
    if(!g.hasData()) {
        return false;
    }
    std::vector<tgt::vec3> points;
    if(!tryExtractPoints(g.getData(), points)) {
        return false;
    }
    if(points.empty()) {
        return false;
    }
    output = std::accumulate(points.begin(), points.end(), tgt::vec3::zero) / static_cast<float>(points.size());
    return true;
}

namespace {

    struct SubGraphNode {
        const VesselGraphNode& node;
        const uint32_t depth;
        const VGNodeID new_id;

        SubGraphNode(const VesselGraphNode& node, uint32_t depth, VGNodeID new_id)
            : node(node)
            , depth(depth)
            , new_id(new_id)
        {
        }
    };

} // anymous namespace

static std::unique_ptr<VesselGraph> extractSubgraph(const VesselGraph& input, const tgt::vec3& starting_point, uint32_t max_edge_distance, bool keepBounds) {
    std::unique_ptr<VesselGraph> output(keepBounds ? new VesselGraph(input.getBounds()) : new VesselGraph());
    auto nodes = input.getNodes();
    if(nodes.empty()) {
        return output;
    }
    float min_dist_sq = std::numeric_limits<float>::infinity();
    const VesselGraphNode* starting_node = nullptr;
    for(auto& node : nodes) {
        float current_dist_sq = tgt::distanceSq(starting_point, node.pos_);
        if(current_dist_sq < min_dist_sq) {
            min_dist_sq = current_dist_sq;
            starting_node = &node;
        }
    }
    if(!starting_node) {
        // No nodes in graph
        return output;
    }
    std::unordered_set<const VesselGraphEdge*> added_edges;
    std::unordered_map<const VesselGraphNode*, VGNodeID> added_nodes;

    std::queue<SubGraphNode> node_queue;
    VGNodeID starting_node_id = output->insertNode(*starting_node);
    added_nodes.insert({ starting_node, starting_node_id });

    node_queue.push(SubGraphNode(*starting_node, 0, starting_node_id));

    while(!node_queue.empty()) {
        const VesselGraphNode& node = node_queue.front().node;
        uint32_t depth = node_queue.front().depth;
        VGNodeID node_id = node_queue.front().new_id;
        node_queue.pop();

        if(depth >= max_edge_distance) {
            continue;
        }

        for(auto& edge: node.getEdges()) {
            if(added_edges.find(&edge.get()) == added_edges.end()) {
                const VesselGraphNode& n1 = edge.get().getNode1();
                const VesselGraphNode& n2 = edge.get().getNode2();

                const VesselGraphNode& new_node = &node == &n1 ? n2 : n1;
                VGNodeID new_node_id;

                auto added_node = added_nodes.find(&new_node);
                if(added_node != added_nodes.end()) {
                    new_node_id = added_node->second;
                } else {
                    new_node_id = output->insertNode(new_node);
                    added_nodes.insert({ &new_node, new_node_id });

                    uint32_t new_depth = depth+1;
                    node_queue.push(SubGraphNode(new_node, new_depth, new_node_id));
                }
                std::vector<VesselSkeletonVoxel> new_voxels(edge.get().getVoxels().begin(), edge.get().getVoxels().end());
                output->insertEdge(node_id, new_node_id, std::move(new_voxels));

                added_edges.insert(&edge.get());
            }
        }
    }

    return output;
}

void SubGraphExtractor::process() {
    const VesselGraph* input = inport_.getData();
    if(!input) {
        outport_.setData(nullptr);
        return;
    }
    if(!enabled_.get()) {
        // Pass input through
        outport_.setData(input, false);
        return;
    }

    tgt::vec3 startingPoint;
    if(!tryExtractCenterPoint(startingPoint_, startingPoint)) {
        LWARNING("No starting point!");
        outport_.setData(nullptr);
        return;
    }

    auto output = extractSubgraph(*input, startingPoint, maxEdgeDistance_.get(), keepBounds_.get());

    outport_.setData(output.release());
}

} // namespace voreen
