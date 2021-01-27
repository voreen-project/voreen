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

#include "vesselgraphrefinement.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/diskarraystorage.h"

#include "tgt/assert.h"

#include <boost/variant.hpp>

#include <unordered_map>
#include <unordered_set>

namespace voreen {

std::unique_ptr<VesselGraph> VesselGraphRefinement::removeEndEdgesRecursively(const VesselGraph& input, RemovableEdgeCheck isRemovableEdge, size_t maxIterations) {
    size_t num_removed_this_it = 0;
    auto output = removeEndEdges(input, isRemovableEdge, num_removed_this_it);
    for(size_t i=1; num_removed_this_it > 0 && i<maxIterations; ++i) {
        num_removed_this_it = 0;
        output = removeEndEdges(*output, isRemovableEdge, num_removed_this_it);
    }
    return output;
}

static std::vector<const VesselGraphEdge*> findDeletable(const VesselGraphNode& node, VesselGraphRefinement::RemovableEdgeCheck isRemovableEdge) {
    std::vector<const VesselGraphEdge*> potentially_deletable;
    const auto& edges = node.getEdges();
    if(edges.size() <= 2) {
        if(node.isAtSampleBorder_ && edges.size() == 1) {
            const auto& the_edge = edges[0].get();
            const auto& the_other_node = the_edge.getOtherNode(node);

            if(!the_edge.isLoop()
                    && !the_other_node.isAtSampleBorder_
                    && the_other_node.isEndNode()
                    && isRemovableEdge(the_edge)) {
                // The edge reaches from the sample border up to a single node in the volume that is not connected to anything else.
                // => It may be deletable
                //
                // In any other cases we cannot delete any edges without changing the topology.
                potentially_deletable.push_back(&the_edge);
            }
        }
        return potentially_deletable;
    }

    for(auto edge : edges) {
        if(isRemovableEdge(edge.get()) && edge.get().isEndStanding()) {
            potentially_deletable.push_back(&edge.get());
        }
    }
    size_t num_kept = edges.size() - potentially_deletable.size();
    if(num_kept >= 2) {
        return potentially_deletable; // We are deleting all candidates, as we retain >= 2 nodes anyways
    } if(num_kept == 1) {
        float highest_length = -std::numeric_limits<float>::infinity();
        const VesselGraphEdge* longest = nullptr;
        for(auto edge : potentially_deletable) {
            float l = edge->getRelativeBulgeSize();
            //tgtAssert(l >= 0, "Invalid bulge size");
            if(l > highest_length) {
                highest_length = l;
                longest = edge;
            }
        }
        std::vector<const VesselGraphEdge*> deletable;
        for(auto e: potentially_deletable) {
            if(e != longest) {
                deletable.push_back(e);
            }
        }
        tgtAssert(deletable.size() == potentially_deletable.size() - 1, "Invalid number of edges deleted"); // We are deleting exactly all except one candidate
        return deletable;
    } else { // => num_kept == 0
        float highest_length = -std::numeric_limits<float>::infinity();
        float second_highest_length = -std::numeric_limits<float>::infinity();
        const VesselGraphEdge* longest = nullptr;
        const VesselGraphEdge* second_longest = nullptr;
        for(auto edge : potentially_deletable) {
            float l = edge->getRelativeBulgeSize();
            //tgtAssert(l >= 0, "Invalid bulge size");
            if(l > highest_length) {
                second_highest_length = highest_length;
                second_longest = longest;
                highest_length = l;
                longest = edge;
            } else if (l > second_highest_length) {
                second_highest_length = l;
                second_longest = edge;
            }
        }
        std::vector<const VesselGraphEdge*> deletable;
        for(auto e : potentially_deletable) {
            if(e != longest && e != second_longest) {
                deletable.push_back(e);
            }
        }
        tgtAssert(deletable.size() == potentially_deletable.size() - 2, "Invalid number of edges deleted"); // We are deleting exactly all except two candidates edges
        return deletable;
    }
}

std::unique_ptr<VesselGraph> VesselGraphRefinement::removeEndEdges(const VesselGraph& input, RemovableEdgeCheck isRemovableEdge, size_t& num_removed_edges) {
    VesselGraphBuilder output(input.getBounds());
    num_removed_edges = 0;

    DiskArrayStorage<bool> edges_deleted_map(VoreenApplication::app()->getUniqueTmpFilePath(".newnodeids"));
    for(const VesselGraphEdge& _ : input.getEdges()) {
        edges_deleted_map.storeElement(false);
    }
    for(const VesselGraphNode& node : input.getNodes()) {
        for(const VesselGraphEdge* edge : findDeletable(node, isRemovableEdge)) {
            edges_deleted_map[edge->getID().raw()] = true;
        }
    }

    // Maps node ids in the old graph to node ids in the new graph
    DiskArrayStorage<VGNodeID> new_node_ids(VoreenApplication::app()->getUniqueTmpFilePath(".newnodeids"));
    for(const VesselGraphNode& _ : input.getNodes()) {
        new_node_ids.storeElement(VGNodeID::INVALID);
    }
    for(const VesselGraphEdge& edge : input.getEdges()) {
        if(edges_deleted_map[edge.getID().raw()]) {
            ++num_removed_edges;
            continue;
        }

        const VesselGraphNode& node1 = edge.getNode1();
        const VesselGraphNode& node2 = edge.getNode2();

        VGNodeID old_node1_id = node1.getID();
        VGNodeID old_node2_id = node2.getID();
        VGNodeID new_node1_id = VGNodeID::INVALID;
        VGNodeID new_node2_id = VGNodeID::INVALID;
        if(new_node_ids[old_node1_id.raw()].isValid()) {
            new_node1_id = new_node_ids[old_node1_id.raw()];
        } else {
            new_node1_id = output.insertNode(node1);
            new_node_ids[old_node1_id.raw()] = new_node1_id;
        }
        if(new_node_ids[old_node2_id.raw()].isValid()) {
            new_node2_id = new_node_ids[old_node2_id.raw()];
        } else {
            new_node2_id = output.insertNode(node2);
            new_node_ids[old_node2_id.raw()] = new_node2_id;
        }
        output.insertEdge(new_node1_id, new_node2_id, edge.getVoxels(), edge.getUUID());
    }
    // Insert all singular nodes as they cannot have been inserted before (via any edges), but should not be removed from the graph
    for(const VesselGraphNode& node : input.getNodes()) {
        if(node.getDegree() == 0) {
            output.insertNode(node);
        }
    }
    return removeDregree2Nodes(*std::move(output).finalize());
}

struct FutureNode {
    static const VGNodeID ID_NOT_ASSIGNED;
    FutureNode(const VesselGraphNode* node)
        : nodes_()
        , edges_()
        , id_(ID_NOT_ASSIGNED)
        , parent_(nullptr)
    {
        nodes_.insert(node);
    }
    //TODO: path compression optimization?
    void merge(FutureNode& other) {
        FutureNode* thisRoot = getRootNode();
        FutureNode* otherRoot = other.getRootNode();
        tgtAssert(thisRoot->parent_ == nullptr, "found non-root");
        tgtAssert(otherRoot->parent_ == nullptr, "found non-root");

        if(thisRoot == otherRoot) {
            // Nodes have already been merged
            return;
        }
        for(auto& node : otherRoot->nodes_) {
            thisRoot->nodes_.insert(node);
        }
        for(auto& edge : otherRoot->edges_) {
            thisRoot->edges_.insert(edge);
        }
        otherRoot->parent_ = thisRoot;
    }
    FutureNode* getRootNode() {
        if(parent_) {
            return parent_->getRootNode();
        } else {
            return this;
        }
    }
    std::unordered_set<const VesselGraphNode*> nodes_;
    std::unordered_set<const VesselGraphEdge*> edges_;
    VGNodeID id_; // Id in the future graph, will only be set when creating the graph
    FutureNode* parent_; // Points to other node if this has been merged with another node
};
const VGNodeID FutureNode::ID_NOT_ASSIGNED = -1;

inline bool isPartOfLoop(const VesselGraphEdge& edge) {
    const VesselGraphNode& node1 = edge.getNode1();
    const VesselGraphNode& node2 = edge.getNode2();
    for(const auto& sibling_edge : node1.getEdges()) {
        if(edge.getID() != sibling_edge.get().getID() //We do not want to check edge itself
                && (sibling_edge.get().getNodeID1() == node2.getID()
                    || sibling_edge.get().getNodeID2() == node2.getID())) {
            // Either node1 == node2 or another edge connects node1 or node2 => loop
            return true;
        }
    }
    return false;
}

std::unique_ptr<VesselGraph> VesselGraphRefinement::removeAllEdges(const VesselGraph& input, RemovableEdgeCheck isRemovableEdge) {

    // -----------------
    // Filtering and joining of nodes
    // -----------------

    std::unordered_map<const VesselGraphNode*, FutureNode*> new_nodes;
    std::vector<FutureNode> future_nodes;
    future_nodes.reserve(input.getNodes().size()); // Avoid pointer invalidation
    for(const VesselGraphNode& node : input.getNodes()) {
        future_nodes.emplace_back(FutureNode(&node));
        new_nodes.insert({ &node, &future_nodes.back()} );
    }

    std::vector<const VesselGraphEdge*> new_edges;
    for(const VesselGraphEdge& edge : input.getEdges()) {
        if(!isRemovableEdge(edge)) {
            new_edges.push_back(&edge);
        } else {
            const VesselGraphNode& node1 = edge.getNode1();
            const VesselGraphNode& node2 = edge.getNode2();
            if(node1.isEndNode() && !node2.isEndNode()) {
                // Only node1 edge is an end node: Remove the whole branch
                new_nodes.erase(&node1);
                continue;
            }
            if(!node1.isEndNode() && node2.isEndNode()) {
                // Only node2 is an end node: Remove the whole branch
                new_nodes.erase(&node2);
                continue;
            }
            tgtAssert(node1.isEndNode() == node2.isEndNode(), "Invalid node configuration");
            { // node1.isEndNode() == node2.isEndNode()
                if(isPartOfLoop(edge) || node1.isAtSampleBorder_ || node2.isAtSampleBorder_) {
                    // We do not want to change the topology of the graph and therefore cannot break off loops
                    // However, this check is relatively costly, so we only perform it now.
                    new_edges.push_back(&edge);
                    continue;
                }
                // Check for loop or multi edge configuration
                // Both nodes are (not) an end Node: Merge branch into a single node.
                FutureNode* future_node1 = new_nodes[&node1];
                FutureNode* future_node2 = new_nodes[&node2];
                future_node1->merge(*future_node2);
                future_node1->edges_.insert(&edge);
                new_nodes[&node2] = future_node1; //path compression optimization
            }
        }
    }

    // -----------------
    // Insertion
    // -----------------

    VesselGraphBuilder builder(input.getBounds());
    // Create the actual future nodes (aka nodes in the new graph) and insert them
    for(const auto& node_pair : new_nodes) {
        FutureNode* future_node = node_pair.second->getRootNode();
        tgtAssert(future_node, "no future node root");
        if(future_node->id_ == FutureNode::ID_NOT_ASSIGNED) {
            std::vector<tgt::vec3> voxels;
            std::vector<tgt::vec3> criticalVoxels;
            for(const auto& old_node : future_node->nodes_) {
                for(const auto& voxel: old_node->voxels_) {
                    voxels.push_back(voxel);
                }
            }
            for(const auto& old_edge : future_node->edges_) {
                for(const auto& voxel: old_edge->getVoxels()) {
                    voxels.push_back(voxel.pos_);
                }
            }
            tgt::vec3 new_pos = tgt::vec3::zero;
            for(const auto& voxel : voxels) {
                new_pos += voxel;
            }
            new_pos /= static_cast<float>(voxels.size());
            future_node->id_ = builder.insertNode(new_pos, std::move(voxels), 0.0f, false /* TODO: should be always not at sample border, right? */);
        }
    }

    // Insert all remaining edges into the new graph using the ids from previously created future nodes
    for(const auto& new_edge : new_edges) {
        VGNodeID future_node_id1 = new_nodes[&new_edge->getNode1()]->getRootNode()->id_;
        VGNodeID future_node_id2 = new_nodes[&new_edge->getNode2()]->getRootNode()->id_;

        tgtAssert(future_node_id1 != FutureNode::ID_NOT_ASSIGNED, "node1 without assigned id")
        tgtAssert(future_node_id2 != FutureNode::ID_NOT_ASSIGNED, "node2 without assigned id")

        std::vector<VesselSkeletonVoxel> voxels;
        //TODO: copy voxels from old edge and insert via ids
        for(auto& voxel: new_edge->getVoxels()) {
            voxels.push_back(voxel);
        }

        //TODO: check order of voxels and nodes do they have to be swapped?
        builder.insertEdge(future_node_id1, future_node_id2, std::move(voxels), new_edge->getUUID());
    }

    return removeDregree2Nodes(*std::move(builder).finalize());
}

struct FutureEdge {
    FutureEdge(const VesselGraphEdge& edge)
        : edges_(&edge)
        , begin_(&edge.getNode1())
        , end_(&edge.getNode2())
        , parent_(nullptr)
    {
    }

    void merge(FutureEdge& other) {
        FutureEdge* thisRoot = getRootNode();
        FutureEdge* otherRoot = other.getRootNode();
        tgtAssert(thisRoot->parent_ == nullptr, "found non-root");
        tgtAssert(otherRoot->parent_ == nullptr, "found non-root");

        if(thisRoot == otherRoot) {
            // Edges have already been merged
            return;
        }

        // Change both edges to be of multiple edge type:
        if(thisRoot->edges_.type() == typeid(const VesselGraphEdge*)) {
            const VesselGraphEdge* edge = boost::get<const VesselGraphEdge*>(thisRoot->edges_);
            std::deque<const VesselGraphEdge*> edges;
            edges.push_back(edge);
            thisRoot->edges_ = EdgesVariant(std::move(edges));
        }
        if(otherRoot->edges_.type() == typeid(const VesselGraphEdge*)) {
            const VesselGraphEdge* edge = boost::get<const VesselGraphEdge*>(otherRoot->edges_);
            std::deque<const VesselGraphEdge*> edges;
            edges.push_back(edge);
            otherRoot->edges_ = EdgesVariant(std::move(edges));
        }

        tgtAssert(edges_.type() == typeid(std::deque<const VesselGraphEdge*>), "Invalid multi edge type");
        tgtAssert(other.edges_.type() == typeid(std::deque<const VesselGraphEdge*>), "Invalid multi edge type");

        auto& thisq = boost::get<std::deque<const VesselGraphEdge*>>(thisRoot->edges_);
        auto& otherq = boost::get<std::deque<const VesselGraphEdge*>>(otherRoot->edges_);

        // This is to preserve the invariant of ordering in edges_
        if(thisRoot->end_ == otherRoot->begin_ && thisRoot->end_->getDegree() == 2){
            thisRoot->end_ = otherRoot->end_;
            for(auto edge = otherq.begin(); edge != otherq.end(); ++edge) {
                thisq.push_back(*edge);
            }
        } else if(thisRoot->end_ == otherRoot->end_ && thisRoot->end_->getDegree() == 2) {
            thisRoot->end_ = otherRoot->begin_;
            for(auto edge = otherq.rbegin(); edge != otherq.rend(); ++edge) {
                thisq.push_back(*edge);
            }
        } else if(thisRoot->begin_ == otherRoot->begin_ && thisRoot->begin_->getDegree() == 2) {
            thisRoot->begin_ = otherRoot->end_;
            for(auto edge = otherq.begin(); edge != otherq.end(); ++edge) {
                thisq.push_front(*edge);
            }
        } else if(thisRoot->begin_ == otherRoot->end_ && thisRoot->begin_->getDegree() == 2) {
            thisRoot->begin_ = otherRoot->begin_;
            for(auto edge = otherq.rbegin(); edge != otherq.rend(); ++edge) {
                thisq.push_front(*edge);
            }
        } else {
            tgtAssert(false, "non connected edges");
        }
        otherRoot->parent_ = thisRoot;

        tgtAssert(thisRoot->begin_ == &thisq.front()->getNode1() || thisRoot->begin_ == &thisq.front()->getNode2(), "merge failed");
        tgtAssert(thisRoot->end_ == &thisq.back()->getNode1() || thisRoot->end_ == &thisq.back()->getNode2(), "merge failed");
    }

    std::vector<VesselSkeletonVoxel> collectVoxels() {
        std::vector<VesselSkeletonVoxel> output;
        const VesselGraphNode* current_start_voxel = begin_;
        if(edges_.type() == typeid(const VesselGraphEdge*)) {
            auto& thisedge = boost::get<const VesselGraphEdge*>(edges_);
            output.insert(output.end(), thisedge->getVoxels().begin(), thisedge->getVoxels().end());
        } else {
            tgtAssert(edges_.type() == typeid(std::deque<const VesselGraphEdge*>), "Invalid edges type");
            auto& thisq = boost::get<std::deque<const VesselGraphEdge*>>(edges_);
            tgtAssert(thisq.size() >= 2, "Single edge in queue");

            const VesselGraphEdge* firstEdge = thisq.front();
            const VesselGraphEdge* lastEdge = thisq.back();

            // Generate fallback data;
            float minDistToSurfaceSum = 0.0;
            float avgDistToSurfaceSum = 0.0;
            float maxDistToSurfaceSum = 0.0;
            size_t totalVoxels = 0.0;
            for(const VesselGraphEdge* edge : thisq) {
                for(const auto& voxel : edge->getVoxels()) {
                    if(voxel.hasValidData()) {
                        minDistToSurfaceSum += voxel.minDistToSurface_;
                        avgDistToSurfaceSum += voxel.avgDistToSurface_;
                        maxDistToSurfaceSum += voxel.maxDistToSurface_;
                        totalVoxels += 1;
                    }
                }
            }
            float fallbackMinDistToSurface = minDistToSurfaceSum / totalVoxels;
            float fallbackAvgDistToSurface = avgDistToSurfaceSum / totalVoxels;
            float fallbackMaxDistToSurface = maxDistToSurfaceSum / totalVoxels;

            auto pushVoxel = [&output,
                 firstEdge,
                 lastEdge,
                 fallbackMinDistToSurface,
                 fallbackAvgDistToSurface,
                 fallbackMaxDistToSurface]
                     (const VesselGraphEdge& edge, size_t voxelIndex) {

                // We want all voxels to be outer if they are between the start of the first and the end of the last edge
                bool isFutureOuter =
                        (&edge == firstEdge && (voxelIndex >= edge.outerPathBeginIndex_ || edge.isEndStanding()))
                     || (&edge == lastEdge && (voxelIndex < edge.outerPathEndIndex_ || edge.isEndStanding()))
                     || (&edge != firstEdge && &edge != lastEdge);

                const VesselSkeletonVoxel& voxel = edge.getVoxels()[voxelIndex];
                if(isFutureOuter && voxel.isInner()) {
                    VesselSkeletonVoxel newVoxel(voxel); //copy
                    newVoxel.nearOtherEdge_ = false;

                    if(!voxel.hasValidData()) {
                        // We need to fix up the properties, as this voxel is supposed to be an outer voxel
                        // Best bet are the average values of the edge itself.
                        newVoxel.numSurfaceVoxels_ = 1;
                        if(edge.hasValidData()) {
                            newVoxel.minDistToSurface_ = edge.getMinRadiusAvg();
                            newVoxel.avgDistToSurface_ = edge.getAvgRadiusAvg();
                            newVoxel.maxDistToSurface_ = edge.getMaxRadiusAvg();
                        } else {
                            newVoxel.minDistToSurface_ = fallbackMinDistToSurface;
                            newVoxel.avgDistToSurface_ = fallbackAvgDistToSurface;
                            newVoxel.maxDistToSurface_ = fallbackMaxDistToSurface;
                        }
                    }

                    output.push_back(newVoxel);
                } else {
                    output.push_back(voxel);
                }
            };

            for(const VesselGraphEdge* edge : thisq) {
                // We know (and use the fact that): edge.getVoxels() are ordered from edge.getNode1() to edge.getNode2()
                if(current_start_voxel == &edge->getNode1()) {
                    for(size_t i=0; i < edge->getVoxels().size(); ++i) {
                        pushVoxel(*edge, i);
                    }
                    current_start_voxel = &edge->getNode2();
                } else {
                    tgtAssert(current_start_voxel == &edge->getNode2(), "edges_ are not ordered!");
                    for(size_t i = edge->getVoxels().size(); i!=0; --i) {
                        pushVoxel(*edge, i-1);
                    }
                    current_start_voxel = &edge->getNode1();
                }
            }
        }
        return output;
    }


    FutureEdge* getRootNode() {
        if(parent_) {
            return parent_->getRootNode();
        } else {
            return this;
        }
    }

    VesselGraphEdgeUUID getUUID() const {
        if(edges_.type() == typeid(const VesselGraphEdge*)) {
            auto& thisedge = boost::get<const VesselGraphEdge*>(edges_);
            return thisedge->getUUID();
        } else {
            tgtAssert(edges_.type() == typeid(std::deque<const VesselGraphEdge*>), "Invalid edges type");
            return VoreenApplication::app()->generateUUID();
        }
    }

    // Ordered so that they form a line
    // (but skeleton lines within can be in any direction!
    // See collectVoxels()!
    //
    // We use variant here because in the vast majority of cases edges are note merged and we save on heap memory this way
    typedef boost::variant<const VesselGraphEdge*, std::deque<const VesselGraphEdge*>> EdgesVariant;
    EdgesVariant edges_;




    const VesselGraphNode* begin_; //Points to the node that is not shared between any edges and is part of edges_.front()
    const VesselGraphNode* end_; //Points to the node that is not shared between any edges and is part of edges_.back()
    FutureEdge* parent_; // Points to other node if this has been merged with another node
};

std::unique_ptr<VesselGraph> VesselGraphRefinement::removeDregree2Nodes(const VesselGraph& input) {
    VesselGraphBuilder builder(input.getBounds());

    // 1. build futureedge for all edges, create hashmap edge -> futureedge
    DiskArrayStorage<FutureEdge> future_edges(VoreenApplication::app()->getUniqueTmpFilePath(".futureedges"));

    for(const auto& edge : input.getEdges()) {
        auto insertPos = future_edges.storeElement(FutureEdge(edge));
        tgtAssert(insertPos == edge.getID().raw(), "edges not properly ordered");
    }
    // 2. find deletable nodes, gather non deleted nodes
    //    also, merge edges via hashmap
    //
    //    insert all current nodes with either new node id or invalid id for deleted nodes
    DiskArrayStorage<VGNodeID> new_node_ids(VoreenApplication::app()->getUniqueTmpFilePath(".futureedges"));
    for(const auto& node : input.getNodes()) {
        size_t insertPos;
        if(node.getDegree() == 2) {
            auto node_edges = node.getEdges();
            tgtAssert(node_edges.size() == 2, "Invalid number of edges");
            const VesselGraphEdge& current_edge_1 = node_edges.at(0);
            const VesselGraphEdge& current_edge_2 = node_edges.at(1);
            FutureEdge& future_edge_1 = future_edges[current_edge_1.getID().raw()];
            FutureEdge& future_edge_2 = future_edges[current_edge_2.getID().raw()];

            future_edge_1.merge(future_edge_2);

            // Node will not be inserted: push invalid id
            insertPos = new_node_ids.storeElement(VGNodeID::INVALID);
        } else {
            VGNodeID node_id = builder.insertNode(node);
            insertPos = new_node_ids.storeElement(node_id);
        }
        tgtAssert(insertPos == node.getID().raw(), "nodes not properly ordered");
    }

    // 3. insert future edges using ends of gathered edges and map generated
    //    - the new VesselSkeletonVoxels can only be gathered from the edges for now,
    //      because they are guaranteed to form a skeletal line.
    for(auto& future_edge : future_edges.asArray()) {
        if(&future_edge != future_edge.getRootNode()) {
            continue; //Not a root component
        }
        VGNodeID node_id_begin = VGNodeID::INVALID;
        VGNodeID node_id_end = VGNodeID::INVALID;
        if(!new_node_ids[future_edge.begin_->getID().raw()].isValid()) {
            // we have to have detected a circle: add a single node for the edge:
            tgtAssert(future_edge.begin_ && future_edge.begin_ == future_edge.end_, "Non circle without nodes");
            const VesselGraphNode* circle_node = future_edge.begin_;
            VGNodeID circle_node_id = builder.insertNode(*circle_node);
            node_id_begin = circle_node_id;
            node_id_end = circle_node_id;
        } else {
            node_id_begin = new_node_ids[future_edge.begin_->getID().raw()];
            node_id_end = new_node_ids[future_edge.end_->getID().raw()];
        }
        tgtAssert(node_id_begin.isValid(), "Invalid node_id_begin");
        tgtAssert(node_id_end.isValid(), "Invalid node_id_begin");

        std::vector<VesselSkeletonVoxel> voxels = future_edge.collectVoxels();
        /*VGEdgeID inserted_edge =*/ builder.insertEdge(node_id_begin, node_id_end, std::move(voxels), future_edge.getUUID());
        //tgtAssert(builder.getEdge(inserted_edge).getLength() > 0, "Invalid edge length");
    }

    return std::move(builder).finalize();
}
} // namespace voreen
