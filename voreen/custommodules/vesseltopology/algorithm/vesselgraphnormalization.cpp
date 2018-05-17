/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2016 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "vesselgraphnormalization.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/assert.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace voreen {

std::unique_ptr<VesselGraph> VesselGraphNormalization::removeEndEdgesRecursively(const VesselGraph& input, RemovableEdgeCheck isRemovableEdge) {
    size_t num_removed_this_it = 0;
    auto output = removeEndEdges(input, isRemovableEdge, num_removed_this_it);
    while(num_removed_this_it > 0) {
        num_removed_this_it = 0;
        output = removeEndEdges(*output, isRemovableEdge, num_removed_this_it);
    }
    return output;
}

static std::vector<const VesselGraphEdge*> findDeletable(const VesselGraphNode& node, VesselGraphNormalization::RemovableEdgeCheck isRemovableEdge) {
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
        float highest_length = 0;
        const VesselGraphEdge* longest = nullptr;
        for(auto edge : potentially_deletable) {
            float l = edge->getLength();
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
        float highest_length = 0;
        float second_highest_length = 0;
        const VesselGraphEdge* longest = nullptr;
        const VesselGraphEdge* second_longest = nullptr;
        for(auto edge : potentially_deletable) {
            float l = edge->getLength();
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

std::unique_ptr<VesselGraph> VesselGraphNormalization::removeEndEdges(const VesselGraph& input, RemovableEdgeCheck isRemovableEdge, size_t& num_removed_edges) {
    std::unique_ptr<VesselGraph> output(new VesselGraph(input.getBounds()));
    num_removed_edges = 0;

    std::unordered_set<const VesselGraphEdge*> edges_to_delete;
    for(const VesselGraphNode& node : input.getNodes()) {
        for(const VesselGraphEdge* edge : findDeletable(node, isRemovableEdge)) {
            edges_to_delete.insert(edge);
        }
    }

    std::unordered_map<const VesselGraphNode*, size_t> new_node_ids;
    for(const VesselGraphEdge& edge : input.getEdges()) {
        if(edges_to_delete.count(&edge) > 0) {
            ++num_removed_edges;
            continue;
        }

        const VesselGraphNode& node1 = edge.getNode1();
        const VesselGraphNode& node2 = edge.getNode2();

        size_t new_node1_id = -1;
        size_t new_node2_id = -1;
        if(new_node_ids.count(&node1) > 0) {
            new_node1_id = new_node_ids[&node1];
        } else {
            new_node1_id = output->insertNode(node1);
            new_node_ids[&node1] = new_node1_id;
        }
        if(new_node_ids.count(&node2) > 0) {
            new_node2_id = new_node_ids[&node2];
        } else {
            new_node2_id = output->insertNode(node2);
            new_node_ids[&node2] = new_node2_id;
        }
        std::vector<VesselSkeletonVoxel> new_voxels(edge.getVoxels());
        output->insertEdge(new_node1_id, new_node2_id, std::move(new_voxels), edge.getUUID());
    }
    for(const VesselGraphNode& node : input.getNodes()) {
        if(node.getDegree() == 0) {
            output->insertNode(node);
        }
    }
    return removeDregree2Nodes(*output);
}

struct FutureNode {
    static const size_t ID_NOT_ASSIGNED = -1;
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
    size_t id_; // Id in the future graph, will only be set when creating the graph
    FutureNode* parent_; // Points to other node if this has been merged with another node
};

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

std::unique_ptr<VesselGraph> VesselGraphNormalization::removeAllEdges(const VesselGraph& input, RemovableEdgeCheck isRemovableEdge) {

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

    std::unique_ptr<VesselGraph> output(new VesselGraph(input.getBounds()));
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
            future_node->id_ = output->insertNode(new_pos, std::move(voxels), false /* TODO: should be always not at sample border, right? */);
        }
    }

    // Insert all remaining edges into the new graph using the ids from previously created future nodes
    for(const auto& new_edge : new_edges) {
        size_t future_node_id1 = new_nodes[&new_edge->getNode1()]->getRootNode()->id_;
        size_t future_node_id2 = new_nodes[&new_edge->getNode2()]->getRootNode()->id_;

        tgtAssert(future_node_id1 != FutureNode::ID_NOT_ASSIGNED, "node1 without assigned id")
        tgtAssert(future_node_id2 != FutureNode::ID_NOT_ASSIGNED, "node2 without assigned id")

        std::vector<VesselSkeletonVoxel> voxels;
        //TODO: copy voxels from old edge and insert via ids
        for(auto& voxel: new_edge->getVoxels()) {
            voxels.push_back(voxel);
        }

        //TODO: check order of voxels and nodes do they have to be swapped?
        output->insertEdge(future_node_id1, future_node_id2, std::move(voxels), new_edge->getUUID());
    }

    return removeDregree2Nodes(*output);
}

struct FutureEdge {
    FutureEdge(const VesselGraphEdge& edge)
        : edges_()
        , begin_(&edge.getNode1())
        , end_(&edge.getNode2())
        , parent_(nullptr)
    {
        edges_.push_back(&edge);
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

        // This is to preserve the invariant of ordering in edges_
        if(thisRoot->end_ == otherRoot->begin_ && thisRoot->end_->getDegree() == 2){
            thisRoot->end_ = otherRoot->end_;
            for(auto edge = otherRoot->edges_.begin(); edge != otherRoot->edges_.end(); ++edge) {
                thisRoot->edges_.push_back(*edge);
            }
        } else if(thisRoot->end_ == otherRoot->end_ && thisRoot->end_->getDegree() == 2) {
            thisRoot->end_ = otherRoot->begin_;
            for(auto edge = otherRoot->edges_.rbegin(); edge != otherRoot->edges_.rend(); ++edge) {
                thisRoot->edges_.push_back(*edge);
            }
        } else if(thisRoot->begin_ == otherRoot->begin_ && thisRoot->begin_->getDegree() == 2) {
            thisRoot->begin_ = otherRoot->end_;
            for(auto edge = otherRoot->edges_.begin(); edge != otherRoot->edges_.end(); ++edge) {
                thisRoot->edges_.push_front(*edge);
            }
        } else if(thisRoot->begin_ == otherRoot->end_ && thisRoot->begin_->getDegree() == 2) {
            thisRoot->begin_ = otherRoot->begin_;
            for(auto edge = otherRoot->edges_.rbegin(); edge != otherRoot->edges_.rend(); ++edge) {
                thisRoot->edges_.push_front(*edge);
            }
        } else {
            tgtAssert(false, "non connected edges");
        }
        otherRoot->parent_ = thisRoot;

        tgtAssert(thisRoot->begin_ == &thisRoot->edges_.front()->getNode1() || thisRoot->begin_ == &thisRoot->edges_.front()->getNode2(), "merge failed");
        tgtAssert(thisRoot->end_ == &thisRoot->edges_.back()->getNode1() || thisRoot->end_ == &thisRoot->edges_.back()->getNode2(), "merge failed");
    }

    std::vector<VesselSkeletonVoxel> collectVoxels() {
        std::vector<VesselSkeletonVoxel> output;
        const VesselGraphNode* current_start_voxel = begin_;
        for(const VesselGraphEdge* edge : edges_) {
            // assumption: edge.getVoxels() are ordered from edge.getNode1() to edge.getNode2()
            // TODO: check that assumption
            if(current_start_voxel == &edge->getNode1()) {
                output.insert(output.end(), edge->getVoxels().begin(), edge->getVoxels().end());
                current_start_voxel = &edge->getNode2();
            } else {
                tgtAssert(current_start_voxel == &edge->getNode2(), "edges_ are not ordered!");
                output.insert(output.end(), edge->getVoxels().rbegin(), edge->getVoxels().rend());
                current_start_voxel = &edge->getNode1();
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
        if(edges_.size() == 1) {
            return edges_[0]->getUUID();
        } else {
            return VoreenApplication::app()->generateUUID();
        }
    }

    std::deque<const VesselGraphEdge*> edges_; // Ordered so that they form a line
                                               // (but skeleton lines within can be in any direction!
                                               // See collectVoxels()!

    const VesselGraphNode* begin_; //Points to the node that is not shared between any edges and is part of edges_.front()
    const VesselGraphNode* end_; //Points to the node that is not shared between any edges and is part of edges_.back()
    FutureEdge* parent_; // Points to other node if this has been merged with another node
};

std::unique_ptr<VesselGraph> VesselGraphNormalization::removeDregree2Nodes(const VesselGraph& input) {
    std::unique_ptr<VesselGraph> output(new VesselGraph(input.getBounds()));

    // 1. build futureedge for all edges, create hashmap edge -> futureedge
    std::vector<FutureEdge> future_edges;
    future_edges.reserve(input.getEdges().size()); // preallocate to avoid pointer invalidation
    std::unordered_map<const VesselGraphEdge*, FutureEdge*> future_edge_map;
    for(const auto& edge : input.getEdges()) {
        future_edges.emplace_back(edge);
        future_edge_map.insert({&edge, &future_edges.back()});
    }
    // 2. find deletable nodes, gather non deleted nodes
    //    also, merge edges via hashmap
    //
    //    insert non-deleted nodes, generate map: old node pointer -> new id
    std::unordered_map<const VesselGraphNode*, size_t> new_node_ids;
    for(const auto& node : input.getNodes()) {
        if(node.getDegree() == 2) {
            auto node_edges = node.getEdges(); //TODO: go with ids if performance is a problem
            tgtAssert(node_edges.size() == 2, "Invalid number of edges");
            const VesselGraphEdge& current_edge_1 = node.getEdges().at(0);
            const VesselGraphEdge& current_edge_2 = node.getEdges().at(1);
            FutureEdge* future_edge_1 = future_edge_map.at(&current_edge_1);
            FutureEdge* future_edge_2 = future_edge_map.at(&current_edge_2);
            tgtAssert(future_edge_1, "no future_edge_1");
            tgtAssert(future_edge_2, "no future_edge_2");

            future_edge_1->merge(*future_edge_2);
        } else {
            size_t node_id = output->insertNode(node);
            new_node_ids.insert({&node, node_id});
        }
    }

    // 3. insert future edges using ends of gathered edges and map generated
    //    - the new VesselSkeletonVoxels can only be gathered from the edges for now,
    //      because they are guaranteed to form a skeletal line.
    for(auto& future_edge : future_edges) {
        if(&future_edge != future_edge.getRootNode()) {
            continue; //Not a root component
        }
        const size_t NODE_NOT_FOUND = -1;
        size_t node_id_begin = NODE_NOT_FOUND;
        size_t node_id_end = NODE_NOT_FOUND;
        if(new_node_ids.count(future_edge.begin_) == 0) {
            // we have to have detected a circle: add a single node for the edge:
            tgtAssert(future_edge.begin_ && future_edge.begin_ == future_edge.end_, "Non circle without nodes");
            const VesselGraphNode* circle_node = future_edge.begin_;
            size_t circle_node_id = output->insertNode(*circle_node);
            node_id_begin = circle_node_id;
            node_id_end = circle_node_id;
        } else {
            try{
                node_id_begin = new_node_ids.at(future_edge.begin_);
            } catch(...) {
                tgtAssert(false, "node lookup failed for begin");
            }
            try {
                node_id_end = new_node_ids.at(future_edge.end_);
            } catch(...) {
                tgtAssert(false, "node lookup failed for begin");
            }
        }
        tgtAssert(node_id_begin != NODE_NOT_FOUND, "Invalid node_id_begin");
        tgtAssert(node_id_end != NODE_NOT_FOUND, "Invalid node_id_begin");

        std::vector<VesselSkeletonVoxel> voxels = future_edge.collectVoxels();
        size_t inserted_edge = output->insertEdge(node_id_begin, node_id_end, std::move(voxels), future_edge.getUUID());
        tgtAssert(output->getEdge(inserted_edge).getLength() > 0, "Invalid edge length");
    }


    //TODO check order of voxels in edges when merging!
    return output;
}
} // namespace voreen
