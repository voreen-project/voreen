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

#include "vesselgraphperturbation.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "../algorithm/vesselgraphrefinement.h"

#include <unordered_set>
#include <set>
#include <unordered_map>


namespace {
    typedef std::mt19937 random_engine_type;
}
namespace voreen {

const std::string VesselGraphPerturbation::loggerCat_("voreen.vesselgraphperturbation");


VesselGraphPerturbation::VesselGraphPerturbation()
    : Processor()
    , inport_(Port::INPORT, "graph.input", "Graph Input", false, Processor::INVALID_RESULT)
    , outport_(Port::OUTPORT, "graph.output", "Perturbed Graph Output", false, Processor::VALID)
    , enabled_("enabled", "Enabled", true)
    , perturbationMethod_("perturbationMethod", "Perturbation Method")
    , perturbationAmount_("perturbationAmount", "Perturbation Amount", 0.5, 0.0, 1.0)
    , usePredeterminedSeed_("usePredeterminedSeed", "Use Predetermined Seed", false)
    , predeterminedSeed_("predeterminedSeed", "Seed", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , randomDevice()
{

    addPort(inport_);
    addPort(outport_);

    addProperty(enabled_);

    addProperty(perturbationMethod_);
        perturbationMethod_.addOption("add_edges", "Add Edges", ADD_EDGES);
        perturbationMethod_.addOption("split_nodes", "Split Nodes", SPLIT_NODES);
        perturbationMethod_.addOption("subdivide_edges", "Subdivide Edges", SUBDIVIDE_EDGES);
        perturbationMethod_.addOption("split_edges", "Split Edges", SPLIT_EDGES);
        perturbationMethod_.addOption("move_nodes", "Move Nodes", MOVE_NODES);
        perturbationMethod_.addOption("change_properties", "Change Properties", CHANGE_PROPERTIES);
        perturbationMethod_.addOption("combined", "Combined", COMBINED);
        perturbationMethod_.selectByValue(ADD_EDGES);

    addProperty(perturbationAmount_);

    addProperty(usePredeterminedSeed_);
        ON_CHANGE_LAMBDA(usePredeterminedSeed_, [this] {
                    predeterminedSeed_.setVisibleFlag(usePredeterminedSeed_.get());
                });
    addProperty(predeterminedSeed_);
}

VesselGraphPerturbation::~VesselGraphPerturbation() {
}

typedef std::normal_distribution<float> normal_distribution;
typedef std::lognormal_distribution<float> lognormal_distribution;

static normal_distribution getGraphEdgePropertyDestribution(const VesselGraph& graph, std::function<float(const VesselGraphEdge&)> property) {
    float mean, stddev;
    graph.getEdgePropertyStats(property, mean, stddev);
    tgtAssert(std::isnormal(mean), "Invalid mean");
    tgtAssert(std::isnormal(stddev), "Invalid stddev");
    graph.getEdgePropertyStats(property, mean, stddev);
    return normal_distribution(mean, stddev);
}

static float make_positive(float val) {
    return std::max(std::numeric_limits<float>::epsilon(), val);
}

static float make_non_negative(float val) {
    return std::max(0.0f, val);
}

static float clamp_01(float val) {
    return std::min(1.0f, std::max(0.0f, val));
}

struct EdgeGenerator {
    normal_distribution distr_distance_;
    normal_distribution distr_length_;
    normal_distribution distr_straightness_;
    normal_distribution distr_minRadiusAvg_;
    normal_distribution distr_minRadiusStdDeviation_;
    normal_distribution distr_maxRadiusAvg_;
    normal_distribution distr_maxRadiusStdDeviation_;
    normal_distribution distr_avgRadiusAvg_;
    normal_distribution distr_avgRadiusStdDeviation_;
    normal_distribution distr_roundnessAvg_;
    normal_distribution distr_roundnessStdDeviation_;

    EdgeGenerator(const VesselGraph& templateGraph)
        : distr_distance_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getDistance)))
        , distr_length_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getLength)))
        , distr_straightness_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getStraightness)))
        , distr_minRadiusAvg_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getMinRadiusAvg)))
        , distr_minRadiusStdDeviation_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getMinRadiusStdDeviation)))
        , distr_maxRadiusAvg_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getMaxRadiusAvg)))
        , distr_maxRadiusStdDeviation_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getMaxRadiusStdDeviation)))
        , distr_avgRadiusAvg_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getAvgRadiusAvg)))
        , distr_avgRadiusStdDeviation_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getAvgRadiusStdDeviation)))
        , distr_roundnessAvg_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getRoundnessAvg)))
        , distr_roundnessStdDeviation_(getGraphEdgePropertyDestribution(templateGraph, std::mem_fn(&VesselGraphEdge::getRoundnessStdDeviation)))
    {
    }

    VesselGraphEdgePathProperties generateEdgePathProperties(random_engine_type random_engine, float distance) {
        VesselGraphEdgePathProperties output;
        //TODO: We might want to do something more sophisticated here and look at correlations between properties.
        //      Might be some amount of work...
        if(distance > 0) {
            output.length_ = distance/distr_straightness_(random_engine);
        } else {
            output.length_ = make_positive(distr_length_(random_engine)); // Account for loops
        }
        tgtAssert(output.length_ >= 0, "Invalid length");
        std::vector<float> radii;
        radii.push_back(make_positive(distr_minRadiusAvg_(random_engine)));
        radii.push_back(make_positive(distr_avgRadiusAvg_(random_engine)));
        radii.push_back(make_positive(distr_maxRadiusAvg_(random_engine)));
        std::sort(radii.begin(), radii.end());
        output.minRadiusAvg_ = radii.at(0);
        output.minRadiusStdDeviation_ = make_non_negative(distr_minRadiusStdDeviation_(random_engine));
        output.avgRadiusAvg_ = radii.at(1);
        output.avgRadiusStdDeviation_ = make_non_negative(distr_avgRadiusStdDeviation_(random_engine));
        output.maxRadiusAvg_ = radii.at(2);
        output.maxRadiusStdDeviation_ = make_non_negative(distr_maxRadiusStdDeviation_(random_engine));
        output.roundnessAvg_ = clamp_01(distr_roundnessAvg_(random_engine));
        output.roundnessStdDeviation_ = make_non_negative(distr_roundnessStdDeviation_(random_engine));

        output.volume_ = output.length_*output.avgRadiusAvg_*output.avgRadiusAvg_;
        return output;
    }
};

static VesselGraphEdgePathProperties modifyEdgePathProperties(random_engine_type& random_engine, float amount, const VesselGraphEdgePathProperties& base) {
    VesselGraphEdgePathProperties output;
    // Generate the edge that subdivides the first
    auto dist = lognormal_distribution(0, amount);

    auto transform = [&random_engine, &dist] (float val) {
        if(val == VesselGraphEdgePathProperties::INVALID_DATA) {
            return val;
        } else {
            float new_val = dist(random_engine)*val;
            tgtAssert(new_val >= 0, "Negative property");
            return new_val;
        }
    };

    output.length_ = transform(base.length_);
    tgtAssert(output.length_ >= 0, "Invalid length");

    output.minRadiusAvg_ = transform(base.minRadiusAvg_);
    output.minRadiusStdDeviation_ = transform(base.minRadiusStdDeviation_);
    output.maxRadiusAvg_ = transform(base.maxRadiusAvg_);
    output.maxRadiusStdDeviation_ = transform(base.maxRadiusStdDeviation_);
    output.avgRadiusAvg_ = transform(base.avgRadiusAvg_);
    output.avgRadiusStdDeviation_ = transform(base.avgRadiusStdDeviation_);
    output.roundnessAvg_ = transform(base.roundnessAvg_);
    output.roundnessStdDeviation_ = transform(base.roundnessStdDeviation_);

    output.volume_ = output.length_*output.avgRadiusAvg_*output.avgRadiusAvg_;

    return output;
}

template<class T>
static std::vector<T> draw_n(std::vector<T>&& input, size_t n, random_engine_type& random_engine) {
    tgtAssert(n <= input.size(), "Invalid number of elements to draw");

    std::vector<T> output;
    for(size_t i=0; i<n; ++i) {
        auto& drawn = input.at(std::uniform_int_distribution<size_t>(0, input.size()-1)(random_engine));
        output.push_back(drawn);
        // Removal via the clausing method:
        std::swap(drawn, input.back());
        input.pop_back();
    }
    return output;
}

static std::vector<const VesselGraphEdge*> draw_edges(const std::vector<const VesselGraphEdge*>& edges, float percentage, random_engine_type random_engine) {
    std::vector<const VesselGraphEdge*> remainig_edge_refs;
    for(const auto& edge : edges) {
        remainig_edge_refs.push_back(edge);
    }
    return draw_n(std::move(remainig_edge_refs), std::floor(percentage*edges.size()), random_engine);
}

static std::vector<const VesselGraphEdge*> select_edges_by_uuid(const VesselGraph& graph, const std::set<VesselGraphEdgeUUID>& uuids) {
    std::vector<const VesselGraphEdge*> result;
    for(const auto& edge : graph.getEdges()) {
        if(uuids.count(edge.getUUID()) > 0) {
            result.push_back(&edge);
        }
    }
    std::sort(result.begin(), result.end(), [] (const VesselGraphEdge* e1, const VesselGraphEdge* e2) {
            return e1->getUUID() < e2->getUUID();
            });
    return result;
}

static std::vector<const VesselGraphNode*> select_nodes_by_uuid(const VesselGraph& graph, const std::set<VesselGraphNodeUUID>& uuids) {
    std::vector<const VesselGraphNode*> result;
    for(const auto& node : graph.getNodes()) {
        if(uuids.count(node.getUUID()) > 0) {
            result.push_back(&node);
        }
    }
    std::sort(result.begin(), result.end(), [] (const VesselGraphNode* e1, const VesselGraphNode* e2) {
            return e1->getID() < e2->getID();
            });
    return result;
}

static std::set<VesselGraphEdgeUUID> get_edge_uuids(const VesselGraph& graph) {
    std::set<VesselGraphEdgeUUID> result;
    for(const auto& edge : graph.getEdges()) {
        result.insert(edge.getUUID());
    }
    return result;
}

static std::set<VesselGraphEdgeUUID> get_node_uuids(const VesselGraph& graph) {
    std::set<VesselGraphNodeUUID> result;
    for(const auto& node : graph.getNodes()) {
        result.insert(node.getUUID());
    }
    return result;
}

static std::vector<const VesselGraphNode*> draw_nodes_if(const std::vector<const VesselGraphNode*>& nodes, float percentage, std::function<bool(const VesselGraphNode&)> viable, random_engine_type random_engine) {

    std::vector <const VesselGraphNode*> viable_nodes;
    for(auto& node: nodes) {
        if(viable(*node)) {
            viable_nodes.push_back(node);
        }
    }
    return draw_n(std::move(viable_nodes), std::floor(percentage*viable_nodes.size()), random_engine);
}

static tgt::vec3 generate_random_direction(random_engine_type& random_engine) {
    std::uniform_real_distribution<float> direction_distr(-1, 1);
    tgt::vec3 new_direction;
    do {
        new_direction = tgt::vec3(
                direction_distr(random_engine),
                direction_distr(random_engine),
                direction_distr(random_engine)
                );
    } while (tgt::lengthSq(new_direction) > 1 || new_direction == tgt::vec3::zero);
    return tgt::normalize(new_direction);
}


static std::unique_ptr<VesselGraph> addEdgesToNodesWithUUIDs(const VesselGraph& input, const std::set<VesselGraphNodeUUID>& uuids, float amount, random_engine_type random_engine) {
    EdgeGenerator edgeGenerator(input);
    //Find pairs of nodes that are both connected to one other node:
    std::vector<std::pair<const VesselGraphNode*, const VesselGraphNode*>> pairs;
    for(auto& node_ptr : select_nodes_by_uuid(input, uuids)) {
        const VesselGraphNode& node = *node_ptr;
        for(auto& b1: node.getEdges()) {
            for(auto& b2: node.getEdges()) {
                if(&b1 == &b2 || b1.get().isLoop() || b2.get().isLoop()) {
                    continue;
                }
                const VesselGraphNode& node1 = b1.get().getOtherNode(node);
                const VesselGraphNode& node2 = b2.get().getOtherNode(node);
                if(&node1 == &node2) {
                    continue;
                }
                pairs.push_back(std::make_pair(&node1, &node2));
            }
        }
    }
    //This may happen in theory, but hopefully not in our case. We want to be able to draw for up to amount=1.0f
    //tgtAssert(pairs.size() > input.getEdges().size(), "Less insert candidates than edges!");

    auto forked_re = random_engine;
    auto selected_pairs = draw_n(std::move(pairs), std::min(pairs.size(), static_cast<size_t>(std::floor(amount*input.getEdges().size()))), forked_re);

    VesselGraphBuilder builder;
    for(auto pair: selected_pairs) {
        const VesselGraphNode& node1 = *pair.first;
        const VesselGraphNode& node2 = *pair.second;

        float distance = tgt::distance(node1.pos_, node2.pos_);
        builder.insertEdge(node1.getID(), node2.getID(), edgeGenerator.generateEdgePathProperties(random_engine, distance));
    }
    return std::move(builder).finalize();
}
static std::unique_ptr<VesselGraph> addEdges(const VesselGraph& input, float amount, random_engine_type random_engine) {
    return addEdgesToNodesWithUUIDs(input, get_node_uuids(input), amount, random_engine);
}
struct NodeSplitEdgeConfiguration {
    VGNodeID node1_id;
    VGNodeID node2_id;
    NodeSplitEdgeConfiguration(const VesselGraphEdge& edge)
        : node1_id(edge.getNodeID1())
        , node2_id(edge.getNodeID2())
    {
    }
};

static std::unique_ptr<VesselGraph> splitNodesWithUUIDs(const VesselGraph& input, const std::set<VesselGraphNodeUUID>& uuids, float amount, random_engine_type random_engine) {

    auto nodes_to_split = draw_nodes_if(select_nodes_by_uuid(input, uuids), amount, [] (const VesselGraphNode& n) {
            return n.getDegree() > 2;
            }, random_engine);

    std::unordered_set<const VesselGraphNode*> nodes_to_split_lookup(nodes_to_split.begin(), nodes_to_split.end());

    VesselGraphBuilder builder;
    for(auto& node : input.getNodes()) {
        builder.insertNode(node);
    }

    std::unordered_map<const VesselGraphEdge*, NodeSplitEdgeConfiguration> split_off_edges;
    for(auto node_ptr : nodes_to_split) {
        const auto& node = *node_ptr;
        //TODO: we should consider to "nudge" the node a little bit, if we want to support splitting with
        //other configurations than n=2+s (see below).
        VGNodeID split_node_id = builder.insertNode(node);

        auto edges_to_split_off = draw_n(node.getEdges(), 2, random_engine);

        for(auto edge_ref: edges_to_split_off) {
            const auto& edge = edge_ref.get();
            if(split_off_edges.count(&edge) == 0) {
                split_off_edges.insert(std::make_pair(&edge, NodeSplitEdgeConfiguration(edge)));
            }
            NodeSplitEdgeConfiguration& edgeconf = split_off_edges.at(&edge);
            if(edge.isLoop()) {
                if(std::uniform_int_distribution<int>(0,1)(random_engine) == 0) {
                    // The edge is a loop, so getNodeID1() == getNodeID2().
                    // Therefore, we just choose one of the "ends" of the edge to try to replace it with the newly
                    // generated node id, and if that fails we choose the other.
                    if(edgeconf.node1_id != split_node_id) {
                        edgeconf.node1_id = split_node_id;
                    } else {
                        tgtAssert(edgeconf.node2_id != split_node_id, "something is wrong with node splitting (loop)");
                        edgeconf.node2_id = split_node_id;
                    }
                } else {
                    if(edgeconf.node2_id != split_node_id) {
                        edgeconf.node2_id = split_node_id;
                    } else {
                        tgtAssert(edgeconf.node1_id != split_node_id, "something is wrong with node splitting (loop)");
                        edgeconf.node1_id = split_node_id;
                    }
                }
            } else {
                // The edge is NOT a loop, so we know getNodeID1() != getNodeID2().
                // Thus, the edge is only connected to this node once, and we search for the node id that
                // matches the id of the node to be split.
                if(edgeconf.node1_id == node.getID()) {
                    edgeconf.node1_id = split_node_id;
                } else {
                    tgtAssert(edgeconf.node2_id == node.getID(), "something is wrong with node splitting");
                    edgeconf.node2_id = split_node_id;
                }
            }
        }
    }

    for(auto& edge : input.getEdges()) {
        VGNodeID node1_id, node2_id;
        if(split_off_edges.count(&edge) > 0) {
            NodeSplitEdgeConfiguration& edgeconf = split_off_edges.at(&edge);
            node1_id = edgeconf.node1_id;
            node2_id = edgeconf.node2_id;
        } else {
            node1_id = edge.getNodeID1();
            node2_id = edge.getNodeID2();
        }
        tgtAssert(edge.getLength() > 0, "Invalid edge length");
        builder.insertEdge(node1_id, node2_id, edge.getVoxels(), edge.getUUID());
    }

    return VesselGraphRefinement::removeDregree2Nodes(*std::move(builder).finalize());
}
static std::unique_ptr<VesselGraph> splitNodes(const VesselGraph& input, float amount, random_engine_type random_engine) {
    return splitNodesWithUUIDs(input, get_node_uuids(input), amount, random_engine);
}

struct ContinuousPathPos {
    int index;
    float interNodePos;
    const VesselGraphEdge& edge;

    //TODO: We're not really sampling splits evenly across the length of the edge, as the
    //      distance between nodes is not constant.
    ContinuousPathPos(const VesselGraphEdge& edge, random_engine_type& random_engine)
        : index(std::uniform_int_distribution<int>(0, edge.getVoxels().size()-2)(random_engine))
        , interNodePos(std::uniform_real_distribution<float>(
                0 + std::numeric_limits<float>::epsilon(),    // avoid zero length edges
                1 - std::numeric_limits<float>::epsilon() // also for 1
                )(random_engine))
        , edge(edge)
    {
        tgtAssert(edge.getVoxels().size() >= 2, "Need at least two voxels to split");
        tgtAssert(0.0 < interNodePos && interNodePos < 1.0f, "Invalid interNodePos");
    }

    tgt::vec3 splitLocation() const {
        const auto& path(edge.getVoxels());
        int path_size = path.size();
        tgtAssert(0 <= index && index < path_size, "invalid index");
        //tgt::vec3 split_pos_base_l = (index == -1)          ? edge.getNode1().pos_ : path.at(index).pos_;
        //tgt::vec3 split_pos_base_r = (index+1 == path_size) ? edge.getNode2().pos_ : path.at(index+1).pos_;
        tgt::vec3 split_pos_base_l = path.at(index).pos_;
        tgt::vec3 split_pos_base_r = path.at(index+1).pos_;

        tgtAssert(split_pos_base_l != split_pos_base_r, "Same split pos base");

        return split_pos_base_l*(1-interNodePos) + split_pos_base_r*(interNodePos); //mix
    }

    std::vector<VesselSkeletonVoxel> splitPathLeft() const {
        const auto& path = edge.getVoxels();
        return std::vector<VesselSkeletonVoxel>(path.begin(), path.begin() + index + 1);
    }

    std::vector<VesselSkeletonVoxel> splitPathRight() const {
        const auto& path = edge.getVoxels();
        return std::vector<VesselSkeletonVoxel>(path.begin() + index + 1, path.end());
    }

    bool operator<(const ContinuousPathPos& other) const {
        tgtAssert(&edge == &other.edge, "Comparing paths from different edges");
        if(index != other.index) {
            return index < other.index;
        } else {
            return interNodePos < other.interNodePos;
        }
    }
};


static std::unique_ptr<VesselGraph> subdivideEdgesWithUUIDs(const VesselGraph& input, const std::set<VesselGraphEdgeUUID>& allowed_uuids, float amount, random_engine_type random_engine) {
    std::vector<const VesselGraphEdge*> input_edges;
    for(const auto& edge : select_edges_by_uuid(input, allowed_uuids)) {
        if(edge->getVoxels().size() >= 2) {
            input_edges.push_back(edge);
        }
    }
    std::vector<const VesselGraphEdge*> edges_to_subdivide = draw_edges(input_edges, amount, random_engine);
    std::unordered_set<const VesselGraphEdge*> edges_to_subdivide_lookup(edges_to_subdivide.begin(), edges_to_subdivide.end());

    VesselGraphBuilder builder;
    for(auto& node : input.getNodes()) {
        builder.insertNode(node);
    }

    // First: Insert all non-subdivided edges
    for(auto& edge : input.getEdges()) {
        if(edges_to_subdivide_lookup.count(&edge) == 0) {
            // Note: We ware using the fact that we are iterating over nodes (and inserted them) in order of ID!
            //       The same applies below.
            builder.insertEdge(edge.getNodeID1(), edge.getNodeID2(), edge, edge.getUUID());
        }
    }

    // Now: Subdivide selected edges and insert resulting edges
    for(auto edge_ptr : edges_to_subdivide) {
        auto& edge = *edge_ptr;
        const auto& path(edge.getVoxels());

        ContinuousPathPos split(edge, random_engine);
        tgt::vec3 split_pos = split.splitLocation();

        VGNodeID split_node_id = builder.insertNode(split_pos, std::vector<tgt::vec3>(), 0.0f, false);

        // First path part
        auto split_left = split.splitPathLeft();
        if(std::count_if(split_left.begin(), split_left.end(), [] (const VesselSkeletonVoxel& v) { return v.hasValidData(); }) > 0) {
            builder.insertEdge(edge.getNodeID1(), split_node_id, std::move(split_left));
        } else {
            VesselGraphEdgePathProperties properties = edge.getPathProperties();
            // Here (and below for the second path) we can really just guess that most of the properties
            // are the same as in the rest of the edge/vessel segment
            properties.length_ = tgt::distance(edge.getNode1().pos_, split_pos);
            properties.volume_ = properties.length_*edge.getAvgCrossSection();
            builder.insertEdge(edge.getNodeID1(), split_node_id, properties);
            tgtAssert(properties.length_ >= 0, "Invalid length");
        }

        // Second path part
        auto split_right = split.splitPathRight();
        if(std::count_if(split_left.begin(), split_left.end(), [] (const VesselSkeletonVoxel& v) { return v.hasValidData(); }) > 0) {
            builder.insertEdge(split_node_id, edge.getNodeID2(), std::move(split_right));
        } else {
            VesselGraphEdgePathProperties properties = edge.getPathProperties();
            // See above
            properties.length_ = tgt::distance(split_pos, edge.getNode2().pos_);
            properties.volume_ = properties.length_*edge.getAvgCrossSection();
            builder.insertEdge(split_node_id, edge.getNodeID2(), properties);
            tgtAssert(properties.length_ >= 0, "Invalid length");
        }

        // Generate the edge that subdivides the first
        float new_length = std::abs(normal_distribution(0, edge.getDistance()/2)(random_engine))
                            + std::numeric_limits<float>::epsilon(); // has to be > 0

        tgt::vec3 new_direction = generate_random_direction(random_engine);

        tgt::vec3 new_node_pos = split_pos + new_direction*new_length;
        VGNodeID new_node_id = builder.insertNode(new_node_pos, std::vector<tgt::vec3>(), 0.0f, false);

        VesselGraphEdgePathProperties new_edge_properties = edge.getPathProperties();
        // Again, we assume that radius, roundness, etc. do not change
        // TODO: Maybe we want to randomize these using EdgeGenerator?
        new_edge_properties.length_ = new_length;
        tgtAssert(new_edge_properties.length_ >= 0, "Invalid length");
        new_edge_properties.volume_ = new_edge_properties.length_*edge.getAvgCrossSection();
        //tgtAssert(new_edge_properties.hasValidData(), "invalid new edge data");
        builder.insertEdge(split_node_id, new_node_id, new_edge_properties);
    }
    return std::move(builder).finalize();
}

static std::unique_ptr<VesselGraph> subdivideEdges(const VesselGraph& input, float amount, random_engine_type random_engine) {
    return subdivideEdgesWithUUIDs(input, get_edge_uuids(input), amount, random_engine);
}

static std::unique_ptr<VesselGraph> splitEdgesWithUUIDs(const VesselGraph& input, const std::set<VesselGraphEdgeUUID>& allowed_uuids, float amount, random_engine_type random_engine) {
    std::vector<const VesselGraphEdge*> input_edges;
    for(const auto& edge : select_edges_by_uuid(input, allowed_uuids)) {
        if(edge->getVoxels().size() >= 2) {
            input_edges.push_back(edge);
        }
    }
    std::vector<const VesselGraphEdge*> edges_to_split = draw_edges(input_edges, amount, random_engine);
    std::unordered_set<const VesselGraphEdge*> edges_to_split_lookup(edges_to_split.begin(), edges_to_split.end());

    VesselGraphBuilder builder;
    for(auto& node : input.getNodes()) {
        builder.insertNode(node);
    }

    // First: Insert all non-split edges
    for(auto& edge : input.getEdges()) {
        if(edges_to_split_lookup.count(&edge) == 0) {
            // Note: We ware using the fact that we are iterating over nodes (and inserted them) in order of ID!
            //       The same applies below.
            builder.insertEdge(edge.getNodeID1(), edge.getNodeID2(), edge, edge.getUUID());
        }
    }

    // Now: Subdivide selected edges and insert resulting edges
    for(auto edge_ptr : edges_to_split) {
        auto& edge = *edge_ptr;
        const auto& path(edge.getVoxels());

        ContinuousPathPos split1(edge, random_engine);
        ContinuousPathPos split2(edge, random_engine);

        ContinuousPathPos& split_l = (split1 < split2) ? split1: split2;
        ContinuousPathPos& split_r = (split1 < split2) ? split2: split1;

        {
            // First path part
            tgt::vec3 split_pos_l = split_l.splitLocation();
            VGNodeID split_node_id_l = builder.insertNode(split_pos_l, std::vector<tgt::vec3>(), 0.0f, false);

            auto path_left = split_l.splitPathLeft();
            if(path_left.empty()) {
                VesselGraphEdgePathProperties properties = edge.getPathProperties();
                // Here (and below for the second path) we can really just guess that most of the properties
                // are the same as in the rest of the edge/vessel segment
                properties.length_ = tgt::distance(edge.getNode1().pos_, split_pos_l);
                properties.volume_ = properties.length_*edge.getAvgCrossSection();
                builder.insertEdge(edge.getNodeID1(), split_node_id_l, properties);
                tgtAssert(properties.length_ >= 0, "Invalid length");
            } else {
                builder.insertEdge(edge.getNodeID1(), split_node_id_l, std::move(path_left));
            }
        }

        {
            // Second path part
            tgt::vec3 split_pos_r = split_r.splitLocation();
            VGNodeID split_node_id_r = builder.insertNode(split_pos_r, std::vector<tgt::vec3>(), 0.0f, false);

            auto path_right = split_r.splitPathRight();
            if(path_right.empty()) {
                VesselGraphEdgePathProperties properties = edge.getPathProperties();
                // See above
                properties.length_ = tgt::distance(split_pos_r, edge.getNode2().pos_);
                properties.volume_ = properties.length_*edge.getAvgCrossSection();
                builder.insertEdge(split_node_id_r, edge.getNodeID2(), properties);
                tgtAssert(properties.length_ >= 0, "Invalid length");
            } else {
                builder.insertEdge(split_node_id_r, edge.getNodeID2(), std::move(path_right));
            }
        }
    }
    return std::move(builder).finalize();
}
static std::unique_ptr<VesselGraph> splitEdges(const VesselGraph& input, float amount, random_engine_type random_engine) {
    return splitEdgesWithUUIDs(input, get_edge_uuids(input), amount, random_engine);
}

static std::unique_ptr<VesselGraph> moveNodes(const VesselGraph& input, float amount, random_engine_type random_engine) {
    float length_sum = 0;
    for(auto& edge : input.getEdges()) {
        length_sum += edge.getLength();
    }
    float length_avg = length_sum / input.getEdges().size();

    float perturbation = 2*amount*length_avg;

    VesselGraphBuilder builder;
    for(auto& node : input.getNodes()) {
        // Generate the edge that subdivides the first
        float shift_length = std::abs(normal_distribution(0, perturbation/2)(random_engine));

        tgt::vec3 shift = shift_length * generate_random_direction(random_engine);

        std::vector<tgt::vec3> voxels;

        for(auto& v : node.voxels_) {
            voxels.push_back(v + shift);
        }

        builder.insertNode(node.pos_+ shift, std::move(voxels), node.getRadius(), node.isAtSampleBorder_, node.getUUID());
    }

    for(auto& edge : input.getEdges()) {
        builder.insertEdge(edge.getNodeID1(), edge.getNodeID2(), edge, edge.getUUID());
    }
    return std::move(builder).finalize();
}

static std::unique_ptr<VesselGraph> changeProperties(const VesselGraph& input, float amount, random_engine_type random_engine) {
    float length_sum = 0;
    for(auto& edge : input.getEdges()) {
        length_sum += edge.getLength();
    }
    float length_avg = length_sum / input.getEdges().size();

    float perturbation = 2*amount*length_avg;

    VesselGraphBuilder builder;
    for(auto& node : input.getNodes()) {
        builder.insertNode(node);
    }

    for(auto& edge : input.getEdges()) {
        auto properties = modifyEdgePathProperties(random_engine, amount, edge.getPathProperties());
        builder.insertEdge(edge.getNodeID1(), edge.getNodeID2(), properties, edge.getUUID());
    }
    return std::move(builder).finalize();
}

static void assertValidData(const VesselGraph& graph) {
    for(auto& edge : graph.getEdges()) {
        tgtAssert(edge.getLength() > 0, "Invalid edge length");
    }

    for(auto& edge : graph.getEdges()) {
        if(edge.hasValidData()) {
            return;
        }
    }
    tgtAssert(false, "No edge has valid data");
}

static void assertIDsPresent(const VesselGraph& graph, const std::set<VesselGraphEdgeUUID>& edge_uuids, const std::set<VesselGraphNodeUUID>& node_uuids) {
    for(auto& edge : graph.getEdges()) {
        if(edge_uuids.count(edge.getUUID()) > 0) {
            goto edges_valid;
        }
    }
    tgtAssert(false, "No edge known uuid");
edges_valid:

    for(auto& node : graph.getNodes()) {
        if(node_uuids.count(node.getUUID()) > 0) {
            goto nodes_valid;
        }
    }
    tgtAssert(false, "No node known uuid");
nodes_valid:
    ;
}

static std::unique_ptr<VesselGraph> perturbCombined(const VesselGraph& input, float amount, random_engine_type random_engine) {
    assertValidData(input);
    auto original_edge_uuids = get_edge_uuids(input);
    auto original_node_uuids = get_node_uuids(input);

    auto g1 = moveNodes(input, amount, random_engine);
    assertValidData(*g1);
    assertIDsPresent(*g1, original_edge_uuids, original_node_uuids);

    auto g2 = splitNodesWithUUIDs(*g1, original_node_uuids, amount, random_engine);
    assertValidData(*g2);
    assertIDsPresent(*g2, original_edge_uuids, original_node_uuids);

    auto g3 = addEdgesToNodesWithUUIDs(*g2, original_node_uuids, amount, random_engine);
    assertValidData(*g3);
    assertIDsPresent(*g3, original_edge_uuids, original_node_uuids);

    auto g4 = subdivideEdgesWithUUIDs(*g3, original_edge_uuids, amount/2 /* we have to divide the edges up 50/50 between the 4th and 5th perturbation*/, random_engine);
    assertValidData(*g4);
    assertIDsPresent(*g4, original_edge_uuids, original_node_uuids);

    auto g5 = splitEdgesWithUUIDs(*g4, original_edge_uuids, amount, random_engine);
    assertValidData(*g5);
    assertIDsPresent(*g5, original_edge_uuids, original_node_uuids);

    auto g6 = changeProperties(*g5, amount, random_engine);
    assertValidData(*g6);

    return g6;
}

void VesselGraphPerturbation::process() {
    const VesselGraph* input = inport_.getData();
    if(!input) {
        outport_.setData(nullptr);
        return;
    }
    if(!enabled_.get()) {
        outport_.setData(input, false);
        return;
    }
    random_engine_type randomEngine(randomDevice());
    if(usePredeterminedSeed_.get()) {
        randomEngine.seed(predeterminedSeed_.get());
    }

    std::unique_ptr<VesselGraph> output(nullptr);
    switch(perturbationMethod_.getValue()) {
        case ADD_EDGES: // See NetMets: software for quantifying [...] network segmentation, figure 4 (a)
            output = addEdges(*input, perturbationAmount_.get(), randomEngine);
            break;
        case SPLIT_NODES: // ... (b)
            output = splitNodes(*input, perturbationAmount_.get(), randomEngine);
            break;
        case SUBDIVIDE_EDGES: // ... (c)
            output = subdivideEdges(*input, perturbationAmount_.get(), randomEngine);
            break;
        case SPLIT_EDGES: // not shown, cut an edge in half, adding two nodes: o-------o -> o--o o--o
            output = splitEdges(*input, perturbationAmount_.get(), randomEngine);
            break;
        case MOVE_NODES:
            output = moveNodes(*input, perturbationAmount_.get(), randomEngine);
            break;
        case CHANGE_PROPERTIES:
            output = changeProperties(*input, perturbationAmount_.get(), randomEngine);
            break;
        case COMBINED:
            output = perturbCombined(*input, perturbationAmount_.get(), randomEngine);
            break;
    }
    tgtAssert(output, "No output graph generated");

    outport_.setData(output.release());
}
} // namespace voreen
