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

#include "vesselgraph.h"
#include <functional>
#include <numeric>
#include "voreen/core/voreenapplication.h"

#include "tgt/logmanager.h"

namespace voreen {

// VesselSkeletonVoxel -------------------------------------------------------------------------
VesselSkeletonVoxel::VesselSkeletonVoxel(const tgt::vec3& pos, float minDistToSurface, float maxDistToSurface, float avgDistToSurface, uint32_t numSurfaceVoxels, float volume)
    : pos_(pos)
    , minDistToSurface_(minDistToSurface)
    , maxDistToSurface_(maxDistToSurface)
    , avgDistToSurface_(avgDistToSurface)
    , numSurfaceVoxels_(numSurfaceVoxels)
    , volume_(volume)
{
}

VesselSkeletonVoxel::VesselSkeletonVoxel()
    : pos_(tgt::vec3::zero)
    , minDistToSurface_(0)
    , maxDistToSurface_(0)
    , avgDistToSurface_(0)
    , numSurfaceVoxels_(0)
    , volume_(0)
{
}

bool VesselSkeletonVoxel::hasValidData() const {
    return std::isfinite(minDistToSurface_);
}

float VesselSkeletonVoxel::roundness() const {
    if(maxDistToSurface_ != 0) {
        return minDistToSurface_ / maxDistToSurface_;
    } else {
        return 1;
    }
}
void VesselSkeletonVoxelSerializable::serialize(Serializer& s) const {
    s.serialize("pos", inner_.pos_);
    s.serialize("minDistToSurface", inner_.minDistToSurface_);
    s.serialize("maxDistToSurface", inner_.maxDistToSurface_);
    s.serialize("avgDistToSurface", inner_.avgDistToSurface_);
    s.serialize("numSurfaceVoxels", inner_.numSurfaceVoxels_);
    s.serialize("volume", inner_.volume_);
}
void VesselSkeletonVoxelSerializable::deserialize(Deserializer& s) {
    s.deserialize("pos", inner_.pos_);
    s.deserialize("minDistToSurface", inner_.minDistToSurface_);
    s.deserialize("maxDistToSurface", inner_.maxDistToSurface_);
    s.deserialize("avgDistToSurface", inner_.avgDistToSurface_);
    s.deserialize("numSurfaceVoxels", inner_.numSurfaceVoxels_);
    s.deserialize("volume", inner_.volume_);
}
VesselSkeletonVoxelSerializable::VesselSkeletonVoxelSerializable()
    : inner_()
{
}
VesselSkeletonVoxelSerializable::VesselSkeletonVoxelSerializable(VesselSkeletonVoxel val)
    : inner_(val)
{
}


// VesselGraphNode -------------------------------------------------------------------------
VesselGraphNode::VesselGraphNode(VesselGraph& graph, VGNodeID id, const tgt::vec3& position, std::vector<tgt::vec3> voxels, float radius, bool isAtSampleBorder, VesselGraphNodeUUID uuid)
    : id_(id)
    , uuid_(uuid)
    , edges_()
    , pos_(position)
    , voxels_(voxels)
    , isAtSampleBorder_(isAtSampleBorder)
    , radius_(radius)
    , graph_(&graph)
{
}
VesselGraphNode::VesselGraphNode(VesselGraphNode&& other)
    : id_(other.id_)
    , uuid_(std::move(other.uuid_))
    , edges_(std::move(other.edges_))
    , pos_(other.pos_)
    , voxels_(std::move(other.voxels_))
    , isAtSampleBorder_(other.isAtSampleBorder_)
    , radius_(other.radius_)
    , graph_(other.graph_)
{
    other.graph_ = nullptr;
}
VesselGraphNode& VesselGraphNode::operator=(VesselGraphNode&& other)
{
    this->~VesselGraphNode();
    new(this) VesselGraphNode(std::move(other));
    return *this;
}

// ONLY used for deserialization
VesselGraphNode::VesselGraphNode()
    : id_(-1)
    , uuid_()
    , edges_()
    , pos_()
    , voxels_()
    , isAtSampleBorder_(false)
    , radius_(std::numeric_limits<float>::quiet_NaN())
    , graph_(nullptr)
{
}

int VesselGraphNode::getDegree() const {
    return edges_.size();
}

bool VesselGraphNode::isEndNode() const {
    return getDegree() <= 1 && !isAtSampleBorder_;
}
std::vector<std::reference_wrapper<const VesselGraphEdge>> VesselGraphNode::getEdges() const {
    std::vector<std::reference_wrapper<const VesselGraphEdge>> edgerefs;
    for(VGEdgeID i : edges_) {
        edgerefs.push_back(graph_->getEdge(i));
    }
    return edgerefs;
}

std::vector<const VesselGraphEdge*> VesselGraphNode::getEdgesAsPtrs() const {
    std::vector<const VesselGraphEdge*> edgerefs;
    for(VGEdgeID i : edges_) {
        edgerefs.push_back(&graph_->getEdge(i));
    }
    return edgerefs;
}

std::vector<const VesselGraphNode*> VesselGraphNode::getNeighbors() const {
    std::vector<const VesselGraphNode*> neighbor_refs;
    for(VGEdgeID i : edges_) {
        auto& edge = graph_->getEdge(i);
        if(!edge.isLoop()) {
            neighbor_refs.push_back(&edge.getOtherNode(*this));
        }
    }
    return neighbor_refs;
}

VesselGraphNodeUUID VesselGraphNode::getUUID() const {
    return uuid_;
}

VGNodeID VesselGraphNode::getID() const {
    return id_;
}
float VesselGraphNode::getRadius() const {
    return radius_;
}
void VesselGraphNode::serialize(Serializer& s) const {
    s.serialize("id", id_.raw());
    std::vector<uint32_t> edges;
    for(auto edge : edges_) {
        edges.push_back(edge.raw());
    }
    s.serialize("edges", edges);
    s.serialize("pos", pos_);
    s.serialize("voxels_", voxels_);
    s.serialize("radius", radius_);
    s.serialize("isAtSampleBorder", isAtSampleBorder_);
}
void VesselGraphNode::deserialize(Deserializer& s) {
    uint32_t id;
    s.deserialize("id", id);
    id_ = id;
    std::vector<uint32_t> edges;
    s.deserialize("edges", edges);
    for(auto edge : edges) {
        edges_.push_back(edge);
    }
    s.deserialize("pos", pos_);
    s.deserialize("voxels_", voxels_);
    s.deserialize("radius", radius_);
    s.deserialize("isAtSampleBorder", isAtSampleBorder_);
}

// VesselGraphEdge -------------------------------------------------------------------------
template<class T, class E>
static void statisticalAnalysis(const T& voxels, std::function<float(const E&)> getVoxelValue, float& /*out*/ avg, float& /*out*/ stddev) {

    float sum = 0;
    int numValidVoxels = 0;
    for(auto& voxel : voxels) {
        if(voxel.hasValidData()) {
            float val = getVoxelValue(voxel);
            sum += val;
            ++numValidVoxels;
        }
    }
    if(numValidVoxels == 0) {
        return;
    }
    avg = sum/numValidVoxels;
    tgtAssert(std::isnormal(avg) || avg == 0, "Unnormal avg");

    float varSum = 0;
    for(auto& voxel : voxels) {
        if(voxel.hasValidData()) {
            float val = getVoxelValue(voxel) - avg;
            varSum += val*val;
        }
    }
    float variance = varSum/numValidVoxels;
    tgtAssert(variance >= 0, "negative variance");
    stddev = std::sqrt(variance);
    tgtAssert(std::isnormal(stddev) || stddev == 0, "Unnormal stddev");

    tgtAssert(!std::isnan(stddev), "Unnormal stddev");
    tgtAssert(!std::isnan(avg), "Unnormal avg");
}

const float VesselGraphEdgePathProperties::INVALID_DATA = -1;
VesselGraphEdgePathProperties::VesselGraphEdgePathProperties()
    : length_(INVALID_DATA)
    , volume_(INVALID_DATA)
    , minRadiusAvg_(INVALID_DATA)
    , minRadiusStdDeviation_(INVALID_DATA)
    , maxRadiusAvg_(INVALID_DATA)
    , maxRadiusStdDeviation_(INVALID_DATA)
    , maxRadiusMax_(INVALID_DATA)
    , avgRadiusAvg_(INVALID_DATA)
    , avgRadiusStdDeviation_(INVALID_DATA)
    , roundnessAvg_(INVALID_DATA)
    , roundnessStdDeviation_(INVALID_DATA)
{
}

bool VesselGraphEdgePathProperties::hasValidData() const {
    bool isValid = minRadiusAvg_ != INVALID_DATA;
    // If minRadiusAvg_ is invalid, all other skeleton voxel based properties should be invalid, as well.
    tgtAssert(isValid ^ (minRadiusStdDeviation_ == INVALID_DATA), "Partially invalid data");
    tgtAssert(isValid ^ (avgRadiusAvg_ == INVALID_DATA), "Partially invalid data");
    tgtAssert(isValid ^ (avgRadiusStdDeviation_ == INVALID_DATA), "Partially invalid data");
    tgtAssert(isValid ^ (maxRadiusAvg_ == INVALID_DATA), "Partially invalid data");
    tgtAssert(isValid ^ (maxRadiusStdDeviation_ == INVALID_DATA), "Partially invalid data");
    tgtAssert(isValid ^ (maxRadiusMax_ == INVALID_DATA), "Partially invalid data");
    tgtAssert(isValid ^ (roundnessAvg_ == INVALID_DATA), "Partially invalid data");
    tgtAssert(isValid ^ (roundnessStdDeviation_ == INVALID_DATA), "Partially invalid data");
    return isValid;
}


VesselGraphEdgePathProperties VesselGraphEdgePathProperties::fromPath(const VesselGraphNode& begin, const VesselGraphNode& end, const DiskArray<VesselSkeletonVoxel>& path) {
    VesselGraphEdgePathProperties output;

    // Compute length
    output.length_ = 0;
    if(path.empty()) {
        output.length_ = tgt::distance(begin.pos_, end.pos_);
    } else {
        for(size_t i=0; i < path.size()-1; ++i) {
            output.length_ += tgt::distance(path[i].pos_, path[i+1].pos_);
        }
        output.length_ += tgt::distance(begin.pos_, path.front().pos_);
        output.length_ += tgt::distance(path.back().pos_, end.pos_);
    }

    // Compute volume
    if(path.size() > 0) {
        output.volume_ = 0;
        for(const VesselSkeletonVoxel& voxel : path) {
            output.volume_ += voxel.volume_;
        }
    }

    // Compute min radius vals
    statisticalAnalysis<DiskArray<VesselSkeletonVoxel>, VesselSkeletonVoxel>(path, [] (const VesselSkeletonVoxel& v) {
            return v.minDistToSurface_;
            },
            output.minRadiusAvg_, output.minRadiusStdDeviation_);

    // Compute avg radius vals
    statisticalAnalysis<DiskArray<VesselSkeletonVoxel>, VesselSkeletonVoxel>(path, [] (const VesselSkeletonVoxel& v) {
            return v.avgDistToSurface_;
            },
            output.avgRadiusAvg_, output.avgRadiusStdDeviation_);

    // Compute max radius vals
    statisticalAnalysis<DiskArray<VesselSkeletonVoxel>, VesselSkeletonVoxel>(path, [] (const VesselSkeletonVoxel& v) {
            return v.maxDistToSurface_;
            },
            output.maxRadiusAvg_, output.maxRadiusStdDeviation_);
    auto maybe_max = std::max_element(path.begin(), path.end(), [](const VesselSkeletonVoxel& v1, const VesselSkeletonVoxel& v2) {
            return v1.maxDistToSurface_ < v2.maxDistToSurface_;
            });
    output.maxRadiusMax_ = maybe_max!=path.end() && maybe_max->hasValidData() ? maybe_max->maxDistToSurface_ : INVALID_DATA;

    // Compute roundness vals
    statisticalAnalysis<DiskArray<VesselSkeletonVoxel>, VesselSkeletonVoxel>(path, [] (const VesselSkeletonVoxel& v) {
            return v.roundness();
            },
            output.roundnessAvg_, output.roundnessStdDeviation_);

    tgtAssert(!std::isnan(output.length_), "Invalid length");
    return output;
}
void VesselGraphEdgePathProperties::serialize(Serializer& s) const {
    s.serialize("length", length_);
    s.serialize("volume", volume_);
    s.serialize("minRadiusAvg", minRadiusAvg_);
    s.serialize("minRadiusStdDeviation", minRadiusStdDeviation_);
    s.serialize("maxRadiusAvg", maxRadiusAvg_);
    s.serialize("maxRadiusStdDeviation", maxRadiusStdDeviation_);
    s.serialize("avgRadiusAvg", avgRadiusAvg_);
    s.serialize("avgRadiusStdDeviation", avgRadiusStdDeviation_);
    s.serialize("roundnessAvg", roundnessAvg_);
    s.serialize("roundnessStdDeviation", roundnessStdDeviation_);
}
void VesselGraphEdgePathProperties::deserialize(Deserializer& s) {
    s.deserialize("length", length_);
    s.deserialize("volume", volume_);
    s.deserialize("minRadiusAvg", minRadiusAvg_);
    s.deserialize("minRadiusStdDeviation", minRadiusStdDeviation_);
    s.deserialize("maxRadiusAvg", maxRadiusAvg_);
    s.deserialize("maxRadiusStdDeviation", maxRadiusStdDeviation_);
    s.deserialize("avgRadiusAvg", avgRadiusAvg_);
    s.deserialize("avgRadiusStdDeviation", avgRadiusStdDeviation_);
    s.deserialize("roundnessAvg", roundnessAvg_);
    s.deserialize("roundnessStdDeviation", roundnessStdDeviation_);
}

VesselGraphEdge::VesselGraphEdge(VesselGraph& graph, VGEdgeID id, VGNodeID node1ID, VGNodeID node2ID, DiskArray<VesselSkeletonVoxel>&& voxels, VesselGraphEdgeUUID uuid)
    : graph_(&graph)
    , id_(id)
    , node1_(node1ID)
    , node2_(node2ID)
    , distance_(std::numeric_limits<float>::quiet_NaN())
    , pathProps_() //Invalid until initialized
    , voxels_(std::move(voxels))
    , uuid_(uuid)
{
    VesselGraphNode& node1 = getNode1();
    VesselGraphNode& node2 = getNode2();
    // Compute distance
    distance_ = tgt::distance(node1.pos_, node2.pos_);

    // Compute path properties
    pathProps_ = VesselGraphEdgePathProperties::fromPath(node1, node2, voxels_);

    tgtAssert(!std::isnan(getLength()), "Invalid length");
}

VesselGraphEdge::VesselGraphEdge(VesselGraph& graph, VGEdgeID id, VGNodeID node1ID, VGNodeID node2ID, VesselGraphEdgePathProperties pathProps, VesselGraphEdgeUUID uuid)
    : graph_(&graph)
    , id_(id)
    , node1_(node1ID)
    , node2_(node2ID)
    , distance_(std::numeric_limits<float>::quiet_NaN())
    , pathProps_(pathProps)
    , voxels_()
    , uuid_(uuid)
{
    VesselGraphNode& node1 = getNode1();
    VesselGraphNode& node2 = getNode2();
    // Compute distance
    distance_ = tgt::distance(node1.pos_, node2.pos_);

    tgtAssert(!std::isnan(getLength()), "Invalid length");
}

VesselGraphEdge::VesselGraphEdge(VesselGraphEdge&& other)
    : graph_(other.graph_)
    , id_(other.id_)
    , node1_(other.node1_)
    , node2_(other.node2_)
    , distance_(other.distance_)
    , pathProps_(other.pathProps_)
    , voxels_(std::move(other.voxels_))
    , uuid_(std::move(other.uuid_))
{
}

void VesselGraphEdge::operator=(VesselGraphEdge&& other)
{
    graph_ = other.graph_;
    id_ = other.id_;
    node1_ = other.node1_;
    node2_ = other.node2_;
    distance_ = other.distance_;
    pathProps_ = other.pathProps_;
    voxels_ = std::move(other.voxels_);

    other.graph_ = nullptr;
}

VesselGraphEdge::VesselGraphEdge()
    : graph_(nullptr)
    , id_(-1)
    , node1_(-1)
    , node2_(-1)
    , distance_(-1)
    , pathProps_()
    , voxels_()
    , uuid_()
{
}

VesselGraphEdgeUUID VesselGraphEdge::getUUID() const {
    return uuid_;
}

const VesselGraphNode& VesselGraphEdge::getNode1() const {
    return graph_->getNode(node1_);
}
const VesselGraphNode& VesselGraphEdge::getNode2() const {
    return graph_->getNode(node2_);
}
const VesselGraphNode& VesselGraphEdge::getOtherNode(const VesselGraphNode& firstNode) const {
    if(firstNode.getID() == node1_) {
        // firstNode is node1
        tgtAssert(&firstNode == &getNode1(), "Passed invalid node to getOtherNode");

        return getNode2();
    } else {
        // firstNode is node2
        tgtAssert(&firstNode == &getNode2(), "Passed invalid node to getOtherNode");

        return getNode1();
    }
}
VesselGraphNode& VesselGraphEdge::getNode1() {
    return graph_->getNode(node1_);
}
VesselGraphNode& VesselGraphEdge::getNode2() {
    return graph_->getNode(node2_);
}
const VesselGraphEdgePathProperties& VesselGraphEdge::getPathProperties() const {
    return pathProps_;
}

VGNodeID VesselGraphEdge::getNodeID1() const {
    return node1_;
}
VGNodeID VesselGraphEdge::getNodeID2() const {
    return node2_;
}

bool VesselGraphEdge::hasValidData() const {
    tgtAssert(distance_ != VesselGraphEdgePathProperties::INVALID_DATA, "Invalid distance, should not happen.");
    return pathProps_.hasValidData();
}

float VesselGraphEdge::getLength() const {
    return pathProps_.length_;
}

float VesselGraphEdge::getDistance() const {
    return distance_;
}

float VesselGraphEdge::getCurveness() const {
    if(getDistance() == 0) {
        // This can happen in loops... Symbolicly set it to curveness to 0.
        return 0.0f;
    } else {
        float curveness =  getLength()/getDistance();
        return curveness;
    }
}
float VesselGraphEdge::getStraightness() const {
    if(getLength() == 0) {
        // This should really not happen
        tgtAssert(false, "Length = 0");
        return 0.0f;
    } else {
        return getDistance()/getLength();
    }
}

float VesselGraphEdge::getVolume() const {
    return pathProps_.volume_;
}

float VesselGraphEdge::getAvgCrossSection() const {
    tgtAssert(pathProps_.length_ > 0, "Invalid length");
    return pathProps_.volume_/pathProps_.length_;
}

float VesselGraphEdge::getMinRadiusAvg() const {
    return pathProps_.minRadiusAvg_;
}

float VesselGraphEdge::getMinRadiusStdDeviation() const {
    return pathProps_.minRadiusStdDeviation_;
}

float VesselGraphEdge::getAvgRadiusAvg() const {
    return pathProps_.avgRadiusAvg_;
}

float VesselGraphEdge::getAvgRadiusStdDeviation() const {
    return pathProps_.avgRadiusStdDeviation_;
}

float VesselGraphEdge::getMaxRadiusAvg() const {
    return pathProps_.maxRadiusAvg_;
}

float VesselGraphEdge::getMaxRadiusMax() const {
    return pathProps_.maxRadiusMax_;
}

float VesselGraphEdge::getMaxRadiusStdDeviation() const {
    return pathProps_.maxRadiusStdDeviation_;
}

float VesselGraphEdge::getRoundnessAvg() const {
    return pathProps_.roundnessAvg_;
}

float VesselGraphEdge::getRoundnessStdDeviation() const {
    return pathProps_.roundnessStdDeviation_;
}

float VesselGraphEdge::getElongation() const {
    float diameter = 2*pathProps_.avgRadiusAvg_;
    if(diameter > 0) {
        return pathProps_.length_ / diameter;
    } else {
        //tgtAssert(false, "Invalid radius");
        return 1; // Return neutral elongation
    }
}

float VesselGraphEdge::getEffectiveLength() const {
    return std::max(0.0f, getLength() - getNode1().getRadius() - getNode2().getRadius());
}

float VesselGraphEdge::getRelativeBulgeSize() const {
    //float edge_length_contributed_length = std::max(getNode1().estimatedRadius(), getNode2().estimatedRadius());
    //return std::max(0.0f,(getLength() / edge_length_contributed_length) - 1);
    if(hasValidData()) {
        if(getNode1().isEndNode() == getNode2().isEndNode()) {
            return -1;
        }
        //const VesselGraphNode* endNode;
        const VesselGraphNode* branchNode;
        if(getNode1().isEndNode()) {
            //endNode = &getNode1();
            branchNode = &getNode2();
        } else {
            //endNode = &getNode2();
            branchNode = &getNode1();
        }
        float parentVesselRadius = getMaxRadiusMax();
        if(branchNode->getRadius() > 0) {
            parentVesselRadius = std::min(parentVesselRadius, branchNode->getRadius());
        }
        tgtAssert(parentVesselRadius > 0, "Invalid parent vessel radius");
        return getLength() / parentVesselRadius;
    } else {
        return -1;
    }
}

bool VesselGraphEdge::isLoop() const {
    return getNodeID1() == getNodeID2();
}

const DiskArray<VesselSkeletonVoxel>& VesselGraphEdge::getVoxels() const {
    return voxels_;
}
VGEdgeID VesselGraphEdge::getID() const {
    return id_;
}
bool VesselGraphEdge::isEndStanding() const {
    const VesselGraphNode& node1 = getNode1();
    const VesselGraphNode& node2 = getNode2();
    return !node1.isAtSampleBorder_ && !node2.isAtSampleBorder_ && node1.isEndNode() != node2.isEndNode();
}
size_t VesselGraphEdge::getNumValidVoxels() const {
    return std::count_if(voxels_.begin(), voxels_.end(), [] (const VesselSkeletonVoxel& v) {
            return v.hasValidData();
            });
}

void VesselGraphEdge::serialize(Serializer& s) const {
    s.serialize("id", id_.raw());
    s.serialize("node1", node1_.raw());
    s.serialize("node2", node2_.raw());
    s.serialize("distance", distance_);

    if(voxels_.empty()) {
        s.serialize("pathProperties", pathProps_);
    } else {
        std::vector<VesselSkeletonVoxelSerializable> voxels;
        voxels.reserve(voxels.size());
        for(const auto& voxel : voxels_) {
            voxels.emplace_back(voxel);
        }
        s.serialize("skeletonVoxels", voxels);
    }
}
void VesselGraphEdge::deserialize(Deserializer& s) {
    uint32_t id, n1, n2;
    s.deserialize("id", id);
    s.deserialize("node1", n1);
    s.deserialize("node2", n2);
    id_ = id;
    node1_ = n1;
    node2_ = n2;
    s.deserialize("distance", distance_);

    bool noSkeletonVoxelsTag = false;
    try {
        std::vector<VesselSkeletonVoxelSerializable> voxels(voxels_.begin(), voxels_.end());
        s.deserialize("skeletonVoxels", voxels);
        auto builder = graph_->voxelStorage_->build();
        for(const auto& voxel : voxels) {
            builder.push(voxel.inner_);
        }
        voxels_ = std::move(builder).finalize();
    } catch (SerializationException s) {
        noSkeletonVoxelsTag = true;
    }
    if(noSkeletonVoxelsTag || voxels_.empty()) {
        s.deserialize("pathProperties", pathProps_);
    }
}

void VesselGraphEdge::updatePathPropertiesFromVoxels() {
    // NOTE:
    // Yes, this is pretty ugly, but we cannot set the graph pointer in the deserialization and we cannot
    // therefore get the nodes in the deserialization. This method will therefore be called after the
    // "normal" deserialization of edges, and only after all nodes have been added to the graph.
    //
    if(!voxels_.empty()) {
        pathProps_ = VesselGraphEdgePathProperties::fromPath(getNode1(), getNode2(), voxels_);
    }
}


// VesselGraph -------------------------------------------------------------------------

VesselGraph::VesselGraph(const tgt::Bounds& bounds)
    : nodes_()
    , edges_()
    , bounds_(bounds)
    , voxelStorage_(new DiskArrayStorage<VesselSkeletonVoxel>(VoreenApplication::app()->getUniqueTmpFilePath(".vgvoxels")))
{
}

VesselGraph::VesselGraph()
    : nodes_()
    , edges_()
    , bounds_()
    , voxelStorage_(new DiskArrayStorage<VesselSkeletonVoxel>(VoreenApplication::app()->getUniqueTmpFilePath(".vgvoxels")))
{
}

VesselGraph::VesselGraph(VesselGraph&& other)
    : nodes_(std::move(other.nodes_))
    , edges_(std::move(other.edges_))
    , bounds_(other.bounds_)
    , voxelStorage_(std::move(other.voxelStorage_))
{
    // Fix up references:
    for(auto& node: nodes_) {
        node.graph_ = this;
    }

    for(auto& edge: edges_) {
        edge.graph_ = this;
    }
}

VesselGraph VesselGraph::clone() const {

    VesselGraph res(getBounds());
    for(auto& node: nodes_) {
        std::vector<tgt::vec3> voxels(node.voxels_);
        res.insertNode(node.pos_, std::move(voxels), node.getRadius(), node.isAtSampleBorder_, node.getUUID());
    }

    for(auto& edge: edges_) {
        res.insertEdge(edge.getNodeID1(), edge.getNodeID2(), edge, edge.getUUID());
    }
    return res;
}

const VesselGraphNode& VesselGraph::getNode(VGNodeID i) const {
    return nodes_.at(i.raw());
}

const VesselGraphEdge& VesselGraph::getEdge(VGEdgeID i) const {
    return edges_.at(i.raw());
}

VesselGraphNode& VesselGraph::getNode(VGNodeID i) {
    return nodes_.at(i.raw());
}

VesselGraphEdge& VesselGraph::getEdge(VGEdgeID i) {
    return edges_.at(i.raw());
}
VGNodeID VesselGraph::insertNode(const VesselGraphNode& base) {
    std::vector<tgt::vec3> new_voxels(base.voxels_);
    return insertNode(base.pos_, std::move(new_voxels), base.getRadius(), base.isAtSampleBorder_, base.getUUID());
}

VGNodeID VesselGraph::insertNode(const tgt::vec3& position, const std::vector<tgt::vec3>&& voxels, float radius, bool isAtSampleBorder) {
    return insertNode(position, std::move(voxels), radius, isAtSampleBorder, VoreenApplication::app()->generateUUID());
}
VGNodeID VesselGraph::insertNode(const tgt::vec3& position, const std::vector<tgt::vec3>&& voxels, float radius, bool isAtSampleBorder, VesselGraphNodeUUID uuid) {
    size_t edgeID = nodes_.size();
    bounds_.addPoint(position);
    //TODO: add voxels as well?
    nodes_.emplace_back(*this, edgeID, position, std::move(voxels), radius, isAtSampleBorder, uuid);
    return edgeID;
}


VGEdgeID VesselGraph::insertEdge(VGNodeID node1, VGNodeID node2, const DiskArray<VesselSkeletonVoxel>& voxels) {
    return insertEdge(node1, node2, voxels, VoreenApplication::app()->generateUUID());
}
VGEdgeID VesselGraph::insertEdge(VGNodeID node1, VGNodeID node2, const DiskArray<VesselSkeletonVoxel>& voxels, VesselGraphEdgeUUID uuid) {
    tgtAssert(node1 < nodes_.size(), "Edge references nonexistent node");
    tgtAssert(node2 < nodes_.size(), "Edge references nonexistent node");
    VesselGraphNode& n1 = nodes_.at(node1.raw());
    VesselGraphNode& n2 = nodes_.at(node2.raw());

    VGEdgeID edgeID = edges_.size();
    edges_.emplace_back(*this, edgeID, node1, node2, voxelStorage_->store(voxels), uuid);

    n1.edges_.push_back(edgeID);
    n2.edges_.push_back(edgeID);
    return edgeID;
}
VGEdgeID VesselGraph::insertEdge(VGNodeID node1, VGNodeID node2, const std::vector<VesselSkeletonVoxel>& voxels) {
    return insertEdge(node1, node2, voxels, VoreenApplication::app()->generateUUID());
}
VGEdgeID VesselGraph::insertEdge(VGNodeID node1, VGNodeID node2, const std::vector<VesselSkeletonVoxel>& voxels, VesselGraphEdgeUUID uuid) {
    tgtAssert(node1 < nodes_.size(), "Edge references nonexistent node");
    tgtAssert(node2 < nodes_.size(), "Edge references nonexistent node");
    VesselGraphNode& n1 = nodes_.at(node1.raw());
    VesselGraphNode& n2 = nodes_.at(node2.raw());

    VGEdgeID edgeID = edges_.size();
    edges_.emplace_back(*this, edgeID, node1, node2, voxelStorage_->store(voxels), uuid);

    n1.edges_.push_back(edgeID);
    n2.edges_.push_back(edgeID);
    return edgeID;
}

VGEdgeID VesselGraph::insertEdge(VGNodeID node1, VGNodeID node2, const VesselGraphEdge& path_definition, VesselGraphEdgeUUID uuid) {
    return insertEdge(node1, node2, path_definition.voxels_, uuid);
}

VGEdgeID VesselGraph::insertEdge(VGNodeID node1, VGNodeID node2, const VesselGraphEdge& path_definition) {
    return insertEdge(node1, node2, path_definition, VoreenApplication::app()->generateUUID());
}

VGEdgeID VesselGraph::insertEdge(VGNodeID node1, VGNodeID node2, VesselGraphEdgePathProperties pathProperties) {
    return insertEdge(node1, node2, pathProperties, VoreenApplication::app()->generateUUID());
}
VGEdgeID VesselGraph::insertEdge(VGNodeID node1, VGNodeID node2, VesselGraphEdgePathProperties pathProperties, VesselGraphEdgeUUID uuid) {
    tgtAssert(node1 < nodes_.size(), "Edge references nonexistent node");
    tgtAssert(node2 < nodes_.size(), "Edge references nonexistent node");
    VesselGraphNode& n1 = nodes_.at(node1.raw());
    VesselGraphNode& n2 = nodes_.at(node2.raw());

    VGEdgeID edgeID = edges_.size();
    edges_.emplace_back(*this, edgeID, node1, node2, pathProperties, uuid);

    n1.edges_.push_back(edgeID);
    n2.edges_.push_back(edgeID);
    return edgeID;
}

const std::vector<VesselGraphNode>& VesselGraph::getNodes() const {
    return nodes_;
}

const std::vector<VesselGraphEdge>& VesselGraph::getEdges() const {
    return edges_;
}

std::vector<VesselGraphNode>& VesselGraph::getNodes() {
    return nodes_;
}

std::vector<VesselGraphEdge>& VesselGraph::getEdges() {
    return edges_;
}

void VesselGraph::getEdgePropertyStats(std::function<float(const VesselGraphEdge&)> f, float& /*out*/ mean, float& /*out*/stddev) const {
    statisticalAnalysis(edges_, f, mean, stddev);
}

const tgt::Bounds& VesselGraph::getBounds() const {
    return bounds_;
}

void VesselGraph::serialize(Serializer& s) const {
    s.serialize("nodes", nodes_);
    s.serialize("edges", edges_);
    s.serialize("bounds", bounds_);
}
void VesselGraph::deserialize(Deserializer& s) {
    // Yeah, this is quite ugly... But we cannot set the graph pointer in Serializable::deserialize()
    s.deserialize("nodes", nodes_);
    for(auto& node : nodes_) {
        node.graph_ = this;
    }
    s.deserialize("edges", edges_);
    for(auto& edge : edges_) {
        edge.graph_ = this;
        edge.updatePathPropertiesFromVoxels();
    }
    s.deserialize("bounds", bounds_);
}

} // namespace voreen
