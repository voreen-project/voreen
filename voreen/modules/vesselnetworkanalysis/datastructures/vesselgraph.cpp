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

#include "vesselgraph.h"
#include <functional>
#include <numeric>

#include "voreen/core/voreenapplication.h"

#include "tgt/logmanager.h"

namespace voreen {

// VesselSkeletonVoxel -------------------------------------------------------------------------
VesselSkeletonVoxel::VesselSkeletonVoxel(const tgt::vec3& pos, float minDistToSurface, float maxDistToSurface, float avgDistToSurface, uint32_t numSurfaceVoxels, float volume, bool nearOtherEdge)
    : pos_(pos)
    , minDistToSurface_(minDistToSurface)
    , maxDistToSurface_(maxDistToSurface)
    , avgDistToSurface_(avgDistToSurface)
    , numSurfaceVoxels_(numSurfaceVoxels)
    , volume_(volume)
    , nearOtherEdge_(nearOtherEdge)
{
    tgtAssert(minDistToSurface_ >= 0, "invalid mindist");
    tgtAssert(avgDistToSurface_ >= 0, "invalid avgdist");
    tgtAssert(maxDistToSurface_ >= 0, "invalid maxdist");
    tgtAssert(volume_ >= 0, "invalid volume");
}

VesselSkeletonVoxel::VesselSkeletonVoxel()
    : pos_(tgt::vec3::zero)
    , minDistToSurface_(std::numeric_limits<float>::infinity())
    , maxDistToSurface_(0)
    , avgDistToSurface_(0)
    , numSurfaceVoxels_(0)
    , volume_(0)
    , nearOtherEdge_(false)
{
}

bool VesselSkeletonVoxel::hasValidData() const {
    return numSurfaceVoxels_ > 0;
}

bool VesselSkeletonVoxel::isInner() const {
    return nearOtherEdge_;
}

bool VesselSkeletonVoxel::isOuter() const {
    return !isInner();
}

float VesselSkeletonVoxel::roundness() const {
    if(maxDistToSurface_ != 0) {
        return minDistToSurface_ / maxDistToSurface_;
    } else {
        return 1;
    }
}

template<typename V>
static void serializeJsonValue(SerializationOutputStream& stream, const V& value) {
    stream << value;
}
template<>
void serializeJsonValue(SerializationOutputStream& stream, const tgt::vec3& value) {
    stream << "[" << value.x << "," << value.y << "," << value.z << "]";
}
template<>
void serializeJsonValue(SerializationOutputStream& stream, const bool& value) {
    if(value) {
        stream << "true";
    } else {
        stream << "false";
    }
}
template<>
void serializeJsonValue(SerializationOutputStream& stream, const VGEdgeID& value) {
    serializeJsonValue(stream, value.raw());
}
template<>
void serializeJsonValue(SerializationOutputStream& stream, const VGNodeID& value) {
    serializeJsonValue(stream, value.raw());
}
template<>
void serializeJsonValue(SerializationOutputStream& stream, const VesselSkeletonVoxel& value) {
    value.serializeToJson(stream);
}
template<>
void serializeJsonValue(SerializationOutputStream& stream, const VesselGraphNode& value) {
    value.serializeToJson(stream);
}
void serializeJsonValue(SerializationOutputStream& stream, const VesselGraphEdge& value) {
    value.serializeToJson(stream);
}
void serializeJsonValue(SerializationOutputStream& stream, const VesselGraphEdgePathProperties& value) {
    value.serializeToJson(stream);
}
void serializeJsonValue(SerializationOutputStream& stream, const float& value) {
    stream.precision(std::numeric_limits<float>::max_digits10);
    if(std::isnan(value)) {
        stream << "NaN";
    } else if (!std::isfinite(value)) {
        if(value < 0) {
            stream << "-Infinity";
        } else {
            stream << "Infinity";
        }
    } else {
        stream << value;
    }
}
void serializeJsonValue(SerializationOutputStream& stream, const double& value) {
    stream.precision(std::numeric_limits<double>::max_digits10);
    if(std::isnan(value)) {
        stream << "NaN";
    } else if (!std::isfinite(value)) {
        if(value < 0) {
            stream << "-Infinity";
        } else {
            stream << "Infinity";
        }
    } else {
        stream << value;
    }
}

static void serializeJsonMemberName(SerializationOutputStream& stream, const std::string& name) {
    stream << "\"" << name << "\":";
}
template<typename V>
static void serializeJsonMember(SerializationOutputStream& stream, const std::string& name, const V& value) {
    serializeJsonMemberName(stream, name);
    serializeJsonValue(stream, value);
}

template<typename V>
static void serializeJsonMemberIterable(SerializationOutputStream& s, const std::string& name, const V& iterable) {
    s << "\"" << name << "\":[";
    bool first = true;
    for(const auto& elm : iterable) {
        if(first) {
            first = false;
        } else {
            s << ",";
        }
        serializeJsonValue(s, elm);
    }
    s << "]";
}

void VesselSkeletonVoxel::serializeToJson(SerializationOutputStream& s) const {
    s << "{";
    serializeJsonMember(s, "pos", pos_); s << ",";
    serializeJsonMember(s, "minDistToSurface", minDistToSurface_); s << ",";
    serializeJsonMember(s, "maxDistToSurface", maxDistToSurface_); s << ",";
    serializeJsonMember(s, "avgDistToSurface", avgDistToSurface_); s << ",";
    serializeJsonMember(s, "numSurfaceVoxels", numSurfaceVoxels_); s << ",";
    serializeJsonMember(s, "volume", volume_); s << ",";
    serializeJsonMember(s, "nearOtherEdge", nearOtherEdge_);
    s << "}";
}

void VesselSkeletonVoxelSerializable::serialize(Serializer& s) const {
    s.serialize("pos", inner_.pos_);
    s.serialize("minDistToSurface", inner_.minDistToSurface_);
    s.serialize("maxDistToSurface", inner_.maxDistToSurface_);
    s.serialize("avgDistToSurface", inner_.avgDistToSurface_);
    s.serialize("numSurfaceVoxels", inner_.numSurfaceVoxels_);
    s.serialize("volume", inner_.volume_);
    s.serialize("nearOtherEdge", inner_.nearOtherEdge_);
}
void VesselSkeletonVoxelSerializable::deserialize(Deserializer& s) {
    s.deserialize("pos", inner_.pos_);
    s.deserialize("minDistToSurface", inner_.minDistToSurface_);
    s.deserialize("maxDistToSurface", inner_.maxDistToSurface_);
    s.deserialize("avgDistToSurface", inner_.avgDistToSurface_);
    s.deserialize("numSurfaceVoxels", inner_.numSurfaceVoxels_);
    s.deserialize("volume", inner_.volume_);
    s.deserialize("nearOtherEdge", inner_.nearOtherEdge_);
}
VesselSkeletonVoxelSerializable::VesselSkeletonVoxelSerializable()
    : inner_()
{
}
VesselSkeletonVoxelSerializable::VesselSkeletonVoxelSerializable(VesselSkeletonVoxel val)
    : inner_(val)
{
}

// VG*IDs ----------------------------------------------------------------------------------
const VGNodeID VGNodeID::INVALID = VGNodeID(-1);
const VGEdgeID VGEdgeID::INVALID = VGEdgeID(-1);


// VesselGraphNode -------------------------------------------------------------------------
VesselGraphNode::VesselGraphNode(VesselGraph& graph, VGNodeID id, const tgt::vec3& position, DiskArray<tgt::vec3>&& voxels, float radius, bool isAtSampleBorder, VesselGraphNodeUUID uuid)
    : id_(id)
    , uuid_(uuid)
    , pos_(position)
    , voxels_(std::move(voxels))
    , isAtSampleBorder_(isAtSampleBorder)
    , radius_(radius)
    , graph_(&graph)
    , edges_(*graph.nodeEdgeIdStorage_)
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

void VesselGraphNode::serializeToJson(SerializationOutputStream& s) const {
    s << "{";
    serializeJsonMember(s, "id", id_.raw()); s << ",";
    serializeJsonMemberIterable(s, "edges", edges_); s << ",";
    serializeJsonMember(s, "pos", pos_); s << ",";

    serializeJsonMemberIterable(s, "voxels_", voxels_); s << ",";
    serializeJsonMember(s, "radius", radius_); s << ",";
    serializeJsonMember(s, "isAtSampleBorder", isAtSampleBorder_);
    s << "}";
}

VesselGraphNodeSerializable::VesselGraphNodeSerializable(const VesselGraphNode& node)
    : inner_(node)
{
}

void VesselGraphNodeSerializable::serialize(Serializer& s) const {
    s.serialize("id", inner_.id_.raw());
    std::vector<uint32_t> edges;
    for(auto edge : inner_.edges_) {
        edges.push_back(edge.raw());
    }
    s.serialize("edges", edges);
    s.serialize("pos", inner_.pos_);

    std::vector<tgt::vec3> voxels;
    for(auto voxel : inner_.voxels_) {
        voxels.push_back(voxel);
    }
    s.serialize("voxels_", voxels);
    s.serialize("radius", inner_.radius_);
    s.serialize("isAtSampleBorder", inner_.isAtSampleBorder_);
}
void VesselGraphNodeSerializable::deserialize(Deserializer&) {
    tgtAssert(false, "Cannot deserialize VesselGraphNodeSerializable");
}
void VesselGraphNodeDeserializable::serialize(Serializer&) const {
    tgtAssert(false, "Cannot serialize VesselGraphNodeDeserializable");
}
void VesselGraphNodeDeserializable::deserialize(Deserializer& s) {
    uint32_t id;
    s.deserialize("id", id);
    id_ = id;
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
            tgtAssert(val >= 0, "Invalid val");
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
    , innerLengthNode1_(INVALID_DATA)
    , innerLengthNode2_(INVALID_DATA)
    , tipRadiusNode1_(INVALID_DATA)
    , tipRadiusNode2_(INVALID_DATA)
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


VesselGraphEdgePathProperties VesselGraphEdgePathProperties::fromPath(const VesselGraphNode& begin, const VesselGraphNode& end, const DiskArray<VesselSkeletonVoxel>& path, size_t outerPathBeginIndex, size_t outerPathEndIndex) {
    VesselGraphEdgePathProperties output;

    // Compute length
    output.length_ = 0;
    float distFromLastOuter;
    float distToFirstOuter;

    if(path.empty()) {
        float dist = tgt::distance(begin.pos_, end.pos_);
        output.length_ = dist;
        distFromLastOuter = 0.0;
        distToFirstOuter = 0.0;

        output.tipRadiusNode1_ = 0.0;
        output.tipRadiusNode2_ = 0.0;
    } else {
        float beginDist = tgt::distance(begin.pos_, path.front().pos_);
        output.length_ += beginDist;
        distFromLastOuter = 0.0;
        if(outerPathBeginIndex != 0) {
            // If the outer path does not begin with an outer voxel, we assume the node itself to be outside
            distToFirstOuter = beginDist;
        } else {
            // Otherwise it must be inside
            distToFirstOuter = 0.0;
        }

        if(outerPathEndIndex == 0) {
            // Path is completely inside, including first node
            distFromLastOuter += beginDist;
        }

        for(size_t i=0; i < path.size()-1; ++i) {
            float dist = tgt::distance(path[i].pos_, path[i+1].pos_);
            output.length_ += dist;

            if(i < outerPathBeginIndex) {
                distToFirstOuter += dist;
            }
            if(i >= outerPathEndIndex) {
                distFromLastOuter += dist;
            }
        }

        float endDist = tgt::distance(path.back().pos_, end.pos_);
        output.length_ += endDist;
        if(outerPathEndIndex < path.size()) {
            // If the path does not end with an outer voxel, it is assumed to be inside.
            distFromLastOuter += endDist;
        }
        if(outerPathBeginIndex == path.size()) {
            // Path is completely inside, including first node
            distToFirstOuter += endDist;
        }

        auto& frontVoxel = path.front();
        auto& backVoxel = path.back();
        output.tipRadiusNode1_ = frontVoxel.hasValidData() ? frontVoxel.minDistToSurface_ : 0.0f;
        output.tipRadiusNode2_ = backVoxel.hasValidData() ? backVoxel.minDistToSurface_ : 0.0f;
    }

    output.innerLengthNode1_ = distToFirstOuter;
    output.innerLengthNode2_ = distFromLastOuter;

    // Compute volume
    if(path.size() > 0) {
        output.volume_ = 0;
        for(const VesselSkeletonVoxel& voxel : path) {
            output.volume_ += voxel.volume_;
        }
    }

    //auto outerPath = path.slice(outerPathBeginIndex, outerPathEndIndex);

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
void VesselGraphEdgePathProperties::serializeToJson(SerializationOutputStream& s) const {
    s << "{";
    serializeJsonMember(s, "length", length_); s << ",";
    serializeJsonMember(s, "volume", volume_); s << ",";
    serializeJsonMember(s, "minRadiusAvg", minRadiusAvg_); s << ",";
    serializeJsonMember(s, "minRadiusStdDeviation", minRadiusStdDeviation_); s << ",";
    serializeJsonMember(s, "maxRadiusAvg", maxRadiusAvg_); s << ",";
    serializeJsonMember(s, "maxRadiusStdDeviation", maxRadiusStdDeviation_); s << ",";
    serializeJsonMember(s, "avgRadiusAvg", avgRadiusAvg_); s << ",";
    serializeJsonMember(s, "avgRadiusStdDeviation", avgRadiusStdDeviation_); s << ",";
    serializeJsonMember(s, "roundnessAvg", roundnessAvg_); s << ",";
    serializeJsonMember(s, "roundnessStdDeviation", roundnessStdDeviation_); s << ",";
    serializeJsonMember(s, "innerLengthNode1", innerLengthNode1_); s << ",";
    serializeJsonMember(s, "innerLengthNode2", innerLengthNode2_); s << ",";
    serializeJsonMember(s, "tipRadiusNode1_", tipRadiusNode1_); s << ",";
    serializeJsonMember(s, "tipRadiusNode2_", tipRadiusNode2_);
    s << "}";
}

VesselGraphEdgePathPropertiesSerializable::VesselGraphEdgePathPropertiesSerializable()
    : inner_()
{
}
void VesselGraphEdgePathPropertiesSerializable::serialize(Serializer& s) const {
    s.serialize("length", inner_.length_);
    s.serialize("volume", inner_.volume_);
    s.serialize("minRadiusAvg", inner_.minRadiusAvg_);
    s.serialize("minRadiusStdDeviation", inner_.minRadiusStdDeviation_);
    s.serialize("maxRadiusAvg", inner_.maxRadiusAvg_);
    s.serialize("maxRadiusStdDeviation", inner_.maxRadiusStdDeviation_);
    s.serialize("avgRadiusAvg", inner_.avgRadiusAvg_);
    s.serialize("avgRadiusStdDeviation", inner_.avgRadiusStdDeviation_);
    s.serialize("roundnessAvg", inner_.roundnessAvg_);
    s.serialize("roundnessStdDeviation", inner_.roundnessStdDeviation_);
    s.serialize("innerLengthNode1", inner_.innerLengthNode1_);
    s.serialize("innerLengthNode2", inner_.innerLengthNode2_);
    s.serialize("tipRadiusNode1_", inner_.tipRadiusNode1_);
    s.serialize("tipRadiusNode2_", inner_.tipRadiusNode2_);
}
void VesselGraphEdgePathPropertiesSerializable::deserialize(Deserializer& s) {
    s.deserialize("length", inner_.length_);
    s.deserialize("volume", inner_.volume_);
    s.deserialize("minRadiusAvg", inner_.minRadiusAvg_);
    s.deserialize("minRadiusStdDeviation", inner_.minRadiusStdDeviation_);
    s.deserialize("maxRadiusAvg", inner_.maxRadiusAvg_);
    s.deserialize("maxRadiusStdDeviation", inner_.maxRadiusStdDeviation_);
    s.deserialize("avgRadiusAvg", inner_.avgRadiusAvg_);
    s.deserialize("avgRadiusStdDeviation", inner_.avgRadiusStdDeviation_);
    s.deserialize("roundnessAvg", inner_.roundnessAvg_);
    s.deserialize("roundnessStdDeviation", inner_.roundnessStdDeviation_);
    s.deserialize("tipRadiusNode1_", inner_.tipRadiusNode1_);
    s.deserialize("tipRadiusNode2_", inner_.tipRadiusNode2_);
}

VesselGraphEdge::VesselGraphEdge(VesselGraph& graph, VGEdgeID id, VGNodeID node1ID, VGNodeID node2ID, DiskArray<VesselSkeletonVoxel>&& voxels, VesselGraphEdgeUUID uuid)
    : graph_(&graph)
    , id_(id)
    , node1_(node1ID)
    , node2_(node2ID)
    , distance_(std::numeric_limits<float>::quiet_NaN())
    , pathProps_() //Invalid until finalized
    , voxels_(std::move(voxels))
    , uuid_(uuid)
{
}
void VesselGraphEdge::finalizeConstruction() {
    VesselGraphNode& node1 = getNode1();
    VesselGraphNode& node2 = getNode2();
    // Compute distance
    distance_ = tgt::distance(node1.pos_, node2.pos_);

    size_t firstInner = voxels_.size();
    size_t lastInner = -1;
    size_t longestRunBegin = 0;
    size_t longestRunEnd = 0;
    bool inRun = false;
    size_t currentRunBegin = 0;
    for(size_t i=0; i < voxels_.size(); ++i) {
        if(voxels_[i].isInner()) {
            if(firstInner == voxels_.size()) {
                firstInner = i;
            }
            lastInner = i;
            if(inRun) {
                size_t currentRunEnd = i;
                size_t len = currentRunEnd - currentRunBegin;
                if(len > longestRunEnd-longestRunBegin) {
                    longestRunBegin = currentRunBegin;
                    longestRunEnd = currentRunEnd;
                }
                inRun = false;
            }
        } else {
            if(!inRun) {
                currentRunBegin = i;
                inRun = true;
            }
        }
    }
    if(inRun) {
        size_t currentRunEnd = voxels_.size();
        size_t len = currentRunEnd - currentRunBegin;
        if(len > longestRunEnd-longestRunBegin) {
            longestRunBegin = currentRunBegin;
            longestRunEnd = currentRunEnd;
        }
    }

    bool allInner = firstInner == 0 && lastInner == voxels_.size()-1;

    // If we do not have any outer voxels, at least try to even out the innerLengths of both nodes.
    if(longestRunBegin == longestRunEnd) {
        longestRunBegin = voxels_.size()/2;
        longestRunEnd = longestRunBegin;
    }

    if(getNode2().isEndNode()) {
        // At least one node is a leaf: The first inner voxel (starting from the node) marks the end of the outer path
        outerPathBeginIndex_ = lastInner+1;
    } else {
        if(getNode1().isEndNode()) {
            if(allInner) {
                // Special case for end standing edges with all inner voxels
                outerPathBeginIndex_ = lastInner+1;
            } else {
                // At least one node is a leaf: the outer path begins directly at that node
                outerPathBeginIndex_ = 0;
            }
        } else {
            // Both nodes are branching points: use the longest run of outer voxels
            outerPathBeginIndex_ = longestRunBegin;
        }
    }
    if(getNode1().isEndNode()) {
        // At least one node is a leaf: The first inner voxel (starting from the node) marks the end of the outer path
        outerPathEndIndex_ = firstInner;
    } else {
        if(getNode2().isEndNode()) {
            if(allInner) {
                // Special case for end standing edges with all inner voxels
                outerPathEndIndex_ = firstInner;
            } else {
                // At least one node is a leaf: the outer path begins directly at that node
                outerPathEndIndex_ = voxels_.size();
            }
        } else {
            // Both nodes are branching points: use the longest run of outer voxels
            outerPathEndIndex_ = longestRunEnd;
        }
    }

    // Compute path properties
    pathProps_ = VesselGraphEdgePathProperties::fromPath(node1, node2, voxels_, outerPathBeginIndex_, outerPathEndIndex_);

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
    , outerPathBeginIndex_(other.outerPathBeginIndex_)
    , outerPathEndIndex_(other.outerPathEndIndex_)
{
    other.graph_ = nullptr;
}

VesselGraphEdge& VesselGraphEdge::operator=(VesselGraphEdge&& other)
{
    if(this != &other) {
        this->~VesselGraphEdge();
        new(this) VesselGraphEdge(std::move(other));
    }
    return *this;
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
        return std::numeric_limits<float>::infinity();
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
    float len = std::max(0.0f, pathProps_.length_ - pathProps_.innerLengthNode1_ - pathProps_.innerLengthNode2_);
    if(diameter > 0) {
        return len / diameter;
    } else {
        //tgtAssert(false, "Invalid radius");
        return 1; // Return neutral elongation
    }
}

float VesselGraphEdge::getRelativeBulgeSize() const {
    //float edge_length_contributed_length = std::max(getNode1().estimatedRadius(), getNode2().estimatedRadius());
    //return std::max(0.0f,(getLength() / edge_length_contributed_length) - 1);
    if(hasValidData()) {
        if(getNode1().isEndNode() == getNode2().isEndNode()) {
            return -1;
        }
        float innerLength;
        float tipRadius;
        if(getNode1().isEndNode()) {
            innerLength = pathProps_.innerLengthNode2_;
            tipRadius = pathProps_.tipRadiusNode1_;
        } else {
            innerLength = pathProps_.innerLengthNode1_;
            tipRadius = pathProps_.tipRadiusNode2_;
        }
        float radius = getAvgRadiusAvg();
        float outerLength = getLength() - innerLength + tipRadius;
        tgtAssert(radius > 0, "Invalid vessel radius");
        tgtAssert(outerLength >= 0, "Invalid outer length");
        return outerLength / radius;
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
DiskArray<VesselSkeletonVoxel> VesselGraphEdge::getOuterVoxels() const {
    if(outerPathBeginIndex_ > outerPathEndIndex_) {
        return voxels_.slice(outerPathBeginIndex_, outerPathBeginIndex_);
    } else {
        return voxels_.slice(outerPathBeginIndex_, outerPathEndIndex_);
    }
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

void VesselGraphEdge::serializeToJson(SerializationOutputStream& s) const {
    s << "{";
    serializeJsonMember(s, "id", id_.raw()); s << ",";
    serializeJsonMember(s, "node1", node1_.raw()); s << ",";
    serializeJsonMember(s, "node2", node2_.raw()); s << ",";
    serializeJsonMember(s, "distance", distance_); s << ",";

    if(voxels_.empty()) {
        serializeJsonMember(s, "pathProperties", pathProps_);
    } else {
        serializeJsonMemberIterable(s, "skeletonVoxels", voxels_);
    }
    s << "}";
}

VesselGraphEdgeSerializable::VesselGraphEdgeSerializable(const VesselGraphEdge& e)
    : inner_(e)
{
}
void VesselGraphEdgeSerializable::serialize(Serializer& s) const {
    s.serialize("id", inner_.id_.raw());
    s.serialize("node1", inner_.node1_.raw());
    s.serialize("node2", inner_.node2_.raw());
    s.serialize("distance", inner_.distance_);

    if(inner_.voxels_.empty()) {
        VesselGraphEdgePathPropertiesSerializable props;
        props.inner_ = inner_.pathProps_;
        s.serialize("pathProperties", props);
    } else {
        std::vector<VesselSkeletonVoxelSerializable> voxels;
        voxels.reserve(voxels.size());
        for(const auto& voxel : inner_.voxels_) {
            voxels.emplace_back(voxel);
        }
        s.serialize("skeletonVoxels", voxels);
    }
}
void VesselGraphEdgeSerializable::deserialize(Deserializer&) {
    tgtAssert(false, "Cannot deserialize VesselGraphEdgeSerializable");
}
void VesselGraphEdgeDeserializable::serialize(Serializer&) const {
    tgtAssert(false, "Cannot serialize VesselGraphEdgeDeserializable");
}
void VesselGraphEdgeDeserializable::deserialize(Deserializer& s) {
    uint32_t id, n1, n2;
    s.deserialize("id", id);
    s.deserialize("node1", n1);
    s.deserialize("node2", n2);
    id_ = id;
    node1_ = n1;
    node2_ = n2;

    try {
        std::vector<VesselSkeletonVoxelSerializable> voxels;
        s.deserialize("skeletonVoxels", voxels);

        for(const auto& voxel : voxels) {
            voxels_.push_back(voxel.inner_);
        }
    } catch (SerializationException& s) {
    }
    if(voxels_.empty()) {
        VesselGraphEdgePathPropertiesSerializable props;
        s.deserialize("pathProperties", props);
        pathProps_ = props.inner_;
    }
}

// VesselGraph -------------------------------------------------------------------------

VesselGraph::VesselGraph(const tgt::Bounds& bounds)
    : nodes_(new DiskArrayStorage<VesselGraphNode>(VoreenApplication::app()->getUniqueTmpFilePath(".vgnodes")))
    , edges_(new DiskArrayStorage<VesselGraphEdge>(VoreenApplication::app()->getUniqueTmpFilePath(".vgedges")))
    , nodeEdgeIdStorage_(new DiskArrayBackedList<VGEdgeID>::Storage(VoreenApplication::app()->getUniqueTmpFilePath(".vgnodeedgerefs")))
    , edgeVoxelStorage_(new DiskArrayStorage<VesselSkeletonVoxel>(VoreenApplication::app()->getUniqueTmpFilePath(".vgedgevoxels")))
    , nodeVoxelStorage_(new DiskArrayStorage<tgt::vec3>(VoreenApplication::app()->getUniqueTmpFilePath(".vgnodevoxels")))
    , bounds_(bounds)
{
}

VesselGraph::VesselGraph()
    : VesselGraph(tgt::Bounds())
{
}

VesselGraph::VesselGraph(VesselGraph&& other)
    : nodes_(std::move(other.nodes_))
    , edges_(std::move(other.edges_))
    , nodeEdgeIdStorage_(std::move(other.nodeEdgeIdStorage_))
    , bounds_(other.bounds_)
    , edgeVoxelStorage_(std::move(other.edgeVoxelStorage_))
    , nodeVoxelStorage_(std::move(other.nodeVoxelStorage_))
{
    // Fix up references:
    for(auto& node: nodes_->asArray()) {
        node.graph_ = this;
    }

    for(auto& edge: edges_->asArray()) {
        edge.graph_ = this;
    }
}

VesselGraph& VesselGraph::operator=(VesselGraph&& other) {
    if(this != &other) {
        this->~VesselGraph();
        new(this) VesselGraph(std::move(other));
    }
    return *this;
}

VesselGraph VesselGraph::clone() const {

    VesselGraphBuilder builder(getBounds());
    for(auto& node: nodes_->asArray()) {
        builder.insertNode(node);
    }

    for(auto& edge: edges_->asArray()) {
        builder.insertEdge(edge.getNodeID1(), edge.getNodeID2(), edge, edge.getUUID());
    }
    auto ptr = std::move(builder).finalize();
    return VesselGraph(std::move(*ptr));
}

const VesselGraphNode& VesselGraph::getNode(VGNodeID i) const {
    return (*nodes_)[i.raw()];
}

const VesselGraphEdge& VesselGraph::getEdge(VGEdgeID i) const {
    return (*edges_)[i.raw()];
}

VesselGraphNode& VesselGraph::getNode(VGNodeID i) {
    return (*nodes_)[i.raw()];
}

VesselGraphEdge& VesselGraph::getEdge(VGEdgeID i) {
    return (*edges_)[i.raw()];
}

DiskArray<VesselGraphNode> VesselGraph::getNodes() const {
    return nodes_->asArray();
}

DiskArray<VesselGraphEdge> VesselGraph::getEdges() const {
    return edges_->asArray();
}

DiskArray<VesselGraphNode> VesselGraph::getNodes() {
    return nodes_->asArray();
}

DiskArray<VesselGraphEdge> VesselGraph::getEdges() {
    return edges_->asArray();
}

void VesselGraph::getEdgePropertyStats(std::function<float(const VesselGraphEdge&)> f, float& /*out*/ mean, float& /*out*/stddev) const {
    statisticalAnalysis(getEdges(), f, mean, stddev);
}

const tgt::Bounds& VesselGraph::getBounds() const {
    return bounds_;
}

void VesselGraph::serializeToJson(SerializationOutputStream& s) const {
    s << "{\"graph\":{";
    serializeJsonMemberIterable(s, "nodes", getNodes()); s << ",";
    serializeJsonMemberIterable(s, "edges", getEdges()); s << ",";
    s << "\"bounds\":{";
    serializeJsonMember(s, "First", bounds_.getLLF()); s << ",";
    serializeJsonMember(s, "Second", bounds_.getURB());
    s << "}}}";
}

void VesselGraph::serialize(Serializer& s) const {
    std::vector<VesselGraphNodeSerializable> nodes;
    for(const auto& node : getNodes()) {
        nodes.emplace_back(node);
    }
    s.serialize("nodes", nodes);


    std::vector<VesselGraphEdgeSerializable> edges;
    for(const auto& edge : getEdges()) {
        edges.emplace_back(edge);
    }
    s.serialize("edges", edges);

    s.serialize("bounds", bounds_);
}
void VesselGraph::deserialize(Deserializer& s) {
    VesselGraphBuilder builder;

    std::vector<VesselGraphNodeDeserializable> nodes;
    s.deserialize("nodes", nodes);
    for(const auto& node : nodes) {
        auto index = builder.insertNode(node.pos_, node.voxels_, node.radius_, node.isAtSampleBorder_);
        tgtAssert(index == node.id_, "Invalid node id after deserialization");
    }

    std::vector<VesselGraphEdgeDeserializable> edges;
    s.deserialize("edges", edges, "");
    for(const auto& edge : edges) {
        VGEdgeID index;
        if(edge.voxels_.empty()) {
            index = builder.insertEdge(edge.node1_, edge.node2_, edge.pathProps_);
        } else {
            index = builder.insertEdge(edge.node1_, edge.node2_, edge.voxels_);
        }
        tgtAssert(index == edge.id_, "Invalid node id after deserialization");
    }
    auto graph = std::move(builder).finalize();
    *this = std::move(*graph);

    s.deserialize("bounds", bounds_);
}

VesselGraphBuilder::VesselGraphBuilder(const tgt::Bounds& bounds)
    : graph_(new VesselGraph(bounds))
{
}

VesselGraphBuilder::VesselGraphBuilder()
    : graph_(new VesselGraph())
{
}

VGNodeID VesselGraphBuilder::insertNode(const VesselGraphNode& base) {
    return insertNode(base.pos_, base.voxels_, base.getRadius(), base.isAtSampleBorder_, base.getUUID());
}

VGNodeID VesselGraphBuilder::insertNode(const tgt::vec3& position, const DiskArray<tgt::vec3>& voxels, float radius, bool isAtSampleBorder) {
    return insertNode(position, voxels, radius, isAtSampleBorder, VoreenApplication::app()->generateUUID());
}
VGNodeID VesselGraphBuilder::insertNode(const tgt::vec3& position, const DiskArray<tgt::vec3>& voxels, float radius, bool isAtSampleBorder, VesselGraphNodeUUID uuid) {
    size_t nodeID = graph_->nodes_->size();
    graph_->bounds_.addPoint(position);
    graph_->nodes_->storeElement(VesselGraphNode(*graph_, nodeID, position, graph_->nodeVoxelStorage_->store(voxels), radius, isAtSampleBorder, uuid));
    return nodeID;
}
VGNodeID VesselGraphBuilder::insertNode(const tgt::vec3& position, const std::vector<tgt::vec3>& voxels, float radius, bool isAtSampleBorder) {
    return insertNode(position, voxels, radius, isAtSampleBorder, VoreenApplication::app()->generateUUID());
}
VGNodeID VesselGraphBuilder::insertNode(const tgt::vec3& position, const std::vector<tgt::vec3>& voxels, float radius, bool isAtSampleBorder, VesselGraphNodeUUID uuid) {
    size_t nodeID = graph_->nodes_->size();
    graph_->bounds_.addPoint(position);
    graph_->nodes_->storeElement(VesselGraphNode(*graph_, nodeID, position, graph_->nodeVoxelStorage_->store(voxels), radius, isAtSampleBorder, uuid));
    return nodeID;
}


VGEdgeID VesselGraphBuilder::insertEdge(VGNodeID node1, VGNodeID node2, const DiskArray<VesselSkeletonVoxel>& voxels) {
    return insertEdge(node1, node2, voxels, VoreenApplication::app()->generateUUID());
}
VGEdgeID VesselGraphBuilder::insertEdge(VGNodeID node1, VGNodeID node2, const DiskArray<VesselSkeletonVoxel>& voxels, VesselGraphEdgeUUID uuid) {
    tgtAssert(node1 < graph_->nodes_->size(), "Edge references nonexistent node");
    tgtAssert(node2 < graph_->nodes_->size(), "Edge references nonexistent node");
    VesselGraphNode& n1 = graph_->getNode(node1);
    VesselGraphNode& n2 = graph_->getNode(node2);

    VGEdgeID edgeID = graph_->edges_->size();
    graph_->edges_->storeElement(VesselGraphEdge(*graph_, edgeID, node1, node2, graph_->edgeVoxelStorage_->store(voxels), uuid));

    n1.edges_.push(edgeID);
    n2.edges_.push(edgeID);
    return edgeID;
}
VGEdgeID VesselGraphBuilder::insertEdge(VGNodeID node1, VGNodeID node2, const std::vector<VesselSkeletonVoxel>& voxels) {
    return insertEdge(node1, node2, voxels, VoreenApplication::app()->generateUUID());
}
VGEdgeID VesselGraphBuilder::insertEdge(VGNodeID node1, VGNodeID node2, const std::vector<VesselSkeletonVoxel>& voxels, VesselGraphEdgeUUID uuid) {
    tgtAssert(node1 < graph_->nodes_->size(), "Edge references nonexistent node");
    tgtAssert(node2 < graph_->nodes_->size(), "Edge references nonexistent node");
    VesselGraphNode& n1 = graph_->getNode(node1);
    VesselGraphNode& n2 = graph_->getNode(node2);

    VGEdgeID edgeID = graph_->edges_->size();
    graph_->edges_->storeElement(VesselGraphEdge(*graph_, edgeID, node1, node2, graph_->edgeVoxelStorage_->store(voxels), uuid));

    n1.edges_.push(edgeID);
    n2.edges_.push(edgeID);
    return edgeID;
}

VGEdgeID VesselGraphBuilder::insertEdge(VGNodeID node1, VGNodeID node2, const VesselGraphEdge& path_definition, VesselGraphEdgeUUID uuid) {
    return insertEdge(node1, node2, path_definition.voxels_, uuid);
}

VGEdgeID VesselGraphBuilder::insertEdge(VGNodeID node1, VGNodeID node2, const VesselGraphEdge& path_definition) {
    return insertEdge(node1, node2, path_definition, VoreenApplication::app()->generateUUID());
}

VGEdgeID VesselGraphBuilder::insertEdge(VGNodeID node1, VGNodeID node2, VesselGraphEdgePathProperties pathProperties) {
    return insertEdge(node1, node2, pathProperties, VoreenApplication::app()->generateUUID());
}
VGEdgeID VesselGraphBuilder::insertEdge(VGNodeID node1, VGNodeID node2, VesselGraphEdgePathProperties pathProperties, VesselGraphEdgeUUID uuid) {
    tgtAssert(node1 < graph_->nodes_->size(), "Edge references nonexistent node");
    tgtAssert(node2 < graph_->nodes_->size(), "Edge references nonexistent node");
    VesselGraphNode& n1 = graph_->getNode(node1);
    VesselGraphNode& n2 = graph_->getNode(node2);

    VGEdgeID edgeID = graph_->edges_->size();
    graph_->edges_->storeElement(VesselGraphEdge(*graph_, edgeID, node1, node2, pathProperties, uuid));

    n1.edges_.push(edgeID);
    n2.edges_.push(edgeID);
    return edgeID;
}

const VesselGraphNode& VesselGraphBuilder::getNode(VGNodeID i) const {
    return graph_->getNode(i);
}
std::unique_ptr<VesselGraph> VesselGraphBuilder::finalize() && {
    for(auto& edge : graph_->getEdges()) {
        edge.finalizeConstruction();
    }
    return std::move(graph_);
}

} // namespace voreen
