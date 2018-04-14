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
#ifndef VRN_VESSELGRAPH_H
#define VRN_VESSELGRAPH_H

#include <vector>
#include <deque>

#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/bounds.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"

#include <functional>

namespace voreen {

/*
 * Classes for storage of graph extracted from segmentations of vessels.
 */

//Forward declarations...
struct VesselGraphNode;
struct VesselGraphEdge;
class VesselGraph;

// A single voxel in a branch in the vessel graph
struct VesselSkeletonVoxel : public Serializable {
    VesselSkeletonVoxel(const tgt::vec3& pos, float minDistToSurface, float maxDistToSurface, float avgDistToSurface, size_t numSurfaceVoxels, float volume);

    tgt::vec3 pos_;
    float minDistToSurface_;
    float maxDistToSurface_;
    float avgDistToSurface_;
    size_t numSurfaceVoxels_; // Surface voxels that belong to this skeleton voxel
    float volume_;           // Volume of the segmentation that belong to this
                             // skeleton voxel.
                             // float, because an object voxel might be closest
                             // to multiple skeleton voxels => volume will be split.

    // Compute the roundness, i.e. minDistToSurface_/maxDistToSurface_
    float roundness() const;
    // Determine whether this voxel has valid data, i.e., has any associated
    // surface voxels that were used to compute the distances and other data.
    bool hasValidData() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:
    // Only for deserialization. you should probably not use this.
    friend class XmlDeserializer;
    friend class JsonDeserializer;
    VesselSkeletonVoxel();
};

// A node within the vessel graph.
// It stores its position, references to edges as well as all voxels that define this node.
struct VesselGraphNode : public Serializable {
    VesselGraphNode(VesselGraph& graph, size_t id, const tgt::vec3& position, std::vector<tgt::vec3> voxels, bool isAtSampleBorder);

    VesselGraphNode(VesselGraphNode&&);
    void operator=(VesselGraphNode&&);

    std::vector<std::reference_wrapper<const VesselGraphEdge>> getEdges() const;
    std::vector<const VesselGraphNode*> getNeighbors() const;
    int getDegree() const;
    bool isEndNode() const;
    size_t getID() const;

    float estimatedRadius() const;

    size_t id_;
    std::vector<size_t> edges_;
    tgt::vec3 pos_;
    std::vector<tgt::vec3> voxels_;
    bool isAtSampleBorder_;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:
    VesselGraph* graph_; // Will never be null (except briefly during deserialization)

private:
    // Disable copy constructors:
    VesselGraphNode(const VesselGraphNode&);
    void operator=(const VesselGraphNode&);

private:
    // Only for deserialization. you should probably not use this.
    friend class XmlDeserializer;
    friend class JsonDeserializer;
    friend class VesselGraph;
    VesselGraphNode();
};

struct VesselGraphEdgePathProperties : public Serializable {
    const static float INVALID_DATA;

    VesselGraphEdgePathProperties();

    float length_;
    float volume_;
    float minRadiusAvg_;
    float minRadiusStdDeviation_;
    float maxRadiusAvg_;
    float maxRadiusStdDeviation_;
    float avgRadiusAvg_;
    float avgRadiusStdDeviation_;
    float roundnessAvg_;
    float roundnessStdDeviation_;

    static VesselGraphEdgePathProperties fromPath(const VesselGraphNode& begin, const VesselGraphNode& end, const std::vector<VesselSkeletonVoxel>& path);
    bool hasValidData() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
};

// An edge within a vessel graph.
// It stores references to its nodes, properties of the associated branch of the vessel network, as well as the medial
// axis of the branch.
struct VesselGraphEdge : public Serializable {
    // Construct an edge implicitly. All properties will be calculated from the voxels
    VesselGraphEdge(VesselGraph& graph, size_t id, size_t node1ID, size_t node2ID, const std::vector<VesselSkeletonVoxel>&& voxels);
    // Construct an edge excplitily, all properties must be given. The path will be a straight line between two nodes
    VesselGraphEdge(VesselGraph& graph, size_t id, size_t node1ID, size_t node2ID, VesselGraphEdgePathProperties pathProps);

    // Move constructor
    VesselGraphEdge(VesselGraphEdge&& other);
    void operator=(VesselGraphEdge&& other);

    /**
     * Indicates whether or not this edge has any valid edge voxels and therefore
     * has valid values for volume etc.
     */
    bool hasValidData() const;

    float getLength() const;
    float getDistance() const;
    float getCurveness() const;
    float getStraightness() const;
    float getVolume() const;
    float getAvgCrossSection() const;
    float getMinRadiusAvg() const;
    float getMinRadiusStdDeviation() const;
    float getMaxRadiusAvg() const;
    float getMaxRadiusStdDeviation() const;
    float getAvgRadiusAvg() const;
    float getAvgRadiusStdDeviation() const;
    float getRoundnessAvg() const;
    float getRoundnessStdDeviation() const;
    const std::vector<VesselSkeletonVoxel>& getVoxels() const;

    float getElongation() const;
    float getEffectiveLength() const;
    float getRelativeBulgeSize() const;

    bool isLoop() const;

    const VesselGraphNode& getNode1() const;
    const VesselGraphNode& getNode2() const;

    // Get the node that this edge connects the given node to.
    // Note: It is assumed that firstNode is indeed a node of this edge!
    // Note: If this edge is a cycle, firstNode is returned.
    const VesselGraphNode& getOtherNode(const VesselGraphNode& firstNode) const;

    VesselGraphNode& getNode1();
    VesselGraphNode& getNode2();

    const VesselGraphEdgePathProperties& getPathProperties() const;

    size_t getNodeID1() const;
    size_t getNodeID2() const;


    size_t getID() const;

    bool isEndStanding() const;
    size_t getNumValidVoxels() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:
    size_t id_;
    size_t node1_;
    size_t node2_;

    float distance_; //cached in order to avoid expensive lookups and distance calculation
    VesselGraphEdgePathProperties pathProps_;

    //NOTE: the path in voxels (geometrically) starts at node1_ and ends in node2_.
    std::vector<VesselSkeletonVoxel> voxels_;

private:
    VesselGraph* graph_; // Will never be null (except briefly during deserialization)

private:
    // Disable copy constructors:
    VesselGraphEdge(const VesselGraphEdge&);
    void operator=(const VesselGraphEdge&);

private:
    // Only for deserialization. you should probably not use this.
    friend class XmlDeserializer;
    friend class JsonDeserializer;
    friend class VesselGraph;
    VesselGraphEdge();
    void updatePathPropertiesFromVoxels();
};

// The vessel graph itself. it stores nodes as well as edges.
// References between nodes and edges are stored within the substrucutres.
//
// To avoid pointer/reference invalidation nodes and edges can only be added to the graph,
// but not removed.
class VesselGraph : public Serializable {
public:
    // Create a graph with predetermined bounds
    VesselGraph(const tgt::Bounds& bounds);
    // Create a graph with undefined bounds
    VesselGraph();
    // Create a deep copy from another graph
    VesselGraph(const VesselGraph& original);

    virtual ~VesselGraph() {}

    // Insert a new node into the graph and clone it from the given existing node
    size_t insertNode(const VesselGraphNode& base);
    // Insert a new node into the graph and construct it from the given parameters
    size_t insertNode(const tgt::vec3& position, const std::vector<tgt::vec3>&& voxels, bool isAtSampleBorder);

    // Insert a new edge by deriving its properties from the provided path
    size_t insertEdge(size_t node1, size_t node2, const std::vector<VesselSkeletonVoxel>&& path);
    // Insert a new edge by deriving its properties from the path of the provided edge
    size_t insertEdge(size_t node1, size_t node2, const VesselGraphEdge& path_definition);
    // Insert a new edge by providing predetermined pathProperties. The path itself will be a straight line between two nodes.
    size_t insertEdge(size_t node1, size_t node2, VesselGraphEdgePathProperties pathProperties);

    const VesselGraphNode& getNode(size_t i) const;
    const VesselGraphEdge& getEdge(size_t i) const;
    VesselGraphNode& getNode(size_t i);
    VesselGraphEdge& getEdge(size_t i);

    const std::vector<VesselGraphNode>& getNodes() const;
    const std::vector<VesselGraphEdge>& getEdges() const;
    std::vector<VesselGraphNode>& getNodes();
    std::vector<VesselGraphEdge>& getEdges();

    void getEdgePropertyStats(std::function<float(const VesselGraphEdge&)>, float& /*out*/ mean, float& /*out*/stddev) const;

    // Get a bounding box that encompasses all nodes within the graph.
    const tgt::Bounds& getBounds() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:
    // Note: We can use vector here because we do not store pointers to other edges/nodes in nodes/edges, but only indices
    std::vector<VesselGraphNode> nodes_;
    std::vector<VesselGraphEdge> edges_;

    tgt::Bounds bounds_;
};

} // namespace voreen

#endif // VRN_VESSELGRAPH_H
