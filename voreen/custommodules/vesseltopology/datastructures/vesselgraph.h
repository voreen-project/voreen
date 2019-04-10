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

#ifndef VRN_VESSELGRAPH_H
#define VRN_VESSELGRAPH_H

#include <vector>

#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/bounds.h"

#include "../datastructures/diskarraystorage.h"

#ifndef VRN_VESSELTOPOLOGY_MINIMAL_VESSELGRAPH
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"
#endif

#include <boost/uuid/uuid.hpp>
#include <functional>

namespace voreen {

/*
 * Classes for storage of graph extracted from segmentations of vessels.
 */

//Forward declarations...
struct VesselGraphNode;
struct VesselGraphEdge;
class VesselGraph;


typedef boost::uuids::uuid VesselGraphEdgeUUID;
typedef boost::uuids::uuid VesselGraphNodeUUID;

// A single voxel in a branch in the vessel graph
struct VesselSkeletonVoxel {
    VesselSkeletonVoxel(const tgt::vec3& pos, float minDistToSurface, float maxDistToSurface, float avgDistToSurface, uint32_t numSurfaceVoxels, float volume, bool nearOtherEdge);

    tgt::vec3 pos_;
    float minDistToSurface_;
    float maxDistToSurface_;
    float avgDistToSurface_;
    uint32_t numSurfaceVoxels_; // Surface voxels that belong to this skeleton voxel
    float volume_;           // Volume of the segmentation that belong to this
                             // skeleton voxel.
                             // float, because an object voxel might be closest
                             // to multiple skeleton voxels => volume will be split.
                             //

    // Whether or not this skeleton voxel has associated surface points that are adjacent
    // to voxels of other edges.
    bool nearOtherEdge_;

    // Compute the roundness, i.e. minDistToSurface_/maxDistToSurface_
    float roundness() const;
    // Determine whether this voxel has valid data, i.e., has any associated
    // surface voxels that were used to compute the distances and other data.
    bool hasValidData() const;

    // Is a voxel contained in an intersection
    bool isInner() const;
    // .. or not
    bool isOuter() const;
private:
    friend struct VesselSkeletonVoxelSerializable;
    VesselSkeletonVoxel();
};

#ifndef VRN_VESSELTOPOLOGY_MINIMAL_VESSELGRAPH
struct VesselSkeletonVoxelSerializable : public Serializable {
    VesselSkeletonVoxelSerializable(VesselSkeletonVoxel);
    VesselSkeletonVoxel inner_;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
private:
    // Only for deserialization. you should probably not use this.
    friend class Deserializer;
    VesselSkeletonVoxelSerializable();
};
#endif

class VGNodeID {
    uint32_t internal_;
public:
    const static VGNodeID INVALID;
    VGNodeID()
        : internal_((uint32_t)-1)
    {
    }

    VGNodeID(uint32_t i)
        : internal_(i)
    {
    };
    inline bool operator==(const VGNodeID& other) const {
        return internal_ == other.internal_;
    }
    inline bool operator!=(const VGNodeID& other) const {
        return internal_ != other.internal_;
    }
    inline bool operator<(const VGNodeID& other) const {
        return internal_ < other.internal_;
    }
    inline uint32_t raw() const {
        return internal_;
    }
    inline bool isValid() const {
        return internal_ != (uint32_t)-1;
    }
};

class VGEdgeID {
    uint32_t internal_;
public:
    const static VGEdgeID INVALID;
    VGEdgeID()
        : internal_((uint32_t)-1)
    {
    }
    VGEdgeID(uint32_t i)
        : internal_(i)
    {
    };
    inline bool operator==(const VGEdgeID& other) const {
        return internal_ == other.internal_;
    }
    inline bool operator!=(const VGEdgeID& other) const {
        return internal_ != other.internal_;
    }
    inline bool operator<(const VGEdgeID& other) const {
        return internal_ < other.internal_;
    }
    inline uint32_t raw() const {
        return internal_;
    }
    inline bool isValid() const {
        return internal_ != (uint32_t)-1;
    }
};

// A node within the vessel graph.
// It stores its position, references to edges as well as all voxels that define this node.
struct VesselGraphNode {
    VesselGraphNode(VesselGraph& graph, VGNodeID id, const tgt::vec3& position, DiskArray<tgt::vec3>&& voxels, float radius, bool isAtSampleBorder, VesselGraphNodeUUID uuid);

    VesselGraphNode(VesselGraphNode&&);
    VesselGraphNode& operator=(VesselGraphNode&&);

    std::vector<std::reference_wrapper<const VesselGraphEdge>> getEdges() const;
    std::vector<const VesselGraphEdge*> getEdgesAsPtrs() const; //All pointers guaranteed to be non-null
    std::vector<const VesselGraphNode*> getNeighbors() const;
    int getDegree() const;
    bool isEndNode() const;
    VGNodeID getID() const;
    float getRadius() const;
    VesselGraphNodeUUID getUUID() const;

    VGNodeID id_;
    VesselGraphNodeUUID uuid_;
    DiskArrayBackedList<VGEdgeID> edges_;
    tgt::vec3 pos_;
    DiskArray<tgt::vec3> voxels_;
    bool isAtSampleBorder_;
    float radius_;

private:
    VesselGraph* graph_; // Will never be null

private:
    // Disable copy constructors:
    VesselGraphNode(const VesselGraphNode&) = delete;
    void operator=(const VesselGraphNode&) = delete;

private:
    friend class VesselGraph;
    friend class VesselGraphBuilder;
    // Only for deserialization. you should probably not use this.
    friend struct VesselGraphNodeDeserializable;
};

#ifndef VRN_VESSELTOPOLOGY_MINIMAL_VESSELGRAPH
struct VesselGraphNodeSerializable : public Serializable {
    VesselGraphNodeSerializable(const VesselGraphNode&);
    const VesselGraphNode& inner_;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
private:
    friend class Deserializer;
};
#endif
#ifndef VRN_VESSELTOPOLOGY_MINIMAL_VESSELGRAPH
struct VesselGraphNodeDeserializable : public Serializable {
    VGNodeID id_;
    tgt::vec3 pos_;
    std::vector<tgt::vec3> voxels_;
    float radius_;
    bool isAtSampleBorder_;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
};
#endif

struct VesselGraphEdgePathProperties {
    const static float INVALID_DATA;

    VesselGraphEdgePathProperties();

    float length_;
    float volume_;
    float minRadiusAvg_;
    float minRadiusStdDeviation_;
    float maxRadiusAvg_;
    float maxRadiusStdDeviation_;
    float maxRadiusMax_;
    float avgRadiusAvg_;
    float avgRadiusStdDeviation_;
    float roundnessAvg_;
    float roundnessStdDeviation_;

    float innerLengthNode1_;
    float innerLengthNode2_;
    float tipRadiusNode1_;
    float tipRadiusNode2_;

    static VesselGraphEdgePathProperties fromPath(const VesselGraphNode& begin, const VesselGraphNode& end, const DiskArray<VesselSkeletonVoxel>& path, size_t outerPathBeginIndex, size_t outerPathEndIndex);
    bool hasValidData() const;
};

#ifndef VRN_VESSELTOPOLOGY_MINIMAL_VESSELGRAPH
struct VesselGraphEdgePathPropertiesSerializable : public Serializable {
    VesselGraphEdgePathProperties inner_;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    VesselGraphEdgePathPropertiesSerializable();
};
#endif

// An edge within a vessel graph.
// It stores references to its nodes, properties of the associated branch of the vessel network, as well as the medial
// axis of the branch.
struct VesselGraphEdge {
    // Construct an edge implicitly. All properties will be calculated from the voxels
    VesselGraphEdge(VesselGraph& graph, VGEdgeID id, VGNodeID node1ID, VGNodeID node2ID, DiskArray<VesselSkeletonVoxel>&& voxels, VesselGraphEdgeUUID uuid);

    // Construct an edge excplitily, all properties must be given. The path will be a straight line between two nodes
    VesselGraphEdge(VesselGraph& graph, VGEdgeID id, VGNodeID node1ID, VGNodeID node2ID, VesselGraphEdgePathProperties pathProps, VesselGraphEdgeUUID uuid);

    // Move constructor
    VesselGraphEdge(VesselGraphEdge&& other);
    VesselGraphEdge& operator=(VesselGraphEdge&& other);

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
    float getMaxRadiusMax() const;
    float getAvgRadiusAvg() const;
    float getAvgRadiusStdDeviation() const;
    float getRoundnessAvg() const;
    float getRoundnessStdDeviation() const;
    const DiskArray<VesselSkeletonVoxel>& getVoxels() const;

    // Only valid for leaf branches!
    DiskArray<VesselSkeletonVoxel> getOuterVoxels() const;

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

    VGNodeID getNodeID1() const;
    VGNodeID getNodeID2() const;


    // Returns the identifier for edges within the graph
    VGEdgeID getID() const;

    // Returns the globally unique identifier
    VesselGraphEdgeUUID getUUID() const;

    bool isEndStanding() const;
    size_t getNumValidVoxels() const;

    size_t outerPathBeginIndex_;
    size_t outerPathEndIndex_;

private:
    VGEdgeID id_; //within the graph
    VGNodeID node1_;
    VGNodeID node2_;

    float distance_; //cached in order to avoid expensive lookups and distance calculation
    VesselGraphEdgePathProperties pathProps_;

    VesselGraphEdgeUUID uuid_; //globally

    //NOTE: the path in voxels (geometrically) starts at node1_ and ends in node2_.
    DiskArray<VesselSkeletonVoxel> voxels_;

private:
    VesselGraph* graph_; // Will never be null

private:
    // Disable copy constructors:
    VesselGraphEdge(const VesselGraphEdge&) = delete;
    void operator=(const VesselGraphEdge&) = delete;

    // Compute remaining properties that my depend on other properties of the graph.
    void finalizeConstruction();

private:
    // Only for deserialization. you should probably not use this.
    friend class VesselGraph;
    friend class VesselGraphBuilder;
    friend struct VesselGraphEdgeSerializable;
    friend struct VesselGraphEdgeDeserializable;
    VesselGraphEdge();
};

#ifndef VRN_VESSELTOPOLOGY_MINIMAL_VESSELGRAPH
struct VesselGraphEdgeSerializable : public Serializable {
    VesselGraphEdgeSerializable(const VesselGraphEdge&);
    const VesselGraphEdge& inner_;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
private:
    // Only for deserialization. you should probably not use this.
    friend class Deserializer;
    //VesselGraphEdgeSerializable();
};
#endif
#ifndef VRN_VESSELTOPOLOGY_MINIMAL_VESSELGRAPH
struct VesselGraphEdgeDeserializable : public Serializable {
    VGEdgeID id_;
    VGNodeID node1_;
    VGNodeID node2_;

    std::vector<VesselSkeletonVoxel> voxels_;
    VesselGraphEdgePathProperties pathProps_;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
};
#endif

// The vessel graph itself. it stores nodes as well as edges.
// References between nodes and edges are stored within the substrucutres.
//
// To avoid pointer/reference invalidation nodes and edges can only be added to the graph,
// but not removed.
#ifndef VRN_VESSELTOPOLOGY_MINIMAL_VESSELGRAPH
class VesselGraph : public Serializable {
#else
class VesselGraph {
#endif
public:
    // Move data from one graph to another
    VesselGraph(VesselGraph&& original);
    VesselGraph& operator=(VesselGraph&& other);

    // Create a deep copy from another graph
    VesselGraph clone() const;

    virtual ~VesselGraph() {}

    const VesselGraphNode& getNode(VGNodeID i) const;
    const VesselGraphEdge& getEdge(VGEdgeID i) const;
    VesselGraphNode& getNode(VGNodeID i);
    VesselGraphEdge& getEdge(VGEdgeID i);

    DiskArray<VesselGraphNode> getNodes() const;
    DiskArray<VesselGraphEdge> getEdges() const;
    DiskArray<VesselGraphNode> getNodes();
    DiskArray<VesselGraphEdge> getEdges();

    void getEdgePropertyStats(std::function<float(const VesselGraphEdge&)>, float& /*out*/ mean, float& /*out*/stddev) const;

    // Get a bounding box that encompasses all nodes within the graph.
    const tgt::Bounds& getBounds() const;

#ifndef VRN_VESSELTOPOLOGY_MINIMAL_VESSELGRAPH
    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
#endif

private:
    friend class VesselGraphBuilder;
    // Create a graph with predetermined bounds
    VesselGraph(const tgt::Bounds& bounds);
    // Create a graph with undefined bounds
    VesselGraph();

    // Note: No need to worry about pointer invalidation here because we do not store pointers
    // to other edges/nodes in nodes/edges, but only indices
    //
    // Note: We store the diskarrays in unique pointers because DiskArrayStorage cannot be moved.
    // The pointers are however guaranteed to be != null at all times.
    std::unique_ptr<DiskArrayStorage<VesselGraphNode>> nodes_; //never null
    std::unique_ptr<DiskArrayStorage<VesselGraphEdge>> edges_; //never null

    // Storage for the lists in which the nodes store their connected edges
    std::unique_ptr<DiskArrayBackedList<VGEdgeID>::Storage> nodeEdgeIdStorage_; //never null

    friend struct VesselGraphEdge;
    friend struct VesselGraphNode;
    friend struct VesselGraphEdgeDeserializable;
    std::unique_ptr<DiskArrayStorage<VesselSkeletonVoxel>> edgeVoxelStorage_; //never null
    std::unique_ptr<DiskArrayStorage<tgt::vec3>> nodeVoxelStorage_; //never null

    tgt::Bounds bounds_;
};

class VesselGraphBuilder {
public:
    // Create a graph with predetermined bounds
    VesselGraphBuilder(const tgt::Bounds& bounds);
    // Create a graph with undefined bounds
    VesselGraphBuilder();

    std::unique_ptr<VesselGraph> finalize() &&;

    const VesselGraphNode& getNode(VGNodeID i) const;

    // Insert a new node into the graph and clone it from the given existing node
    VGNodeID insertNode(const VesselGraphNode& base);
    // Insert a new node into the graph and construct it from the given parameters
    VGNodeID insertNode(const tgt::vec3& position, const DiskArray<tgt::vec3>& voxels, float radius, bool isAtSampleBorder);
    VGNodeID insertNode(const tgt::vec3& position, const DiskArray<tgt::vec3>& voxels, float radius, bool isAtSampleBorder, VesselGraphNodeUUID uuid);
    VGNodeID insertNode(const tgt::vec3& position, const std::vector<tgt::vec3>& voxels, float radius, bool isAtSampleBorder);
    VGNodeID insertNode(const tgt::vec3& position, const std::vector<tgt::vec3>& voxels, float radius, bool isAtSampleBorder, VesselGraphNodeUUID uuid);

    // Insert a new edge by deriving its properties from the provided path
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, const DiskArray<VesselSkeletonVoxel>& path);
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, const DiskArray<VesselSkeletonVoxel>& path, VesselGraphEdgeUUID uuid);
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, const std::vector<VesselSkeletonVoxel>& path);
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, const std::vector<VesselSkeletonVoxel>& path, VesselGraphEdgeUUID uuid);
    // Insert a new edge by deriving its properties from the path of the provided edge
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, const VesselGraphEdge& path_definition);
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, const VesselGraphEdge& path_definition, VesselGraphEdgeUUID uuid);
    // Insert a new edge by providing predetermined pathProperties. The path itself will be a straight line between two nodes.
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, VesselGraphEdgePathProperties pathProperties);
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, VesselGraphEdgePathProperties pathProperties, VesselGraphEdgeUUID uuid);

private:
    std::unique_ptr<VesselGraph> graph_;
};

} // namespace voreen

#endif // VRN_VESSELGRAPH_H
