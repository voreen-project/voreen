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

#pragma once
#include <vector>
#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/memory.h"

#include "voreen/core/datastructures/diskarraystorage.h"
#include "voreen/core/datastructures/volume/volume.h"

#include "modules/bigdataimageprocessing/volumefiltering/slicereader.h"
#include "modules/bigdataimageprocessing/datastructures/lz4slicevolume.h"

#include "../algorithm/idvolume.h"
#include "../datastructures/kdtree.h"
#include "vesselgraph.h"

namespace voreen {

struct ProtoVesselGraphEdgeElement {
    typedef float CoordType;
    const tgt::Vector3<CoordType>& getPos() const {
        return pos_;
    }

    ProtoVesselGraphEdgeElement(const tgt::vec3 pos, uint64_t voxelIndex)
        : pos_(pos)
        , padding_(0.0)
        , voxelIndex_(voxelIndex)
    {
    }
    ProtoVesselGraphEdgeElement(const ProtoVesselGraphEdgeElement& other)
        : pos_(other.pos_)
        , padding_(0.0)
        , voxelIndex_(other.voxelIndex_)
    {
    }
    ProtoVesselGraphEdgeElement& operator=(const ProtoVesselGraphEdgeElement& other) {
        pos_ = other.pos_;
        voxelIndex_ = other.voxelIndex_;
        return *this;
    }
    tgt::vec3 pos_;
    float padding_; // Avoid invalid write with uninitialized memory warnings in valgrind
    uint64_t voxelIndex_;
};

struct ProtoVesselGraph;
struct ProtoVesselGraphNode;

struct ProtoVesselGraphEdge {
    typedef static_kdtree::Tree<ProtoVesselGraphEdgeElement, static_kdtree::SharedNodeStorage<ProtoVesselGraphEdgeElement>> ElementTree;
    ProtoVesselGraphEdge(const tgt::mat4& toRWMatrix, VGEdgeID id, const ProtoVesselGraphNode& node1, const ProtoVesselGraphNode& node2, const DiskArray<tgt::svec3>& voxels, ProtoVesselGraph& graph);
    static_kdtree::SearchNearestResultSet<ProtoVesselGraphEdgeElement> findClosestVoxelIndex(tgt::vec3) const;
    DiskArray<tgt::vec3>& voxels() {
        return voxelsRw_;
    }
    const DiskArray<tgt::vec3>& voxels() const {
        return voxelsRw_;
    }

    VGEdgeID id_;
    VGNodeID node1_;
    VGNodeID node2_;
    DiskArray<tgt::svec3> voxels_; //Voxel position of the original centerline
    DiskArray<tgt::vec3> voxelsRw_; //Voxel position transformed into real world
    DiskArray<tgt::vec3> voxelsRwSmooth_; //Smoothed centerline in real world coordinates
    ElementTree tree_;
};

struct ProtoVesselGraphNode {

    ProtoVesselGraphNode(VGNodeID id, DiskArray<tgt::svec3>&& voxels, bool atSampleBorder, ProtoVesselGraph& graph);

    VGNodeID id_;
    DiskArray<tgt::svec3> voxels_;
    tgt::vec3 voxelPos_;
    bool atSampleBorder_;
    DiskArrayBackedList<VGEdgeID> edges_;
};

struct ProtoVesselGraph {
    typedef static_kdtree::SharedMemoryTreeBuilder<ProtoVesselGraphEdgeElement> TreeBuilder;
    ProtoVesselGraph(tgt::mat4 toRWMatrix);

    VGNodeID insertNode(const std::vector<tgt::svec3>& voxels, bool atSampleBorder);
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, const DiskArray<tgt::svec3>& voxels);

    std::unique_ptr<VesselGraph> createVesselGraph(const LZ4SliceVolume<uint32_t>& ccaVolume, const boost::optional<LZ4SliceVolume<uint8_t>>& sampleMask, ProgressReporter& progress);

    DiskArrayStorage<ProtoVesselGraphNode> nodes_;
    DiskArrayStorage<ProtoVesselGraphEdge> edges_;

    DiskArrayStorage<tgt::svec3> voxelStorage_;
    DiskArrayStorage<tgt::vec3> rwvoxelStorage_;

    // Storage for the lists in which the nodes store their connected edges
    DiskArrayBackedList<VGEdgeID>::Storage nodeEdgeIdStorage_;


    tgt::mat4 toRWMatrix_;
    TreeBuilder treeBuilder_;
};

}
