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
#include "../datastructures/kdtree.h"
#include "../datastructures/diskarraystorage.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "custommodules/bigdataimageprocessing/volumefiltering/slicereader.h"
#include "custommodules/bigdataimageprocessing/datastructures/lz4slicevolume.h"
#include "../algorithm/idvolume.h"
#include "vesselgraph.h"

namespace voreen {

struct ProtoVesselGraphEdgeElement {
    typedef float CoordType;
    const tgt::Vector3<CoordType>& getPos() const {
        return pos_;
    }

    ProtoVesselGraphEdgeElement(const tgt::vec3 pos, uint64_t voxelIndex)
        : pos_(pos)
        , voxelIndex_(voxelIndex)
    {
    }
    ProtoVesselGraphEdgeElement(const ProtoVesselGraphEdgeElement& other)
        : pos_(other.pos_)
        , voxelIndex_(other.voxelIndex_)
    {
    }
    ProtoVesselGraphEdgeElement& operator=(const ProtoVesselGraphEdgeElement& other) {
        pos_ = other.pos_;
        voxelIndex_ = other.voxelIndex_;
        return *this;
    }
    tgt::vec3 pos_;
    uint64_t voxelIndex_;
};

struct ProtoVesselGraph;

struct ProtoVesselGraphEdge {
    typedef static_kdtree::Tree<ProtoVesselGraphEdgeElement, static_kdtree::SharedNodeStorage<ProtoVesselGraphEdgeElement>> ElementTree;
    ProtoVesselGraphEdge(const tgt::mat4& toRWMatrix, VGEdgeID id, VGNodeID node1, VGNodeID node2, DiskArray<tgt::svec3> voxels, ProtoVesselGraph& graph);
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
    DiskArray<tgt::svec3> voxels_;
    DiskArray<tgt::vec3> voxelsRw_;
    ElementTree tree_;
};

struct ProtoVesselGraphNode {

    ProtoVesselGraphNode(VGNodeID id, std::vector<tgt::svec3>&& voxels, bool atSampleBorder);

    VGNodeID id_;
    std::vector<tgt::svec3> voxels_;
    tgt::vec3 voxelPos_;
    bool atSampleBorder_;
    std::vector<VGEdgeID> edges_;
};

struct BranchIdVolumeReader;

struct ProtoVesselGraph {
    typedef static_kdtree::SharedMemoryTreeBuilder<ProtoVesselGraphEdgeElement> TreeBuilder;
    ProtoVesselGraph(tgt::mat4 toRWMatrix);

    VGNodeID insertNode(std::vector<tgt::svec3>&& voxels, bool atSampleBorder);
    VGEdgeID insertEdge(VGNodeID node1, VGNodeID node2, DiskArray<tgt::svec3> voxels);

    std::unique_ptr<VesselGraph> createVesselGraph(BranchIdVolumeReader& segmentedVolumeReader, const boost::optional<LZ4SliceVolume<uint8_t>>& sampleMask, ProgressReporter& progress);

    std::vector<ProtoVesselGraphNode> nodes_;
    std::vector<ProtoVesselGraphEdge> edges_;

    DiskArrayStorage<tgt::svec3> voxelStorage_;
    DiskArrayStorage<tgt::vec3> rwvoxelStorage_;

    tgt::mat4 toRWMatrix_;
    TreeBuilder treeBuilder_;
};

struct BranchIdVolumeReader {
    BranchIdVolumeReader(const LZ4SliceVolume<uint32_t>& ccaVolume)
        : branchIdReader_(ccaVolume)
    {
        branchIdReader_.seek(-1);
    }

    VGEdgeID getEdgeId(const tgt::ivec3& xypos) const {
        auto vox = branchIdReader_.getVoxel(xypos);
        tgtAssert(vox, "Invalid voxel positions");
        return *vox;
    }

    void advance() {
        branchIdReader_.advance();
    }
    bool isValidEdgeId(VGEdgeID id) const {
        return id.raw() != IdVolume::UNLABELED_FOREGROUND_VALUE && id.raw() != IdVolume::BACKGROUND_VALUE;
    }

    bool isObject(const tgt::ivec3& pos) const {
        auto voxel = branchIdReader_.getVoxel(pos);
        return voxel && *voxel != IdVolume::BACKGROUND_VALUE;
    }

    tgt::vec3 getSpacing() const {
        return branchIdReader_.getVolume().getMetaData().getSpacing();
    }

    tgt::mat4 getVoxelToWorldMatrix() const {
        return branchIdReader_.getVolume().getMetaData().getVoxelToWorldMatrix();
    }

    tgt::svec3 getDimensions() const {
        return branchIdReader_.getVolume().getDimensions();
    }


    LZ4SliceVolumeReader<uint32_t, 1> branchIdReader_;
};

}
