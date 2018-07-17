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

#pragma once
#include <vector>
#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/memory.h"
#include "../util/kdtreebuilder.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "custommodules/bigdataimageprocessing/volumefiltering/slicereader.h"
#include "custommodules/bigdataimageprocessing/datastructures/lz4slicevolume.h"
#include "../algorithm/idvolume.h"
#include "vesselgraph.h"

namespace voreen {

struct ProtoVesselGraphEdgeVoxel {
    typedef tgt::vec3 VoxelType;
    VoxelType rwpos_;
    ProtoVesselGraphEdgeVoxel(VoxelType v)
        : rwpos_(v)
    {
    }

    inline typename VoxelType::ElemType at(int dim) const {
        return rwpos_[dim];
    }
};

struct ProtoVoxelGraphEdgeVoxelSet {
    typedef uint64_t IndexType;
    typedef float DistanceType;

    std::vector<uint64_t> closestIDs_;
    float currentDist_;

    ProtoVoxelGraphEdgeVoxelSet()
        : closestIDs_()
        , currentDist_(std::numeric_limits<float>::infinity())
    {
    }

    inline bool found() {
        return currentDist_ < std::numeric_limits<float>::infinity();
    }

    inline void init() {
        clear();
    }

    inline void clear() {
        closestIDs_.clear();
    }

    inline size_t size() const {
        return closestIDs_.size();
    }

    inline bool full() const {
        //Not sure what to do here...
        //This is analogous to the implementation of the max radius set of nanoflann
        return true;
    }

    inline DistanceType worstDist() const {
        return currentDist_*1.001;
    }

    /**
     * Find the worst result (furtherest neighbor) without copying or sorting
     * Pre-conditions: size() > 0
     */
    std::pair<IndexType,DistanceType> worst_item() const {
        if (closestIDs_.empty()) throw std::runtime_error("Cannot invoke RadiusResultSet::worst_item() on an empty list of results.");
        uint64_t id = closestIDs_[0];
        return std::make_pair<IndexType,DistanceType>(std::move(id), worstDist());
    }

    /**
     * Called during search to add an element matching the criteria.
     * @return true if the search should be continued, false if the results are sufficient
     */
    inline bool addPoint(float dist, size_t index)
    {
        if(dist <= currentDist_) {
            if(dist < currentDist_) {
                closestIDs_.clear();
                currentDist_ = dist;
            }
            closestIDs_.push_back(index);
        }
        return true;
    }
};

struct EdgeVoxelFinder {
    typedef KDTreeBuilder<ProtoVesselGraphEdgeVoxel> Builder;
    typedef nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<typename Builder::coord_t, Builder>,
		Builder,
		3 /* dim */
		> Index;

    EdgeVoxelFinder(Builder&& builder)
        : storage_(std::move(builder))
        , index_(3 /*dim */, storage_, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* recommended as a sensible value by library creator */))
    {
        index_.buildIndex();
    }

    EdgeVoxelFinder(EdgeVoxelFinder&& other)
        : EdgeVoxelFinder(std::move(other.storage_))
    {
    }


    /**
     * Find the closest skeleton voxel(s)
     */
    std::vector<uint64_t> findClosest(tgt::vec3 rwPos) const {
        nanoflann::SearchParams params;
        params.sorted = false;

        ProtoVoxelGraphEdgeVoxelSet set;
        index_.radiusSearchCustomCallback(rwPos.elem, set, params);

        return std::move(set.closestIDs_);
    }

    Builder storage_;
    Index index_;
};

struct ProtoVesselGraphEdge {
    ProtoVesselGraphEdge(const tgt::mat4& toRWMatrix, uint64_t id, uint64_t node1, uint64_t node2, std::vector<tgt::svec3>&& voxels);
    std::vector<uint64_t> findClosestVoxelIndex(tgt::vec3) const;
    std::vector<ProtoVesselGraphEdgeVoxel>& voxels() {
        return rwvoxels_.storage_.points();
    }
    const std::vector<ProtoVesselGraphEdgeVoxel>& voxels() const {
        return rwvoxels_.storage_.points();
    }

    uint64_t id_;
    uint64_t node1_;
    uint64_t node2_;
    std::vector<tgt::svec3> voxels_;
    EdgeVoxelFinder rwvoxels_;
};

struct ProtoVesselGraphNode {

    ProtoVesselGraphNode(uint64_t id, std::vector<tgt::svec3>&& voxels, bool atSampleBorder);

    uint64_t id_;
    std::vector<tgt::svec3> voxels_;
    tgt::vec3 voxelPos_;
    bool atSampleBorder_;
    std::vector<uint64_t> edges_;
};

struct BranchIdVolumeReader;

struct ProtoVesselGraph {
    ProtoVesselGraph(tgt::mat4 toRWMatrix);

    uint64_t insertNode(std::vector<tgt::svec3>&& voxels, bool atSampleBorder);
    uint64_t insertEdge(size_t node1, size_t node2, std::vector<tgt::svec3>&& voxels);

    std::unique_ptr<VesselGraph> createVesselGraph(BranchIdVolumeReader& segmentedVolumeReader, const boost::optional<LZ4SliceVolume<uint8_t>>& sampleMask, ProgressReporter& progress);

    std::vector<ProtoVesselGraphNode> nodes_;
    std::vector<ProtoVesselGraphEdge> edges_;
    tgt::mat4 toRWMatrix_;

};

struct BranchIdVolumeReader {
    BranchIdVolumeReader(const LZ4SliceVolume<uint32_t>& ccaVolume)
        : branchIdReader_(ccaVolume)
    {
        branchIdReader_.seek(-1);
    }

    uint64_t getEdgeId(const tgt::ivec3& xypos) const {
        auto vox = branchIdReader_.getVoxel(xypos);
        tgtAssert(vox, "Invalid voxel positions");
        return *vox;
    }

    void advance() {
        branchIdReader_.advance();
    }
    bool isValidEdgeId(uint64_t id) const {
        uint32_t id32 = static_cast<uint32_t>(id);
        return id32 != IdVolume::UNLABELED_FOREGROUND_VALUE && id32 != IdVolume::BACKGROUND_VALUE;
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
