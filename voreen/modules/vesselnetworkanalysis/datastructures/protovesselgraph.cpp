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

#include <memory>

#include "protovesselgraph.h"
#include "../util/tasktimelogger.h"
#include "voreen/core/io/progressreporter.h"
#include "../algorithm/idvolume.h"

#include "tgt/assert.h"

namespace {
inline static float surfaceDistanceSq(const tgt::vec3& skeletonVoxel, const tgt::vec3 surfaceVoxel, const tgt::vec3 spacing) {
    //return tgt::distanceSq(skeletonVoxel, surfaceVoxel);
    float avgToSurfaceDist = (spacing.x + spacing.y + spacing.z)/(3 /*average of three components */ * 2 /*surface is somewhere in the middle*/);
    float dist = tgt::distance(skeletonVoxel, surfaceVoxel) + avgToSurfaceDist;
    return dist*dist;
    //tgt::vec3 fromTo = surfaceVoxel - skeletonVoxel;
    //if(tgt::lengthSq(fromTo) == 0) {
    //    float expectedDegenerateDist = (spacing.x + spacing.y + spacing.z)/(3 /*average of three components */ * 2 /*surface is somewhere in the middle*/);
    //    return expectedDegenerateDist*expectedDegenerateDist;
    //}
    //tgt::vec3 fromToNormalized = tgt::normalize(fromTo);
    //fromTo += fromToNormalized*spacing*0.5f;
    //return tgt::lengthSq(fromTo);
}
}

namespace voreen {

ProtoVesselGraphNode::ProtoVesselGraphNode(VGNodeID id, DiskArray<tgt::svec3>&& voxels, bool atSampleBorder, ProtoVesselGraph& graph)
    : id_(id)
    , voxels_(std::move(voxels))
    , atSampleBorder_(atSampleBorder)
    , edges_(graph.nodeEdgeIdStorage_)
{
    tgtAssert(!voxels_.empty(), "No voxels");
    tgt::svec3 voxelSum = std::accumulate(voxels_.begin(), voxels_.end(), tgt::svec3::zero);
    voxelPos_ = tgt::vec3(voxelSum)/static_cast<float>(voxels_.size());
}

static DiskArray<tgt::vec3> transformVoxels(const tgt::mat4& toRWMatrix, const DiskArray<tgt::svec3>& voxels, DiskArrayStorage<tgt::vec3>& storage) {
    std::vector<tgt::vec3> vec; //TODO: somehow avoid copy here?
    for(const auto& v: voxels) {
        vec.push_back(toRWMatrix.transform(tgt::vec3(v)));
    }
    return storage.store(vec);
}
static DiskArray<tgt::vec3> smoothVoxels(const DiskArray<tgt::vec3>& voxels, tgt::vec3 node1PosRw, tgt::vec3 node2PosRw, DiskArrayStorage<tgt::vec3>& storage) {
    std::vector<tgt::vec3> controlPoints;
    controlPoints.push_back(node1PosRw);
    for(const auto& v: voxels) {
        controlPoints.push_back(v);
    }
    controlPoints.push_back(node2PosRw);


    std::vector<tgt::vec3> vec;
    int numPoints = voxels.size();
    for(int i=1; i<=numPoints; ++i) {
        // No smoothing
        //tgt::vec3 result = controlPoints[i];

        // Local Bezier curve of degree 2 (i.e., mix the neighbors):
        //tgt::vec3 result = 0.5f*controlPoints[i];
        //result += 0.25f*controlPoints[i-1];
        //result += 0.25f*controlPoints[i+1];

        // Local Bezier curve of degree 4 (i.e., mix the neighbors):
        tgt::vec3 result = 6.0f*controlPoints[i];
        result += 4.0f*controlPoints[i-1];
        result += 4.0f*controlPoints[i+1];
        result += controlPoints[std::max(i-2, 0)];
        result += controlPoints[std::min(i+2, numPoints+1)];
        result /= 16.0f;

        vec.push_back(result);
    }

    return storage.store(vec);
}

static ProtoVesselGraphEdge::ElementTree buildTree(const DiskArray<tgt::vec3>& voxels, static_kdtree::SharedMemoryTreeBuilder<ProtoVesselGraphEdgeElement>& treeBuilder) {
    voreen::static_kdtree::ElementArrayBuilder<ProtoVesselGraphEdgeElement> builder(VoreenApplication::app()->getUniqueTmpFilePath(".kdtreestorage"));
    uint64_t i = 0;
    for(const auto& v: voxels) {
        builder.push(ProtoVesselGraphEdgeElement(v, i));
        ++i;
    }
    return treeBuilder.buildTree(std::move(builder));
}

ProtoVesselGraphEdge::ProtoVesselGraphEdge(const tgt::mat4& toRWMatrix, VGEdgeID id, const ProtoVesselGraphNode& node1, const ProtoVesselGraphNode& node2, const DiskArray<tgt::svec3>& voxels, ProtoVesselGraph& graph)
    : id_(id)
    , node1_(node1.id_)
    , node2_(node2.id_)
    , voxels_(graph.voxelStorage_.store(voxels))
    , voxelsRw_(transformVoxels(toRWMatrix, voxels_, graph.rwvoxelStorage_))
    , voxelsRwSmooth_(smoothVoxels(voxelsRw_, toRWMatrix.transform(node1.voxelPos_), toRWMatrix.transform(node2.voxelPos_), graph.rwvoxelStorage_))
    , tree_(buildTree(voxelsRwSmooth_, graph.treeBuilder_))
{
}
static_kdtree::SearchNearestResultSet<ProtoVesselGraphEdgeElement> ProtoVesselGraphEdge::findClosestVoxelIndex(tgt::vec3 v) const {
    auto resultSet = tree_.findAllNearest(v);
    tgtAssert(!resultSet.elements_.empty(), "Should be: No voxels => no id!");
    return resultSet;
}

VGNodeID ProtoVesselGraph::insertNode(const std::vector<tgt::svec3>& voxels, bool atSampleBorder) {
    size_t id = nodes_.size();
    size_t storageId = nodes_.storeElement(ProtoVesselGraphNode(id, voxelStorage_.store(voxels), atSampleBorder, *this));
    tgtAssert(id == storageId, "invalid storage id");
    return id;
}
VGEdgeID ProtoVesselGraph::insertEdge(VGNodeID node1, VGNodeID node2, const DiskArray<tgt::svec3>& voxels) {
    tgtAssert(node1 < nodes_.size(), "Edge references nonexistent node");
    tgtAssert(node2 < nodes_.size(), "Edge references nonexistent node");
    ProtoVesselGraphNode& n1 = nodes_[node1.raw()];
    ProtoVesselGraphNode& n2 = nodes_[node2.raw()];

    VGEdgeID edgeID = edges_.size();
    //edges_.push_back(ProtoVesselGraphEdge(toRWMatrix_, edgeID, node1, node2, std::move(voxels)));
    size_t storageId = edges_.storeElement(ProtoVesselGraphEdge(toRWMatrix_, edgeID, n1, n2, voxels, *this));
    tgtAssert(edgeID.raw() == storageId, "invalid storage id");

    n1.edges_.push(edgeID);
    n2.edges_.push(edgeID);
    return edgeID;
}

ProtoVesselGraph::ProtoVesselGraph(tgt::mat4 toRWMatrix)
    : nodes_(VoreenApplication::app()->getUniqueTmpFilePath(".protonodes"))
    , edges_(VoreenApplication::app()->getUniqueTmpFilePath(".protoedges"))
    , toRWMatrix_(toRWMatrix)
    , treeBuilder_(VoreenApplication::app()->getUniqueTmpFilePath(".kdtrees"))
    , voxelStorage_(VoreenApplication::app()->getUniqueTmpFilePath(".voxelstorage"))
    , rwvoxelStorage_(VoreenApplication::app()->getUniqueTmpFilePath(".rwvoxelstorage"))
    , nodeEdgeIdStorage_(VoreenApplication::app()->getUniqueTmpFilePath(".nodeedgeidrefs"))
{
}

struct ProtoNodeRef {
    const ProtoVesselGraphNode& node_;
    float radius_;
    tgt::vec3 rwPos_;

    ProtoNodeRef(const ProtoVesselGraphNode& node, const tgt::mat4& voxelToRw)
        : node_(node)
        , radius_(0)
        , rwPos_(voxelToRw.transform(node.voxelPos_))
    {
    }
};

const int MAX_CHUNK_SIZE = 8;

typedef OverlappingSlabReader<uint32_t, MAX_CHUNK_SIZE, 1> BranchIdVolumeReader;

static bool isValidEdgeId(VGEdgeID id) {
    return id.raw() != IdVolume::UNLABELED_FOREGROUND_VALUE && id.raw() != IdVolume::BACKGROUND_VALUE;
}

static bool isObject(const BranchIdVolumeReader& r, const tgt::ivec3& pos) {
    auto voxel = r.getVoxel(pos);
    return voxel && *voxel != IdVolume::BACKGROUND_VALUE;
}

static VGEdgeID getEdgeId(const BranchIdVolumeReader& r, const tgt::ivec3& pos) {
    auto voxel = r.getVoxel(pos);
    tgtAssert(voxel, "Invalid voxel position");
    return *voxel;
}

static boost::optional<VGEdgeID> findNearOtherEdgeID(const BranchIdVolumeReader& segmentedVolumeReader, const tgt::ivec3& dimensions, const tgt::ivec3& currentPos, VGEdgeID centerID) {
    tgt::ivec3 minPos = tgt::max(tgt::ivec3::zero, currentPos - tgt::ivec3(1, 1, 1));
    tgt::ivec3 maxPos = tgt::min(dimensions, currentPos + tgt::ivec3(2,2,2));
    for(int z = minPos.z; z < maxPos.z; ++z) {
        for(int y = minPos.y; y < maxPos.y; ++y) {
            for(int x = minPos.x; x < maxPos.x; ++x) {
                tgt::ivec3 pos(x,y,z);
                if (pos == currentPos) {
                    continue;
                }
                VGEdgeID id = getEdgeId(segmentedVolumeReader, pos);
                if(id != centerID && isValidEdgeId(id)) {
                    return id;
                }
            }
        }
    }
    return boost::none;
}

std::unique_ptr<VesselGraph> ProtoVesselGraph::createVesselGraph(const LZ4SliceVolume<uint32_t>& ccaVolume, const boost::optional<LZ4SliceVolume<uint8_t>>& sampleMask, ProgressReporter& progress) {
    TaskTimeLogger _("Extract edge features", tgt::Info);


    const tgt::svec3 dimensions = ccaVolume.getDimensions();
    const tgt::vec3 spacing = ccaVolume.getMetaData().getSpacing();
    const tgt::mat4 toRWMatrix = ccaVolume.getMetaData().getVoxelToWorldMatrix();

    BranchIdVolumeReader segmentedVolumeReader(ccaVolume);

    tgtAssert(!sampleMask || sampleMask->getDimensions() == dimensions, "Invalid segmentation volume dimensions");

    DiskArrayStorage<ProtoNodeRef> nodeRefs(VoreenApplication::app()->getUniqueTmpFilePath(".protonoderefs"));
    for(auto& node : nodes_.asArray()) {
        nodeRefs.storeElement(ProtoNodeRef(node, toRWMatrix));
    }

    // Precreate VoxelSkeletonLists
    DiskArrayStorage<VesselSkeletonVoxel> tmpSkeletonVoxelListStorage(VoreenApplication::app()->getUniqueTmpFilePath(".voxellists"));

    DiskArrayStorage<DiskArray<VesselSkeletonVoxel>> skeletonVoxelLists(VoreenApplication::app()->getUniqueTmpFilePath(".skelvoxlists"));

    for(const auto& edge : edges_.asArray()) {
        auto builder = tmpSkeletonVoxelListStorage.build();
        for(const auto& voxel : edge.voxelsRwSmooth_) {
            builder.push(VesselSkeletonVoxel(voxel, std::numeric_limits<float>::infinity(), 0, 0, 0, 0, false));
        }
        skeletonVoxelLists.storeElement(DiskArray<VesselSkeletonVoxel>(std::move(builder).finalize()));
    }

    auto calc_num_chunks = [] (int dim, int size) {
        if((dim % size) == 0) {
            return dim/size;
        } else {
            return dim/size+1;
        }
    };

    auto calc_chunk_size = [&] (int chunk, int num_chunks, int dimensions) {
        int size = (chunk != num_chunks - 1 || (dimensions % MAX_CHUNK_SIZE == 0)) ? MAX_CHUNK_SIZE : dimensions % MAX_CHUNK_SIZE;
        tgtAssert(size > 0, "Invalid chunk size");
        return size;
    };

    const int num_chunks_x = calc_num_chunks(dimensions.x, MAX_CHUNK_SIZE);
    const int num_chunks_y = calc_num_chunks(dimensions.y, MAX_CHUNK_SIZE);
    const int num_chunks_z = calc_num_chunks(dimensions.z, MAX_CHUNK_SIZE);

    tgt::ivec3 idim = dimensions;
    for(int chunk_z = 0; chunk_z < num_chunks_z; ++chunk_z) {
        progress.setProgress(static_cast<float>(chunk_z)/num_chunks_z);

        const int chunk_size_z = calc_chunk_size(chunk_z, num_chunks_z, dimensions.z);

        int slab_begin = chunk_z*MAX_CHUNK_SIZE;
        int slab_end = slab_begin + chunk_size_z;

        boost::optional<VolumeAtomic<uint8_t>> sampleMaskSlab = sampleMask ? boost::optional<VolumeAtomic<uint8_t>>(sampleMask->loadSlab(slab_begin, slab_end)) : boost::none;

        segmentedVolumeReader.advance();
        tgtAssert(slab_begin == segmentedVolumeReader.currentSlabStart_, "Invalid start");

        for(int chunk_y = 0; chunk_y < num_chunks_y; ++chunk_y) {
            const int chunk_size_y = calc_chunk_size(chunk_y, num_chunks_y, dimensions.y);

            for(int chunk_x = 0; chunk_x < num_chunks_x; ++chunk_x) {
                const int chunk_size_x = calc_chunk_size(chunk_x, num_chunks_x, dimensions.x);
                for(int dz = 0; dz<chunk_size_z; ++dz) {
                    int z = chunk_z*MAX_CHUNK_SIZE + dz;

                    for(int dy = 0; dy<chunk_size_y; ++dy) {
                        int y = chunk_y*MAX_CHUNK_SIZE + dy;

                        for(int dx = 0; dx<chunk_size_x; ++dx) {
                            int x = chunk_x*MAX_CHUNK_SIZE + dx;

                            tgt::ivec3 ipos(x, y, z);
                            VGEdgeID id = getEdgeId(segmentedVolumeReader, ipos);
                            if(!isValidEdgeId(id)) {
                                continue;
                            }

                            // Ignore values outside of sample (may happen due to inaccuracies in mask generation)
                            if(sampleMaskSlab && sampleMaskSlab->voxel(x, y, 0) == 0) {
                                continue;
                            }

                            tgt::vec3 rwVoxel = toRWMatrix.transform(tgt::vec3(x, y, z));
                            auto neared_result = edges_[id.raw()].findClosestVoxelIndex(rwVoxel);

                            float volume = tgt::hmul(spacing) / neared_result.elements_.size();
                            for(auto element : neared_result.elements_) {
                                VesselSkeletonVoxel& voxel = skeletonVoxelLists[id.raw()].at(element->voxelIndex_);

                                voxel.volume_ += volume;
                            }

                            // Determine if it is a surface voxel:
                            // Only consider those that have a background 6-neighbor
                            bool isSurfaceVoxel =
                                !( (x == 0        || isObject(segmentedVolumeReader, tgt::ivec3(x-1, y  , z  )))
                                && (x == idim.x-1 || isObject(segmentedVolumeReader, tgt::ivec3(x+1, y  , z  )))
                                && (y == 0        || isObject(segmentedVolumeReader, tgt::ivec3(x  , y-1, z  )))
                                && (y == idim.y-1 || isObject(segmentedVolumeReader, tgt::ivec3(x  , y+1, z  )))
                                && (z == 0        || isObject(segmentedVolumeReader, tgt::ivec3(x  , y  , z-1)))
                                && (z == idim.z-1 || isObject(segmentedVolumeReader, tgt::ivec3(x  , y  , z+1))));

                            // The rest is only relevant for surface voxels, so...
                            if(!isSurfaceVoxel) {
                                continue;
                            }

                            auto otherNearEdgeID = findNearOtherEdgeID(segmentedVolumeReader, idim, ipos, id);

                            for(auto element : neared_result.elements_) {
                                VesselSkeletonVoxel& voxel = skeletonVoxelLists[id.raw()].at(element->voxelIndex_);

                                float dist = std::sqrt(surfaceDistanceSq(voxel.pos_, rwVoxel, spacing));
                                tgtAssert(dist >= 0, "Invalid voxel distance");

                                if(dist > voxel.maxDistToSurface_) {
                                    voxel.maxDistToSurface_ = dist;
                                }
                                if(dist < voxel.minDistToSurface_) {
                                    voxel.minDistToSurface_ = dist;
                                }
                                voxel.avgDistToSurface_ = (voxel.avgDistToSurface_*voxel.numSurfaceVoxels_ + dist)/(voxel.numSurfaceVoxels_+1);
                                ++voxel.numSurfaceVoxels_;

                                if(otherNearEdgeID) {
                                    voxel.nearOtherEdge_ = true;
                                }
                            }

                            // Determine radius for nodes
                            if(otherNearEdgeID) {
                                auto& edge1 = edges_[id.raw()];
                                auto& edge2 = edges_[otherNearEdgeID->raw()];
                                std::vector<VGNodeID> candidates;
                                if(edge1.node1_ == edge2.node1_) { candidates.push_back(edge1.node1_); }
                                if(edge1.node1_ == edge2.node2_) { candidates.push_back(edge1.node1_); }
                                if(edge1.node2_ == edge2.node1_) { candidates.push_back(edge1.node2_); }
                                if(edge1.node2_ == edge2.node2_) { candidates.push_back(edge1.node2_); }

                                VGNodeID best = -1;
                                float bestDistSq = std::numeric_limits<float>::infinity();
                                for(VGNodeID candidate : candidates) {
                                    float distSq = tgt::distanceSq(nodeRefs[candidate.raw()].rwPos_, rwVoxel);
                                    if(distSq < bestDistSq) {
                                        bestDistSq = distSq;
                                        best = candidate;
                                    }
                                }

                                if(best != -1) {
                                    ProtoNodeRef& node = nodeRefs[best.raw()];
                                    node.radius_ = std::max(node.radius_, std::sqrt(bestDistSq));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    progress.setProgress(1.0f);

    VesselGraphBuilder builder;

    // Create nodes
    for(const ProtoNodeRef& nodeRef : nodeRefs.asArray()) {
        const ProtoVesselGraphNode& node = nodeRef.node_;
        std::vector<tgt::vec3> rwBrachVoxels;
        for(const auto& p : node.voxels_) {
            tgt::vec3 rwpos = toRWMatrix.transform(tgt::vec3(p));
            rwBrachVoxels.push_back(rwpos);
        }
        tgtAssert(!std::isnan(tgt::hmul(nodeRef.rwPos_)), "Invalid pos");
        VGNodeID id = builder.insertNode(nodeRef.rwPos_, std::move(rwBrachVoxels), nodeRef.radius_, node.atSampleBorder_);
        tgtAssert(id == node.id_, "ID mismatch");
    }


    // Create edges from branches
    for(size_t i = 0; i < skeletonVoxelLists.size(); ++i) {
        auto& skeletonVoxelList = skeletonVoxelLists[i];
        auto& protoEdge = edges_[i];

        builder.insertEdge(protoEdge.node1_, protoEdge.node2_, std::move(skeletonVoxelList));
    }

    return std::move(builder).finalize();
}

}
