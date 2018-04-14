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

ProtoVesselGraphNode::ProtoVesselGraphNode(uint64_t id, std::vector<tgt::svec3>&& voxels, bool atSampleBorder)
    : id_(id)
    , voxels_(std::move(voxels))
    , atSampleBorder_(atSampleBorder)
{
}

static EdgeVoxelFinder::Builder transformVoxels(const tgt::mat4& toRWMatrix, const std::vector<tgt::svec3>& voxels) {
    EdgeVoxelFinder::Builder builder;
    for(const auto& v: voxels) {
        builder.push(ProtoVesselGraphEdgeVoxel(toRWMatrix.transform(tgt::vec3(v))));
    }
    return builder;
}

ProtoVesselGraphEdge::ProtoVesselGraphEdge(const tgt::mat4& toRWMatrix, uint64_t id, uint64_t node1, uint64_t node2, std::vector<tgt::svec3>&& voxels)
    : id_(id)
    , node1_(node1)
    , node2_(node2)
    , voxels_(std::move(voxels))
    , rwvoxels_(transformVoxels(toRWMatrix, voxels_))
{
}
std::vector<uint64_t> ProtoVesselGraphEdge::findClosestVoxelIndex(tgt::vec3 v) const {
    auto resultSet = rwvoxels_.findClosest(v);
    tgtAssert(!resultSet.empty(), "Should be: No voxels => no id!");
    return resultSet;
}

uint64_t ProtoVesselGraph::insertNode(std::vector<tgt::svec3>&& voxels, bool atSampleBorder) {
    size_t id = nodes_.size();
    nodes_.emplace_back(id, std::move(voxels), atSampleBorder);
    return id;
}
uint64_t ProtoVesselGraph::insertEdge(size_t node1, size_t node2, std::vector<tgt::svec3>&& voxels) {
    tgtAssert(node1 < nodes_.size(), "Edge references nonexistent node");
    tgtAssert(node2 < nodes_.size(), "Edge references nonexistent node");
    ProtoVesselGraphNode& n1 = nodes_.at(node1);
    ProtoVesselGraphNode& n2 = nodes_.at(node2);

    uint64_t edgeID = edges_.size();
    //edges_.push_back(ProtoVesselGraphEdge(toRWMatrix_, edgeID, node1, node2, std::move(voxels)));
    edges_.emplace_back(toRWMatrix_, edgeID, node1, node2, std::move(voxels));

    n1.edges_.push_back(edgeID);
    n2.edges_.push_back(edgeID);
    return edgeID;
}

ProtoVesselGraph::ProtoVesselGraph(tgt::mat4 toRWMatrix)
    : nodes_()
    , edges_()
    , toRWMatrix_(toRWMatrix)
{
}

std::unique_ptr<VesselGraph> ProtoVesselGraph::createVesselGraph(BranchIdVolumeReader& segmentedVolumeReader, const VolumeBase* sampleMask, ProgressReporter& progress) {
    TaskTimeLogger _("Extract edge features", tgt::Info);

    const tgt::vec3 spacing = segmentedVolumeReader.getSpacing();
    const tgt::svec3 dimensions = segmentedVolumeReader.getDimensions();
    //std::unique_ptr<VesselGraph> graph(new VesselGraph(skeleton.getBoundingBox().getBoundingBox()));
    std::unique_ptr<VesselGraph> graph(new VesselGraph());
    const tgt::mat4 toRWMatrix = segmentedVolumeReader.getVoxelToWorldMatrix();

    tgtAssert(segmentedVolumeReader.getDimensions() == dimensions, "Invalid segmentation volume dimensions");

    // For now get RAM representation for volume mask.
    // One may consider switching to a slice based approach later, but that
    // may get complicated (or would have to be employed during CCA).
    const VolumeRAM* sampleMaskRAM = sampleMask ? sampleMask->getRepresentation<VolumeRAM>() : nullptr;
    tgtAssert(!sampleMaskRAM || sampleMaskRAM->getDimensions() == dimensions, "Invalid sampleMask dimensions");

    // Create nodes
    for(const auto& node : nodes_) {
        size_t numVoxels = 0;
        tgt::vec3 sum = tgt::vec3::zero;
        std::vector<tgt::vec3> rwBrachVoxels;
        for(const auto& p : node.voxels_) {
            tgt::vec3 rwpos = toRWMatrix.transform(tgt::vec3(p));
            sum += rwpos;
            rwBrachVoxels.push_back(rwpos);
            ++numVoxels;
        }
        size_t id = graph->insertNode(sum/static_cast<float>(numVoxels), std::move(rwBrachVoxels), node.atSampleBorder_);
        tgtAssert(id == node.id_, "ID mismatch");
    }

    // Precreate VoxelSkeletonLists
    std::vector<std::vector<VesselSkeletonVoxel>> skeletonVoxelsLists;
    skeletonVoxelsLists.reserve(edges_.size());
    for(const auto& edge : edges_) {
        skeletonVoxelsLists.emplace_back();
        auto& skeletonVoxelsList = skeletonVoxelsLists.back();
        skeletonVoxelsList.reserve(edge.voxels().size());
        for(const auto& voxel : edge.voxels()) {
            skeletonVoxelsList.emplace_back(voxel.rwpos_, std::numeric_limits<float>::infinity(), 0, 0, 0, 0);
        }
    }

    const float criticalVoxelDistDiff = 1.001f*tgt::length(spacing);

    for(int z = 0; z < dimensions.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dimensions.z);
        segmentedVolumeReader.advance();

        for(int y = 0; y < dimensions.y; ++y) {
            for(int x = 0; x < dimensions.x; ++x) {

                tgt::ivec3 ipos(x, y, z);
                uint64_t id = segmentedVolumeReader.getEdgeId(ipos.xy());
                if(!segmentedVolumeReader.isValidEdgeId(id)) {
                    continue;
                }

                // Ignore values outside of sample (may happen due to inaccuracies in mask generation)
                if(sampleMaskRAM && sampleMaskRAM->getVoxelNormalized(ipos) == 0) {
                    continue;
                }

                tgt::vec3 rwVoxel = toRWMatrix.transform(tgt::vec3(x, y, z));
                std::vector<uint64_t> voxelids = edges_.at(id).findClosestVoxelIndex(rwVoxel);

                float volume = tgt::hmul(spacing) / voxelids.size();
                for(uint64_t voxelid : voxelids) {
                    VesselSkeletonVoxel& voxel = skeletonVoxelsLists.at(id).at(voxelid);

                    voxel.volume_ += volume;
                }

                // Determine if it is a surface voxel:
                // Only consider those that have a background 6-neighbor
                bool isSurfaceVoxel =
                        !( (x == 0              || segmentedVolumeReader.isObject(tgt::ivec3(x-1, y  , z  )))
                        && (x == dimensions.x-1 || segmentedVolumeReader.isObject(tgt::ivec3(x+1, y  , z  )))
                        && (y == 0              || segmentedVolumeReader.isObject(tgt::ivec3(x  , y-1, z  )))
                        && (y == dimensions.y-1 || segmentedVolumeReader.isObject(tgt::ivec3(x  , y+1, z  )))
                        && (z == 0              || segmentedVolumeReader.isObject(tgt::ivec3(x  , y  , z-1)))
                        && (z == dimensions.z-1 || segmentedVolumeReader.isObject(tgt::ivec3(x  , y  , z+1))));

                // The rest is only relevant for surface voxels, so...
                if(!isSurfaceVoxel) {
                    continue;
                }

                for(uint64_t voxelid : voxelids) {
                    VesselSkeletonVoxel& voxel = skeletonVoxelsLists.at(id).at(voxelid);

                    float dist = std::sqrt(surfaceDistanceSq(voxel.pos_, rwVoxel, spacing));

                    if(dist > voxel.maxDistToSurface_) {
                        voxel.maxDistToSurface_ = dist;
                    }
                    if(dist < voxel.minDistToSurface_) {
                        voxel.minDistToSurface_ = dist;
                    }
                    voxel.avgDistToSurface_ = (voxel.avgDistToSurface_*voxel.numSurfaceVoxels_ + dist)/(voxel.numSurfaceVoxels_+1);
                    ++voxel.numSurfaceVoxels_;
                }
            }
        }
    }
    progress.setProgress(1.0f);

    // Create edges from branches
    for(size_t i = 0; i < skeletonVoxelsLists.size(); ++i) {
        auto& skeletonVoxelsList = skeletonVoxelsLists.at(i);
        auto& protoEdge = edges_.at(i);

        graph->insertEdge(protoEdge.node1_, protoEdge.node2_, std::move(skeletonVoxelsList));
    }

    return graph;
}
const uint64_t BranchIdVolumeReader::INVALID_EDGE_ID = -1;

}
