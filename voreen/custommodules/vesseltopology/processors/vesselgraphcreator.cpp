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

#include "vesselgraphcreator.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumecontainer.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"
#include "modules/hdf5/io/hdf5filevolume.h"

#include "../algorithm/streaminggraphcreation.h"
#include "../algorithm/volumemask.h"
#include "../algorithm/idvolume.h"
#include "../algorithm/vesselgraphnormalization.h"
#include "../algorithm/boundshierarchy.h"

#include "custommodules/bigdataimageprocessing/processors/connectedcomponentanalysis.h"

#include "tgt/bounds.h"
#include "tgt/filesystem.h"

#include <vector>
#include <memory>
#include <unordered_set>

namespace voreen {

const std::string VesselGraphCreator::loggerCat_("voreen.vesseltopology.vesselgraphcreator");

VesselGraphCreator::VesselGraphCreator()
    : segmentedVolumeInport_(Port::INPORT, "vesselgraphcreator.segmentedVolume.inport", "Segmentation Volume")
    , sampleMaskInport_(Port::INPORT, "vesselgraphcreator.samplemask.inport", "Sample Mask (optional)")
    , fixedForegroundPointInport_(Port::INPORT, "vesselgraphcreator.fixedForegroundPointInport", "Fixed Foreground Points", false, Processor::INVALID_RESULT)
    , graphOutport_(Port::OUTPORT, "vesselgraphcreator_graph.outport", "Graph", false, Processor::VALID)
    , nodeOutport_(Port::OUTPORT, "vesselgraphcreator_node.outport", "Nodes Voxels", false, Processor::VALID)
    , edgeOutport_(Port::OUTPORT, "vesselgraphcreator_edge.outport", "Edges Voxels", false, Processor::VALID)
    , generatedSkeletonsOutport_(Port::OUTPORT, "generatedSkeletons.outport", "Generated Skeletons", false, Processor::VALID)
    , numRefinementIterations_("numRefinementIterations", "Refinement Iterations", 0, 0, 10, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_ADVANCED)
    , minVoxelLength_("minVoxelLength", "Min Voxel Length", 0, 0, 50)
    , minElongation_("minElongation", "Minimum Elongation", 0, 0, 5)
    , minBulgeSize_("minBulgeSize", "Minimum Bulge Size", 0, 0, 10)
    , binarizationThresholdSegmentation_("binarizationThresholdSegmentation", "Binarization Threshold (Segmentation)", 0.5f, 0.0f, 1.0f, Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
{
    addPort(segmentedVolumeInport_);
    addPort(sampleMaskInport_);
    addPort(fixedForegroundPointInport_);
    addPort(graphOutport_);
    addPort(nodeOutport_);
    addPort(edgeOutport_);
    addPort(generatedSkeletonsOutport_);

        addProperty(numRefinementIterations_);
            numRefinementIterations_.setGroupID("refinement");
        addProperty(minVoxelLength_);
            minVoxelLength_.setGroupID("refinement");
        addProperty(minElongation_);
            minElongation_.setGroupID("refinement");
        addProperty(minBulgeSize_);
            minBulgeSize_.setGroupID("refinement");
    setPropertyGroupGuiName("refinement", "Refinement");

        addProperty(binarizationThresholdSegmentation_);
            binarizationThresholdSegmentation_.setGroupID("binarization");
    setPropertyGroupGuiName("binarization", "Binarization");
}

VesselGraphCreator::~VesselGraphCreator() {
}
VoreenSerializableObject* VesselGraphCreator::create() const {
    return new VesselGraphCreator();
}

bool VesselGraphCreator::isReady() const {
    return isInitialized() && segmentedVolumeInport_.isReady() &&
        (graphOutport_.isReady() || nodeOutport_.isReady() || edgeOutport_.isReady() || generatedSkeletonsOutport_.isReady());
}

struct SparseBinaryVolume {
    size_t linearPos(const tgt::svec3& pos) {
        return pos.x + dimensions_.x*(pos.y + dimensions_.y*pos.z);
    }
    SparseBinaryVolume(const std::vector<tgt::svec3>& voxels, const tgt::svec3& dimensions)
        : linearVoxelPositions_()
        , dimensions_(dimensions)
    {
        for(auto& voxel : voxels) {
            linearVoxelPositions_.emplace_back(linearPos(voxel));
        }
        std::sort(linearVoxelPositions_.begin(), linearVoxelPositions_.end());
    }

    bool isObject(const tgt::svec3& pos) {
        size_t posLinear = linearPos(pos);
        return std::binary_search(linearVoxelPositions_.begin(), linearVoxelPositions_.end(), posLinear);
    }
    std::vector<size_t> linearVoxelPositions_;
    tgt::svec3 dimensions_;
};

struct EndNodeVoxelExtractor {
    static void addVoxels(std::vector<tgt::svec3>& voxels, const VesselGraphNode& node, const tgt::mat4& rwToVoxel) {
        if(node.isEndNode()) {
            voxels.push_back(tgt::round(transform(rwToVoxel, node.pos_)));
        }
    }
};

struct GraphNodeVoxelReader {
    template<class S>
    static GraphNodeVoxelReader create(const VesselGraph& graph, const tgt::mat4& rwToVoxel, const tgt::svec3& dimensions) {
        std::vector<tgt::svec3> voxels;

        for(auto& node : graph.getNodes()) {
            S::addVoxels(voxels, node, rwToVoxel);
        }
        return GraphNodeVoxelReader(SparseBinaryVolume(voxels, dimensions));
    }

    GraphNodeVoxelReader(SparseBinaryVolume&& vol)
        : nodeVoxels_(std::move(vol))
        , z_(-1)
    {
    }
    void advance() {
        ++z_;
    }
    bool isForeground(size_t x, size_t y) {
        return nodeVoxels_.isObject(tgt::svec3(x, y, static_cast<size_t>(z_)));
    }
    SparseBinaryVolume nodeVoxels_;
    int z_;
};

struct EdgeVoxelRef {
    typedef tgt::vec3 VoxelType;

    VoxelType rwPos; // Cached for better performance
    const ProtoVesselGraphEdge& edge;

    inline typename VoxelType::ElemType at(int dim) const {
        return rwPos[dim];
    }

    EdgeVoxelRef(VoxelType rwPos, const ProtoVesselGraphEdge& edge)
        : rwPos(rwPos)
        , edge(edge)
    {
    }
};

struct EdgeVoxelRefResultSetBase {
    typedef size_t IndexType;
    typedef float DistanceType;

    const std::vector<EdgeVoxelRef>& storage_;
    std::vector<size_t> closestIDs_;

    EdgeVoxelRefResultSetBase(const std::vector<EdgeVoxelRef>& storage)
        : storage_(storage)
        , closestIDs_()
    {
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

    std::vector<const EdgeVoxelRef*> getVoxels() const {
        std::vector<const EdgeVoxelRef*> result;
        for(auto id : closestIDs_) {
            result.push_back(&getEdgeVoxelRef(id));
        }
        return result;
    }

    const EdgeVoxelRef& getEdgeVoxelRef(size_t id) const {
        return storage_.at(id);
    }
};

struct SkeletonVoxelsWithinAdaptor : public EdgeVoxelRefResultSetBase {
    float maxdistSq_;
    float worstDistSoFar_;
    size_t worstIndexSoFar_;

    SkeletonVoxelsWithinAdaptor(const std::vector<EdgeVoxelRef>& storage, float maxdistSq)
        : EdgeVoxelRefResultSetBase(storage)
        , maxdistSq_(maxdistSq)
        , worstDistSoFar_(-1)
        , worstIndexSoFar_(-1)
    {
    }

    /**
     * Called during search to add an element matching the criteria.
     * @return true if the search should be continued, false if the results are sufficient
     */
    inline bool addPoint(float dist, size_t index)
    {
        if(dist <= maxdistSq_) {
            closestIDs_.push_back(index);
            if(worstDistSoFar_ < dist) {
                worstDistSoFar_ = dist;
                worstIndexSoFar_ = index;
            }
        }
        return true;
    }

    inline DistanceType worstDist() const {
        return maxdistSq_;
    }

    /**
     * Find the worst result (furtherest neighbor) without copying or sorting
     * Pre-conditions: size() > 0
     */
    std::pair<IndexType,DistanceType> worst_item() const
    {
        if (closestIDs_.empty()) throw std::runtime_error("Cannot invoke RadiusResultSet::worst_item() on an empty list of results.");
        return std::pair<IndexType,DistanceType>(worstIndexSoFar_, worstDistSoFar_);
    }
};

struct ClosestSkeletonVoxelsAdaptor : public EdgeVoxelRefResultSetBase {
    float currentDist_;

    ClosestSkeletonVoxelsAdaptor(const std::vector<EdgeVoxelRef>& storage)
        : EdgeVoxelRefResultSetBase(storage)
        , currentDist_(std::numeric_limits<float>::infinity())
    {
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
            }
            currentDist_ = dist;
            closestIDs_.push_back(index);
        }
        return true;
    }

    inline DistanceType worstDist() const {
        return currentDist_*1.001;
    }

    /**
     * Find the worst result (furtherest neighbor) without copying or sorting
     * Pre-conditions: size() > 0
     */
    std::pair<IndexType,DistanceType> worst_item() const
    {
        if (closestIDs_.empty()) throw std::runtime_error("Cannot invoke RadiusResultSet::worst_item() on an empty list of results.");
        size_t id = closestIDs_[0];
        return std::make_pair<IndexType,DistanceType>(std::move(id), worstDist());
    }
};

static std::pair<std::unique_ptr<HDF5FileVolume>, std::unique_ptr<HDF5FileVolume>> splitSegmentationCriticalVoxels(const std::string& tmpPathNoCritical, const std::string& tmpPathOnlyCritical, const ProtoVesselGraph& graph, SkeletonClassReader&& skeletonClassReader, const VolumeBase& segmentation, const float segmentationBinarizationThreshold, ProgressReporter& progress) {
    TaskTimeLogger _("Split segmentation wrt critical voxels", tgt::Info);
    const auto rwToVoxel = segmentation.getWorldToVoxelMatrix();
    const auto voxelToRw = segmentation.getVoxelToWorldMatrix();
    const auto dimensions = segmentation.getDimensions();
    const float criticalVoxelDistDiff = 1.001f*tgt::length(segmentation.getSpacing());

    const std::string format = "uint8";
    const tgt::svec3 slicedim(dimensions.xy(), 1);

    std::unique_ptr<HDF5FileVolume> onlyCriticalVoxels = HDF5FileVolume::createVolume(tmpPathOnlyCritical, HDF5VolumeWriter::VOLUME_DATASET_NAME, format, dimensions, 1, true, 1, slicedim, false);
    std::unique_ptr<HDF5FileVolume> noCriticalVoxels = HDF5FileVolume::createVolume(tmpPathNoCritical, HDF5VolumeWriter::VOLUME_DATASET_NAME, format, dimensions, 1, true, 1, slicedim, false);

    onlyCriticalVoxels->writeOffset(segmentation.getOffset());
    noCriticalVoxels->writeOffset(segmentation.getOffset());

    onlyCriticalVoxels->writeSpacing(segmentation.getSpacing());
    noCriticalVoxels->writeSpacing(segmentation.getSpacing());

    onlyCriticalVoxels->writePhysicalToWorldTransformation(segmentation.getPhysicalToWorldMatrix());
    noCriticalVoxels->writePhysicalToWorldTransformation(segmentation.getPhysicalToWorldMatrix());

    KDTreeBuilder<EdgeVoxelRef> finderBuilder;

    for(auto& edge : graph.edges_) {
        for(auto& rwvoxel : edge.voxels()) {
            finderBuilder.push(EdgeVoxelRef(rwvoxel.rwpos_, edge));
        }
    }

    KDTreeVoxelFinder<EdgeVoxelRef> finder(std::move(finderBuilder));

    for(size_t z = 0; z<dimensions.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dimensions.z);
        std::unique_ptr<const VolumeRAM> slice(segmentation.getSlice(z));
        std::unique_ptr<VolumeRAM> outputSliceOnlyCV(VolumeFactory().create(format, slicedim));
        std::unique_ptr<VolumeRAM> outputSliceNoCV(VolumeFactory().create(format, slicedim));

        skeletonClassReader.advance();
        for(size_t y = 0; y<dimensions.y; ++y) {
            for(size_t x = 0; x<dimensions.x; ++x) {
                const tgt::svec3 p(x,y,0);
                if(slice->getVoxelNormalized(p) >= segmentationBinarizationThreshold) {
                    tgt::ivec3 ipos(x,y,z);
                    tgt::vec3 rwpos = transform(voxelToRw, ipos);

                    bool isCritical;

                    if(finder.storage_.points().empty()) {
                        isCritical = false;
                    } else {
                        ClosestSkeletonVoxelsAdaptor closestVoxels(finder.storage_.points());
                        finder.findClosest(rwpos, closestVoxels);
                        tgtAssert(closestVoxels.size() > 0, "No voxels");

                        float closestDistSq = closestVoxels.currentDist_;

                        float criticalDistance = std::sqrt(closestDistSq) + criticalVoxelDistDiff;
                        float criticalDistanceSq = criticalDistance*criticalDistance;
                        SkeletonVoxelsWithinAdaptor criticalVoxelsRes(finder.storage_.points(), criticalDistanceSq);
                        finder.findClosest(rwpos, criticalVoxelsRes);
                        tgtAssert(closestVoxels.size() <= criticalVoxelsRes.size(), "Found fewer voxels than before");

                        std::set<const ProtoVesselGraphEdge*> criticalEdges;
                        for(const auto& vox : criticalVoxelsRes.getVoxels()) {
                            criticalEdges.insert(&vox->edge);
                        }
                        isCritical = criticalEdges.size() >= 2;
                    }


                    /*
                    // The voxel can only be critical if any of the edges share a node:
                    std::set<size_t> nodeIDs;
                    bool foundCriticalSharingANode = false;
                    for(const auto& edge : criticalEdges) {
                        if(nodeIDs.count(edge->getNodeID1()) > 0) {
                            foundCriticalSharingANode = true;
                            break;
                        }
                        nodeIDs.insert(edge->getNodeID1());

                        if(edge->getNodeID1() != edge->getNodeID2() && nodeIDs.count(edge->getNodeID2()) > 0) {
                            foundCriticalSharingANode = true;
                            break;
                        }
                        nodeIDs.insert(edge->getNodeID2());
                    }

                    bool isCritical = foundCriticalSharingANode;
                    */


                    if(skeletonClassReader.getClass(tgt::ivec2(x, y)) != 2 && isCritical) {
                        outputSliceNoCV->setVoxelNormalized(0, p);
                        outputSliceOnlyCV->setVoxelNormalized(1, p);
                    } else {
                        outputSliceNoCV->setVoxelNormalized(1, p);
                        outputSliceOnlyCV->setVoxelNormalized(0, p);
                    }
                } else {
                    outputSliceNoCV->setVoxelNormalized(0, p);
                    outputSliceOnlyCV->setVoxelNormalized(0, p);
                }
            }
        }
        noCriticalVoxels->writeSlices(outputSliceNoCV.get(), z);
        onlyCriticalVoxels->writeSlices(outputSliceOnlyCV.get(), z);
    }

    progress.setProgress(1.0f);
    return std::make_pair(std::move(noCriticalVoxels), std::move(onlyCriticalVoxels));
}

class NoopMetadata {
public:
    NoopMetadata() {}
    NoopMetadata(tgt::svec2 yzPo, size_t lowerBound, size_t upperBound) {}
    ~NoopMetadata() {}
    NoopMetadata& operator+=(const NoopMetadata& rhs) {
        return *this;
    }
};

static void addFixedForegroundPointsToMask(const std::vector<tgt::vec3>& points, VolumeMask& mask) {
    for(auto& p: points) {
        mask.set(tgt::svec3(tgt::iround(p)), VolumeMask::FIXED_OBJECT);
    }
}

static std::unique_ptr<HDF5FileVolume> createCCAVolume(const VolumeBase& input, const std::string& tmpVolumePath, size_t& numComponents, ProgressReporter& progress) {
    TaskTimeLogger _("Create CCA volume", tgt::Info);

    StreamingComponents<0, NoopMetadata> sc;

    std::function<bool(const VolumeRAM* vol, tgt::svec3 pos)> isOne = [](const VolumeRAM* slice, tgt::svec3 pos) {
        return slice->getVoxelNormalized(pos) > 0.5f;
    };
    auto writeMetaData = [] (uint32_t, const NoopMetadata&) {};

    auto componentConstraintTest = [] (const NoopMetadata& metaData) {
        return true;
    };

    const std::string baseType = "uint32";
    const std::string volumeFilePath = tmpVolumePath;
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
    const tgt::svec3 dim = input.getDimensions();
    const int deflateLevel = 1;
    const tgt::svec3 chunkSize(dim.xy(), 1);

    std::unique_ptr<HDF5FileVolume> outputVolume = HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, input.getDimensions(), 1, true, deflateLevel, chunkSize, false);

    auto stats = sc.cca(input, *outputVolume, writeMetaData, isOne, true, componentConstraintTest, progress);
    numComponents = stats.numComponents;

    return outputVolume;
}

static uint64_t toLinearPos(const tgt::svec3& pos, const tgt::svec3& dimensions) {
    return pos.x + dimensions.x*(pos.y + dimensions.y*pos.z);
}

struct IdVolumeInitializer {
    uint32_t id_;
    IdVolumeStorageInitializer storage_;
    SurfaceBuilder surface_;
    size_t numUnlabeledForegroundVoxels_;
    size_t numTotalVoxels_;

    tgt::svec3 offset_;
    tgt::svec3 size_;

    IdVolumeInitializer(uint32_t id, tgt::svec3 offset, tgt::svec3 size)
        : id_(id)
        , storage_(VoreenApplication::app()->getUniqueTmpFilePath(".raw"))
        , surface_()
        , numUnlabeledForegroundVoxels_(0)
        , numTotalVoxels_(0)
        , offset_(offset)
        , size_(size)
    {
    }

    IdVolumeInitializer(IdVolumeInitializer&& old)
        : id_(old.id_)
        , storage_(std::move(old.storage_))
        , surface_(std::move(old.surface_))
        , numUnlabeledForegroundVoxels_(old.numUnlabeledForegroundVoxels_)
        , numTotalVoxels_(old.numTotalVoxels_)
        , offset_(old.offset_)
        , size_(old.size_)
    {
    }

    tgt::SBounds getBounds() const {
        tgtAssert(tgt::hmul(size_) > 0, "Invalid size");
        return tgt::SBounds(offset_, offset_+size_-tgt::svec3::one);
    }

    void pushLabel(IdVolume::Value v) {
        storage_.push(v);
        ++numTotalVoxels_;
    }

    void pushSurfaceVoxel(const tgt::svec3& p) {
        surface_.push(toLinearPos(p - offset_, size_));
    }
};

static void initializeIdVolumes(HDF5FileVolume& branchIds, const HDF5FileVolume& holeIds, std::vector<IdVolumeInitializer>& initializers, ProgressReporter& progress) {
    TaskTimeLogger _("Initialize IdVolumes", tgt::Info);
    IdVolumeInitializationReader initializationReader(branchIds, holeIds);

    tgt::svec3 dim = branchIds.getDimensions();
    tgtAssert(holeIds.getDimensions() == dim, "Volume dimension mismatch");

    std::vector<std::pair<IdVolumeInitializer*, tgt::SBounds>> initializer_finder_init;
    for(auto& initializer: initializers) {
        initializer_finder_init.emplace_back(&initializer, initializer.getBounds());
    }
    BoundsHierarchy<size_t, IdVolumeInitializer*> initializer_finder(std::move(initializer_finder_init));

    for(size_t z=0; z < dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);
        initializationReader.advance();

        for(size_t y=0; y < dim.y; ++y) {
            for(size_t x=0; x < dim.x; ++x) {
                tgt::svec3 pos(x,y,z);

                for(IdVolumeInitializer* initializer : initializer_finder.findBounds(pos)) {
                    uint32_t label;

                    if(initializationReader.isObject(pos, initializer->id_)) {
                        if(initializationReader.isLabeled(pos)) {
                            label = initializationReader.getBranchId(pos);

                            bool isLabelSurfaceVoxel =
                                !( (x == 0       || initializationReader.isLabeled(tgt::ivec3(x-1, y  , z  )))
                                && (x == dim.x-1 || initializationReader.isLabeled(tgt::ivec3(x+1, y  , z  )))
                                && (y == 0       || initializationReader.isLabeled(tgt::ivec3(x  , y-1, z  )))
                                && (y == dim.y-1 || initializationReader.isLabeled(tgt::ivec3(x  , y+1, z  )))
                                && (z == 0       || initializationReader.isLabeled(tgt::ivec3(x  , y  , z-1)))
                                && (z == dim.z-1 || initializationReader.isLabeled(tgt::ivec3(x  , y  , z+1))));

                            if(isLabelSurfaceVoxel) {
                                initializer->pushSurfaceVoxel(pos);
                            }
                        } else {
                            label = IdVolume::UNLABELED_FOREGROUND_VALUE;
                            initializer->numUnlabeledForegroundVoxels_++;
                        }
                    } else {
                        label = IdVolume::BACKGROUND_VALUE;
                    }
                    initializer->pushLabel(label);
                }
            }
        }
    }
    progress.setProgress(1.0f);
}

struct IdVolumeFinalizer {
    uint32_t id_;
    IdVolume vol_;

    tgt::svec3 offset_;
    tgt::svec3 size_;

    IdVolumeFinalizer(IdVolumeInitializer&& volume, ProgressReporter& progress)
        : id_(volume.id_)
        , vol_(std::move(volume.storage_), SurfaceBuilder::finalize(std::move(volume.surface_)), volume.size_, volume.numUnlabeledForegroundVoxels_)
        , offset_(volume.offset_)
        , size_(volume.size_)
    {
        tgtAssert(volume.numTotalVoxels_ == tgt::hmul(size_), "Incomplete IdVolumeInitializer");
        vol_.floodFromLabels(progress, std::numeric_limits<size_t>::max());
    }

    bool containsPoint(const tgt::svec3& p) const {
        return
            offset_.x <= p.x && p.x < size_.x+offset_.x
         && offset_.y <= p.y && p.y < size_.y+offset_.y
         && offset_.z <= p.z && p.z < size_.z+offset_.z;
    }
    tgt::SBounds getBounds() const {
        tgtAssert(tgt::hmul(size_) > 0, "Invalid size");
        return tgt::SBounds(offset_, offset_+size_-tgt::svec3::one);
    }

    IdVolume::Value getValue(const tgt::svec3& p) const {
        tgtAssert(containsPoint(p), "Invalid p");
        return vol_.data_->get(p - offset_);
    }
};

static void finalizeIdVolumes(HDF5FileVolume& branchIds, const HDF5FileVolume& holeIds, std::vector<IdVolumeFinalizer>& finalizers, ProgressReporter& progress) {
    TaskTimeLogger _("Initialize IdVolumes", tgt::Info);
    tgtAssert(branchIds.getDimensions() == holeIds.getDimensions(), "Invalid dimensions");

    const tgt::svec3 dim = branchIds.getDimensions();

    std::vector<std::pair<IdVolumeFinalizer*, tgt::SBounds>> finalizer_finder_init;
    for(auto& finalizer: finalizers) {
        finalizer_finder_init.emplace_back(&finalizer, finalizer.getBounds());
    }
    BoundsHierarchy<size_t, IdVolumeFinalizer*> finalizer_finder(std::move(finalizer_finder_init));

    for(size_t z=0; z < dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);

        std::unique_ptr<VolumeRAM_UInt32> branchSlice(dynamic_cast<VolumeRAM_UInt32*>(branchIds.loadSlices(z, z)));
        std::unique_ptr<VolumeRAM_UInt32> holeSlice(dynamic_cast<VolumeRAM_UInt32*>(holeIds.loadSlices(z, z)));
        tgtAssert(branchSlice, "Invalid volume format");
        tgtAssert(holeSlice, "Invalid volume format");

        for(size_t y=0; y < dim.y; ++y) {
            for(size_t x=0; x < dim.x; ++x) {
                tgt::svec3 pos(x,y,z);
                for(IdVolumeFinalizer* finalizer : finalizer_finder.findBounds(pos)) {
                    if(holeSlice->voxel(x,y,0) == finalizer->id_) {
                        uint32_t floodedId = finalizer->getValue(pos);
                        tgtAssert(floodedId != IdVolume::BACKGROUND_VALUE, "flooded label is background");
                        if(floodedId != IdVolume::UNLABELED_FOREGROUND_VALUE) {
                            branchSlice->voxel(x,y,0) = floodedId;
                        }
                    }
                }
            }
        }
        branchIds.writeSlices(branchSlice.get(), z);
    }
    progress.setProgress(1.0f);
}

struct UnfinishedRegions {
    std::unique_ptr<HDF5FileVolume> holeIds;
    std::map<uint32_t, tgt::SBounds> regions;

    /*
    void floodRegion(uint32_t regionId, HDF5FileVolume& branchIds, ProgressReporter& progress) {
        auto bounds = regions.at(regionId);

        tgt::svec3 llf = bounds.getLLF();
        tgt::svec3 urb = bounds.getURB();

        // Expand region in order to provide initial labels for flooding
        llf.x = llf.x>0 ? llf.x-1 : 0;
        llf.y = llf.y>0 ? llf.y-1 : 0;
        llf.z = llf.z>0 ? llf.z-1 : 0;

        // urb for TemplateBounds<size_t> is just the maximum position, but we need the first index that is outside the box
        //                               |
        //                              \|/
        urb.x = std::min<size_t>(urb.x+1+1, branchIds.getDimensions().x);
        urb.y = std::min<size_t>(urb.y+1+1, branchIds.getDimensions().y);
        urb.z = std::min<size_t>(urb.z+1+1, branchIds.getDimensions().z);

        tgt::svec3 offset = llf;
        tgt::svec3 dim = urb - llf;

        IdVolumeInitializationReader reader(branchIds, *holeIds, offset, dim, regionId);
        IdVolume idVol(reader);

        idVol.floodFromLabels(progress, 1000);

        IdVolumeReader idReader(idVol);

        for(size_t z=0; z < dim.z; ++z) {
            idReader.advance();
            tgt::svec3 sliceOffset(offset.xy(), offset.z+z);
            tgt::svec3 sliceDim(dim.xy(), 1);

            std::unique_ptr<VolumeRAM_UInt32> branchSlice(dynamic_cast<VolumeRAM_UInt32*>(branchIds.loadBrick(sliceOffset, sliceDim)));
            std::unique_ptr<VolumeRAM_UInt32> holeSlice(dynamic_cast<VolumeRAM_UInt32*>(holeIds->loadBrick(sliceOffset, sliceDim)));
            tgtAssert(branchSlice, "Invalid volume format");
            tgtAssert(holeSlice, "Invalid volume format");

            for(size_t y=0; y < dim.y; ++y) {
                for(size_t x=0; x < dim.x; ++x) {
                    if(holeSlice->voxel(x,y,0) == regionId) {
                        uint32_t floodedId = idReader.getId(tgt::svec3(x,y,z));
                        tgtAssert(floodedId != IdVolume::BACKGROUND_VALUE, "flooded label is background");
                        if(floodedId != IdVolume::UNLABELED_FOREGROUND_VALUE) {
                            branchSlice->voxel(x,y,0) = floodedId;
                        }
                    }
                }
            }
            branchIds.writeBrick(branchSlice.get(), sliceOffset);
        }
    }

    void floodAllRegions(HDF5FileVolume& branchIds, ProgressReporter& progress) {
        TaskTimeLogger _("Flood remaining unlabled regions", tgt::Info);
        size_t i=0;
        for(auto& kv : regions) {
            SubtaskProgressReporter sp(progress, tgt::vec2(i, i+1)/static_cast<float>(regions.size()));
            floodRegion(kv.first, branchIds, sp);
        }
    }
    */

    void floodAllRegions(HDF5FileVolume& branchIds, ProgressReporter& progress) {
        TaskTimeLogger _("Flood remaining unlabled regions", tgt::Info);

        SubtaskProgressReporterCollection<3> subtaskReporters(progress);

        std::vector<IdVolumeInitializer> initializers;
        for(auto& kv: regions) {

            uint32_t id = kv.first;
            auto bounds = kv.second;

            tgt::svec3 llf = bounds.getLLF();
            tgt::svec3 urb = bounds.getURB();

            // Expand region in order to provide initial labels for flooding
            llf.x = llf.x>0 ? llf.x-1 : 0;
            llf.y = llf.y>0 ? llf.y-1 : 0;
            llf.z = llf.z>0 ? llf.z-1 : 0;

            // urb for TemplateBounds<size_t> is just the maximum position, but we need the first index that is outside the box
            //                               |
            //                              \|/
            urb.x = std::min<size_t>(urb.x+1+1, branchIds.getDimensions().x);
            urb.y = std::min<size_t>(urb.y+1+1, branchIds.getDimensions().y);
            urb.z = std::min<size_t>(urb.z+1+1, branchIds.getDimensions().z);

            tgt::svec3 offset = llf;
            tgt::svec3 dim = urb - llf;

            initializers.emplace_back(id, offset, dim);
        }

        initializeIdVolumes(branchIds, *holeIds, initializers, subtaskReporters.get<0>());

        std::vector<IdVolumeFinalizer> finalizers;
        size_t i=0;
        for(auto& initializer: initializers) {
            SubtaskProgressReporter sp(subtaskReporters.get<1>(), tgt::vec2(i, i+1)/static_cast<float>(regions.size()));
            finalizers.emplace_back(std::move(initializer), sp);
            ++i;
        }

        finalizeIdVolumes(branchIds, *holeIds, finalizers, subtaskReporters.get<2>());
    }
};

static UnfinishedRegions collectUnfinishedRegions(const VolumeBase& input, const std::string& tmpVolumePath, ProgressReporter& progress) {
    TaskTimeLogger _("Collect unlabled regions", tgt::Info);

    std::map<uint32_t, tgt::SBounds> regions;

    StreamingComponents<0, CCANodeMetaData> sc;

    std::function<bool(const VolumeRAM* vol, tgt::svec3 pos)> isOne = [](const VolumeRAM* slice, tgt::svec3 pos) {
        return slice->getVoxelNormalized(pos) > 0.5f;
    };
    auto writeMetaData = [&regions] (uint32_t id, const CCANodeMetaData& metadata) {
        regions.emplace(id, metadata.bounds_);
    };

    auto componentConstraintTest = [] (const CCANodeMetaData&) {
        return true;
    };

    const std::string baseType = "uint32";
    const std::string volumeFilePath = tmpVolumePath;
    const std::string volumeLocation = HDF5VolumeWriter::VOLUME_DATASET_NAME;
    const tgt::svec3 dim = input.getDimensions();
    const int deflateLevel = 1;
    const tgt::svec3 chunkSize(dim.xy(), 1);

    std::unique_ptr<HDF5FileVolume> outputVolume = HDF5FileVolume::createVolume(volumeFilePath, volumeLocation, baseType, input.getDimensions(), 1, true, deflateLevel, chunkSize, false);

    sc.cca(input, *outputVolume, writeMetaData, isOne, true, componentConstraintTest, progress);

    return UnfinishedRegions {
        std::move(outputVolume), regions
    };
}

static void addClippedOffRegionsToCritical(HDF5FileVolume& criticalVoxels, const HDF5FileVolume& ccaVolume, const ProtoVesselGraph& graph, size_t numComponents, const VolumeBase& segmentation, float segmentationBinarizationThreshold, ProgressReporter& progress) {
    TaskTimeLogger _("Add clipped off regions to critical voxels", tgt::Info);

    BranchIdVolumeReader ccaReader(ccaVolume, graph, numComponents, segmentation, segmentationBinarizationThreshold);
    auto dim = criticalVoxels.getDimensions();

    for(size_t z = 0; z < dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);

        std::unique_ptr<VolumeRAM> slice(criticalVoxels.loadSlices(z,z));
        tgtAssert(slice, "Invalid volume format");

        ccaReader.advance();

        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                tgt::ivec3 p(x,y,z);

                if(ccaReader.isObject(p) && !ccaReader.isValidEdgeId(ccaReader.getEdgeId(p.xy()))) {
                    slice->setVoxelNormalized(1.0f, x, y, 0);
                }
            }
        }

        criticalVoxels.writeSlices(slice.get(), z);
    }
    progress.setProgress(1.0f);
}

std::unique_ptr<VesselGraph> createGraphFromMask(VesselGraphCreatorInput& input, VolumeMask&& skeleton, VolumeList& generatedSkeletons, ProgressReporter& progress) {
    SubtaskProgressReporterCollection<8> subtaskReporters(progress);

    // Create new protograph
    std::unique_ptr<MetaDataCollector> mdc = cca(NeighborCountVoxelClassifier(skeleton), subtaskReporters.get<0>());
    std::unique_ptr<ProtoVesselGraph> protograph = mdc->createProtoVesselGraph(input.segmentation.getDimensions(), input.segmentation.getVoxelToWorldMatrix(), input.sampleMask, subtaskReporters.get<1>());


    // Create an assignment edge <-> vessel component
    //  1. Split vessel components at nodes edges
    const std::string segNoCriticalTmpPath = VoreenApplication::app()->getUniqueTmpFilePath(".h5");
    const std::string segOnlyCriticalTmpPath = VoreenApplication::app()->getUniqueTmpFilePath(".h5");

    auto segmentationWithAndWithoutCriticalVoxels = splitSegmentationCriticalVoxels(segNoCriticalTmpPath, segOnlyCriticalTmpPath, *protograph, SkeletonClassReader(skeleton), input.segmentation, input.binarizationThresholdSegmentationNormalized, subtaskReporters.get<2>());
    std::unique_ptr<HDF5FileVolume> segNoCriticalVoxelsHDF5 = std::move(segmentationWithAndWithoutCriticalVoxels.first);
    std::unique_ptr<HDF5FileVolume> segOnlyCriticalVoxelsHDF5 = std::move(segmentationWithAndWithoutCriticalVoxels.second);

    segNoCriticalVoxelsHDF5.reset(nullptr); //drop hdf5 vol

    //  2. Perform a cca to distinguish components
    size_t numComponents;
    std::unique_ptr<VolumeBase> segNoCriticalVoxels(HDF5VolumeReader().read(segNoCriticalTmpPath)->at(0));
    const std::string ccaNoCriticalTmpPath = VoreenApplication::app()->getUniqueTmpFilePath(".h5");
    std::unique_ptr<HDF5FileVolume> branchIdSegmentation = createCCAVolume(*segNoCriticalVoxels, ccaNoCriticalTmpPath, numComponents, subtaskReporters.get<3>());

    // 4. Find unlabeled regions
    //
    // Add cut off unlabeled regions to the critical voxels (i.e., unlabeled regions) to allow flooding from valid edge label regions later
    addClippedOffRegionsToCritical(*segOnlyCriticalVoxelsHDF5, *branchIdSegmentation, *protograph, numComponents, input.segmentation, input.binarizationThresholdSegmentationNormalized, subtaskReporters.get<4>());
    segOnlyCriticalVoxelsHDF5.reset(nullptr); //drop hdf5 vol
    std::unique_ptr<VolumeBase> segOnlyCriticalVoxels(HDF5VolumeReader().read(segOnlyCriticalTmpPath)->at(0)); //reopen hdf5vol as VolumeBase
    const std::string ccaOnlyCriticalTmpPath = VoreenApplication::app()->getUniqueTmpFilePath(".h5");
    auto unfinishedRegions = collectUnfinishedRegions(*segOnlyCriticalVoxels, ccaOnlyCriticalTmpPath, subtaskReporters.get<5>());

    //  5. Flood remaining unlabeled regions locally
    unfinishedRegions.floodAllRegions(*branchIdSegmentation, subtaskReporters.get<6>());

    //  6. Create a reader that reads edge ids from the underlying segmentation
    BranchIdVolumeReader ccaReader(*branchIdSegmentation, *protograph, numComponents, input.segmentation, input.binarizationThresholdSegmentationNormalized);

    //  7. Create better VesselGraph from protograph and and the edge-id-segmentation
    auto output = protograph->createVesselGraph(ccaReader, input.sampleMask, subtaskReporters.get<7>());


#ifdef VRN_DEBUG
    generatedSkeletons.add(HDF5VolumeReader().read(ccaOnlyCriticalTmpPath)->at(0));
    generatedSkeletons.add(HDF5VolumeReader().read(ccaNoCriticalTmpPath)->at(0));
#endif

    // Clean up files
    tgt::FileSystem::deleteFile(segNoCriticalTmpPath);
    tgt::FileSystem::deleteFile(segOnlyCriticalTmpPath);
    tgt::FileSystem::deleteFile(ccaNoCriticalTmpPath);
    tgt::FileSystem::deleteFile(ccaOnlyCriticalTmpPath);

    return output;
}

std::unique_ptr<VesselGraph> refineVesselGraph(VesselGraphCreatorInput& input, const VesselGraph& prevGraph, VolumeList& generatedSkeletons, ProgressReporter& progress) {
    TaskTimeLogger _("Refinement iteration", tgt::Info);

    SubtaskProgressReporterCollection<3> subtaskReporters(progress);
    progress.setProgress(0.0f);

    // Refine previous graph
    auto isEdgeDeletable = [&input] (const VesselGraphEdge& edge) {
        return !edge.hasValidData() || edge.getVoxels().size() < input.minVoxelLength || edge.getElongation() < input.minElongation || edge.getRelativeBulgeSize() < input.minBulgeSize;
    };
    std::unique_ptr<VesselGraph> normalizedGraph = VesselGraphNormalization::removeEndEdgesRecursively(prevGraph, isEdgeDeletable);

    tgt::mat4 rwToVoxel = input.segmentation.getWorldToVoxelMatrix();

    // Create new voxelmask
    auto fixedForegroundMaskOnlyEndvoxels = GraphNodeVoxelReader::create<EndNodeVoxelExtractor>(*normalizedGraph, rwToVoxel, input.segmentation.getDimensions());
    VolumeMask mask(input.segmentation, input.sampleMask, fixedForegroundMaskOnlyEndvoxels, input.binarizationThresholdSegmentationNormalized, subtaskReporters.get<0>());
    addFixedForegroundPointsToMask(input.fixedForegroundPoints, mask);

    // Skeletonize new mask
    mask.skeletonize<VolumeMask::IMPROVED_NO_LINE_PRESERVATION>(std::numeric_limits<size_t>::max(), subtaskReporters.get<1>());

    return createGraphFromMask(input, std::move(mask), generatedSkeletons, subtaskReporters.get<2>());
}

std::unique_ptr<VesselGraph> createInitialVesselGraph(VesselGraphCreatorInput& input, VolumeList& generatedSkeletons, ProgressReporter& progress) {
    TaskTimeLogger _("Create initial VesselGraph", tgt::Info);
    SubtaskProgressReporterCollection<3> subtaskReporters(progress);
    progress.setProgress(0.0f);

    VolumeMask mask(input.segmentation, input.sampleMask, NoFixedForeground(), input.binarizationThresholdSegmentationNormalized, subtaskReporters.get<0>());
    addFixedForegroundPointsToMask(input.fixedForegroundPoints, mask);
    mask.skeletonize<VolumeMask::IMPROVED>(std::numeric_limits<size_t>::max(), subtaskReporters.get<1>());

    return createGraphFromMask(input, std::move(mask), generatedSkeletons, subtaskReporters.get<2>());
}

static bool tryExtractPoints(const Geometry* geometry, std::vector<tgt::vec3>& output, const tgt::mat4& rwToVoxel) {
    if(!geometry) {
        return false;
    }
    auto transformation = rwToVoxel*geometry->getTransformationMatrix();
    auto transform = [&transformation] (const tgt::vec3& v) {
        return (transformation*tgt::vec4(v, 1)).xyz();
    };
    if(const PointSegmentListGeometryVec3* seedList = dynamic_cast<const PointSegmentListGeometryVec3*>(geometry)) {
        for(int i = 0; i < seedList->getNumSegments(); i++) {
            auto segment = seedList->getSegment(i);
            for(const tgt::vec3& v: segment) {
                output.push_back(transform(v));
            }
        }
    } else if(const PointListGeometryVec3* seeds = dynamic_cast<const PointListGeometryVec3*>(geometry)) {
        for(const tgt::vec3& v: *seeds) {
            output.push_back(transform(v));
        }
    } else {
        return false;
    }
    return true;
}

VesselGraphCreatorInput VesselGraphCreator::prepareComputeInput() {
    if(segmentedVolumeInport_.hasChanged()) {
        clearOutports();
    }
    auto segmentation = segmentedVolumeInport_.getThreadSafeData();
    auto sampleMask = sampleMaskInport_.getThreadSafeData(); //May be null and still valid!
    if(!segmentation) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }
    if(sampleMask && sampleMask->getDimensions() != segmentation->getDimensions()) {
        throw InvalidInputException("Dimensions of sampleMask and segmented volume differ", InvalidInputException::S_ERROR);
    }


    std::vector<tgt::vec3> fixedForegroundPoints;
    tryExtractPoints(fixedForegroundPointInport_.getThreadSafeData(), fixedForegroundPoints, segmentation->getWorldToVoxelMatrix());

    const float binarizationThresholdSegmentationNormalized = segmentation->getRealWorldMapping().realWorldToNormalized(binarizationThresholdSegmentation_.get());

    generatedSkeletonsOutport_.clear();
    return VesselGraphCreatorInput(
        *segmentation,
        sampleMask,
        std::move(fixedForegroundPoints),
        binarizationThresholdSegmentationNormalized,
        numRefinementIterations_.get(),
        minVoxelLength_.get(),
        minElongation_.get(),
        minBulgeSize_.get()
    );
}

VesselGraphCreatorOutput VesselGraphCreator::compute(VesselGraphCreatorInput input, ProgressReporter& progressReporter) const {
    TaskTimeLogger _("Extract VesselGraph (total)", tgt::Info);
    float progressPerIteration = 1.0f/(input.numRefinementIterations+1);
    SubtaskProgressReporter initialProgress(progressReporter, tgt::vec2(0, progressPerIteration));
    std::unique_ptr<VolumeList> generatedSkeletons(new VolumeContainer());
    std::unique_ptr<VesselGraph> graph = createInitialVesselGraph(input, *generatedSkeletons, initialProgress);
    for(int i=0; i < input.numRefinementIterations; ++i) {
        SubtaskProgressReporter refinementProgress(progressReporter, progressPerIteration*tgt::vec2(i+1, i+2));
        try {
            graph = refineVesselGraph(input, *graph, *generatedSkeletons, refinementProgress);
        } catch(tgt::IOException e) {
            LERROR("Could not create hdf5 output volume for branchIdSegmentation.");
            std::cout << e.what() << std::endl;
            break;
        }
    }
    return VesselGraphCreatorOutput {
        std::move(graph),
        std::move(generatedSkeletons)
    };
}

void VesselGraphCreator::processComputeOutput(VesselGraphCreatorOutput output) {
    const VesselGraph* graph = output.graph.release();
    graphOutport_.setData(graph);
    if(graph) {
        PointListGeometryVec3* nodeGeom = new PointListGeometryVec3();
        for(auto& node : graph->getNodes()) {
            nodeGeom->addPoint(node.pos_);
        }
        nodeOutport_.setData(nodeGeom);

        PointSegmentListGeometryVec3* edgeGeom = new PointSegmentListGeometryVec3();
        for(auto& edge : graph->getEdges()) {
            std::vector<tgt::vec3> segment(edge.getVoxels().size());
            std::transform(edge.getVoxels().begin(), edge.getVoxels().end(), segment.begin(), [] (const VesselSkeletonVoxel& voxel) {
                    return voxel.pos_;
                    });
            edgeGeom->addSegment(segment);
        }
        edgeOutport_.setData(edgeGeom);
    } else {
        nodeOutport_.setData(nullptr);
        edgeOutport_.setData(nullptr);
    }
    generatedSkeletonsOutport_.setData(output.generatedSkeletons.release());
}
void VesselGraphCreator::adjustPropertiesToInput() {
    const VolumeBase* segmentation = segmentedVolumeInport_.getData();
    if(segmentation) {
        if(!segmentation->hasDerivedData<VolumeMinMax>()) {
            LWARNING("Calculating VolumeMinMax. This may take a while...");
        }
        const VolumeMinMax* mm = segmentation->getDerivedData<VolumeMinMax>();
        bool minMaxWereEqual = binarizationThresholdSegmentation_.getMinValue() == binarizationThresholdSegmentation_.getMaxValue();
        binarizationThresholdSegmentation_.setMinValue(mm->getMin());
        binarizationThresholdSegmentation_.setMaxValue(mm->getMax());
        binarizationThresholdSegmentation_.adaptDecimalsToRange(3);

        if(minMaxWereEqual) {
            binarizationThresholdSegmentation_.set((binarizationThresholdSegmentation_.getMinValue() + binarizationThresholdSegmentation_.getMaxValue())*0.5f);
        }
    }
}

} // namespace voreen
