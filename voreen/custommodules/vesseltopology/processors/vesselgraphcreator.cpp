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
#include "voreen/core/utils/stringutils.h"

#include "../algorithm/streaminggraphcreation.h"
#include "../algorithm/volumemask.h"
#include "../algorithm/idvolume.h"
#include "../algorithm/vesselgraphnormalization.h"
#include "../algorithm/boundshierarchy.h"

#include "custommodules/bigdataimageprocessing/datastructures/lz4slicevolume.h"
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
    , generatedVolumesOutport_(Port::OUTPORT, "generatedSkeletons.outport", "Generated Volumes", false, Processor::VALID)
    , generatedGraphsOutport_(Port::OUTPORT, "generatedGraphs.outport", "Generated Graphs", false, Processor::VALID)
    , numRefinementIterations_("numRefinementIterations", "Refinement Iterations", 0, 0, 100, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_ADVANCED)
    , minVoxelLength_("minVoxelLength", "Min Voxel Length", 0, 0, 50, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
    , minElongation_("minElongation", "Minimum Elongation", 0, 0, 5, Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_DEBUG)
    , minBulgeSize_("minBulgeSize", "Minimum Bulge Size", 0, 0, 10, Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_DEFAULT)
    , saveDebugData_("saveDebugData", "Generate Debug Data", false, Processor::INVALID_RESULT, Property::LOD_DEBUG)
    , binarizationThresholdSegmentation_("binarizationThresholdSegmentation", "Binarization Threshold (Segmentation)", 0.5f, 0.0f, 1.0f, Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
    , tmpStorageSizeInfo_("tmpStorageSizeInfo", "Required temporary storage (estimated)", "")
{
    addPort(segmentedVolumeInport_);
    addPort(sampleMaskInport_);
    addPort(fixedForegroundPointInport_);
    addPort(graphOutport_);
    addPort(nodeOutport_);
    addPort(edgeOutport_);
    addPort(generatedVolumesOutport_);
    addPort(generatedGraphsOutport_);

        addProperty(numRefinementIterations_);
            numRefinementIterations_.setGroupID("refinement");
        addProperty(minVoxelLength_);
            minVoxelLength_.setGroupID("refinement");
        addProperty(minElongation_);
            minElongation_.setGroupID("refinement");
        addProperty(minBulgeSize_);
            minBulgeSize_.setGroupID("refinement");
        addProperty(saveDebugData_);
            saveDebugData_.setGroupID("refinement");
    setPropertyGroupGuiName("refinement", "Refinement");

        addProperty(binarizationThresholdSegmentation_);
            binarizationThresholdSegmentation_.setGroupID("binarization");
    setPropertyGroupGuiName("binarization", "Binarization");

        addProperty(tmpStorageSizeInfo_);
            tmpStorageSizeInfo_.setReadOnlyFlag(true);
            tmpStorageSizeInfo_.setGroupID("info");
    setPropertyGroupGuiName("info", "Information");
}

VesselGraphCreator::~VesselGraphCreator() {
}
VoreenSerializableObject* VesselGraphCreator::create() const {
    return new VesselGraphCreator();
}

bool VesselGraphCreator::isReady() const {
    return isInitialized() && segmentedVolumeInport_.isReady() &&
        (graphOutport_.isReady() || nodeOutport_.isReady() || edgeOutport_.isReady() || generatedVolumesOutport_.isReady() || generatedGraphsOutport_.isReady());
}

struct VesselGraphCreatorProcessedInput {
    VesselGraphCreatorProcessedInput(VesselGraphCreatorInput& input, ProgressReporter& progress)
        : segmentation(binarizeVolume(input.segmentation, input.binarizationThresholdSegmentationNormalized,
                    input.sampleMask ? SubtaskProgressReporter(progress, tgt::vec2(0,0.5)) : SubtaskProgressReporter(progress, tgt::vec2(0,1))))
        , sampleMask(input.sampleMask
                ? boost::optional<LZ4SliceVolume<uint8_t>>(binarizeVolume(*input.sampleMask, input.binarizationThresholdSegmentationNormalized, SubtaskProgressReporter(progress, tgt::vec2(0,0.5))))
                : boost::none)
        , fixedForegroundPoints(std::move(input.fixedForegroundPoints))
        , numRefinementIterations(input.numRefinementIterations)
        , minVoxelLength(input.minVoxelLength)
        , minElongation(input.minElongation)
        , minBulgeSize(input.minBulgeSize)
        , saveDebugData(input.saveDebugData)
    {
    }

    LZ4SliceVolume<uint8_t> segmentation;
    boost::optional<LZ4SliceVolume<uint8_t>> sampleMask;

    std::vector<tgt::vec3> fixedForegroundPoints;
    int numRefinementIterations;
    int minVoxelLength;
    float minElongation;
    float minBulgeSize;
    bool saveDebugData;
};

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

static std::pair<LZ4SliceVolume<uint8_t>, LZ4SliceVolume<uint8_t>> splitSegmentationCriticalVoxels(const std::string& tmpPathNoCritical, const std::string& tmpPathOnlyCritical, const ProtoVesselGraph& graph, SkeletonClassReader&& skeletonClassReader, const VesselGraphCreatorProcessedInput& input, ProgressReporter& progress) {
    TaskTimeLogger _("Split segmentation wrt critical voxels", tgt::Info);

    const auto voxelToRw = input.segmentation.getMetaData().getVoxelToWorldMatrix();
    const auto dimensions = input.segmentation.getDimensions();
    const float criticalVoxelDistDiff = 1.001f*tgt::length(input.segmentation.getMetaData().getSpacing());

    const std::string format = "uint8";
    const tgt::svec3 slicedim(dimensions.xy(), 1);

    KDTreeBuilder<EdgeVoxelRef> finderBuilder;

    for(auto& edge : graph.edges_) {
        for(auto& rwvoxel : edge.voxels()) {
            finderBuilder.push(EdgeVoxelRef(rwvoxel.rwpos_, edge));
        }
    }

    KDTreeVoxelFinder<EdgeVoxelRef> finder(std::move(finderBuilder));

    LZ4SliceVolumeBuilder<uint8_t> onlyCriticalVoxels(tmpPathOnlyCritical, input.segmentation.getMetaData());
    LZ4SliceVolumeBuilder<uint8_t> noCriticalVoxels(tmpPathNoCritical, input.segmentation.getMetaData());

    for(size_t z = 0; z<dimensions.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dimensions.z);
        auto slice = input.segmentation.loadSlice(z);

        auto outputSliceOnlyCV = onlyCriticalVoxels.getNextWritableSlice();
        auto outputSliceNoCV = noCriticalVoxels.getNextWritableSlice();

        skeletonClassReader.advance();
        for(size_t y = 0; y<dimensions.y; ++y) {
            for(size_t x = 0; x<dimensions.x; ++x) {
                const tgt::svec3 p(x,y,0);
                if(slice.voxel(p) > 0) {
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
                        outputSliceNoCV->voxel(p) = 0;
                        outputSliceOnlyCV->voxel(p) = 1;
                    } else {
                        outputSliceNoCV->voxel(p) = 1;
                        outputSliceOnlyCV->voxel(p) = 0;
                    }
                } else {
                    outputSliceNoCV->voxel(p) = 0;
                    outputSliceOnlyCV->voxel(p) = 0;
                }
            }
        }
    }

    progress.setProgress(1.0f);
    return std::make_pair(
            std::move(noCriticalVoxels).finalize(),
            std::move(onlyCriticalVoxels).finalize());
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

struct LZ4SliceVolumeCCAInputWrapper {
    LZ4SliceVolumeCCAInputWrapper(const LZ4SliceVolume<uint8_t>& vol)
        : volume_(vol)
    {
    }

    tgt::svec3 getDimensions() const {
        return volume_.getDimensions();
    }
    tgt::vec3 getSpacing() const {
        return tgt::vec3::one;
    }
    tgt::vec3 getOffset() const {
        return tgt::vec3::zero;
    }
    tgt::mat4 getPhysicalToWorldMatrix() const {
        return tgt::mat4::identity;
    }
    std::string getBaseType() const {
        return "uint32";
    }
    VolumeRAM* getSlice(size_t z) const {
        return new VolumeAtomic<uint8_t>(volume_.loadSlice(z));
    }
    const LZ4SliceVolume<uint8_t>& volume_;
};
struct LZ4SliceVolumeCCABuilderWrapper {
    LZ4SliceVolumeCCABuilderWrapper(LZ4SliceVolumeBuilder<uint32_t>& builder)
        : builder_(builder)
    {
    }

    tgt::svec3 getDimensions() const {
        return builder_.getDimensions();
    }
    std::string getBaseType() const {
        return "uint32";
    }
    void writeSlices(VolumeAtomic<uint32_t>* slice, size_t _z) const {
        tgtAssert(slice, "No slice");
        builder_.pushSlice(*slice);
    }
    void writeSlices(VolumeRAM* slice, size_t _z) const {
        tgtAssert(false, "Invalid slice");
    }
    void writeSpacing(const tgt::svec3&) {
        //ignored for now
    }
    void writeOffset(const tgt::svec3&) {
        //ignored for now
    }
    void writePhysicalToWorldTransformation(const tgt::mat4&) {
        //ignored for now
    }
    void writeRealWorldMapping(const RealWorldMapping&) {
        //ignored for now
    }
    void writeVolumeMinMax(const VolumeMinMax*) {
        //ignored for now
    }
    LZ4SliceVolumeBuilder<uint32_t>& builder_;
};

static LZ4SliceVolume<uint32_t> createCCAVolume(const LZ4SliceVolume<uint8_t>& input, const std::string& tmpVolumePath, size_t& numComponents, ProgressReporter& progress) {
    TaskTimeLogger _("Create CCA volume", tgt::Info);

    StreamingComponents<0, NoopMetadata> sc;

    std::function<bool(const VolumeRAM* vol, tgt::svec3 pos)> isOne = [](const VolumeRAM* slice, tgt::svec3 pos) {
        return slice->getVoxelNormalized(pos) > 0.0f;
    };
    auto writeMetaData = [] (uint32_t, const NoopMetadata&) {};

    auto componentConstraintTest = [] (const NoopMetadata& metaData) {
        return true;
    };

    LZ4SliceVolumeBuilder<uint32_t> outputVolumeBuilder(tmpVolumePath, input.getMetaData().withRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint32_t>()));
    LZ4SliceVolumeCCABuilderWrapper outputWrapper(outputVolumeBuilder);

    auto stats = sc.cca(LZ4SliceVolumeCCAInputWrapper(input), outputWrapper, writeMetaData, isOne, true, componentConstraintTest, progress);
    numComponents = stats.numComponents;

    return std::move(outputVolumeBuilder).finalize();
}

static uint64_t toLinearPos(const tgt::svec3& pos, const tgt::svec3& dimensions) {
    return pos.x + dimensions.x*(pos.y + dimensions.y*pos.z);
}

struct IdVolumeInitializer {
    uint32_t id_;
    LZ4SliceVolumeVoxelBuilder<IdVolume::Value> storage_;
    SurfaceBuilder surface_;
    size_t numUnlabeledForegroundVoxels_;
    size_t numTotalVoxels_;

    tgt::svec3 offset_;
    tgt::svec3 size_;

    IdVolumeInitializer(uint32_t id, tgt::svec3 offset, tgt::svec3 size)
        : id_(id)
        , storage_(VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION), LZ4SliceVolumeMetadata(size))
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
        storage_.pushVoxel(v);
        ++numTotalVoxels_;
    }

    void pushSurfaceVoxel(const tgt::svec3& p) {
        surface_.push(toLinearPos(p - offset_, size_));
    }
};

static void initializeIdVolumes(LZ4SliceVolume<uint32_t>& branchIds, const LZ4SliceVolume<uint32_t>& holeIds, std::vector<IdVolumeInitializer>& initializers, ProgressReporter& progress) {
    TaskTimeLogger _("Initialize IdVolumes", tgt::Info);
    IdVolumeInitializationReader initializationReader(branchIds, holeIds);

    tgt::svec3 dim = branchIds.getDimensions();
    tgtAssert(holeIds.getDimensions() == dim, "Volume dimension mismatch");

    std::vector<boost::optional<BoundsHierarchy<size_t, IdVolumeInitializer*>>> initializer_finders;
    for(size_t z=0; z < dim.z; ++z) {
        std::vector<std::pair<IdVolumeInitializer*, tgt::SBounds>> initializer_finder_init;
        for(auto& initializer: initializers) {
            tgt::SBounds full_bounds = initializer.getBounds();
            if(full_bounds.getLLF().z <= z && z <= full_bounds.getURB().z) {
                tgt::SBounds bounds(tgt::svec3(full_bounds.getLLF().xy(), 0), tgt::svec3(full_bounds.getURB().xy(), 0));
                initializer_finder_init.emplace_back(&initializer, bounds);
            }
        }
        if(initializer_finder_init.empty()) {
            initializer_finders.emplace_back(boost::none);
        } else {
            initializer_finders.emplace_back(std::move(initializer_finder_init));
        }
    }

    for(size_t z=0; z < dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);
        initializationReader.advance();

        if(!initializer_finders.at(z)) {
            continue;
        }
        for(size_t y=0; y < dim.y; ++y) {
            for(size_t x=0; x < dim.x; ++x) {
                tgt::svec3 pos(x,y,z);

                tgt::svec3 slicePos(x,y,0);
                for(IdVolumeInitializer* initializer : initializer_finders.at(z)->findBounds(slicePos)) {
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

static LZ4SliceVolume<IdVolume::Value> flood(IdVolumeInitializer&& volume, ProgressReporter& progress) {
    tgtAssert(volume.numTotalVoxels_ == tgt::hmul(volume.size_), "Incomplete IdVolumeInitializer");

    auto lz4vol = std::move(volume.storage_).finalize();

    std::string filename = VoreenApplication::app()->getUniqueTmpFilePath(".raw");
    IdVolumeStorage idStorage(lz4vol, filename);
    IdVolume vol(std::move(idStorage), std::move(volume.surface_).finalize(), volume.numUnlabeledForegroundVoxels_);
    vol.floodFromLabels(progress, std::numeric_limits<size_t>::max());

    std::move(lz4vol).deleteFromDisk();

    LZ4SliceVolumeBuilder<IdVolume::Value> builder(VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION), lz4vol.getMetaData());
    for(size_t z = 0; z<vol.getDimensions().z; ++z) {
        builder.pushSlice(vol.data_->getSlice(z));
    }
    return std::move(builder).finalize();
}

struct IdVolumeFinalizer {
    uint32_t id_;
    LZ4SliceVolume<IdVolume::Value> vol_;
    LZ4SliceVolumeSliceCacher<IdVolume::Value> slices_;

    tgt::svec3 offset_;
    tgt::svec3 size_;

    IdVolumeFinalizer(IdVolumeInitializer&& volume, ProgressReporter& progress)
        : id_(volume.id_)
        , vol_(flood(std::move(volume), progress))
        , slices_(vol_)
        , offset_(volume.offset_)
        , size_(volume.size_)
    {
    }

    IdVolumeFinalizer(const IdVolumeFinalizer& other) = delete;
    IdVolumeFinalizer(IdVolumeFinalizer&& other)
        : id_(other.id_)
        , vol_(std::move(other.vol_))
        , slices_(vol_)
        , offset_(other.offset_)
        , size_(other.size_)
    {
    }

    ~IdVolumeFinalizer()
    {
        std::move(vol_).deleteFromDisk();
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
        tgt::svec3 pos = p - offset_;
        return slices_.getSlice(pos.z).voxel(pos.x, pos.y, 0);
    }
};

static void finalizeIdVolumes(LZ4SliceVolume<uint32_t>& branchIds, const LZ4SliceVolume<uint32_t>& holeIds, std::vector<IdVolumeFinalizer>& finalizers, ProgressReporter& progress) {
    TaskTimeLogger _("Finalize IdVolumes", tgt::Info);
    tgtAssert(branchIds.getDimensions() == holeIds.getDimensions(), "Invalid dimensions");

    const tgt::svec3 dim = branchIds.getDimensions();

    std::vector<boost::optional<BoundsHierarchy<size_t, IdVolumeFinalizer*>>> finalizer_finders;
    for(size_t z=0; z < dim.z; ++z) {
        std::vector<std::pair<IdVolumeFinalizer*, tgt::SBounds>> finalizer_finder_init;
        for(auto& finalizer: finalizers) {
            tgt::SBounds full_bounds = finalizer.getBounds();
            if(full_bounds.getLLF().z <= z && z <= full_bounds.getURB().z) {
                tgt::SBounds bounds(tgt::svec3(full_bounds.getLLF().xy(), 0), tgt::svec3(full_bounds.getURB().xy(), 0));
                finalizer_finder_init.emplace_back(&finalizer, bounds);
            }
        }
        if(finalizer_finder_init.empty()) {
            finalizer_finders.emplace_back(boost::none);
        } else {
            finalizer_finders.emplace_back(std::move(finalizer_finder_init));
        }
    }

    for(size_t z=0; z < dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);

        if(!finalizer_finders.at(z)) {
            continue;
        }

        auto branchSlice = branchIds.getWritableSlice(z);
        auto holeSlice = holeIds.loadSlice(z);

        for(size_t y=0; y < dim.y; ++y) {
            for(size_t x=0; x < dim.x; ++x) {
                tgt::svec3 pos(x,y,z);
                tgt::svec3 slicePos(x,y,0);

                for(IdVolumeFinalizer* finalizer : finalizer_finders.at(z)->findBounds(slicePos)) {
                    if(holeSlice.voxel(x,y,0) == finalizer->id_) {
                        uint32_t floodedId = finalizer->getValue(pos);
                        tgtAssert(floodedId != IdVolume::BACKGROUND_VALUE, "flooded label is background");
                        if(floodedId != IdVolume::UNLABELED_FOREGROUND_VALUE) {
                            branchSlice->voxel(x,y,0) = floodedId;
                        }
                    }
                }
            }
        }
    }
    progress.setProgress(1.0f);
}

struct UnfinishedRegions {
    LZ4SliceVolume<uint32_t> holeIds;
    std::map<uint32_t, tgt::SBounds> regions;

    void floodAllRegions(LZ4SliceVolume<uint32_t>& branchIds, ProgressReporter& progress) {
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

        initializeIdVolumes(branchIds, holeIds, initializers, subtaskReporters.get<0>());

        std::vector<IdVolumeFinalizer> finalizers;
        size_t i=0;
        for(auto& initializer: initializers) {
            SubtaskProgressReporter sp(subtaskReporters.get<1>(), tgt::vec2(i, i+1)/static_cast<float>(regions.size()));
            finalizers.emplace_back(std::move(initializer), sp);
            ++i;
        }

        finalizeIdVolumes(branchIds, holeIds, finalizers, subtaskReporters.get<2>());
    }
};

static UnfinishedRegions collectUnfinishedRegions(const LZ4SliceVolume<uint8_t>& input, const std::string& tmpVolumePath, ProgressReporter& progress) {
    TaskTimeLogger _("Collect unlabled regions", tgt::Info);

    std::map<uint32_t, tgt::SBounds> regions;

    StreamingComponents<0, CCANodeMetaData> sc;

    std::function<bool(const VolumeRAM* vol, tgt::svec3 pos)> isOne = [](const VolumeRAM* slice, tgt::svec3 pos) {
        return slice->getVoxelNormalized(pos) > 0.0f;
    };
    auto writeMetaData = [&regions] (uint32_t id, const CCANodeMetaData& metadata) {
        regions.emplace(id, metadata.bounds_);
    };

    auto componentConstraintTest = [] (const CCANodeMetaData&) {
        return true;
    };

    LZ4SliceVolumeBuilder<uint32_t> outputVolumeBuilder(tmpVolumePath, input.getMetaData().withRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint32_t>()));
    LZ4SliceVolumeCCABuilderWrapper outputWrapper(outputVolumeBuilder);

    sc.cca(LZ4SliceVolumeCCAInputWrapper(input), outputWrapper, writeMetaData, isOne, true, componentConstraintTest, progress);

    return UnfinishedRegions {
        std::move(outputVolumeBuilder).finalize(), regions
    };
}

static void addClippedOffRegionsToCritical(LZ4SliceVolume<uint8_t>& criticalVoxels, const LZ4SliceVolume<uint32_t>& ccaVolume, const ProtoVesselGraph& graph, size_t numComponents, const LZ4SliceVolume<uint8_t>& segmentation, ProgressReporter& progress) {
    TaskTimeLogger _("Add clipped off regions to critical voxels", tgt::Info);

    BranchIdVolumeReader ccaReader(ccaVolume, graph, numComponents, segmentation);
    auto dim = criticalVoxels.getDimensions();

    for(size_t z = 0; z < dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);

        auto slice = criticalVoxels.getWritableSlice(z);

        ccaReader.advance();

        for(size_t y = 0; y < dim.y; ++y) {
            for(size_t x = 0; x < dim.x; ++x) {
                tgt::ivec3 p(x,y,z);

                if(ccaReader.isObject(p) && !ccaReader.isValidEdgeId(ccaReader.getEdgeId(p))) {
                    slice->voxel(x, y, 0) = 1;
                }
            }
        }
    }
    progress.setProgress(1.0f);
}


std::unique_ptr<VesselGraph> createGraphFromMask(VesselGraphCreatorProcessedInput& input, VolumeMask&& skeleton, VolumeList& generatedVolumes, ProgressReporter& progress) {
    SubtaskProgressReporterCollection<8> subtaskReporters(progress);

    // Create new protograph
    std::unique_ptr<ProtoVesselGraph> protograph(nullptr);
    {
        std::unique_ptr<MetaDataCollector> mdc = cca(NeighborCountVoxelClassifier(skeleton), subtaskReporters.get<0>());
        protograph = mdc->createProtoVesselGraph(input.segmentation.getDimensions(), input.segmentation.getMetaData().getVoxelToWorldMatrix(), input.sampleMask, subtaskReporters.get<1>());
    }


    // Create an assignment edge <-> vessel component
    //  1. Split vessel components at nodes edges
    auto segmentationWithAndWithoutCriticalVoxels = splitSegmentationCriticalVoxels(
            VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION),
            VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION),
            *protograph,
            SkeletonClassReader(skeleton),
            input,
            subtaskReporters.get<2>());
    {
        // Drop VolumeMask and thus free the used non-volatile storage space.
        VolumeMask _dump(std::move(skeleton));
    }
    LZ4SliceVolume<uint8_t> segNoCriticalVoxels = std::move(segmentationWithAndWithoutCriticalVoxels.first);
    LZ4SliceVolume<uint8_t> segOnlyCriticalVoxels = std::move(segmentationWithAndWithoutCriticalVoxels.second);

    //  2. Perform a cca to distinguish components
    size_t numComponents;
    LZ4SliceVolume<uint32_t> branchIdSegmentation = createCCAVolume(segNoCriticalVoxels, VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION), numComponents, subtaskReporters.get<3>());
    if(input.saveDebugData) {
        generatedVolumes.add(std::move(segNoCriticalVoxels).toVolume().release());
    } else {
        std::move(segNoCriticalVoxels).deleteFromDisk();
    }

    // 4. Find unlabeled regions
    //
    // Add cut off unlabeled regions to the critical voxels (i.e., unlabeled regions) to allow flooding from valid edge label regions later
    addClippedOffRegionsToCritical(segOnlyCriticalVoxels, branchIdSegmentation, *protograph, numComponents, input.segmentation, subtaskReporters.get<4>());
    auto unfinishedRegions = collectUnfinishedRegions(segOnlyCriticalVoxels, VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION), subtaskReporters.get<5>());
    if(input.saveDebugData) {
        generatedVolumes.add(std::move(segOnlyCriticalVoxels).toVolume().release());
    } else {
        std::move(segOnlyCriticalVoxels).deleteFromDisk();
    }

    //  5. Flood remaining unlabeled regions locally
    unfinishedRegions.floodAllRegions(branchIdSegmentation, subtaskReporters.get<6>());

    if(input.saveDebugData) {
        generatedVolumes.add(std::move(unfinishedRegions.holeIds).toVolume().release());
    } else {
        std::move(unfinishedRegions.holeIds).deleteFromDisk();
    }

    //  6. Create a reader that reads edge ids from the underlying segmentation
    BranchIdVolumeReader ccaReader(branchIdSegmentation, *protograph, numComponents, input.segmentation);

    //  7. Create better VesselGraph from protograph and and the edge-id-segmentation
    auto output = protograph->createVesselGraph(ccaReader, input.sampleMask, subtaskReporters.get<7>());
    if(input.saveDebugData) {
        generatedVolumes.add(std::move(branchIdSegmentation).toVolume().release());
    } else {
        std::move(branchIdSegmentation).deleteFromDisk();
    }

    return output;
}

std::unique_ptr<VesselGraph> refineVesselGraph(VesselGraphCreatorProcessedInput& input, const VesselGraph& prevGraph, VolumeList& generatedVolumes, ProgressReporter& progress) {
    TaskTimeLogger _("Refinement iteration", tgt::Info);

    SubtaskProgressReporterCollection<3> subtaskReporters(progress);
    progress.setProgress(0.0f);

    // Refine previous graph
    auto isEdgeDeletable = [&input] (const VesselGraphEdge& edge) {
        return !edge.hasValidData() || edge.getVoxels().size() < input.minVoxelLength || edge.getElongation() < input.minElongation || edge.getRelativeBulgeSize() < input.minBulgeSize;
    };
    std::unique_ptr<VesselGraph> normalizedGraph = VesselGraphNormalization::removeEndEdgesRecursively(prevGraph, isEdgeDeletable);

    tgt::mat4 rwToVoxel;
    bool inverted = input.segmentation.getMetaData().getVoxelToWorldMatrix().invert(rwToVoxel);
    tgtAssert(inverted, "Matrix inversion failed!");

    // Create new voxelmask
    auto fixedForegroundMaskOnlyEndvoxels = GraphNodeVoxelReader::create<EndNodeVoxelExtractor>(*normalizedGraph, rwToVoxel, input.segmentation.getMetaData().getDimensions());
    VolumeMask mask(input.segmentation, input.sampleMask, fixedForegroundMaskOnlyEndvoxels, subtaskReporters.get<0>());
    addFixedForegroundPointsToMask(input.fixedForegroundPoints, mask);

    // Skeletonize new mask
    mask.skeletonize<VolumeMask::IMPROVED_NO_LINE_PRESERVATION>(std::numeric_limits<size_t>::max(), subtaskReporters.get<1>());

    return createGraphFromMask(input, std::move(mask), generatedVolumes, subtaskReporters.get<2>());
}

std::unique_ptr<VesselGraph> createInitialVesselGraph(VesselGraphCreatorProcessedInput& input, VolumeList& generatedVolumes, ProgressReporter& progress) {
    TaskTimeLogger _("Create initial VesselGraph", tgt::Info);
    SubtaskProgressReporterCollection<3> subtaskReporters(progress);
    progress.setProgress(0.0f);

    VolumeMask mask(input.segmentation, input.sampleMask, NoFixedForeground(), subtaskReporters.get<0>());
    addFixedForegroundPointsToMask(input.fixedForegroundPoints, mask);
    mask.skeletonize<VolumeMask::IMPROVED>(std::numeric_limits<size_t>::max(), subtaskReporters.get<1>());

    return createGraphFromMask(input, std::move(mask), generatedVolumes, subtaskReporters.get<2>());
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

    generatedVolumesOutport_.clear();
    generatedGraphsOutport_.clear();
    return VesselGraphCreatorInput(
        *segmentation,
        sampleMask,
        std::move(fixedForegroundPoints),
        binarizationThresholdSegmentationNormalized,
        numRefinementIterations_.get(),
        minVoxelLength_.get(),
        minElongation_.get(),
        minBulgeSize_.get(),
        saveDebugData_.get()
    );
}

static bool iterationMadeProgress(const VesselGraph& before, const VesselGraph& after) {
    return before.getBounds() != after.getBounds()
        || before.getNodes().size() != after.getNodes().size()
        || before.getEdges().size() != after.getEdges().size();
}

VesselGraphCreatorOutput VesselGraphCreator::compute(VesselGraphCreatorInput input, ProgressReporter& progressReporter) const {
    TaskTimeLogger _("Extract VesselGraph (total)", tgt::Info);

    SubtaskProgressReporterCollection<2> progressCollection(progressReporter, {0.01,0.99});

    VesselGraphCreatorProcessedInput processedInput(input, progressCollection.get<0>());

    float progressPerIteration = 1.0f/(input.numRefinementIterations+1);
    SubtaskProgressReporter initialProgress(progressCollection.get<1>(), tgt::vec2(0, progressPerIteration));

    std::unique_ptr<VolumeList> generatedVolumes(new VolumeContainer());
    std::unique_ptr<std::vector<VesselGraph>> generatedGraphs(new std::vector<VesselGraph>());
    LINFO("Begin graph extraction iteration: " << 0);
    std::unique_ptr<VesselGraph> graph = createInitialVesselGraph(processedInput, *generatedVolumes, initialProgress);
    LINFO("Inital graph: " << graph->getNodes().size() << " Nodes, " << graph->getEdges().size() << " Edges.");

    for(int i=0; i < input.numRefinementIterations; ++i) {
        LINFO("Begin graph extraction iteration: " << (i+1));
        SubtaskProgressReporter refinementProgress(progressCollection.get<1>(), progressPerIteration*tgt::vec2(i+1, i+2));
        try {
            std::unique_ptr<VesselGraph> prev_graph = std::move(graph);
            graph = refineVesselGraph(processedInput, *prev_graph, *generatedVolumes, refinementProgress);
            bool done = !iterationMadeProgress(*prev_graph, *graph);
            LINFO("Iteration " << (i+1) << ": " << graph->getNodes().size() << " Nodes, " << graph->getEdges().size() << " Edges.");
            if(input.saveDebugData) {
                generatedGraphs->push_back(std::move(*prev_graph));
            }
            if(done) {
                LINFO("Refinement reached fixed point after " << (i+1) << " iterations. Aborting early.");
                break;
            }
        } catch(tgt::IOException e) {
            LERROR("IO Exception occured in VesselGraphCreator compute thread");
            std::cout << e.what() << std::endl;
            break;
        }
    }
    std::move(processedInput.segmentation).deleteFromDisk();
    if(processedInput.sampleMask) {
        std::move(*processedInput.sampleMask).deleteFromDisk();
    }

    if(input.saveDebugData) {
        generatedGraphs->push_back(graph->clone());
    }

    return VesselGraphCreatorOutput {
        std::move(graph),
        std::move(generatedVolumes),
        std::move(generatedGraphs)
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
    generatedVolumesOutport_.setData(output.generatedVolumes.release());
    generatedGraphsOutport_.setData(output.generatedGraphs.release());
}
void VesselGraphCreator::adjustPropertiesToInput() {
    const VolumeBase* segmentation = segmentedVolumeInport_.getData();
    if(segmentation) {
        if(!segmentation->hasDerivedData<VolumeMinMax>()) {
            LINFO("Calculating VolumeMinMax. This may take a while...");
        }
        const VolumeMinMax* mm = segmentation->getDerivedData<VolumeMinMax>();
        bool minMaxWereEqual = binarizationThresholdSegmentation_.getMinValue() == binarizationThresholdSegmentation_.getMaxValue();
        binarizationThresholdSegmentation_.setMinValue(mm->getMin());
        binarizationThresholdSegmentation_.setMaxValue(mm->getMax());
        binarizationThresholdSegmentation_.adaptDecimalsToRange(3);

        if(minMaxWereEqual) {
            binarizationThresholdSegmentation_.set((binarizationThresholdSegmentation_.getMinValue() + binarizationThresholdSegmentation_.getMaxValue())*0.5f);
        }

        tmpStorageSizeInfo_.set(formatMemorySize(segmentation->getNumVoxels()/8*2 /*2 bits per voxel */));
    }
}

} // namespace voreen
