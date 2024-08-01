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

#include "vesselgraphcreator.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/utils/stringutils.h"

#include "../algorithm/streaminggraphcreation.h"
#include "../algorithm/volumemask.h"
#include "../algorithm/idvolume.h"
#include "../algorithm/vesselgraphrefinement.h"
#include "../../bigdataimageprocessing/algorithm/intervalwalker.h"
#include "../datastructures/kdtree.h"

#include "modules/bigdataimageprocessing/datastructures/lz4slicevolume.h"
#include "modules/bigdataimageprocessing/processors/connectedcomponentanalysis.h"

#include "tgt/bounds.h"
#include "tgt/filesystem.h"

#include <vector>
#include <memory>
#include <unordered_set>

namespace voreen {

const std::string VesselGraphCreator::loggerCat_("voreen.vesselnetworkanalysis.vesselgraphcreator");

VesselGraphCreator::VesselGraphCreator()
    : segmentedVolumeInport_(Port::INPORT, "vesselgraphcreator.segmentedVolume.inport", "Segmentation Volume")
    , sampleMaskInport_(Port::INPORT, "vesselgraphcreator.samplemask.inport", "Sample Mask (optional)")
    , fixedForegroundPointInport_(Port::INPORT, "vesselgraphcreator.fixedForegroundPointInport", "Fixed Foreground Points", false, Processor::INVALID_RESULT)
    , graphOutport_(Port::OUTPORT, "vesselgraphcreator_graph.outport", "Graph", false, Processor::VALID)
    , generatedVolumesOutport_(Port::OUTPORT, "generatedSkeletons.outport", "Generated Volumes", false, Processor::VALID)
    , generatedGraphsOutport_(Port::OUTPORT, "generatedGraphs.outport", "Generated Graphs", false, Processor::VALID)
    , numRefinementIterations_("numRefinementIterations", "Refinement Iterations", 100, 0, 100, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_ADVANCED)
    , minVoxelLength_("minVoxelLength", "Min Voxel Length", 0, 0, 50, Processor::INVALID_RESULT, IntProperty::STATIC, Property::LOD_DEBUG)
    , minElongation_("minElongation", "Minimum Elongation", 0, 0, 5, Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_DEBUG)
    , minBulgeSize_("minBulgeSize", "Minimum Bulge Size", 1.0, 0, 10, Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_DEFAULT)
    , saveDebugData_("saveDebugData", "Generate Debug Data", false, Processor::INVALID_RESULT, Property::LOD_DEBUG)
    , binarizationThresholdSegmentation_("binarizationThresholdSegmentation", "Binarization Threshold (Segmentation)", 0.5f, std::numeric_limits<float>::lowest(), std::numeric_limits<float>::max(), Processor::INVALID_RESULT, FloatProperty::STATIC, Property::LOD_ADVANCED)
    , tmpStorageSizeInfo_("tmpStorageSizeInfo", "Required temporary storage (estimated)", "")
{
    addPort(segmentedVolumeInport_);
    addPort(sampleMaskInport_);
    addPort(fixedForegroundPointInport_);
    addPort(graphOutport_);
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
        (graphOutport_.isReady() || generatedVolumesOutport_.isReady() || generatedGraphsOutport_.isReady());
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

static void fixEndVoxelsInMask(const VesselGraph& graph, const tgt::mat4& rwToVoxel, VolumeMask& mask) {
    for(auto& node : graph.getNodes()) {
        if(node.isEndNode()) {
            auto p = tgt::iround(transform(rwToVoxel, node.pos_));
            mask.set(tgt::svec3(p), VolumeMaskValue::FIXED_OBJECT);
        }
    }
}

struct EdgeVoxelRef {
    typedef float CoordType;

    tgt::vec3 rwPos;
    float padding_;
    const ProtoVesselGraphEdge* edge;

    inline const tgt::Vector3<CoordType>& getPos() const {
        return rwPos;
    }

    EdgeVoxelRef(tgt::vec3 rwPos, const ProtoVesselGraphEdge& edge)
        : rwPos(rwPos)
        , padding_(0.0)
        , edge(&edge)
    {
    }
    EdgeVoxelRef(const EdgeVoxelRef& other)
        : rwPos(other.rwPos)
        , padding_(0.0)
        , edge(other.edge)
    {
    }
    EdgeVoxelRef& operator=(const EdgeVoxelRef& other) {
        rwPos = other.rwPos;
        edge = other.edge;
        return *this;
    }
};

static LZ4SliceVolume<uint32_t> createClosestIDVolume(const std::string& tmpPath, ProtoVesselGraph& graph, const VesselGraphCreatorProcessedInput& input, ProgressReporter& progress) {
    TaskTimeLogger _("Create closest id volume", tgt::Info);

    const auto voxelToRw = input.segmentation.getMetaData().getVoxelToWorldMatrix();
    const auto dimensions = input.segmentation.getDimensions();

    static_kdtree::ElementArrayBuilder<EdgeVoxelRef> finderBuilder(VoreenApplication::app()->getUniqueTmpFilePath(".kdtreestorage"));

    for(const auto& edge : graph.edges_.asArray()) {
        for(auto& rwvoxel : edge.voxels()) {
            finderBuilder.push(EdgeVoxelRef(rwvoxel, edge));
        }
    }

    static_kdtree::Tree<EdgeVoxelRef> finder(VoreenApplication::app()->getUniqueTmpFilePath(".kdtree"), std::move(finderBuilder));

    LZ4SliceVolumeBuilder<uint32_t> outputIDs(tmpPath, input.segmentation.getMetaData().withRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint32_t>()));

    const size_t max_chunk_size = 8;

    auto calc_num_chunks = [] (size_t dim, size_t size) {
        if((dim % size) == 0) {
            return dim/size;
        } else {
            return dim/size+1;
        }
    };

    auto calc_chunk_size = [&] (size_t chunk, size_t num_chunks, size_t dimensions) {
        size_t size = chunk != num_chunks - 1 || (dimensions % max_chunk_size == 0) ? max_chunk_size : dimensions % max_chunk_size;
        tgtAssert(size > 0, "Invalid chunk size");
        return size;
    };

    const size_t num_chunks_x = calc_num_chunks(dimensions.x, max_chunk_size);
    const size_t num_chunks_y = calc_num_chunks(dimensions.y, max_chunk_size);
    const size_t num_chunks_z = calc_num_chunks(dimensions.z, max_chunk_size);

    for(size_t chunk_z = 0; chunk_z < num_chunks_z; ++chunk_z) {
        progress.setProgress(static_cast<float>(chunk_z)/num_chunks_z);

        const size_t chunk_size_z = calc_chunk_size(chunk_z, num_chunks_z, dimensions.z);
        auto outputSlab = outputIDs.getNextWriteableSlab(chunk_size_z);

        size_t slab_begin = chunk_z*max_chunk_size;
        size_t slab_end = slab_begin + chunk_size_z;
        auto slab = input.segmentation.loadSlab(slab_begin, slab_end);

        for(size_t chunk_y = 0; chunk_y < num_chunks_y; ++chunk_y) {
            const size_t chunk_size_y = calc_chunk_size(chunk_y, num_chunks_y, dimensions.y);

            for(size_t chunk_x = 0; chunk_x < num_chunks_x; ++chunk_x) {
                const size_t chunk_size_x = calc_chunk_size(chunk_x, num_chunks_x, dimensions.x);
                for(size_t dz = 0; dz<chunk_size_z; ++dz) {
                    size_t z = chunk_z*max_chunk_size + dz;

                    for(size_t dy = 0; dy<chunk_size_y; ++dy) {
                        size_t y = chunk_y*max_chunk_size + dy;

                        for(size_t dx = 0; dx<chunk_size_x; ++dx) {
                            size_t x = chunk_x*max_chunk_size + dx;

                            const tgt::svec3 p(x,y,dz);
                            uint32_t label;
                            if(slab.voxel(p) > 0) {
                                tgt::ivec3 ipos(x,y,z);
                                tgt::vec3 rwpos = transform(voxelToRw, ipos);

                                auto result = finder.findAllNearest(rwpos);
                                if(!result.found()) {
                                    label = 0xFEEDBEEF; //Doesn't really matter. only happens in edge cases anyway
                                } else {
                                    auto min_elm = std::min_element(result.elements_.begin(), result.elements_.end(), [] (const EdgeVoxelRef*& e1, const EdgeVoxelRef*& e2) {
                                            auto p1 = e1->getPos();
                                            auto p2 = e2->getPos();
                                            float cmp = p1.x - p2.x;
                                            cmp = (cmp == 0) ? p1.y - p2.y : cmp;
                                            cmp = (cmp == 0) ? p1.z - p2.z : cmp;
                                            return cmp < 0;
                                            });
                                    tgtAssert(min_elm != result.elements_.end(), "Empty result set");
                                    label = (*min_elm)->edge->id_.raw();
                                }
                            } else {
                                label = IdVolume::BACKGROUND_VALUE;
                            }
                            outputSlab->voxel(p) = label;
                        }
                    }
                }
            }
        }
    }

    progress.setProgress(1.0f);
    return std::move(outputIDs).finalize();
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
        mask.set(tgt::svec3(tgt::iround(p)), VolumeMaskValue::FIXED_OBJECT);
    }
}

template<typename T>
struct LZ4SliceVolumeCCAInputWrapper {
    LZ4SliceVolumeCCAInputWrapper(const LZ4SliceVolume<T>& vol)
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
        return new VolumeAtomic<T>(volume_.loadSlice(z));
    }
    const LZ4SliceVolume<T>& volume_;
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

struct CCAUint32Label {
    uint32_t label_;

    CCAUint32Label(uint32_t label)
        : label_(label)
    {
    }
    bool operator==(const CCAUint32Label& other) const {
        return label_ == other.label_;
    }
    bool operator!=(const CCAUint32Label& other) const {
        return label_ != other.label_;
    }
};

static LZ4SliceVolume<uint32_t> createCCAVolume(const LZ4SliceVolume<uint32_t>& input, const std::string& tmpVolumePath, size_t& numComponents, ProgressReporter& progress) {
    TaskTimeLogger _("Create CCA volume", tgt::Info);

    typedef StreamingComponents<0, NoopMetadata, CCAUint32Label, VolumeAtomic<uint32_t>> SC;
    SC sc;

    typename SC::getClassFunc getClass = [](const VolumeAtomic<uint32_t>& slice, tgt::svec3 pos) {
        uint32_t label = slice.voxel(pos);
        return label != IdVolume::BACKGROUND_VALUE ? boost::optional<CCAUint32Label>(label) : boost::none;
    };
    auto writeMetaData = [] (uint32_t, const NoopMetadata&) {};

    auto componentConstraintTest = [] (const NoopMetadata& metaData) {
        return true;
    };

    LZ4SliceVolumeBuilder<uint32_t> outputVolumeBuilder(tmpVolumePath, input.getMetaData().withRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint32_t>()));
    LZ4SliceVolumeCCABuilderWrapper outputWrapper(outputVolumeBuilder);

    auto stats = sc.cca(LZ4SliceVolumeCCAInputWrapper<uint32_t>(input), outputWrapper, writeMetaData, getClass, true, componentConstraintTest, progress);
    numComponents = stats.numComponents;

    return std::move(outputVolumeBuilder).finalize();
}

static uint64_t toLinearPos(const tgt::svec3& pos, const tgt::svec3& dimensions) {
    return pos.x + dimensions.x*(pos.y + dimensions.y*pos.z);
}

struct IdVolumeInitializer {
    uint32_t id_;
    LZ4SliceVolumeVoxelBuilder<IdVolume::Value> storage_;
    NoFileSurfaceBuilder surface_;
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


    std::vector<Interval<size_t, IdVolumeInitializer*>> intervals;
    for(auto& initializer: initializers) {
        tgt::SBounds full_bounds = initializer.getBounds();

        // Here and below: the upper bound of tgt::SBounds is inclusive, but the intervals of
        // intervalwalker have an exclusive upper bound
        intervals.emplace_back(full_bounds.getLLF().z, full_bounds.getURB().z + 1, &initializer);
    }
    IntervalWalker<size_t, IdVolumeInitializer*> zwalker(0, std::move(intervals));

    for(size_t z=0; z < dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);
        initializationReader.advance();

        std::vector<Interval<size_t, IdVolumeInitializer*>> intervals;

        auto zit = zwalker.next();
        for(auto it = zit.next(); it != zit.end(); it = zit.next()) {
            auto& initializer = it->value;
            tgt::SBounds full_bounds = initializer->getBounds();
            intervals.emplace_back(full_bounds.getLLF().y, full_bounds.getURB().y + 1, initializer);
        }
        if(intervals.empty()) {
            continue;
        }
        IntervalWalker<size_t, IdVolumeInitializer*> ywalker(0, std::move(intervals));

        for(size_t y=0; y < dim.y; ++y) {

            std::vector<Interval<size_t, IdVolumeInitializer*>> intervals;

            auto yit = ywalker.next();
            for(auto it = yit.next(); it != yit.end(); it = yit.next()) {
                auto& initializer = it->value;
                tgt::SBounds full_bounds = initializer->getBounds();
                intervals.emplace_back(full_bounds.getLLF().x, full_bounds.getURB().x + 1, initializer);
            }
            if(intervals.empty()) {
                continue;
            }
            IntervalWalker<size_t, IdVolumeInitializer*> xwalker(0, std::move(intervals));

            for(size_t x=0; x < dim.x; ++x) {
                tgt::svec3 pos(x,y,z);

                tgt::svec3 slicePos(x,y,0);
                auto xit = xwalker.next();
                for(auto it = xit.next(); it != xit.end(); it = xit.next()) {
                    auto& initializer = it->value;
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

    std::vector<Interval<size_t, IdVolumeFinalizer*>> intervals;
    for(auto& finalizer: finalizers) {
        tgt::SBounds full_bounds = finalizer.getBounds();
        // Here and below: the upper bound of tgt::SBounds is inclusive, but the intervals of
        // intervalwalker have an exclusive upper bound
        intervals.emplace_back(full_bounds.getLLF().z, full_bounds.getURB().z + 1, &finalizer);
    }
    IntervalWalker<size_t, IdVolumeFinalizer*> zwalker(0, std::move(intervals));

    for(size_t z=0; z < dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);

        std::vector<Interval<size_t, IdVolumeFinalizer*>> intervals;
        auto zit = zwalker.next();
        for(auto it = zit.next(); it != zit.end(); it = zit.next()) {
            auto& finalizer = it->value;
            tgt::SBounds full_bounds = finalizer->getBounds();
            intervals.emplace_back(full_bounds.getLLF().y, full_bounds.getURB().y + 1, finalizer);
        }
        if(intervals.empty()) {
            continue;
        }
        IntervalWalker<size_t, IdVolumeFinalizer*> ywalker(0, std::move(intervals));

        auto branchSlice = branchIds.getWriteableSlice(z);
        auto holeSlice = holeIds.loadSlice(z);

        for(size_t y=0; y < dim.y; ++y) {

            std::vector<Interval<size_t, IdVolumeFinalizer*>> intervals;
            auto yit = ywalker.next();
            for(auto it = yit.next(); it != yit.end(); it = yit.next()) {
                auto& finalizer = it->value;
                tgt::SBounds full_bounds = finalizer->getBounds();
                intervals.emplace_back(full_bounds.getLLF().x, full_bounds.getURB().x + 1, finalizer);
            }
            if(intervals.empty()) {
                continue;
            }
            IntervalWalker<size_t, IdVolumeFinalizer*> xwalker(0, std::move(intervals));

            for(size_t x=0; x < dim.x; ++x) {
                tgt::svec3 pos(x,y,z);
                tgt::svec3 slicePos(x,y,0);

                auto xit = xwalker.next();
                for(auto it = xit.next(); it != xit.end(); it = xit.next()) {
                    auto& finalizer = it->value;
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

static void mapEdgeIds(LZ4SliceVolume<uint32_t>& regions, size_t numComponents, ProtoVesselGraph& graph, ProgressReporter& progress) {
    TaskTimeLogger _("Map CCA to edge IDs", tgt::Info);
    std::vector<uint32_t> ccaToEdgeIdTable(numComponents+1, IdVolume::UNLABELED_FOREGROUND_VALUE);
    ccaToEdgeIdTable[0] = IdVolume::BACKGROUND_VALUE;

    //populate ccaToEdgeIdTable
    {
        std::vector<std::pair<uint32_t, tgt::svec3>> query_positions;
        for(const auto& edge: graph.edges_.asArray()) {
            if(edge.voxels_.empty()) {
                query_positions.push_back(std::make_pair(edge.id_.raw(), tgt::svec3(-1)));
            } else {
                query_positions.push_back(std::make_pair(edge.id_.raw(), edge.voxels_[0]));
            }
        }
        std::sort(query_positions.begin(), query_positions.end(), [] (const std::pair<uint32_t, tgt::svec3>& p1, const std::pair<uint32_t, tgt::svec3>& p2) {
                return p1.second.z < p2.second.z;
                });

        LZ4SliceVolumeReader<uint32_t,0> reader(regions);
        reader.seek(0);
        for(auto& pair : query_positions) {
            tgt::svec3& p = pair.second;
            if(p.z == (size_t)-1) {
                continue;
            }

            tgtAssert(p.z >= (size_t)reader.getCurrentZPos(), "invalid reader pos");
            while(p.z > (size_t)reader.getCurrentZPos()) {
                reader.advance();
            }
            auto ccaindex = reader.getVoxelRelative(p.xy(), 0);
            tgtAssert(ccaindex, "Read invalid voxel");
            tgtAssert(*ccaindex != 0, "sample at background voxel");

            tgtAssert(ccaToEdgeIdTable[*ccaindex] == IdVolume::UNLABELED_FOREGROUND_VALUE, "Duplicate edge mapping");

            ccaToEdgeIdTable[*ccaindex] = pair.first;
        }
    }


    // Map ids in regions inplace
    auto dim = regions.getDimensions();
    for(size_t z=0; z<dim.z; ++z) {
        auto slice = regions.getWriteableSlice(z);
        for(size_t y=0; y<dim.y; ++y) {
            for(size_t x=0; x<dim.x; ++x) {
                tgt::svec3 p(x,y,0);
                uint32_t& vox = slice->voxel(p);
                vox = ccaToEdgeIdTable.at(vox);
            }
        }
    }
}

static UnfinishedRegions collectUnfinishedRegions(const LZ4SliceVolume<uint32_t>& input, const std::string& tmpVolumePath, ProgressReporter& progress) {
    TaskTimeLogger _("Collect unlabled regions", tgt::Info);

    std::map<uint32_t, tgt::SBounds> regions;

    typedef StreamingComponents<0, CCANodeMetaData, CCAVoidLabel, VolumeAtomic<uint32_t>> SC;
    SC sc;

    typename SC::getClassFunc getClass = [](const VolumeAtomic<uint32_t>& slice, tgt::svec3 pos) {
        return slice.voxel(pos) == IdVolume::UNLABELED_FOREGROUND_VALUE ? CCAVoidLabel::some() : boost::none;
    };
    auto writeMetaData = [&regions] (uint32_t id, const CCANodeMetaData& metadata) {
        regions.emplace(id, metadata.bounds_);
    };

    auto componentConstraintTest = [] (const CCANodeMetaData&) {
        return true;
    };

    LZ4SliceVolumeBuilder<uint32_t> outputVolumeBuilder(tmpVolumePath, input.getMetaData().withRealWorldMapping(RealWorldMapping::createDenormalizingMapping<uint32_t>()));
    LZ4SliceVolumeCCABuilderWrapper outputWrapper(outputVolumeBuilder);

    sc.cca(LZ4SliceVolumeCCAInputWrapper<uint32_t>(input), outputWrapper, writeMetaData, getClass, true, componentConstraintTest, progress);

    return UnfinishedRegions {
        std::move(outputVolumeBuilder).finalize(), regions
    };
}

std::unique_ptr<VesselGraph> createGraphFromMask(VesselGraphCreatorProcessedInput& input, VolumeMask&& skeleton, std::vector<std::unique_ptr<VolumeBase>>& generatedVolumes, ProgressReporter& progress) {
    SubtaskProgressReporterCollection<8> subtaskReporters(progress);

    // Create new protograph
    std::unique_ptr<ProtoVesselGraph> protograph(nullptr);
    {
        std::unique_ptr<MetaDataCollector> mdc = cca(NeighborCountVoxelClassifier(skeleton), subtaskReporters.get<0>());
        protograph = mdc->createProtoVesselGraph(input.segmentation.getDimensions(), input.segmentation.getMetaData().getVoxelToWorldMatrix(), input.sampleMask, subtaskReporters.get<1>());
    }

    {
        // Drop VolumeMask and thus free the used non-volatile storage space.
        VolumeMask _dump(std::move(skeleton));
    }


    //  1. Create an assignment edge <-> vessel component
    auto closestIDs = createClosestIDVolume(
            VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION),
            *protograph,
            input,
            subtaskReporters.get<2>());

    //  2. Perform a cca to distinguish components
    size_t numComponents;
    LZ4SliceVolume<uint32_t> branchIdSegmentation = createCCAVolume(closestIDs, VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION), numComponents, subtaskReporters.get<3>());
    if(input.saveDebugData) {
        generatedVolumes.push_back(std::move(closestIDs).toVolume());
    } else {
        std::move(closestIDs).deleteFromDisk();
    }

    //  3. Change branchIdSegmentation inplace and assign proper edge ids from protograph. Unassigned regions get id UNLABELED_REGION_ID;
    mapEdgeIds(branchIdSegmentation, numComponents, *protograph, subtaskReporters.get<4>());

    //  4. Collect bounding boxes from unfinished (i.e., cut-off) regions.
    auto unfinishedRegions = collectUnfinishedRegions(branchIdSegmentation, VoreenApplication::app()->getUniqueTmpFilePath("." + LZ4SliceVolumeBase::FILE_EXTENSION), subtaskReporters.get<5>());

    //  5. Flood remaining unlabeled regions locally
    unfinishedRegions.floodAllRegions(branchIdSegmentation, subtaskReporters.get<6>());

    if(input.saveDebugData) {
        generatedVolumes.push_back(std::move(unfinishedRegions.holeIds).toVolume());
    } else {
        std::move(unfinishedRegions.holeIds).deleteFromDisk();
    }

    //  6. Create better VesselGraph from protograph and the edge-id-segmentation
    auto output = protograph->createVesselGraph(branchIdSegmentation, input.sampleMask, subtaskReporters.get<7>());
    if(input.saveDebugData) {
        generatedVolumes.push_back(std::move(branchIdSegmentation).toVolume());
    } else {
        std::move(branchIdSegmentation).deleteFromDisk();
    }

    return output;
}

std::unique_ptr<VesselGraph> refineVesselGraph(VesselGraphCreatorProcessedInput& input, const VesselGraph& prevGraph, std::vector<std::unique_ptr<VolumeBase>>& generatedVolumes, ProgressReporter& progress) {
    TaskTimeLogger _("Refinement iteration", tgt::Info);

    SubtaskProgressReporterCollection<3> subtaskReporters(progress);
    progress.setProgress(0.0f);

    // Refine previous graph
    auto isEdgeDeletable = [&input] (const VesselGraphEdge& edge) {
        return !edge.hasValidData() || edge.getVoxels().size() < (size_t)input.minVoxelLength || edge.getElongation() < input.minElongation || edge.getRelativeBulgeSize() < input.minBulgeSize;
    };
    std::unique_ptr<VesselGraph> refinedGraph = VesselGraphRefinement::removeEndEdgesRecursively(prevGraph, isEdgeDeletable);

    tgt::mat4 rwToVoxel;
    bool inverted = input.segmentation.getMetaData().getVoxelToWorldMatrix().invert(rwToVoxel);
    tgtAssert(inverted, "Matrix inversion failed!");

    // Create new voxelmask
    VolumeMask mask(input.segmentation, input.sampleMask, subtaskReporters.get<0>());

    fixEndVoxelsInMask(*refinedGraph, rwToVoxel, mask);
    refinedGraph.reset(); // free graph (and associated ressources) early
    addFixedForegroundPointsToMask(input.fixedForegroundPoints, mask);

    // Skeletonize new mask
    mask.skeletonize<VolumeMask::IMPROVED_NO_LINE_PRESERVATION>(std::numeric_limits<size_t>::max(), subtaskReporters.get<1>());

    return createGraphFromMask(input, std::move(mask), generatedVolumes, subtaskReporters.get<2>());
}

std::unique_ptr<VesselGraph> createInitialVesselGraph(VesselGraphCreatorProcessedInput& input, std::vector<std::unique_ptr<VolumeBase>>& generatedVolumes, ProgressReporter& progress) {
    TaskTimeLogger _("Create initial VesselGraph", tgt::Info);
    SubtaskProgressReporterCollection<3> subtaskReporters(progress);
    progress.setProgress(0.0f);

    VolumeMask mask(input.segmentation, input.sampleMask, subtaskReporters.get<0>());
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
    LINFO("Starting graph extraction from volume '" << input.segmentation.getOrigin().getURL() << "' (dimensions: " << input.segmentation.getDimensions() << ")");

    std::vector<std::unique_ptr<VolumeBase>> generatedVolumes;
    std::unique_ptr<std::vector<VesselGraph>> generatedGraphs(new std::vector<VesselGraph>());
    std::unique_ptr<VesselGraph> graph(nullptr);
    try {
        TaskTimeLogger _("Extract VesselGraph (total)", tgt::Info);
        SubtaskProgressReporterCollection<2> progressCollection(progressReporter, {0.01f,0.99f});

        VesselGraphCreatorProcessedInput processedInput(input, progressCollection.get<0>());

        float progressPerIteration = 1.0f/(input.numRefinementIterations+1);
        SubtaskProgressReporter initialProgress(progressCollection.get<1>(), tgt::vec2(0, progressPerIteration));

        graph = createInitialVesselGraph(processedInput, generatedVolumes, initialProgress);
        LINFO("Inital graph: " << graph->getNodes().size() << " Nodes, " << graph->getEdges().size() << " Edges.");

        for(int i=0; i < input.numRefinementIterations; ++i) {
            SubtaskProgressReporter refinementProgress(progressCollection.get<1>(), progressPerIteration*tgt::vec2(i+1, i+2));

            if(input.saveDebugData) {
                generatedGraphs->push_back(graph->clone());
            }
            std::unique_ptr<VesselGraph> next_graph = refineVesselGraph(processedInput, *graph, generatedVolumes, refinementProgress);

            bool done = !iterationMadeProgress(*graph, *next_graph);

            graph = std::move(next_graph);

            LINFO("Iteration " << (i+1) << ": " << graph->getNodes().size() << " Nodes, " << graph->getEdges().size() << " Edges.");
            if(done) {
                LINFO("Refinement reached fixed point after " << (i+1) << " iterations. Aborting early.");
                break;
            }
        }

        if(input.saveDebugData) {
            generatedGraphs->push_back(graph->clone());
        }

        std::move(processedInput.segmentation).deleteFromDisk();
        if(processedInput.sampleMask) {
            std::move(*processedInput.sampleMask).deleteFromDisk();
        }

    // This does not work with the current AsyncComputeProcessor/interrupt
    // model. We will need to provide a proper way for producing intermediate
    // results later.
    //} catch(InterruptionException&) {
    //    // Finish up work and save out results collected so far
    //    LINFO("Graph extraction interrupted.");
    } catch(tgt::IOException& e) {
        LERROR("IO Exception occured in VesselGraphCreator compute thread");
        std::cout << e.what() << std::endl;
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
    lastGeneratedVolumes_ = std::move(output.generatedVolumes);
    VolumeList* volList = new VolumeList();
    for(auto& vol : lastGeneratedVolumes_) {
        volList->add(vol.get());
    }
    generatedVolumesOutport_.setData(volList);
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
