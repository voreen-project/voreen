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

#ifndef VRN_STREAMING_COMPONENTS_VESSELTOPOLOGY_H
#define VRN_STREAMING_COMPONENTS_VESSELTOPOLOGY_H

#include "voreen/core/datastructures/volume/volume.h"
#include "modules/hdf5/io/hdf5filevolume.h"
#include "voreen/core/io/progressreporter.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include "../datastructures/vesselgraph.h"
#include "../datastructures/protovesselgraph.h"
#include "../algorithm/volumemask.h"
#include "../util/tasktimelogger.h"
#include "../util/kdtreebuilder.h"
#include "custommodules/bigdataimageprocessing/volumefiltering/slicereader.h"
#include "custommodules/bigdataimageprocessing/datastructures/lz4slicevolume.h"

#include "tgt/vector.h"
#include "tgt/memory.h"

#include <memory>

/////////////////////////////////////////////////////////////////////
// Helper classes and functions for streaming vessel graph extraction.
// See vesselgraph.h and vesselgraphextractor.h
//////////////////////////////////////////////////////////////////////


namespace voreen {

// Helper functions
inline static tgt::vec3 transform(const tgt::mat4& m, const tgt::vec3& vec) {
    return (m*tgt::vec4(vec, 1)).xyz();
}
inline static tgt::vec3 transform(const tgt::mat4& m, const tgt::svec3& vec) {
    return (m*tgt::vec4(vec, 1)).xyz();
}
inline static tgt::vec3 transform(const tgt::mat4& m, const tgt::ivec3& vec) {
    return (m*tgt::vec4(vec, 1)).xyz();
}
static bool are26Neighbors(const tgt::svec3& v1, const tgt::svec3& v2) {
    return std::abs(static_cast<int>(v1.x) - static_cast<int>(v2.x)) <= 1
        && std::abs(static_cast<int>(v1.y) - static_cast<int>(v2.y)) <= 1
        && std::abs(static_cast<int>(v1.z) - static_cast<int>(v2.z)) <= 1;
}

// RunPosition ------------------------------------------------------------------------------
// The position of a single run (see RunLeaf)
struct RunPosition {
    RunPosition(tgt::svec2 yz, size_t xlow, size_t xhigh);
    tgt::svec3 begin() const;
    tgt::svec3 end() const;

    const size_t y_;
    const size_t z_;
    const size_t xlow_;
    const size_t xhigh_;
};

/// RunTrees -----------------------------------------------------------------------------
// Regular voxels that make up a branch in a vessel graph are stored in tree to allow
// for fast merging of components.
struct RunTree {
    virtual ~RunTree() {}
    //Collect als voxels from this tree
    virtual void collectVoxels(std::vector<tgt::svec3>&, bool inverted = false) const = 0;

    //Signal that the order of regular voxel under this node has been inverted
    void invert();
    bool inverted_ = false;
};

// RunNode: An internal node in a tree of (regular voxel) runs that make up a branch
struct RunNode : public RunTree {
    RunNode(std::unique_ptr<RunTree>&& left, std::unique_ptr<RunTree>&& right);
    void collectVoxels(std::vector<tgt::svec3>&, bool inverted = false) const;
    std::unique_ptr<RunTree> left_;
    std::unique_ptr<RunTree> right_;
};

// RunLeaf: A leaf in a run tree
struct RunLeaf : public RunTree {
    RunLeaf(const RunPosition&);
    void collectVoxels(std::vector<tgt::svec3>&, bool inverted = false) const;
    RunPosition run_;
};

struct SkeletonClassReader {
public:
    SkeletonClassReader(const VolumeMask& skeleton);
    virtual ~SkeletonClassReader() {}

    // Get the class of the specified voxel within the active slice
    uint8_t getClass(const tgt::ivec2& xypos) const;

    void advance();

    const tgt::svec3& getDimensions() const;

private:
    // Direct accessor for easy implementation, but private so that api may be changed
    bool isObject(const tgt::ivec3& xyz) const;

    const VolumeMask& skeleton_;
    int z_;
};

/// NeighborCountVoxelClassifier -----------------------------------------------------------
// A slice reader for binary volumes that can classify voxels as:
//  end voxels (< 2 neighbors)
//  regular voxels (= 2 neighbors)
//  branch voxels (> 2 neighbors)
//
//  See slicereader.h for more details on SliceReaders
class NeighborCountVoxelClassifier : public SkeletonClassReader {
public:
    NeighborCountVoxelClassifier(const VolumeMask& skeleton);
    virtual ~NeighborCountVoxelClassifier();

    // 0: no predefined id
    uint32_t getPredefinedComponentId(const tgt::ivec2& xypos) const;
};

/// NeighborCountAndBranchSegmentationVoxelClassifier --------------------------------------
// A slice reader for binary volumes that can classify voxels as:
//  end voxels (< 2 neighbors)
//  regular voxels (= 2 neighbors)
//  branch voxels (> 2 neighbors or outside of the branchSegmentation)
//
//  See slicereader.h for more details on SliceReaders
class NeighborCountAndBranchSegmentationVoxelClassifier {
public:
    NeighborCountAndBranchSegmentationVoxelClassifier(const VolumeMask& skeleton, const VolumeBase& branchSegmentation);
    virtual ~NeighborCountAndBranchSegmentationVoxelClassifier();

    // Get the class of the specified voxel within the active slice
    uint8_t getClass(const tgt::ivec2& xypos) const;
    const tgt::svec3& getDimensions() const;
    // 0: no predefined id
    uint32_t getPredefinedComponentId(const tgt::ivec2& xypos) const;

    void advance();
private:
    SkeletonClassReader skeletonClassReader_;
    VolumeSliceReader branchSegmentationReader_;
};

/// MetaData -----------------------------------------------------------------------------
// Metadata for streaming graph extraction must:
//  - be constructible from a runposition (i.e., metadata for a new (for now) component)
//  - be constructible from two other metadata types (i.e., merge)
//  - be able to consume the metadata of another merged component

// Metadata for end voxel components
struct EndData {
    EndData(EndData&&, EndData&&);
    EndData(const RunPosition&, uint32_t predeterminedComponentId);
    void consume(EndData&& rhs);

    tgt::svec3 pos_;
};

// Metadata for branch voxel components
struct BranchData {
    BranchData(BranchData&&, BranchData&&);
    BranchData(const RunPosition&, uint32_t predeterminedComponentId);

    BranchData(const BranchData&) = delete;
    BranchData(BranchData&& other);

    void consume(BranchData&& rhs);

    std::unique_ptr<RunTree> voxels_;
};

// Metadata for regular voxel components
struct RegularData {
    RegularData(RegularData&&, RegularData&&);
    RegularData(const RunPosition&, uint32_t predeterminedComponentId);

    RegularData(const RegularData&) = delete;
    RegularData(RegularData&& other);

    void consume(RegularData&& rhs);

    std::unique_ptr<RunTree> voxels_; //Invariant: voxels_ is ordered so that voxels_[i] and voxels_[i+1] are neighbros in the volume!
    tgt::svec3 leftEnd_;
    tgt::svec3 rightEnd_;
    uint32_t predeterminedComponentId_;
};

/// ----------------------------------------------------------------------------
/// Helper classes for fast KDTree based finding of skeleton voxels ------------
/// ----------------------------------------------------------------------------
struct VesselSkeletonVoxelRef {
    VesselSkeletonVoxelRef(tgt::vec3 rwPos, VesselSkeletonVoxel* voxel, uint32_t predeterminedComponentId)
        : rwPos(rwPos)
        , voxel(voxel)
        , predeterminedComponentId_(predeterminedComponentId)
    {
    }

    typedef tgt::vec3 VoxelType;

    tgt::vec3 rwPos; // Cached for better performance
    VesselSkeletonVoxel* voxel; //May be null => not part of skeleton
    uint32_t predeterminedComponentId_; // Id assigned to skeleton voxels from other knowledge

    inline typename VoxelType::ElemType at(int dim) const {
        return rwPos[dim];
    }
};

struct VoxelRefResultSetBase {
    typedef size_t IndexType;
    typedef float DistanceType;

    std::vector<VesselSkeletonVoxelRef>* storage_;
    std::vector<size_t> closestIDs_;
    float currentDist_;

    VoxelRefResultSetBase()
        : storage_(nullptr)
        , closestIDs_()
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

    void setStorageRef(std::vector<VesselSkeletonVoxelRef>& storage) {
        storage_ = &storage;
    }

    inline DistanceType worstDist() const {
        return currentDist_*1.001;
    }

    std::vector<const VesselSkeletonVoxelRef*> getVoxels() const {
        tgtAssert(storage_, "Storage ref not initialized");
        std::vector<const VesselSkeletonVoxelRef*> result;
        for(auto id : closestIDs_) {
            result.push_back(&getSkeletonVoxelRef(id));
        }
        return result;
    }

    std::vector<VesselSkeletonVoxelRef*> getVoxels() {
        tgtAssert(storage_, "Storage ref not initialized");
        std::vector<VesselSkeletonVoxelRef*> result;
        for(auto id : closestIDs_) {
            result.push_back(&getSkeletonVoxelRef(id));
        }
        return result;
    }

    const VesselSkeletonVoxelRef& getSkeletonVoxelRef(size_t id) const {
        tgtAssert(storage_, "Storage ref not initialized");
        return storage_->at(id);
    }

    VesselSkeletonVoxelRef& getSkeletonVoxelRef(size_t id) {
        tgtAssert(storage_, "Storage ref not initialized");
        return storage_->at(id);
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

struct ClosestVoxelsAdaptor : public VoxelRefResultSetBase {
    ClosestVoxelsAdaptor()
        : VoxelRefResultSetBase()
    {
    }

    /**
     * Called during search to add an element matching the criteria.
     * @return true if the search should be continued, false if the results are sufficient
     */
    inline bool addPoint(float dist, size_t index)
    {
        const VesselSkeletonVoxelRef& voxel = getSkeletonVoxelRef(index);
        if(dist <= currentDist_) {
            if(dist < currentDist_) {
                closestIDs_.clear();
            }
            currentDist_ = dist;
            if(voxel.voxel) {
                closestIDs_.push_back(index);
            }
        }
        return true;
    }
};

struct ClosestVoxelsAdaptorWithPredeterminedID : public VoxelRefResultSetBase {
    size_t wantedID_;

    ClosestVoxelsAdaptorWithPredeterminedID(size_t wantedID)
        : VoxelRefResultSetBase()
        , wantedID_(wantedID)
    {
    }

    /**
     * Called during search to add an element matching the criteria.
     * @return true if the search should be continued, false if the results are sufficient
     */
    inline bool addPoint(float dist, size_t index)
    {
        const VesselSkeletonVoxelRef& voxel = getSkeletonVoxelRef(index);
        if(voxel.predeterminedComponentId_ != wantedID_) {
            return true;
        }
        if(dist <= currentDist_) {
            if(dist < currentDist_) {
                closestIDs_.clear();
            }
            currentDist_ = dist;
            if(voxel.voxel) {
                closestIDs_.push_back(index);
            }
        }
        return true;
    }
};

template<typename E>
struct KDTreeVoxelFinder {
    typedef KDTreeBuilder<E> Builder;
    typedef nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<typename Builder::coord_t, Builder>,
		Builder,
		3 /* dim */
		> Index;

    KDTreeVoxelFinder(Builder&& builder)
        : storage_(std::move(builder))
        , index_(3 /*dim */, storage_, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* recommended as a sensible value by library creator */))
    {
        index_.buildIndex();
    }

    /**
     * Find the closest skeleton voxel(s)
     */
    template<class ResultSet>
    void findClosest(tgt::vec3 rwPos, ResultSet& set) {
        nanoflann::SearchParams params;
        params.sorted = false;

        index_.radiusSearchCustomCallback(rwPos.elem, set, params);
    }

    Builder storage_;
    Index index_;
};

/// A reader for the initial graph extraction from a generic input segmentation
/// It will be binarized on-the-fly
struct InitialSegmentationSliceReader : public CachingSliceReader {
    InitialSegmentationSliceReader(const VolumeBase& volume, float normalizedThreshold, int extent = 1)
        : CachingSliceReader(std::unique_ptr<VolumeSliceReader>(new VolumeSliceReader(volume)), extent)
        , normalizedThreshold_(normalizedThreshold)
        , volume_(volume)
    {
        seek(-1);
    }
    bool isObject(const tgt::ivec3& pos) const {
        return getVoxelNormalized(pos) >= normalizedThreshold_;
    }
    bool mayBeBranchObject(const tgt::ivec3& pos) const {
        //return isObject(pos);
        return true;
    }
    ClosestVoxelsAdaptor createVoxelFindingResultSet(const tgt::ivec3& pos) {
        return ClosestVoxelsAdaptor();
    }
    bool shouldFindSkeletonVoxelsToAnySurfaceVoxel() const {
        return true;
    }

    tgt::vec3 getSpacing() const {
        return volume_.getSpacing();
    }

    tgt::mat4 getVoxelToWorldMatrix() const {
        return volume_.getVoxelToWorldMatrix();
    }

    float normalizedThreshold_;
    const VolumeBase& volume_;
};

struct CCASegmentationSliceReader {
    CCASegmentationSliceReader(const VolumeBase& segmentation, const VolumeBase& branchIdVolume, float segmentationBinarizationThreshold)
        : segmentationReader_(tgt::make_unique<VolumeSliceReader>(segmentation), 1)
        , branchIdReader_(tgt::make_unique<VolumeSliceReader>(branchIdVolume), 1)
        , segmentationBinarizationThreshold_(segmentationBinarizationThreshold)
        , branchIdVolume_(branchIdVolume)
    {
        tgtAssert(branchIdVolume.getFormat() == "uint32", "Invalid volume format");
        tgtAssert(branchIdVolume.getDimensions() == branchIdVolume.getDimensions(), "Dimension mismatch");
        segmentationReader_.seek(-1);
        branchIdReader_.seek(-1);
    }
    uint32_t getId(const tgt::ivec3& pos) const {
        tgtAssert(pos.z >= 0 && pos.z < getDimensions().z, "Invalid z pos");
        int dz = pos.z - branchIdReader_.getCurrentZPos();
        auto slice = dynamic_cast<const VolumeRAM_UInt32*>(branchIdReader_.getSlice(dz));
        tgtAssert(slice, "Invalid slice type");
        return slice->voxel(pos.x, pos.y, 0);
    }
    void advance() {
        segmentationReader_.advance();
        branchIdReader_.advance();
    }
    const tgt::svec3& getDimensions() const {
        return segmentationReader_.getDimensions();
    }

    bool isObject(const tgt::ivec3& pos) const {
        return segmentationReader_.getVoxelNormalized(pos) >= segmentationBinarizationThreshold_;
    }
    bool mayBeBranchObject(const tgt::ivec3& pos) const {
        return getId(pos) != 0;
    }
    ClosestVoxelsAdaptorWithPredeterminedID createVoxelFindingResultSet(const tgt::ivec3& pos) {
        return ClosestVoxelsAdaptorWithPredeterminedID(getId(pos));
    }
    bool shouldFindSkeletonVoxelsToAnySurfaceVoxel() const {
        return false;
    }

    tgt::vec3 getSpacing() const {
        return branchIdVolume_.getSpacing();
    }

    tgt::mat4 getVoxelToWorldMatrix() const {
        return branchIdVolume_.getVoxelToWorldMatrix();
    }

    CachingSliceReader segmentationReader_;
    CachingSliceReader branchIdReader_;
    const VolumeBase& branchIdVolume_;
    float segmentationBinarizationThreshold_;
};

/// MetaDataCollector -------------------------------------------------------------------------------------------------------------
struct RegularSequence {
    RegularSequence(std::vector<tgt::svec3> voxels, uint32_t predeterminedComponentId);

    std::vector<tgt::svec3> voxels_;
    uint32_t predeterminedComponentId_;
};
// Helper that is used to collect finalized (end, regular, or branch) voxel components and creates the feature-annotated VesselGraph
struct MetaDataCollector {
    // Collect metadata from finalized voxel components
    void collect(EndData&&);
    void collect(RegularData&&);
    void collect(BranchData&&);

    // Create the vessel graph using the collected metadata and the segmentedVolume.
    // skeleton is currently only used for volume metadata
    // sampleMask may be null
    template<class S>
    std::unique_ptr<VesselGraph> createVesselGraph(S& segmentedVolumeReader, const VolumeBase* sampleMask, ProgressReporter& progress);

    std::unique_ptr<ProtoVesselGraph> createProtoVesselGraph(tgt::svec3 dimensions, const tgt::mat4& toRwMatrix, const boost::optional<LZ4SliceVolume<uint8_t>>& sampleMask, ProgressReporter& progress);

    std::vector<tgt::svec3> endPoints_;
    std::vector<RegularSequence> regularSequences_;
    std::vector<std::vector<tgt::svec3>> branchPoints_;

    static const std::string loggerCat_;
};


/// Streaming component extraction classes -------------------------------------
//
//Forward declaration:
template<typename MetaData>
class RunComposition;

// Defines the interface for a node (i.e., either a leaf (Run) or an internal node (RunComposition)) in the connected component tree
template<typename MetaData>
class Node {
public:
    Node(MetaDataCollector&);
    virtual ~Node();
    virtual MetaData getMetaData() = 0;
    // Add a node to this component
    virtual void addNode(Node* newRoot) = 0;
    // A score for how suitable this node is to be the root of a component tree
    // (used in merges)
    virtual uint32_t getRootAptitude() const = 0;

    // Get the root node of this component.
    // Note: this performs path compression!
    Node<MetaData>* getRootNode();
    // Set the parent of this tree
    void setParent(RunComposition<MetaData>* parent);
    // Has to be called before this node is deleted (used for reference counting)
    void finalize();

protected:
    RunComposition<MetaData>* parent_;
    MetaDataCollector& metaDataCollector_;
};

// An internal node in the component tree
template<typename MetaData>
class RunComposition final : public Node<MetaData> {
public:
    RunComposition(MetaDataCollector&, Node<MetaData>* r1, Node<MetaData>* r2);
    ~RunComposition();
    void addNode(Node<MetaData>* newRoot);
    MetaData getMetaData();
    RunComposition<MetaData>* getRoot();

    // Signal an additional reference to this node
    void ref();
    // Signal the removal of one reference to this node.
    // Node: This function may delete this object!
    void unref();
    uint32_t getRootAptitude() const;
private:
    MetaData metaData_;
    uint32_t refCount_;
};

const uint32_t UNDEFINED_COMPONENT_ID = 0;

// A leaf in the component tree
template<typename MetaData>
class Run final : public Node<MetaData> {
public:
    Run(MetaDataCollector&, tgt::svec2 yzPos_, size_t lowerBound_, size_t upperBound_, uint32_t  predeterminedComponentId = UNDEFINED_COMPONENT_ID);
    ~Run();

    // Try to merge with another run.
    void tryMerge(Run& other);

    MetaData getMetaData();
    void addNode(Node<MetaData>* newRoot);
    uint32_t getRootAptitude() const;

    const RunPosition pos_;
    const uint32_t predeterminedComponentId_;
};

// A row within a binary volume, it stores runs of end, regular, and branch voxels.
class Row {
public:
    // For performance reasons: (Re-)initialize this row from the given slice and row number
    // S is a slice reader that classifies voxels
    template<class S>
    void init(const S& slice, MetaDataCollector& mdc, size_t sliceNum, size_t rowNum);
    Row();
    ~Row() {}

    // Connect runs within this rows with runs in other rows
    void connect(Row& other);
    // Release all runs
    void finalizeRuns();
private:
    std::vector<Run<EndData>> endRuns_; ///< Sorted!
    std::vector<Run<RegularData>> regularRuns_; ///< Sorted!
    std::vector<Run<BranchData>> branchRuns_; ///< Sorted!
};

// Helper class for storage of active rows.
// Only the latest <num rows in slice>+2 rows are kept in memory.
class RowStorage {
public:
    RowStorage(const tgt::svec3& volumeDimensions);
    ~RowStorage();
    // Add a row at the given position
    template<class S>
    void add(const S& slice, MetaDataCollector& mdc, size_t sliceNum, size_t row);
    // Get the newest row
    Row& latest() const;

    // Get a row relative to the newest
    template<int DY, int DZ>
    Row& latest_DP() const;

    // Get a row relative to the newest with the newest row
    template<int DY, int DZ>
    void connectLatestWith();

    // Finalize all active rows and release components
    void finalizeRows();

    // Only for debug purposes
    Row* getRows() const;
private:
    // Get a row via storage position
    Row& get(size_t pos) const;

    const size_t storageSize_;
    const size_t rowsPerSlice_;
    Row* rows_;
    size_t storagePos_;
};

// Perform the connected component analysis and collect metadata for end, regular, and branch voxel components
// S is a slice reader that allows for voxel classification
template<class S>
std::unique_ptr<MetaDataCollector> cca(S&& classSliceReader, ProgressReporter& progress);

/// -----------------------------------------------------------------------------------------------------------------------
/// Implementation of template classes and methods ------------------------------------------------------------------------
/// -----------------------------------------------------------------------------------------------------------------------

/// Connected component related classes ----------------------------------------
static float surfaceDistanceSq(const tgt::vec3& skeletonVoxel, const tgt::vec3 surfaceVoxel, const tgt::vec3 spacing) {
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

static bool isAtSampleBorder(const tgt::svec3& p, const tgt::svec3& volumedim) {
    return p.x == 0             || p.y == 0             || p.z == 0
        || p.x >= volumedim.x-1 || p.y >= volumedim.y-1 || p.z >= volumedim.z-1;
}

struct VoxelKDElement {
    typedef tgt::ivec3 VoxelType;

    VoxelKDElement(tgt::ivec3 pos, uint64_t nodeID)
        : pos_(pos)
        , nodeID_(nodeID)
    {
    }

    inline typename VoxelType::ElemType at(int dim) const {
        return pos_[dim];
    }

    VoxelType pos_;
    uint64_t nodeID_;
};

template<typename E>
class KDTreeNeighborhoodFinder {
public:
    typedef nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<typename KDTreeBuilder<E>::coord_t, KDTreeBuilder<E>>,
		KDTreeBuilder<E>,
		3 /* dim */
		> Index;

    KDTreeNeighborhoodFinder(KDTreeBuilder<E>&& builder)
        : adaptor_(std::move(builder))
        , index_(3 /*dim */, adaptor_, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* recommended as a sensible value by library creator */))
    {
        index_.buildIndex();
    }

    std::vector<const E*> find_26_neighbors(const typename E::VoxelType& pos) {
        auto search_radius = static_cast<typename KDTreeBuilder<E>::coord_t>(4); // (L2 distance)^2 < 4 <=> 26_neighbors (or equal)

        nanoflann::SearchParams params;
        params.sorted = false;

        std::vector<std::pair<size_t, typename KDTreeBuilder<E>::coord_t> > index_results;
        index_.radiusSearch(pos.elem, search_radius, index_results, params);

        std::vector<const E*> results;
        for(const auto& index_result : index_results) {
            const E& p = adaptor_.points()[index_result.first];
            if(p.pos_ != pos) {
                results.push_back(&p);
            }
        }
        return results;
    }

private:
    KDTreeBuilder<E> adaptor_;
    Index index_;
};

#ifdef _MSC_VER
#define RELEASE_ASSERT(assertion, errormsg) \
    if(!(assertion)) { \
        LERROR((errormsg)); \
        __debugbreak(); \
    }
#else
#define RELEASE_ASSERT(assertion, errormsg) \
    if(!(assertion)) { \
        LERROR((errormsg)); \
        asm("int $3"); \
    }
#endif

/// MetaDataCollector ----------------------------------------------------------

template<class S>
std::unique_ptr<VesselGraph> MetaDataCollector::createVesselGraph(S& segmentedVolumeReader, const VolumeBase* sampleMask, ProgressReporter& progress) {
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

    KDTreeBuilder<VoxelKDElement> nodeVoxelTreeBuilder;
    KDTreeBuilder<VesselSkeletonVoxelRef> finderBuilder;
    // Create nodes

    for(const auto& p : endPoints_) {
        std::vector<tgt::vec3> voxels;
        tgt::vec3 rwpos = transform(toRWMatrix, p);
        voxels.push_back(rwpos);
        size_t id = graph->insertNode(rwpos, std::move(voxels), isAtSampleBorder(p, dimensions));
        finderBuilder.push(VesselSkeletonVoxelRef(rwpos, nullptr, UNDEFINED_COMPONENT_ID));
        nodeVoxelTreeBuilder.push(VoxelKDElement(tgt::ivec3(p),id));
    }
    for(const auto& branchPoint : branchPoints_) {
        size_t numVoxels = 0;
        bool junctionAtSampleBorder = false;
        tgt::vec3 sum = tgt::vec3::zero;
        std::vector<tgt::vec3> rwBrachVoxels;
        for(const auto& p : branchPoint) {
            tgt::vec3 rwpos = transform(toRWMatrix, p);
            sum += rwpos;
            rwBrachVoxels.push_back(rwpos);
            ++numVoxels;
            finderBuilder.push(VesselSkeletonVoxelRef(rwpos, nullptr, UNDEFINED_COMPONENT_ID));
            junctionAtSampleBorder |= isAtSampleBorder(p, dimensions);
        }
        size_t id = graph->insertNode(sum/static_cast<float>(numVoxels), std::move(rwBrachVoxels), junctionAtSampleBorder);
        for(const auto& p : branchPoint) {
            nodeVoxelTreeBuilder.push(VoxelKDElement(tgt::ivec3(p),id));
        }
    }

    // Precreate VoxelSkeletonLists
    std::vector<std::vector<VesselSkeletonVoxel>> skeletonVoxelsLists;
    skeletonVoxelsLists.reserve(regularSequences_.size());
    for(const auto& regularSequence : regularSequences_) {
        skeletonVoxelsLists.emplace_back();
        auto& skeletonVoxelsList = skeletonVoxelsLists.back();
        skeletonVoxelsList.reserve(regularSequence.voxels_.size());
        for(const tgt::svec3& voxel : regularSequence.voxels_) {
            tgt::vec3 rwpos = transform(toRWMatrix, voxel);
            skeletonVoxelsList.emplace_back(rwpos, std::numeric_limits<float>::infinity(), 0, 0, 0, 0);
            finderBuilder.push(VesselSkeletonVoxelRef(rwpos, &skeletonVoxelsList.back(), regularSequence.predeterminedComponentId_));
        }
    }

    const float criticalVoxelDistDiff = 1.001f*tgt::length(spacing);
    KDTreeVoxelFinder<VesselSkeletonVoxelRef> kdfinder(std::move(finderBuilder));
    KDTreeNeighborhoodFinder<VoxelKDElement> nodeVoxelTree(std::move(nodeVoxelTreeBuilder));

    for(int z = 0; z < dimensions.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dimensions.z);
        segmentedVolumeReader.advance();

        for(int y = 0; y < dimensions.y; ++y) {
            for(int x = 0; x < dimensions.x; ++x) {

                tgt::ivec3 ipos(x, y, z);

                // Only consider object points
                if(!segmentedVolumeReader.isObject(ipos)) {
                    continue;
                }

                // Ignore values outside of sample (may happen due to inaccuracies in mask generation)
                if(sampleMaskRAM && sampleMaskRAM->getVoxelNormalized(ipos) == 0) {
                    continue;
                }

                bool mayBePartOfBranch = segmentedVolumeReader.mayBeBranchObject(ipos);

                if(mayBePartOfBranch) {
                    tgt::vec3 rwVoxel = transform(toRWMatrix, tgt::vec3(x, y, z));

                    auto resultSet = segmentedVolumeReader.createVoxelFindingResultSet(ipos);
                    resultSet.setStorageRef(kdfinder.storage_.points());
                    kdfinder.findClosest(rwVoxel, resultSet);

                    if(segmentedVolumeReader.shouldFindSkeletonVoxelsToAnySurfaceVoxel() && !resultSet.found()) {
                        RELEASE_ASSERT(false, "Did not find closest voxel");
                    }
                    std::vector<std::reference_wrapper<VesselSkeletonVoxel>> minDistVoxels;
                    for(const auto& ref: resultSet.getVoxels()) {
                        if(ref->voxel) {
                            minDistVoxels.push_back(std::reference_wrapper<VesselSkeletonVoxel>(*ref->voxel));
                        }
                    }

                    if(minDistVoxels.empty()) {
                        // No skeleton voxels, only branch voxels
                        continue;
                    }

                    float volume = tgt::hmul(spacing) / minDistVoxels.size();
                    for(VesselSkeletonVoxel& voxel : minDistVoxels) {
                        voxel.volume_ += volume;
                    }

                    // Determine if it is a surface voxel:
                    // Only consider those that have a background 6-neighbor
                    bool isSurfaceVoxel = !( (x == 0              || segmentedVolumeReader.isObject(tgt::ivec3(x-1, y  , z  )))
                            && (x == dimensions.x-1 || segmentedVolumeReader.isObject(tgt::ivec3(x+1, y  , z  )))
                            && (y == 0              || segmentedVolumeReader.isObject(tgt::ivec3(x  , y-1, z  )))
                            && (y == dimensions.y-1 || segmentedVolumeReader.isObject(tgt::ivec3(x  , y+1, z  )))
                            && (z == 0              || segmentedVolumeReader.isObject(tgt::ivec3(x  , y  , z-1)))
                            && (z == dimensions.z-1 || segmentedVolumeReader.isObject(tgt::ivec3(x  , y  , z+1))));

                    // The rest is only relevant for surface voxels, so...
                    if(!isSurfaceVoxel) {
                        continue;
                    }

                    float dist = std::sqrt(surfaceDistanceSq(minDistVoxels.front().get().pos_, rwVoxel, spacing));

                    for(VesselSkeletonVoxel& voxel : minDistVoxels) {
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
    }
    progress.setProgress(1.0f);

    // Create edges from branches
    for(size_t i = 0; i < skeletonVoxelsLists.size(); ++i) {
        auto& regularSequence = regularSequences_[i].voxels_;
        auto& skeletonVoxelsList = skeletonVoxelsLists[i];

        tgtAssert(!regularSequence.empty(), "Empty sequence");
        tgt::ivec3 leftEnd = regularSequence.front();
        tgt::ivec3 rightEnd = regularSequence.back();
        const uint64_t NO_NODE_FOUND = -1;
        uint64_t leftEndNode = NO_NODE_FOUND;
        uint64_t rightEndNode = NO_NODE_FOUND;

        if(leftEnd == rightEnd) {
            // One voxel long branch
            auto neighbors = nodeVoxelTree.find_26_neighbors(leftEnd);
            RELEASE_ASSERT(neighbors.size() == 2, "Invalid number of neighbors");
            leftEndNode = neighbors.at(0)->nodeID_;
            rightEndNode = neighbors.at(1)->nodeID_;
        } else {
            auto leftNeighbors = nodeVoxelTree.find_26_neighbors(leftEnd);
            auto rightNeighbors = nodeVoxelTree.find_26_neighbors(rightEnd);
            if(leftNeighbors.empty() || rightNeighbors.empty()) {
                RELEASE_ASSERT(are26Neighbors(leftEnd, rightEnd), "non-connected, non-loop sequence");
                RELEASE_ASSERT(leftNeighbors.empty() && rightNeighbors.empty(), "left xor right neighbors are empty!");

                // We found a freestanding loop and have to add an extra node to support it
                std::vector<tgt::vec3> voxels;
                tgt::vec3 leftTransformed = transform(toRWMatrix, leftEnd);
                tgt::vec3 rightTransformed = transform(toRWMatrix, rightEnd);
                voxels.push_back(leftTransformed);
                voxels.push_back(rightTransformed);
                tgt::vec3 center = (leftTransformed + rightTransformed)*0.5f;
                size_t newNode = graph->insertNode(center, std::move(voxels), isAtSampleBorder(tgt::svec3(leftEnd), dimensions) || isAtSampleBorder(tgt::svec3(rightEnd), dimensions));
                leftEndNode = newNode;
                rightEndNode = newNode;
            } else {
                // "Regular" branch
                RELEASE_ASSERT(leftNeighbors.size() == 1, "Invalid number of left end neighbors");
                leftEndNode = leftNeighbors.at(0)->nodeID_;

                RELEASE_ASSERT(rightNeighbors.size() == 1, "Invalid number of right end neighbors");
                rightEndNode = rightNeighbors.at(0)->nodeID_;
            }
        }

        graph->insertEdge(leftEndNode, rightEndNode, std::move(skeletonVoxelsList));
    }

    // Create edges from end voxels directly connected to branchPoints
    size_t endNodeId = 0;
    for(const auto& e : endPoints_) {
        tgt::ivec3 endPointI(e);
        for(const auto& p : nodeVoxelTree.find_26_neighbors(endPointI)) {
            graph->insertEdge(p->nodeID_, endNodeId, std::vector<VesselSkeletonVoxel>());
            break;
        }
        ++endNodeId;
    }

    return graph;
}

/// Node -----------------------------------------------------------------------------------------------------------------------

template<typename MetaData>
Node<MetaData>::Node(MetaDataCollector& mdc)
    : parent_(nullptr)
    , metaDataCollector_(mdc)
{
}
template<typename MetaData>
Node<MetaData>::~Node() {
}

template<typename MetaData>
Node<MetaData>* Node<MetaData>::getRootNode() {
    if(!parent_) {
        return this;
    }
    setParent(parent_->getRoot());
    return parent_;
}

template<typename MetaData>
void Node<MetaData>::setParent(RunComposition<MetaData>* newParent) {
    tgtAssert(newParent, "newParent is null");
    tgtAssert(newParent != this, "newParent=this");
    // First ref the new parent, THEN unref the old in case they are the same
    newParent->ref();
    if(parent_) {
        parent_->unref();
    }
    parent_ = newParent;

}

template<typename MetaData>
void Node<MetaData>::finalize() {
    if(parent_) {
        parent_->unref();
    } else {
        this->metaDataCollector_.collect(std::move(getMetaData()));
    }
}


/// RunComposition -----------------------------------------------------------------------------------------------------------------------

template<typename MetaData>
uint32_t RunComposition<MetaData>::getRootAptitude() const {
    return refCount_;
}


template<typename MetaData>
void RunComposition<MetaData>::addNode(Node<MetaData>* other) {
    tgtAssert(other, "newroot is null");
    tgtAssert(!this->parent_, "Parent not null");
    metaData_.consume(other->getMetaData());
    other->setParent(this);
}


template<typename MetaData>
RunComposition<MetaData>::RunComposition(MetaDataCollector& mdc, Node<MetaData>* r1, Node<MetaData>* r2)
    : Node<MetaData>(mdc)
    , metaData_(r1->getMetaData(), r2->getMetaData())
    , refCount_(0)
{
    r1->setParent(this);
    r2->setParent(this);
}

template<typename MetaData>
RunComposition<MetaData>::~RunComposition() {
}


template<typename MetaData>
MetaData RunComposition<MetaData>::getMetaData() {
    return std::move(metaData_);
}

template<typename MetaData>
RunComposition<MetaData>* RunComposition<MetaData>::getRoot() {
    if(!this->parent_) {
        return this;
    }
    this->setParent(this->parent_->getRoot());
    return this->parent_;
}

template<typename MetaData>
void RunComposition<MetaData>::ref() {
    refCount_++;
}


template<typename MetaData>
void RunComposition<MetaData>::unref() {
    --refCount_;
    if(refCount_ == 0) {
        // Nobody likes me :(
        this->finalize();
        delete this;
    }
}

/// Run -----------------------------------------------------------------------------------------------------------------------


template<typename MetaData>
Run<MetaData>::Run(MetaDataCollector& mdc, tgt::svec2 yzPos, size_t lowerBound, size_t upperBound, uint32_t predeterminedComponentId)
    : Node<MetaData>(mdc)
    , pos_(yzPos, lowerBound, upperBound)
    , predeterminedComponentId_(predeterminedComponentId)
{
}


template<typename MetaData>
Run<MetaData>::~Run() {
}


template<typename MetaData>
void Run<MetaData>::addNode(Node<MetaData>* other) {
    tgtAssert(other, "newroot is null");
    tgtAssert(!this->parent_, "parent is not null");
    // Construct a new root
    RunComposition<MetaData>* newRoot = new RunComposition<MetaData>(this->metaDataCollector_, this, other);
}


template<typename MetaData>
void Run<MetaData>::tryMerge(Run& other) {
    if(pos_.xlow_ > other.pos_.xhigh_ || other.pos_.xlow_ > pos_.xhigh_) {
        return;
    }

    Node<MetaData>* thisRoot = this->getRootNode();
    Node<MetaData>* otherRoot = other.getRootNode();
    if(thisRoot == otherRoot) {
        return;
    }
    if(thisRoot->getRootAptitude() > otherRoot->getRootAptitude()) {
        thisRoot->addNode(otherRoot);
    } else {
        otherRoot->addNode(thisRoot);
    }
}


template<typename MetaData>
MetaData Run<MetaData>::getMetaData() {
    return std::move(MetaData(pos_, predeterminedComponentId_));
}

template<typename MetaData>
uint32_t Run<MetaData>::getRootAptitude() const {
    return 0;
}

/// Row ------------------------------------------------------------------------

template<class S>
void Row::init(const S& slice, MetaDataCollector& mdc, size_t sliceNum, size_t rowNum)
{
    // Finalize previous:
    finalizeRuns();
    endRuns_.clear();
    regularRuns_.clear();
    branchRuns_.clear();

    // Insert new runs
    size_t lastClass = 0;
    size_t runStart = 0;
    size_t rowLength = slice.getDimensions().x;
    tgt::svec2 yzPos(rowNum, sliceNum);
    for(size_t x = 0; x < rowLength; ++x) {
        size_t currentClass = slice.getClass(tgt::ivec2(x, rowNum));
        if(currentClass != lastClass) {
            switch(lastClass) {
                case 1:
                    endRuns_.emplace_back(mdc, yzPos, runStart, x);
                    break;
                case 2:
                    regularRuns_.emplace_back(mdc, yzPos, runStart, x, slice.getPredefinedComponentId(tgt::ivec2(runStart, rowNum)));
                    break;
                case 3:
                    branchRuns_.emplace_back(mdc, yzPos, runStart, x);
                    break;
            }
            runStart = x;
            lastClass = currentClass;
        }
    }
    switch(lastClass) {
        case 1:
            endRuns_.emplace_back(mdc, yzPos, runStart, rowLength);
            break;
        case 2:
            regularRuns_.emplace_back(mdc, yzPos, runStart, rowLength, slice.getPredefinedComponentId(tgt::ivec2(runStart, rowNum)));
            break;
        case 3:
            branchRuns_.emplace_back(mdc, yzPos, runStart, rowLength);
            break;
    }
}

/// RowStorage -----------------------------------------------------------------

template<class S>
void RowStorage::add(const S& slice, MetaDataCollector& mdc, size_t sliceNum, size_t rowNum) {
    storagePos_ = (storagePos_ + 1)%storageSize_;
    rows_[storagePos_].init(slice, mdc, sliceNum, rowNum);
}

template<int DY, int DZ>
Row& RowStorage::latest_DP() const {
    static_assert(-1 <= DY && DY <= 1, "Invalid DY");
    static_assert(-1 <= DZ && DZ <= 0, "Invalid DZ");
    return get(storagePos_ + storageSize_ + DY + rowsPerSlice_ * DZ);
}

template<int DY, int DZ>
void RowStorage::connectLatestWith() {
    latest().connect(latest_DP<DY,DZ>());
}

template<class S>
std::unique_ptr<MetaDataCollector> cca(S&& classSliceReader, ProgressReporter& progress) {
    TaskTimeLogger _("Extract graph components", tgt::Info);
    const tgt::svec3 dim = classSliceReader.getDimensions();

    std::unique_ptr<MetaDataCollector> mdc = std::unique_ptr<MetaDataCollector>(new MetaDataCollector());
    {
        RowStorage rows(dim);
        // First layer
        {
            classSliceReader.advance();
            rows.add(classSliceReader, *mdc, 0, 0);
            for(size_t y = 1; y<dim.y; ++y) {
                // Create new row at z=0
                rows.add(classSliceReader, *mdc, 0, y);

                // merge with row (-1, 0)
                rows.template connectLatestWith<-1, 0>();
            }
        }

        // The rest of the layers
        for(size_t z = 1; z<dim.z; ++z) {
            progress.setProgress(static_cast<float>(z)/dim.z);
            classSliceReader.advance();

            // Create new row at y=0
            rows.add(classSliceReader, *mdc, z, 0);

            // merge with row (0, -1)
            rows.template connectLatestWith< 0,-1>();

            // merge with row ( 1,-1)
            rows.template connectLatestWith< 1,-1>();

            for(size_t y = 1; y<dim.y-1; ++y) {
                // Create new row
                rows.add(classSliceReader, *mdc, z, y);

                // merge with row (-1, 0)
                rows.template connectLatestWith<-1, 0>();

                // merge with row (-1,-1)
                rows.template connectLatestWith<-1,-1>();

                // merge with row ( 0,-1)
                rows.template connectLatestWith< 0,-1>();

                // merge with row ( 1,-1)
                rows.template connectLatestWith< 1,-1>();

            }

            // Last row at y = dim.y-1
            rows.add(classSliceReader, *mdc, z, dim.y-1);

            // merge with row (-1, 0)
            rows.template connectLatestWith<-1, 0>();

            // merge with row (-1,-1)
            rows.template connectLatestWith<-1,-1>();

            // merge with row ( 0,-1)
            rows.template connectLatestWith< 0,-1>();
        }
        // Finalize all left over rows:
        rows.finalizeRows();
    } // RowStorage rows goes out of scope and destroys and finalizes left over runs
    // -> mdc has collected all metadata
    progress.setProgress(1.0f);
    return mdc;
}

} //namespace voreen

#endif //VRN_STREAMING_COMPONENTS_VESSELTOPOLOGY_H
