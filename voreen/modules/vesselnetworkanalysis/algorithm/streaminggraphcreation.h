/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_STREAMING_COMPONENTS_VESSELNETWORKANALYSIS_H
#define VRN_STREAMING_COMPONENTS_VESSELNETWORKANALYSIS_H

#include "voreen/core/io/progressreporter.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/diskarraystorage.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include "modules/bigdataimageprocessing/volumefiltering/slicereader.h"
#include "modules/bigdataimageprocessing/datastructures/lz4slicevolume.h"
#include "modules/hdf5/io/hdf5filevolume.h"

#include "../algorithm/volumemask.h"
#include "../datastructures/vesselgraph.h"
#include "../datastructures/protovesselgraph.h"
#include "../util/tasktimelogger.h"

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
    ~RunNode();
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
    EndData(const RunPosition&);
    void consume(EndData&& rhs);

    tgt::svec3 pos_;
};

// Metadata for branch voxel components
struct BranchData {
    BranchData(BranchData&&, BranchData&&);
    BranchData(const RunPosition&);

    BranchData(const BranchData&) = delete;
    BranchData(BranchData&& other);

    void consume(BranchData&& rhs);

    std::unique_ptr<RunTree> voxels_;
};

// Metadata for regular voxel components
struct RegularData {
    RegularData(RegularData&&, RegularData&&);
    RegularData(const RunPosition&);

    RegularData(const RegularData&) = delete;
    RegularData(RegularData&& other);

    void consume(RegularData&& rhs);

    std::unique_ptr<RunTree> voxels_; //Invariant: voxels_ is ordered so that voxels_[i] and voxels_[i+1] are neighbros in the volume!
    tgt::svec3 leftEnd_;
    tgt::svec3 rightEnd_;
};

/// MetaDataCollector -------------------------------------------------------------------------------------------------------------
// Helper that is used to collect finalized (end, regular, or branch) voxel components and creates the feature-annotated VesselGraph
struct MetaDataCollector {
    MetaDataCollector();

    // Collect metadata from finalized voxel components
    void collect(EndData&&);
    void collect(RegularData&&);
    void collect(BranchData&&);

    // Create the vessel graph using the collected metadata and the segmentedVolume.
    // skeleton is currently only used for volume metadata
    // sampleMask may be null
    std::unique_ptr<ProtoVesselGraph> createProtoVesselGraph(tgt::svec3 dimensions, const tgt::mat4& toRwMatrix, const boost::optional<LZ4SliceVolume<uint8_t>>& sampleMask, ProgressReporter& progress);


    std::vector<tgt::svec3> endPoints_;
    std::vector<std::vector<tgt::svec3>> branchPoints_;

    DiskArrayStorage<tgt::svec3> voxelStorage_;
    std::vector<DiskArray<tgt::svec3>> regularSequences_;

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
    Run(MetaDataCollector&, tgt::svec2 yzPos_, size_t lowerBound_, size_t upperBound_);
    ~Run();

    // Try to merge with another run.
    void tryMerge(Run& other);

    MetaData getMetaData();
    void addNode(Node<MetaData>* newRoot);
    uint32_t getRootAptitude() const;

    const RunPosition pos_;
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
Run<MetaData>::Run(MetaDataCollector& mdc, tgt::svec2 yzPos, size_t lowerBound, size_t upperBound)
    : Node<MetaData>(mdc)
    , pos_(yzPos, lowerBound, upperBound)
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
    return std::move(MetaData(pos_));
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
                    regularRuns_.emplace_back(mdc, yzPos, runStart, x);
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
            regularRuns_.emplace_back(mdc, yzPos, runStart, rowLength);
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

#endif //VRN_STREAMING_COMPONENTS_VESSELNETWORKANALYSIS_H
