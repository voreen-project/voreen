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

#include "streaminggraphcreation.h"
#include "tgt/vector.h"

namespace voreen {

// RunPosition ------------------------------------------------------------------------------

RunPosition::RunPosition(tgt::svec2 yz, size_t xlow, size_t xhigh)
    : y_(yz.x)
    , z_(yz.y)
    , xlow_(xlow)
    , xhigh_(xhigh)
{
}
tgt::svec3 RunPosition::begin() const {
    return tgt::svec3(xlow_, y_, z_);
}
tgt::svec3 RunPosition::end() const {
    return tgt::svec3(xhigh_-1, y_, z_);
}


/// RunTrees -----------------------------------------------------------------------------
void RunTree::invert() {
    inverted_ = ! inverted_;
}

// RunNode -------------------------------------------------------------------------------
RunNode::RunNode(std::unique_ptr<RunTree>&& left, std::unique_ptr<RunTree>&& right)
    : left_(std::move(left))
    , right_(std::move(right))
{
}

void RunNode::collectVoxels(std::vector<tgt::svec3>& vec, bool inverted) const {
    tgtAssert(left_, "left child is null")
    tgtAssert(left_, "right child is null")

    bool finalInverted = inverted^inverted_;
    if(finalInverted) {
        right_->collectVoxels(vec, finalInverted);
        left_->collectVoxels(vec, finalInverted);
    } else {
        left_->collectVoxels(vec, finalInverted);
        right_->collectVoxels(vec, finalInverted);
    }
}
// RunLeaf -------------------------------------------------------------------------------
RunLeaf::RunLeaf(const RunPosition& run)
    : run_(run)
{
}
void RunLeaf::collectVoxels(std::vector<tgt::svec3>& vec, bool inverted) const {
    if(inverted_ ^ inverted) {
        for(int x = run_.xhigh_ - 1; x >= static_cast<int>(run_.xlow_); --x) {
            vec.emplace_back(x, run_.y_, run_.z_);
        }
    } else {
        for(size_t x = run_.xlow_; x < run_.xhigh_; ++x) {
            vec.emplace_back(x, run_.y_, run_.z_);
        }
    }
}

/// SkeletonClassReader --------------------------------------------------------
SkeletonClassReader::SkeletonClassReader(const VolumeMask& skeleton)
    : skeleton_(skeleton)
    , z_(-1)
{
}

bool SkeletonClassReader::isObject(const tgt::ivec3& xyz) const {
    return skeleton_.get(xyz, VolumeMaskValue::BACKGROUND) != VolumeMaskValue::BACKGROUND;
}

uint8_t SkeletonClassReader::getClass(const tgt::ivec2& xypos) const {
    uint8_t numNeighbors = 0;
    const int xmin = xypos.x-1;
    const int xmax = xypos.x+2;
    const int ymin = xypos.y-1;
    const int ymax = xypos.y+2;
    const int zmin = z_-1;
    const int zmax = z_+2;
    if(!isObject(tgt::ivec3(xypos, z_))) {
        return 0;
    }
    for(int z=zmin; z<zmax; ++z) {
        for(int y=ymin; y<ymax; ++y) {
            for(int x=xmin; x<xmax; ++x) {
                if(isObject(tgt::ivec3(x,y,z))) {
                    ++numNeighbors;
                }
            }
        }
    }
    // 0 neighbors => endvoxel = class 1
    return std::max<uint8_t>(1, std::min<uint8_t>(numNeighbors-1 /* minus itself */, 3));
}
void SkeletonClassReader::advance() {
    ++z_;
}
const tgt::svec3& SkeletonClassReader::getDimensions() const {
    return skeleton_.getDimensions();
}

/// NeighborCountVoxelClassifier -----------------------------------------------------------
uint32_t NeighborCountVoxelClassifier::getPredefinedComponentId(const tgt::ivec2& xyz) const {
    return UNDEFINED_COMPONENT_ID;
}

NeighborCountVoxelClassifier::NeighborCountVoxelClassifier(const VolumeMask& skeleton)
    : SkeletonClassReader(skeleton)
{
}

NeighborCountVoxelClassifier::~NeighborCountVoxelClassifier() {
}

/// NeighborCountVoxelClassifier -----------------------------------------------------------

NeighborCountAndBranchSegmentationVoxelClassifier::NeighborCountAndBranchSegmentationVoxelClassifier(const VolumeMask& skeleton, const VolumeBase& branchSegmentation)
    : skeletonClassReader_(skeleton)
    , branchSegmentationReader_(branchSegmentation)
{
    tgtAssert(branchSegmentation.getFormat() == "uint32", "Invalid volume format");
    branchSegmentationReader_.seek(-1);
}

NeighborCountAndBranchSegmentationVoxelClassifier::~NeighborCountAndBranchSegmentationVoxelClassifier() {
}

uint8_t NeighborCountAndBranchSegmentationVoxelClassifier::getClass(const tgt::ivec2& xypos) const {
    uint8_t neighborCountClass = skeletonClassReader_.getClass(xypos);
    return neighborCountClass;
    /*
    if(neighborCountClass == 2 && getPredefinedComponentId(xypos) == 0) {
        return 3; //Forcibly classify regular voxels in branch balls as branch voxels
    } else {
        return neighborCountClass;
    }
    */
}
uint32_t NeighborCountAndBranchSegmentationVoxelClassifier::getPredefinedComponentId(const tgt::ivec2& xypos) const {
    auto slice = dynamic_cast<const VolumeRAM_UInt32*>(branchSegmentationReader_.getCurrentSlice());
    tgtAssert(slice, "No slice");
    return slice->voxel(xypos.x, xypos.y, 0);
}

void NeighborCountAndBranchSegmentationVoxelClassifier::advance() {
    skeletonClassReader_.advance();
    branchSegmentationReader_.advance();
}
const tgt::svec3& NeighborCountAndBranchSegmentationVoxelClassifier::getDimensions() const {
    return skeletonClassReader_.getDimensions();
}

/// MetaData -----------------------------------------------------------------------------
/// EndData ------------------------------------------------------------------------------
EndData::EndData(EndData&& l, EndData&&)
    : pos_(l.pos_)
{
    // This may happen for two mutually connected endvoxels
    // In this case we will just keep one of them.
    //tgtAssert(false, "Connected endvoxels: Composition");
}
EndData::EndData(const RunPosition& p, uint32_t)
    : pos_(p.xlow_, p.y_, p.z_)
{
    // This may happen for two mutually connected endvoxels
    // In this case we will just keep one of them.
    //tgtAssert(p.xlow_ + 1 == p.xhigh_, "Connected endvoxels: Run");
}
void EndData::consume(EndData&&) {
    tgtAssert(false, "Connected endvoxels: Consume");
}
/// BranchData ---------------------------------------------------------------------------
BranchData::BranchData(BranchData&& other)
    : voxels_(other.voxels_.release())
{
}
BranchData::BranchData(BranchData&& b1, BranchData&& b2)
    : voxels_(new RunNode(std::move(b1.voxels_), std::move(b2.voxels_)))
{
}
BranchData::BranchData(const RunPosition& rp, uint32_t)
    : voxels_(new RunLeaf(rp))
{
}
void BranchData::consume(BranchData&& rhs) {
    voxels_ = std::unique_ptr<RunTree>(new RunNode(std::move(voxels_), std::move(rhs.voxels_)));
}
/// RegularData --------------------------------------------------------------------------
RegularData::RegularData(RegularData&& other)
    : voxels_(other.voxels_.release())
    , leftEnd_(other.leftEnd_)
    , rightEnd_(other.rightEnd_)
    , predeterminedComponentId_(other.predeterminedComponentId_)
{
}


RegularData::RegularData(RegularData&& b1, RegularData&& b2)
    : voxels_(nullptr)
    , predeterminedComponentId_(b1.predeterminedComponentId_)
{
    tgtAssert(
               b1.predeterminedComponentId_ == UNDEFINED_COMPONENT_ID
            || b2.predeterminedComponentId_ == UNDEFINED_COMPONENT_ID
            || b1.predeterminedComponentId_ == b2.predeterminedComponentId_, "predeterminedComponentId mismatch");
    if(b2.predeterminedComponentId_ == UNDEFINED_COMPONENT_ID) {
        predeterminedComponentId_ = b1.predeterminedComponentId_;
    } else {
        predeterminedComponentId_ = b2.predeterminedComponentId_;
    }

    if(are26Neighbors(b1.leftEnd_, b2.leftEnd_)) {
        b1.voxels_->invert();
        leftEnd_ = b1.rightEnd_;
        rightEnd_ = b2.rightEnd_;
    } else if(are26Neighbors(b1.leftEnd_, b2.rightEnd_)) {
        b1.voxels_->invert();
        b2.voxels_->invert();
        leftEnd_ = b1.rightEnd_;
        rightEnd_ = b2.leftEnd_;
    } else if(are26Neighbors(b1.rightEnd_, b2.leftEnd_)) {
        leftEnd_ = b1.leftEnd_;
        rightEnd_ = b2.rightEnd_;
    } else if(are26Neighbors(b1.rightEnd_, b2.rightEnd_)) {
        b2.voxels_->invert();
        leftEnd_ = b1.leftEnd_;
        rightEnd_ = b2.leftEnd_;
    } else {
        tgtAssert(false, "invalid connection");
    }
    voxels_ = std::unique_ptr<RunTree>(new RunNode(std::move(b1.voxels_), std::move(b2.voxels_)));
}


RegularData::RegularData(const RunPosition& rp, uint32_t predeterminedComponentId)
    : voxels_(new RunLeaf(rp))
    , leftEnd_(rp.begin())
    , rightEnd_(rp.end())
    , predeterminedComponentId_(predeterminedComponentId)
{
}

void RegularData::consume(RegularData&& rhs) {
    tgtAssert(
               predeterminedComponentId_ == UNDEFINED_COMPONENT_ID
            || rhs.predeterminedComponentId_ == UNDEFINED_COMPONENT_ID
            || predeterminedComponentId_ == rhs.predeterminedComponentId_, "predeterminedComponentId mismatch");
    if(predeterminedComponentId_ == UNDEFINED_COMPONENT_ID) {
        predeterminedComponentId_ = rhs.predeterminedComponentId_;
    }
    if(are26Neighbors(leftEnd_, rhs.leftEnd_)) {
        voxels_->invert();
        leftEnd_ = rightEnd_;
        rightEnd_ = rhs.rightEnd_;
    } else if(are26Neighbors(leftEnd_, rhs.rightEnd_)) {
        voxels_->invert();
        rhs.voxels_->invert();
        leftEnd_ = rightEnd_;
        rightEnd_ = rhs.leftEnd_;
    } else if(are26Neighbors(rightEnd_, rhs.leftEnd_)) {
        rightEnd_ = rhs.rightEnd_;
    } else if(are26Neighbors(rightEnd_, rhs.rightEnd_)) {
        rhs.voxels_->invert();
        rightEnd_ = rhs.leftEnd_;
    } else {
        tgtAssert(false, "invalid connection");
    }
    voxels_ = std::unique_ptr<RunTree>(new RunNode(std::move(voxels_), std::move(rhs.voxels_)));
}

// MetadataExtractor -----------------------------------------------------------------------


RegularSequence::RegularSequence(std::vector<tgt::svec3> voxels, uint32_t predeterminedComponentId)
    : voxels_(std::move(voxels))
    , predeterminedComponentId_(predeterminedComponentId)
{
}
/// MetaDataCollector --------------------------------------------------------------------------------
const std::string MetaDataCollector::loggerCat_("voreen.vesseltopology.metadatacollector");
void MetaDataCollector::collect(EndData&& data) {
    endPoints_.push_back(data.pos_);
}
void MetaDataCollector::collect(RegularData&& data) {
    std::vector<tgt::svec3> v;
    data.voxels_->collectVoxels(v);
    regularSequences_.emplace_back(v, data.predeterminedComponentId_);
}
void MetaDataCollector::collect(BranchData&& data) {
    std::vector<tgt::svec3> v;
    data.voxels_->collectVoxels(v);
    branchPoints_.push_back(v);
}
std::unique_ptr<ProtoVesselGraph> MetaDataCollector::createProtoVesselGraph(tgt::svec3 dimensions, const tgt::mat4& toRwMatrix, const boost::optional<LZ4SliceVolume<uint8_t>>& sampleMask, ProgressReporter& progress) {
    TaskTimeLogger _("Create Protograph", tgt::Info);
    std::unique_ptr<ProtoVesselGraph> graph(new ProtoVesselGraph(toRwMatrix));


    std::vector<std::pair<tgt::svec3, uint64_t>> nodeVoxelMap;

    KDTreeBuilder<VoxelKDElement> nodeVoxelTreeBuilder;
    // Insert nodes
    for(const auto& p : endPoints_) {
        std::vector<tgt::svec3> voxels;
        voxels.push_back(p);
        uint64_t id = graph->insertNode(std::move(voxels), isAtSampleBorder(p, dimensions));
        nodeVoxelTreeBuilder.push(VoxelKDElement(tgt::ivec3(p),id));
        nodeVoxelMap.push_back(std::make_pair(p, id));
    }
    for(const auto& branchPoint : branchPoints_) {
        bool junctionAtSampleBorder = false;
        std::vector<tgt::svec3> voxels;
        for(const auto& p : branchPoint) {
            voxels.push_back(p);
            junctionAtSampleBorder |= isAtSampleBorder(p, dimensions);
        }
        uint64_t id = graph->insertNode(std::move(voxels), junctionAtSampleBorder);
        for(const auto& p : branchPoint) {
            nodeVoxelTreeBuilder.push(VoxelKDElement(tgt::ivec3(p),id));
        }
        for(const auto& p : branchPoint) {
            nodeVoxelMap.push_back(std::make_pair(p, id));
        }
    }

    KDTreeNeighborhoodFinder<VoxelKDElement> nodeVoxelTree(std::move(nodeVoxelTreeBuilder));

    size_t progressCounter = 0;
    for(const auto& regularSequence : regularSequences_) {
        progress.setProgress(static_cast<float>(progressCounter++)/regularSequences_.size());
        tgtAssert(!regularSequence.voxels_.empty(), "Empty sequence");
        tgt::ivec3 leftEnd = regularSequence.voxels_.front();
        tgt::ivec3 rightEnd = regularSequence.voxels_.back();
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
                std::vector<tgt::svec3> voxels;
                voxels.push_back(tgt::svec3(leftEnd));
                voxels.push_back(tgt::svec3(rightEnd));
                uint64_t newNode = graph->insertNode(std::move(voxels), isAtSampleBorder(tgt::svec3(leftEnd), dimensions) || isAtSampleBorder(tgt::svec3(rightEnd), dimensions));

                nodeVoxelMap.push_back(std::make_pair(leftEnd, newNode));
                nodeVoxelMap.push_back(std::make_pair(rightEnd, newNode));

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
        std::vector<tgt::svec3> voxels(regularSequence.voxels_);

        graph->insertEdge(leftEndNode, rightEndNode, std::move(voxels));
    }

    // Correct sample mask border information if sampleMask is present
    if(sampleMask) {
        LZ4SliceVolumeReader<uint8_t, 1> sampleMaskReader(*sampleMask);
        std::sort(nodeVoxelMap.begin(), nodeVoxelMap.end(), [] (const std::pair<tgt::svec3, uint64_t>& n1, const std::pair<tgt::svec3, uint64_t>& n2) {
                return n1.first.z < n2.first.z;
                });
        for(const auto& pair: nodeVoxelMap) {
            const tgt::svec3& p = pair.first;
            const uint64_t id = pair.second;

            sampleMaskReader.seek(p.z);

            for(int dz = -1; dz <= 1; ++dz) {
                for(int dy = -1; dy <= 1; ++dy) {
                    for(int dx = -1; dx <= 1; ++dx) {
                        if(sampleMaskReader.getSlice(p.z+dz)->voxel(tgt::svec3(p.x+dx, p.y+dy, 0)) == 0) {
                            graph->nodes_.at(id).atSampleBorder_ = true;
                            goto sample_mask_search_done;
                        }
                    }
                }
            }
            sample_mask_search_done: ;
        }
    }

    // Create edges from end voxels directly connected to branchPoints
    for(const auto& node : graph->nodes_) {
        if(node.edges_.size() != 0) {
            continue;
        }

        tgt::ivec3 endPointI(node.voxels_.at(0));
        for(const auto& p : nodeVoxelTree.find_26_neighbors(endPointI)) {
            graph->insertEdge(p->nodeID_, node.id_, std::vector<tgt::svec3>());
            break;
        }
    }
    progress.setProgress(1.0f);
    return graph;
}


/// Connected component related classes ----------------------------------------

/// Row ------------------------------------------------------------------------

Row::Row()
    : endRuns_()
    , regularRuns_()
    , branchRuns_()
{
}

template<typename MetaData>
static void connectRunVecs(std::vector<Run<MetaData>>& r1, std::vector<Run<MetaData>>& r2) {
    auto thisRun = r1.begin();
    auto otherRun = r2.begin();
    while(thisRun != r1.end() && otherRun != r2.end()) {
        thisRun->tryMerge(*otherRun);

        // Advance the run that cannot overlap with the follower of the current other
        // If both end on the voxel, we can advance both.
        const size_t thisUpper = thisRun->pos_.xhigh_;
        const size_t otherUpper = otherRun->pos_.xhigh_;
        if(thisUpper <= otherUpper) {
            ++thisRun;
        }
        if(otherUpper <= thisUpper) {
            ++otherRun;
        }
    }
}

void Row::connect(Row& other) {
    connectRunVecs(endRuns_, other.endRuns_);
    connectRunVecs(regularRuns_, other.regularRuns_);
    connectRunVecs(branchRuns_, other.branchRuns_);
}

void Row::finalizeRuns() {
    for(auto run : endRuns_) {
        run.finalize();
    }
    for(auto run : regularRuns_) {
        run.finalize();
    }
    for(auto run : branchRuns_) {
        run.finalize();
    }
}

/// RowStorage -----------------------------------------------------------------

RowStorage::RowStorage(const tgt::svec3& volumeDimensions)
    : storageSize_(volumeDimensions.y + 2)
    //: storageSize_(tgt::hmul(volumeDimensions.yz()))
    , rowsPerSlice_(volumeDimensions.y)
    , rows_(new Row[storageSize_])
    , storagePos_(-1)
{
}

RowStorage::~RowStorage() {
    delete[] rows_;
}

Row& RowStorage::latest() const {
    return rows_[storagePos_];
}

Row* RowStorage::getRows() const {
    return rows_;
}

Row& RowStorage::get(size_t pos) const {
    return rows_[pos%storageSize_];
}

void RowStorage::finalizeRows() {
    for(size_t i = 0; i < storageSize_; ++i) {
        rows_[i].finalizeRuns();
    }
}

} //namespace voreen
