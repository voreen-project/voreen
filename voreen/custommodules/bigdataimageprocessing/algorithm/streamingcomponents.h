/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "voreen/core/datastructures/volume/volume.h"
#include "modules/hdf5/io/hdf5filevolume.h"
#include "modules/hdf5/utils/hdf5utils.h"
#include "voreen/core/io/progressreporter.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include <fstream>


#include "tgt/vector.h"

namespace voreen {
#define SC_TEMPLATE template<int ADJACENCY, typename MetaData, typename ClassID>
#define SC_NS StreamingComponents<ADJACENCY, MetaData, ClassID>

struct CCAVoidLabel {
    CCAVoidLabel() { }
    static boost::optional<CCAVoidLabel> some() {
        return boost::optional<CCAVoidLabel>(CCAVoidLabel());
    }
    bool operator==(const CCAVoidLabel&) const {
        return true;
    }
    bool operator!=(const CCAVoidLabel&) const {
        return false;
    }
};

struct StreamingComponentsStats {
    uint32_t numComponents;
    uint64_t numVoxels;
};

SC_TEMPLATE
class StreamingComponents {

public:
    typedef std::function<boost::optional<ClassID>(const VolumeRAM* vol, tgt::svec3 pos)> getClassFunc;
    typedef std::function<bool(const MetaData&)> componentConstraintTest;
    typedef std::function<void(uint32_t id, const MetaData&)> ComponentCompletionCallback;

    template<typename InputType, typename OutputType>
    StreamingComponentsStats cca(const InputType& input, OutputType& output, ComponentCompletionCallback componentCompletionCallback, getClassFunc getClass, bool applyLabeling, componentConstraintTest meetsComponentConstraints, ProgressReporter& progress) const;
private:
    class RowStorage; // Forward declaration

    template<typename outputBaseType, typename OutputType>
    void writeRowsToStorage(RowStorage& rows, OutputType& output, ComponentCompletionCallback componentCompletionCallback, componentConstraintTest meetsComponentConstraints, bool applyLabeling, uint32_t& idCounter, uint64_t& voxelCounter, ProgressReporter& progress) const;


protected:
    static const std::string loggerCat_;

private:
    class RunComposition;
    class Node {
    public:
        Node();
        virtual ~Node();
        virtual MetaData getMetaData() const = 0;
        virtual void addNode(Node* newRoot) = 0;
        virtual uint32_t getRootAptitude() const = 0;

        typename SC_NS::Node* getRootNode();
        uint32_t assignID(uint32_t& idCounter);
        void setParent(RunComposition* parent);

    protected:
        RunComposition* parent_;
        uint32_t id_;
    };

    class Run;
    class RunComposition final : public Node {
    public:
        RunComposition(Node* r1, Node* r2);
        ~RunComposition() {}
        void addNode(Node* newRoot);
        MetaData getMetaData() const;
        typename SC_NS::RunComposition* getRoot();
        void ref();
        void unref();
        uint32_t getRootAptitude() const;
    private:
        MetaData metaData_;
        uint32_t refCount_;
    };

    class Run final : public Node {
    public:
        Run(tgt::svec2 yzPos_, size_t lowerBound_, size_t upperBound_, ClassID cls);
        ~Run();

        template<int ROW_MH_DIST>
        void tryMerge(Run& other);

        MetaData getMetaData() const;
        void addNode(Node* newRoot);
        uint32_t getRootAptitude() const;

        const tgt::svec2 yzPos_;
        const size_t lowerBound_;
        const size_t upperBound_;
        const ClassID class_;
    };

    class Row {
    public:
        void init(const VolumeRAM* slice, size_t sliceNum, size_t rowNum, getClassFunc getClass);
        Row();
        ~Row() {}

        template<int ROWADJACENCY>
        void connect(Row& other);
        std::vector<Run>& getRuns();
    private:
        std::vector<Run> runs_; ///< Sorted!
    };

    class RowStorage {
    public:
        RowStorage(const tgt::svec3& volumeDimensions, getClassFunc getClass);
        ~RowStorage();
        void add(const VolumeRAM* slice, size_t sliceNum, size_t row);
        Row& latest() const;

        template<int DY, int DZ>
        Row& latest_DP() const;

        template<int DY, int DZ>
        void connectLatestWith();

        // Only for debug purposes
        Row* getRows() const;
    private:
        Row& get(size_t pos) const;

        const size_t storageSize_;
        const size_t rowsPerSlice_;
        const getClassFunc getClass_;
        Row* rows_;
        size_t storagePos_;
    };


};

SC_TEMPLATE
const std::string SC_NS::loggerCat_("voreen.bigdataimageprocessing.streamingcomponents");

SC_TEMPLATE
template<typename outputBaseType, typename OutputType>
void SC_NS::writeRowsToStorage(RowStorage& rows, OutputType& output, ComponentCompletionCallback componentCompletionCallback, componentConstraintTest meetsComponentConstraints, bool applyLabeling, uint32_t& idCounter, uint64_t& voxelCounter, ProgressReporter& progress) const {

    const tgt::svec3 dim = output.getDimensions();
    VolumeAtomic<outputBaseType> slice(tgt::vec3(dim.x, dim.y, 1));

    for(size_t z = 0; z<dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);
        // Initialize slice with 0s
        slice.clear();
        for(size_t y = 0; y<dim.y; ++y) {
            Row& currentRow = rows.getRows()[z*dim.y+y];
            for(auto& run : currentRow.getRuns()) {
                uint32_t id = 0;
                Node* n = run.getRootNode();
                MetaData metaData = n->getMetaData();
                if(meetsComponentConstraints(metaData)) {
                    uint32_t prevIdCounter = idCounter;
                    id = n->assignID(idCounter);

                    // This can be moved into the destructor of the nodes if we ever decide not to hold all rows in memory
                    if(prevIdCounter != idCounter) {
                        componentCompletionCallback(id, metaData);
                    }

                    for(size_t x = run.lowerBound_; x < run.upperBound_; ++x) {
                        slice.voxel(x,y,0) = applyLabeling ? id : 1;
                        ++voxelCounter;
                    }

                }
            }
        }
        output.writeSlices(&slice, z);
    }
}

SC_TEMPLATE
template<typename InputType, typename OutputType>
StreamingComponentsStats SC_NS::cca(const InputType& input, OutputType& output, ComponentCompletionCallback componentCompletionCallback, getClassFunc getClass, bool applyLabeling, componentConstraintTest meetsComponentConstraints, ProgressReporter& progress) const {
    const tgt::svec3 dim = input.getDimensions();
    tgtAssert(input.getDimensions() == output.getDimensions(), "dimensions of input and output differ");
    tgtAssert(tgt::hand(tgt::greaterThan(input.getDimensions(), tgt::svec3::one)), "Degenerated volume dimensions");

    progress.setProgressRange(tgt::vec2(0, 0.5));

    RowStorage rows(dim, getClass);
    // First layer
    {
        std::unique_ptr<const VolumeRAM> activeLayer(input.getSlice(0));
        rows.add(activeLayer.get(), 0, 0);
        for(size_t y = 1; y<dim.y; ++y) {
            // Create new row at z=0
            rows.add(activeLayer.get(), 0, y);

            // merge with row (-1, 0)
            rows.template connectLatestWith<-1, 0>();
        }
    }

    // The rest of the layers
    for(size_t z = 1; z<dim.z; ++z) {
        progress.setProgress(static_cast<float>(z)/dim.z);
        std::unique_ptr<const VolumeRAM> activeLayer(input.getSlice(z));

        // Create new row at y=0
        rows.add(activeLayer.get(), z, 0);

        // merge with row (0, -1)
        rows.template connectLatestWith< 0,-1>();

        // merge with row ( 1,-1)
        rows.template connectLatestWith< 1,-1>();

        for(size_t y = 1; y<dim.y; ++y) {
            // Create new row
            rows.add(activeLayer.get(), z, y);

            // merge with row (-1, 0)
            rows.template connectLatestWith<-1, 0>();

            // merge with row ( 0,-1)
            rows.template connectLatestWith< 0,-1>();

            if(y != dim.y-1) {
                // merge with row ( 1,-1), but only if we are not at the end of the slice
                rows.template connectLatestWith< 1,-1>();
            }

            // merge with row (-1,-1)
            rows.template connectLatestWith<-1,-1>();
        }
    }

    uint32_t idCounter = 1;
    uint64_t voxelCounter = 0;

    progress.setProgressRange(tgt::vec2(0.5, 1));
    if(applyLabeling) {
        tgtAssert(output.getBaseType() == "uint32", "data type mismatch");
        writeRowsToStorage<uint32_t>(rows, output, componentCompletionCallback, meetsComponentConstraints, applyLabeling, idCounter, voxelCounter, progress);
    } else {
        tgtAssert(output.getBaseType() == "uint8", "data type mismatch");
        writeRowsToStorage<uint8_t>(rows, output, componentCompletionCallback, meetsComponentConstraints, applyLabeling, idCounter, voxelCounter, progress);
    }
    progress.setProgress(1.0f);

    const uint32_t numComponents = idCounter - 1;
    const float minValue = voxelCounter < tgt::hmul(dim) ? 0.0f : 1.0f; // Check if there are any background voxels at all
    const float maxValue = voxelCounter > 0 ? (applyLabeling ? numComponents : 1.0f) : 0.0f; // Check if there are any foreground voxels, if there are, check if we applied labeling
    std::unique_ptr<VolumeRAM> helper(VolumeFactory().create(output.getBaseType(), tgt::vec3(1))); // Used to get element range
    const tgt::vec2 range = helper->elementRange(); // Used for normalization

    // Write metadata we can save from input or have determined during creation of the volume
    output.writeSpacing(input.getSpacing());
    output.writeOffset(input.getOffset());
    output.writePhysicalToWorldTransformation(input.getPhysicalToWorldMatrix());
    output.writeRealWorldMapping(RealWorldMapping::createDenormalizingMapping(output.getBaseType()));
    const VolumeMinMax vmm(minValue,  maxValue, (minValue+range.x)/(range.x+range.y), (maxValue+range.x)/(range.x+range.y));
    output.writeVolumeMinMax(&vmm);
    return StreamingComponentsStats {
        numComponents,
        voxelCounter
    };
}

SC_TEMPLATE
SC_NS::Node::Node()
    : parent_(nullptr)
    , id_(0)
{
}
SC_TEMPLATE
SC_NS::Node::~Node() {
    if(parent_) {
        parent_->unref();
    } else {
        // TODO output something
        //std::cout << "Node with volume " << volume_ << std::endl;
    }
}
SC_TEMPLATE
uint32_t SC_NS::RunComposition::getRootAptitude() const {
    return refCount_;
}

SC_TEMPLATE
typename SC_NS::Node* SC_NS::Node::getRootNode() {
    if(!parent_) {
        return this;
    }
    setParent(parent_->getRoot());
    return parent_;
}

SC_TEMPLATE
uint32_t SC_NS::Node::assignID(uint32_t& idCounter) {
    if(id_ == 0) {
        id_ = idCounter++;
    }
    return id_;
}

SC_TEMPLATE
void SC_NS::Node::setParent(typename SC_NS::RunComposition* newParent) {
    tgtAssert(newParent, "newParent is null");
    tgtAssert(newParent != this, "newParent=this");
    // First ref the new parent, THEN unref the old in case they are the same
    RunComposition* prevParent = parent_;
    parent_ = newParent;
    parent_->ref();
    if(prevParent) {
        prevParent->unref();
    }
}

SC_TEMPLATE
void SC_NS::RunComposition::addNode(Node* other) {
    tgtAssert(other, "newroot is null");
    tgtAssert(!this->parent_, "Parent not null");
    metaData_ += other->getMetaData();
    other->setParent(this);
}

SC_TEMPLATE
SC_NS::RunComposition::RunComposition(Node* r1, Node* r2)
    : metaData_() //using default constructor
    , refCount_(0)
{
    //TODO just do this in a metadataconstructor?
    addNode(r1);
    addNode(r2);
}

SC_TEMPLATE
MetaData SC_NS::RunComposition::getMetaData() const {
    return metaData_;
}
SC_TEMPLATE
typename SC_NS::RunComposition* SC_NS::RunComposition::getRoot() {
    if(!this->parent_) {
        return this;
    }
    this->setParent(this->parent_->getRoot());
    return this->parent_;
}
SC_TEMPLATE
void SC_NS::RunComposition::ref() {
    refCount_++;
}

SC_TEMPLATE
void SC_NS::RunComposition::unref() {
    if(!--refCount_) {
        // Nobody likes me :(
        delete this;
    }
}

SC_TEMPLATE
SC_NS::Run::Run(tgt::svec2 yzPos, size_t lowerBound, size_t upperBound, ClassID cls)
    : yzPos_(yzPos)
    , lowerBound_(lowerBound) // first voxel part of run
    , upperBound_(upperBound) // first voxel not part of run
    , class_(cls)
{
}

SC_TEMPLATE
SC_NS::Run::~Run() {
    // Is called in superclass
    /*
    if(parent_) {
        parent_->unref();
    } else {
        //std::cout << "Single row with volume " << (upperBound_-lowerBound_) << std::endl;
    }
    */
}

SC_TEMPLATE
void SC_NS::Run::addNode(Node* other) {
    tgtAssert(other, "newroot is null");
    tgtAssert(!this->parent_, "parent is not null");
    // Construct a new root
    RunComposition* newRoot = new RunComposition(this, other);
}

SC_TEMPLATE
template<int ROW_MH_DIST>
void SC_NS::Run::tryMerge(Run& other) {
    if(class_ != other.class_) {
        return;
    }

    static const int MAX_MH_DIST = 3 - ADJACENCY;

    // The two rows are already too far apart in the yz-dimension
    if(MAX_MH_DIST - ROW_MH_DIST <  0) {
        return;
    }
    // The two rows are almost to far in the xy-dimension apart, so they have to actually overlap in the x dimension
    if(MAX_MH_DIST - ROW_MH_DIST == 0 && (lowerBound_ >= other.upperBound_ || other.lowerBound_ >= upperBound_)) {
        return;
    }
    // The two rows close enough in the yz-dimension, so that they only need to be next to each other in the x dimension
    if(MAX_MH_DIST - ROW_MH_DIST >  0 && (lowerBound_ > other.upperBound_ || other.lowerBound_ > upperBound_)) {
        return;
    }

    Node* thisRoot = this->getRootNode();
    Node* otherRoot = other.getRootNode();
    if(thisRoot == otherRoot) {
        return;
    }
    if(thisRoot->getRootAptitude() > otherRoot->getRootAptitude()) {
        thisRoot->addNode(otherRoot);
    } else {
        otherRoot->addNode(thisRoot);
    }
}

SC_TEMPLATE
MetaData SC_NS::Run::getMetaData() const {
    return MetaData(yzPos_, lowerBound_, upperBound_);
}

SC_TEMPLATE
void SC_NS::Row::init(const VolumeRAM* slice, size_t sliceNum, size_t rowNum, getClassFunc getClass)
{
    // Finalize previous:
    runs_.clear();

    // Insert new runs
    boost::optional<ClassID> prevClass = boost::none;
    size_t runStart = 0;
    size_t rowLength = slice->getDimensions().x;
    tgt::svec2 yzPos(rowNum, sliceNum);
    for(size_t x = 0; x < rowLength; ++x) {
        auto currentClass = getClass(slice, tgt::svec3(x, rowNum, 0));
        if(prevClass != currentClass) {
            if(prevClass) {
                runs_.emplace_back(yzPos, runStart, x, *prevClass);
            }
            if(currentClass) {
                runStart = x;
            }
        }
        prevClass = currentClass;
    }
    if(prevClass) {
        runs_.emplace_back(yzPos, runStart, rowLength, *prevClass);
    }
}

SC_TEMPLATE
uint32_t SC_NS::Run::getRootAptitude() const {
    return 0;
}


SC_TEMPLATE
SC_NS::Row::Row()
    : runs_()
{
}

SC_TEMPLATE
template<int ROW_MH_DIST>
void SC_NS::Row::connect(typename SC_NS::Row& other) {
    auto thisRun = runs_.begin();
    auto otherRun = other.runs_.begin();
    while(thisRun != runs_.end() && otherRun != other.runs_.end()) {
        thisRun->template tryMerge<ROW_MH_DIST>(*otherRun);

        // Advance the run that cannot overlap with the follower of the current other
        // If both end on the voxel, we can advance both.
        const size_t thisUpper = thisRun->upperBound_;
        const size_t otherUpper = otherRun->upperBound_;
        if(thisUpper <= otherUpper) {
            ++thisRun;
        }
        if(otherUpper <= thisUpper) {
            ++otherRun;
        }
    }
}

SC_TEMPLATE
std::vector<typename SC_NS::Run>& SC_NS::Row::getRuns() {
    return runs_;
}

SC_TEMPLATE
SC_NS::RowStorage::RowStorage(const tgt::svec3& volumeDimensions, getClassFunc getClass)
    //: storageSize_(volumeDimensions.y + 2)
    : storageSize_(tgt::hmul(volumeDimensions.yz()))
    , rowsPerSlice_(volumeDimensions.y)
    , rows_(new Row[storageSize_])
    , storagePos_(-1)
    , getClass_(getClass)
{
}

SC_TEMPLATE
SC_NS::RowStorage::~RowStorage() {
    delete[] rows_;
}


SC_TEMPLATE
void SC_NS::RowStorage::add(const VolumeRAM* slice, size_t sliceNum, size_t rowNum) {
    storagePos_ = (storagePos_ + 1)%storageSize_;
    rows_[storagePos_].init(slice, sliceNum, rowNum, getClass_);
}

SC_TEMPLATE
typename SC_NS::Row& SC_NS::RowStorage::latest() const {
    return rows_[storagePos_];
}

SC_TEMPLATE
template<int DY, int DZ>
typename SC_NS::Row& SC_NS::RowStorage::latest_DP() const {
    static_assert(-1 <= DY && DY <= 1, "Invalid DY");
    static_assert(-1 <= DZ && DZ <= 0, "Invalid DZ");
    return get(storagePos_ + storageSize_ + DY + rowsPerSlice_ * DZ);
}

SC_TEMPLATE
template<int DY, int DZ>
void SC_NS::RowStorage::connectLatestWith() {
    latest().template connect<DY*DY+DZ*DZ>(latest_DP<DY,DZ>());
}

SC_TEMPLATE
typename SC_NS::Row* SC_NS::RowStorage::getRows() const {
    return rows_;
}

SC_TEMPLATE
typename SC_NS::Row& SC_NS::RowStorage::get(size_t pos) const {
    return rows_[pos%storageSize_];
}

} //namespace voreen
