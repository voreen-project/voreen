/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "voreen/core/datastructures/volume/volume.h"
#include "modules/hdf5/io/hdf5filevolume.h"
#include "modules/hdf5/utils/hdf5utils.h"
#include "voreen/core/io/progressreporter.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include <fstream>
#include <queue>
#include <unordered_map>
#include <boost/iostreams/device/mapped_file.hpp>


#include "tgt/vector.h"

namespace voreen {
#define SC_TEMPLATE template<int ADJACENCY, typename MetaData, typename ClassID, typename SliceType>
#define SC_NS StreamingComponents<ADJACENCY, MetaData, ClassID, SliceType>

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

class CCARunID {
public:
    CCARunID(uint64_t id)
        : id_(id)
        , written_(false)
        , meetsConstraints_(false)
    {
        tgtAssert((id >> 62) == 0, "Invalid id, must be smaller than 2^63");
    }
    uint64_t id() const {
        return id_;
    }
    void markWritten() {
        written_ = 1;
    }
    bool written() const {
        return written_ > 0;
    }
    void setMeetsConstraints(bool v) {
        meetsConstraints_ = v;
    }
    bool meetsConstraints() const {
        return meetsConstraints_;
    }
private:
    uint64_t id_ : 62;
    uint64_t written_: 1;
    uint64_t meetsConstraints_: 1;
};

template<int ADJACENCY, typename MetaData, typename ClassID, typename SliceType=VolumeRAM>
class StreamingComponents {

private:
    class RootFile;
public:
    typedef std::function<boost::optional<ClassID>(const SliceType& vol, tgt::svec3 pos)> getClassFunc;
    typedef std::function<bool(const MetaData&)> componentConstraintTest;
    typedef std::function<void(uint32_t id, const MetaData&)> ComponentCompletionCallback;
    typedef std::function<bool(const MetaData&, const MetaData&)> ComponentComparator;

    template<typename InputType, typename OutputType>
    StreamingComponentsStats cca(const InputType& input, OutputType& output, ComponentCompletionCallback componentCompletionCallback, getClassFunc getClass, bool applyLabeling, componentConstraintTest meetsComponentConstraints, ProgressReporter& progress, ComponentComparator componentComparator = nullptr) const;
private:
    class RowStorage; // Forward declaration

    template<typename outputBaseType, typename InputType, typename OutputType>
    void writeOutputVolume(const RootFile& rootFile, const InputType& input, OutputType& output, getClassFunc getClass, uint64_t& voxelCounter, ProgressReporter& progress) const;


protected:
    static const std::string loggerCat_;

private:
    struct MergerInfo {
        CCARunID idRoot;
        CCARunID idOther;

        bool notesSingleNonConnectedComponent() {
            return idRoot.id() == idOther.id();
        }
    };
    class MergerFile {
    public:
        MergerFile(const std::string& filename)
            : output_(filename, std::fstream::binary | std::fstream::in | std::fstream::out | std::fstream::trunc)
            , filename_(filename)
            , numMergers_(0)
            , metadata_()
        {
        }
        MergerFile(MergerFile&& other)
            : output_(std::move(other.output_))
            , filename_(other.filename_)
            , numMergers_(other.numMergers_)
            , metadata_(std::move(other.metadata_))
        {
            other.filename_ = "";
            other.numMergers_ = -1;
        }
        ~MergerFile() {
            output_.close();
            tgt::FileSystem::deleteFile(filename_);
        }
        void noteMerger(MergerInfo e, MetaData metadata) {
            uint64_t rootid = e.idRoot.id();
            uint64_t otherid = e.idOther.id();
            metadata_[rootid] = metadata;
            if(!e.notesSingleNonConnectedComponent()) {
                metadata_.erase(otherid);
            }
            output_.write(reinterpret_cast<char*>(&e), sizeof(e));
            tgtAssert(output_.good(), "Write failed");
            ++numMergers_;
        }
        std::fstream& getFile() & {
            output_.flush();
            return output_;
        }
        uint64_t getNumMergers() const {
            return numMergers_;
        }

        const std::unordered_map<uint64_t, MetaData>& getMetadata() const {
            return metadata_;
        }
    private:
        std::fstream output_;
        std::string filename_;
        uint64_t numMergers_;
        // Metadata of all root components, mapped via the runid
        std::unordered_map<uint64_t, MetaData> metadata_;
    };
    class RootFile {
    public:
        RootFile(MergerFile&& mergerFile, const std::string& filename, uint64_t numUsedIds, ComponentCompletionCallback cccallback, ComponentComparator componentComparator);
        ~RootFile();
        uint32_t getRootID(uint64_t) const;
        uint32_t getMaxRootId() const;
        const std::vector<uint32_t>& getFinalIdRemappingTable() const;
    private:
        boost::iostreams::mapped_file file_;
        std::string filename_;
        uint32_t maxRootID_;
        std::vector<uint32_t> finalIdRemappingTable_;
    };

    class RunComposition;
    class Node {
    public:
        Node(CCARunID id);
        virtual ~Node();
        virtual MetaData getMetaData() const = 0;
        virtual void addNode(Node* newRoot) = 0;

        typename SC_NS::Node* getRootNode();
        void setParent(RunComposition* parent);
        RunComposition* getParent() const;
        CCARunID getComponentIdForMerger() {
            CCARunID res = id_;
            id_.markWritten();
            return res;
        }
        uint64_t getId() const {
            return id_.id();
        }

    protected:
        RunComposition* parent_;
        CCARunID id_;
    };

    class Run;
    class RunComposition final : public Node {
    public:
        RunComposition(CCARunID id, Node* r1, Node* r2);
        ~RunComposition() {}
        void addNode(Node* other);
        MetaData getMetaData() const;
        typename SC_NS::RunComposition* getRoot();
        void ref();
        void unref();
    private:
        MetaData metaData_;
        uint32_t refCount_;
    };

    class Run final : public Node {
    public:
        Run(uint64_t id, tgt::svec2 yzPos_, size_t lowerBound_, size_t upperBound_, ClassID cls);
        ~Run();

        // If merge between runs is successfull, return the mergerInfo and the metadata
        template<int ROW_MH_DIST>
        boost::optional<std::pair<MergerInfo, MetaData>> tryMerge(Run& other, const componentConstraintTest& cctest);

        MetaData getMetaData() const;
        void addNode(Node* other);

        const tgt::svec2 yzPos_;
        const size_t lowerBound_;
        const size_t upperBound_;
        const ClassID class_;
    };

    class Row {
    public:
        void finalize(MergerFile& mergerFile, const componentConstraintTest& cctest);
        void init(const SliceType& slice, size_t sliceNum, size_t rowNum, getClassFunc getClass, uint64_t& idCounter);
        Row();
        ~Row() {}

        template<int ROWADJACENCY>
        void connect(typename SC_NS::Row& other, MergerFile& mergerFile, const componentConstraintTest& cctest);
        std::vector<Run>& getRuns();
    private:
        std::vector<Run> runs_; ///< Sorted!
    };

    class RowStorage {
    public:
        RowStorage(const tgt::svec3& volumeDimensions, MergerFile& mergerFile, getClassFunc getClass, componentConstraintTest cctest);
        ~RowStorage();
        void add(const SliceType& slice, size_t sliceNum, size_t row, uint64_t& idCounter);
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
        const componentConstraintTest cctest_;
        Row* rows_;
        MergerFile& mergerFile_;
        size_t storagePos_;
    };


};

SC_TEMPLATE
const std::string SC_NS::loggerCat_("voreen.bigdataimageprocessing.streamingcomponents");

SC_TEMPLATE
template<typename outputBaseType, typename InputType, typename OutputType>
void SC_NS::writeOutputVolume(const RootFile& rootFile, const InputType& input, OutputType& output, getClassFunc getClass, uint64_t& voxelCounter, ProgressReporter& progress) const
{
    const std::vector<uint32_t>& idRemappingTable = rootFile.getFinalIdRemappingTable();
    uint64_t runIdCounter = 1;
    tgt::svec3 dim = output.getDimensions();
    for(size_t z = 0; z<dim.z; ++z) {
        progress.setProgress(static_cast<float>(z) / static_cast<float>(dim.z));
        std::unique_ptr<const SliceType> activeLayer(dynamic_cast<const SliceType*>(input.getSlice(z)));
        VolumeAtomic<outputBaseType> slice(tgt::vec3(dim.x, dim.y, 1));
        slice.clear(); //Default initialize with zeros (background)
        tgtAssert(activeLayer, "No slice or invalid type");
        for(size_t y = 0; y<dim.y; ++y) {
            Row row;
            row.init(*activeLayer.get(), z, y, getClass, runIdCounter);
            auto& runs = row.getRuns();
            auto run = runs.begin();
            for(size_t x = 0; x<dim.x;) {
                if(run == runs.end()) {
                    break;
                }

                uint32_t id;
                tgtAssert(run->upperBound_ > x, "overlapping runs");
                if(run->lowerBound_ <= x) {
                    id = rootFile.getRootID(run->getId());
                    // Perform remapping, if available.
                    if(!idRemappingTable.empty()) {
                        id = idRemappingTable[id];
                    }
                    ++voxelCounter;
                } else {
                    id = 0;
                }
                if(std::is_same<outputBaseType, uint32_t>::value) {
                    slice.voxel(x,y,0) = id;
                } else {
                    tgtAssert((std::is_same<outputBaseType, uint8_t>::value), "invalid output type for binary volume");
                    slice.voxel(x,y,0) = id > 0;
                }

                ++x;
                if(x >= run->upperBound_) {
                    tgtAssert(run->upperBound_ == x, "we lost the correct run somehow");
                    ++run;
                }
            }
        }
        output.writeSlices(&slice, z);
    }
    progress.setProgress(1.0f);
}
SC_TEMPLATE
template<typename InputType, typename OutputType>
StreamingComponentsStats SC_NS::cca(const InputType& input, OutputType& output, ComponentCompletionCallback componentCompletionCallback, getClassFunc getClass, bool applyLabeling, componentConstraintTest meetsComponentConstraints, ProgressReporter& progress, ComponentComparator componentComparator) const {
    const tgt::svec3 dim = input.getDimensions();
    tgtAssert(input.getDimensions() == output.getDimensions(), "dimensions of input and output differ");
    tgtAssert(tgt::hand(tgt::greaterThan(input.getDimensions(), tgt::svec3::one)), "Degenerated volume dimensions");

    progress.setProgressRange(tgt::vec2(0, 0.5));

    uint64_t runIdCounter = 1;
    std::string mergerFileName = VoreenApplication::app()->getUniqueTmpFilePath(".merger");
    MergerFile mergers(mergerFileName);

    {
        RowStorage rows(dim, mergers, getClass, meetsComponentConstraints);
        // First layer
        {
            std::unique_ptr<const SliceType> activeLayer(dynamic_cast<const SliceType*>(input.getSlice(0)));
            tgtAssert(activeLayer, "No slice or invalid slice format");
            rows.add(*activeLayer.get(), 0, 0, runIdCounter);
            for(size_t y = 1; y<dim.y; ++y) {
                // Create new row at z=0
                rows.add(*activeLayer.get(), 0, y, runIdCounter);

                // merge with row (-1, 0)
                rows.template connectLatestWith<-1, 0>();
            }
        }

        // The rest of the layers
        for(size_t z = 1; z<dim.z; ++z) {
            progress.setProgress(static_cast<float>(z)/dim.z);
            std::unique_ptr<const SliceType> activeLayer(dynamic_cast<const SliceType*>(input.getSlice(z)));
            tgtAssert(activeLayer, "No slice or invalid type");

            // Create new row at y=0
            rows.add(*activeLayer.get(), z, 0, runIdCounter);

            // merge with row (0, -1)
            rows.template connectLatestWith< 0,-1>();

            // merge with row ( 1,-1)
            rows.template connectLatestWith< 1,-1>();

            for(size_t y = 1; y<dim.y; ++y) {
                // Create new row
                rows.add(*activeLayer.get(), z, y, runIdCounter);

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
        // RowStorage is destructed and all remaining row info is written.
    }

    RootFile rootFile(std::move(mergers), VoreenApplication::app()->getUniqueTmpFilePath(".root"), runIdCounter-1, componentCompletionCallback, componentComparator);

    uint64_t voxelCounter = 0;
    progress.setProgressRange(tgt::vec2(0.5, 1));

    if(applyLabeling) {
        tgtAssert(output.getBaseType() == "uint32", "data type mismatch");
        writeOutputVolume<uint32_t>(rootFile, input, output, getClass, voxelCounter, progress);
    } else {
        tgtAssert(output.getBaseType() == "uint8", "data type mismatch");
        writeOutputVolume<uint8_t>(rootFile, input, output, getClass, voxelCounter, progress);
    }
    progress.setProgress(1.0f);

    const uint32_t numComponents = rootFile.getMaxRootId();
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

struct NodeWithId {
    const CCARunID id;
    NodeWithId* parent;

    uint32_t finalId;

    NodeWithId(CCARunID id)
        : id(id)
        , parent(nullptr)
        , finalId(0)
    {
    }

    NodeWithId* getRoot() {
        NodeWithId* root = this;
        while(root->parent) {
            root = root->parent;
        }
        return root;
    }

    uint32_t getFinalId(uint32_t& finalIdCounter) {
        if(id.meetsConstraints() && finalId == 0) {
            finalId = finalIdCounter++;
        }
        return finalId;
    }
};
struct NodeWithIdMaxHeapComp {
    bool operator()(NodeWithId* lhs, NodeWithId* rhs) {
        tgtAssert(lhs, "Nullptr in heap");
        tgtAssert(rhs, "Nullptr in heap");
        return lhs->id.id() < rhs->id.id();
    }
};

SC_TEMPLATE
SC_NS::RootFile::RootFile(MergerFile&& in, const std::string& filename, uint64_t numUsedIds, ComponentCompletionCallback cccallback, ComponentComparator componentComparator)
    : file_()
    , filename_(filename)
{
    auto next_finalized_node_id = numUsedIds;

    size_t fileSize = next_finalized_node_id * sizeof(next_finalized_node_id);

    // Ensure that fileSize is > 0 so that we create and do not try to open it.
    fileSize = std::max<size_t>(1, fileSize);

    boost::iostreams::mapped_file_params openParams;
    openParams.path = filename_;
    openParams.mode = std::ios::in | std::ios::out;
    openParams.length = fileSize;
    openParams.new_file_size = fileSize; // Option: Do create the file (overwrite if it does not exist)!

    file_.open(openParams);

    uint32_t* data = reinterpret_cast<uint32_t*>(file_.data());

    MergerFile tmp(std::move(in));
    std::fstream& mergerFile = tmp.getFile();
    const auto& metadataMap = tmp.getMetadata();

    MergerInfo mergerIds{0, 0};
    NodeWithId* nodes[2];

    std::unordered_map<uint64_t, NodeWithId*> hash;
    std::priority_queue<NodeWithId*, std::vector<NodeWithId*>, NodeWithIdMaxHeapComp> prior;

    // An optional id sorting based on meta data is performed.
    // We map the run id to its final id based on a user defined sorting.
    std::function<bool(const uint64_t& a, const uint64_t& b)> compare = std::less<uint64_t>();
    if(componentComparator) {
        compare = [metadataMap, componentComparator](const uint64_t& a, const uint64_t& b) {
            return componentComparator(metadataMap.at(a), metadataMap.at(b));
        };
    }
    std::map<uint64_t, uint32_t, decltype(compare)> sortedFinalIds(compare);

    uint32_t finalIdCounter = 1;
    for(int64_t next_merger_id = tmp.getNumMergers()-1; next_merger_id >= 0; next_merger_id--) {
        mergerFile.seekg(next_merger_id*sizeof(mergerIds));
        mergerFile.read(reinterpret_cast<char*>(&mergerIds), sizeof(mergerIds));
        tgtAssert(mergerFile.good(), "Read error");

        if (mergerIds.notesSingleNonConnectedComponent()) {
            prior.push(new NodeWithId(mergerIds.idRoot));
        } else {
            auto handleId = [&] (const CCARunID& id, NodeWithId*& node) {
                auto possibleNode = hash.find(id.id());
                if (possibleNode != hash.end()) {
                    node = possibleNode->second;
                } else {
                    node = new NodeWithId(id);
                    hash.insert({id.id(), node});
                }
                if (!id.written()) {
                    hash.erase(id.id());
                    prior.push(node);
                }
            };
            tgtAssert(mergerIds.idRoot.id() < mergerIds.idOther.id(), "Invalid merge");
            handleId(mergerIds.idOther, nodes[1]);
            handleId(mergerIds.idRoot, nodes[0]);
            nodes[1]->parent = nodes[0];
        }
        while (next_finalized_node_id > 0 && !prior.empty() && prior.top()->id.id() == next_finalized_node_id) {
            NodeWithId* node = prior.top();
            prior.pop();
            NodeWithId* root = node->getRoot();
            uint32_t finalid = root->getFinalId(finalIdCounter);
            if(root == node && finalid != 0) {
                // Node is a root => a finished component
                uint64_t runid = root->id.id();
                sortedFinalIds[runid] = finalid;
            }
            data[next_finalized_node_id] = finalid;
            delete node;
            next_finalized_node_id--;
        }
    }
    maxRootID_ = finalIdCounter - 1;
    tgtAssert(prior.empty(), "Queue not empty");

    // Create remapping table.
    if(componentComparator) {
        finalIdRemappingTable_.resize(finalIdCounter);
        finalIdRemappingTable_[0] = 0; // Zero does not get remapped.
    }
    uint32_t sortedFinalIdCounter = 1;
    for(auto iter = sortedFinalIds.begin(); iter != sortedFinalIds.end(); iter++) {
        const uint64_t& runid = iter->first;
        const uint32_t& finalid = iter->second;

        // If sorting is enabled, we map the old final Id to the sorted final id.
        if(componentComparator) {
            finalIdRemappingTable_[finalid] = sortedFinalIdCounter;
            sortedFinalIdCounter++;
        }

        const auto& entry = metadataMap.find(runid);
        tgtAssert(entry != metadataMap.end(), "No metadata for runid");
        cccallback(sortedFinalIdCounter, entry->second);
    }
}

SC_TEMPLATE
SC_NS::RootFile::~RootFile() {
    file_.close();
    tgt::FileSystem::deleteFile(filename_);
}

SC_TEMPLATE
uint32_t SC_NS::RootFile::getRootID(uint64_t id) const {
    uint32_t* data = reinterpret_cast<uint32_t*>(file_.data());
    return data[id];
}

SC_TEMPLATE
uint32_t SC_NS::RootFile::getMaxRootId() const {
    return maxRootID_;
}

SC_TEMPLATE
const std::vector<uint32_t>& SC_NS::RootFile::getFinalIdRemappingTable() const {
    return finalIdRemappingTable_;
}

SC_TEMPLATE
SC_NS::Node::Node(CCARunID id)
    : parent_(nullptr)
    , id_(id)
{
}
SC_TEMPLATE
SC_NS::Node::~Node() {
    if(parent_) {
        parent_->unref();
    }
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
typename SC_NS::RunComposition* SC_NS::Node::getParent() const {
    return parent_;
}

SC_TEMPLATE
void SC_NS::RunComposition::addNode(Node* other) {
    tgtAssert(other, "newroot is null");
    tgtAssert(!this->parent_, "Parent not null");
    metaData_ += other->getMetaData();
    other->setParent(this);
}

SC_TEMPLATE
SC_NS::RunComposition::RunComposition(CCARunID id, Node* r1, Node* r2)
    : Node(id)
    , metaData_() //using default constructor
    , refCount_(0)
{
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
SC_NS::Run::Run(uint64_t id, tgt::svec2 yzPos, size_t lowerBound, size_t upperBound, ClassID cls)
    : Node(id)
    , yzPos_(yzPos)
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
    RunComposition* newRoot = new RunComposition(Node::id_, this, other);
}

SC_TEMPLATE
template<int ROW_MH_DIST>
boost::optional<std::pair<typename SC_NS::MergerInfo, MetaData>> SC_NS::Run::tryMerge(Run& other, const componentConstraintTest& cctest) {
    if(class_ != other.class_) {
        return boost::none;
    }

    static const int MAX_MH_DIST = 3 - ADJACENCY;

    // The two rows are already too far apart in the yz-dimension
    if(MAX_MH_DIST - ROW_MH_DIST <  0) {
        return boost::none;
    }
    // The two rows are almost to far in the xy-dimension apart, so they have to actually overlap in the x dimension
    if(MAX_MH_DIST - ROW_MH_DIST == 0 && (lowerBound_ >= other.upperBound_ || other.lowerBound_ >= upperBound_)) {
        return boost::none;
    }
    // The two rows close enough in the yz-dimension, so that they only need to be next to each other in the x dimension
    if(MAX_MH_DIST - ROW_MH_DIST >  0 && (lowerBound_ > other.upperBound_ || other.lowerBound_ > upperBound_)) {
        return boost::none;
    }

    Node* thisRoot = this->getRootNode();
    Node* otherRoot = other.getRootNode();
    if(thisRoot == otherRoot) {
        return boost::none;
    }
    // Note: It is important for the creation of the root file that the component with the lower id is the new root!
    Node* root;
    Node* child;
    if(thisRoot->getId() < otherRoot->getId()) {
        root = thisRoot;
        child = otherRoot;
    } else {
        root = otherRoot;
        child = thisRoot;
    }
    auto rootId = root->getComponentIdForMerger();
    auto childId = child->getComponentIdForMerger();
    root->addNode(child);
    rootId.setMeetsConstraints(cctest(root->getMetaData()));
    return std::make_pair(MergerInfo {
        rootId,
        childId,
    }, root->getMetaData());
}

SC_TEMPLATE
MetaData SC_NS::Run::getMetaData() const {
    return MetaData(yzPos_, lowerBound_, upperBound_);
}

SC_TEMPLATE
void SC_NS::Row::finalize(MergerFile& mergerFile, const componentConstraintTest& cctest)
{
    for(auto& run : runs_) {
        if(run.Node::getParent() == nullptr) {
            // Run is not connected to anything
            CCARunID id = run.getComponentIdForMerger();
            const auto& metadata = run.getMetaData();
            id.setMeetsConstraints(cctest(metadata));
            mergerFile.noteMerger(MergerInfo {
                    id,
                    id
                    }, metadata);
        }
    }
    runs_.clear();
}

SC_TEMPLATE
void SC_NS::Row::init(const SliceType& slice, size_t sliceNum, size_t rowNum, getClassFunc getClass, uint64_t& idCounter)
{
    // Insert new runs
    boost::optional<ClassID> prevClass = boost::none;
    size_t runStart = 0;
    size_t rowLength = slice.getDimensions().x;
    tgt::svec2 yzPos(rowNum, sliceNum);
    for(size_t x = 0; x < rowLength; ++x) {
        auto currentClass = getClass(slice, tgt::svec3(x, rowNum, 0));
        if(prevClass != currentClass) {
            if(prevClass) {
                runs_.emplace_back(idCounter++, yzPos, runStart, x, *prevClass);
            }
            if(currentClass) {
                runStart = x;
            }
        }
        prevClass = currentClass;
    }
    if(prevClass) {
        runs_.emplace_back(idCounter++, yzPos, runStart, rowLength, *prevClass);
    }
}


SC_TEMPLATE
SC_NS::Row::Row()
    : runs_()
{
}

SC_TEMPLATE
template<int ROW_MH_DIST>
void SC_NS::Row::connect(typename SC_NS::Row& other, MergerFile& mergerFile, const componentConstraintTest& cctest) {
    auto thisRun = runs_.begin();
    auto otherRun = other.runs_.begin();
    while(thisRun != runs_.end() && otherRun != other.runs_.end()) {
        auto mergeEvent = thisRun->template tryMerge<ROW_MH_DIST>(*otherRun, cctest);
        if(mergeEvent) {
            mergerFile.noteMerger(mergeEvent->first, mergeEvent->second);
        }

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
SC_NS::RowStorage::RowStorage(const tgt::svec3& volumeDimensions, MergerFile& mergerFile, getClassFunc getClass, componentConstraintTest cctest)
    : storageSize_(volumeDimensions.y + 2)
    //: storageSize_(tgt::hmul(volumeDimensions.yz()))
    , rowsPerSlice_(volumeDimensions.y)
    , rows_(new Row[storageSize_])
    , storagePos_(-1)
    , mergerFile_(mergerFile)
    , getClass_(getClass)
    , cctest_(cctest)
{
}

SC_TEMPLATE
SC_NS::RowStorage::~RowStorage() {
    for(size_t i=0; i<storageSize_; ++i) {
        rows_[i].finalize(mergerFile_, cctest_);
    }
    delete[] rows_;
}


SC_TEMPLATE
void SC_NS::RowStorage::add(const SliceType& slice, size_t sliceNum, size_t rowNum, uint64_t& idCounter) {
    storagePos_ = (storagePos_ + 1)%storageSize_;
    rows_[storagePos_].finalize(mergerFile_, cctest_);
    rows_[storagePos_].init(slice, sliceNum, rowNum, getClass_, idCounter);
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
    latest().template connect<DY*DY+DZ*DZ>(latest_DP<DY,DZ>(), mergerFile_, cctest_);
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
