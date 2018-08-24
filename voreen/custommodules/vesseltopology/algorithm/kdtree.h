#include "tgt/vector.h"
#include <fstream>
#include <boost/iostreams/device/mapped_file.hpp>

namespace voreen {
namespace kdtree {

//struct ExampleElement {
//    typedef float CoordType;
//    const tgt::Vector3<CoordType>& getPos() const;
//};

template<typename Element>
class ElementArrayView;

template<typename Element>
class ElementArray { //Memory mapped version of file created from ElementArrayBuilder
public:
    ElementArray(const std::string& filename, size_t numElements);
    ~ElementArray();
    ElementArrayView<Element> view();
    Element* data();
private:
    boost::iostreams::mapped_file file_;
    std::string filename_;
    size_t numElements_;
};

template<typename Element>
class ElementArrayView {
public:
    ElementArrayView(Element* data, Element* end);
    std::tuple<ElementArrayView<Element>, Element, ElementArrayView<Element>> split(int dimension);
    size_t size() const;

private:
    void sort(int dimension);
    Element* begin_;
    Element* end_;
};


template<typename Element>
class ElementArrayBuilder { //Pushes elements to end of file
public:
    ElementArrayBuilder(const std::string& filename);
    ElementArrayBuilder(ElementArrayBuilder&& other);
    void push(const Element& element); //push to end of file
    ElementArray<Element> finalize() &&;
private:
    std::ofstream file_;
    std::string filename_;
    size_t numElements_;
};

template<typename Element>
struct Node {
    Node(Element elm, size_t left_child, size_t right_child);

    Element elm_;
    size_t left_child_;
    size_t right_child_;
};

template<typename Element>
class NodeStorage { //Memory mapped version of file created from ElementArrayBuilder
public:
    NodeStorage();
    NodeStorage(const std::string& filename, size_t numElements);
    NodeStorage(NodeStorage&& other);
    NodeStorage<Element>& operator=(NodeStorage&& other);

    ~NodeStorage();
    const Node<Element>& operator[](size_t index) const;
private:
    boost::iostreams::mapped_file_source file_;
    std::string filename_;
    size_t numElements_;
};

template<typename Element>
class NodeStorageBuilder {
public:
    NodeStorageBuilder(const std::string& filename);
    NodeStorageBuilder(NodeStorageBuilder&& other);

    size_t push(const Node<Element>& node); //push to end of file
    NodeStorage<Element> finalize() &&;
private:
    std::ofstream file_;
    std::string filename_;
    size_t numNodes_;
};

template<typename Element>
class Tree {
public:
    typedef tgt::Vector3<typename Element::CoordType> PosType;
    Tree(const std::string& storagefilename, ElementArrayBuilder<Element>&& elements);

    std::pair<typename Element::CoordType, const Element&> findNearest(const PosType& pos) const;

    const Node<Element>& root() const;
private:
    NodeStorage<Element> nodes_;
    size_t root_;
};

/// Implementation -------------------------------------------------------------

/// ElementArrayBuilder --------------------------------------------------------
template<typename Element>
ElementArrayBuilder<Element>::ElementArrayBuilder(const std::string& filename)
    : file_(filename, std::ofstream::binary | std::ofstream::out | std::ofstream::trunc)
    , filename_(filename)
    , numElements_(0)
{
}
template<typename Element>
ElementArrayBuilder<Element>::ElementArrayBuilder(ElementArrayBuilder&& other)
    : file_(std::move(other.file_))
    , filename_(std::move(other.filename_))
    , numElements_(other.numElements_)
{
}

template<typename Element>
void ElementArrayBuilder<Element>::push(const Element& elm) {
    file_.write(reinterpret_cast<const char*>(&elm), sizeof(Element));
    ++numElements_;
    tgtAssert(file_.good(), "Write failed");
}

template<typename Element>
ElementArray<Element> ElementArrayBuilder<Element>::finalize() && {
    auto tmp = std::move(*this);
    tmp.file_.flush();
    tmp.file_.close();
    return ElementArray<Element>(tmp.filename_, tmp.numElements_);
}

/// ElementArray ---------------------------------------------------------------
template<typename Element>
ElementArray<Element>::ElementArray(const std::string& filename, size_t numElements)
    : file_()
    , filename_(filename)
    , numElements_(numElements)
{
    size_t fileSize = numElements_ * sizeof(Element);
    boost::iostreams::mapped_file_params openParams;
    openParams.path = filename_;
    openParams.mode = std::ios::in | std::ios::out;
    openParams.length = fileSize;

    file_.open(openParams);
    tgtAssert(file_.is_open(), "File not open");
}
template<typename Element>
ElementArray<Element>::~ElementArray() {
    file_.close();
    //tgt::FileSystem::deleteFile(filename_); //TODO
}

template<typename Element>
ElementArrayView<Element> ElementArray<Element>::view() {
    Element* begin = data();
    Element* end = begin+numElements_;
    return ElementArrayView<Element>(begin, end);
}

template<typename Element>
Element* ElementArray<Element>::data()
{
    return reinterpret_cast<Element*>(file_.data());
}

/// ElementArrayView -----------------------------------------------------------
template<typename Element>
ElementArrayView<Element>::ElementArrayView(Element* begin, Element* end)
    : begin_(begin)
    , end_(end)
{
}

template<typename Element>
std::tuple<ElementArrayView<Element>, Element, ElementArrayView<Element>> ElementArrayView<Element>::split(int dimension)
{
    size_t dist = size();
    tgtAssert(dist > 0, "tried to split empty array view");
    sort(dimension);
    Element* center = begin_ + dist/2;
    return {
        ElementArrayView(begin_, center),
        *center,
        ElementArrayView(center+1, end_)
    };
}

template<typename Element>
size_t ElementArrayView<Element>::size() const {
    tgtAssert(begin_ <= end_, "Invalid begin/end");
    return std::distance(begin_, end_);
}

template<typename Element, int dim>
struct SortElementsInDim {
    bool operator()(const Element& e1, const Element& e2) {
        static_assert(0 <= dim <= 2, "invalid dim");
        return e1.getPos().elem[dim] < e2.getPos().elem[dim];
    }
};

template<typename Element>
void ElementArrayView<Element>::sort(int dimension)
{
    switch(dimension) {
        case 0:
            std::sort(begin_, end_, SortElementsInDim<Element, 0>());
            break;
        case 1:
            std::sort(begin_, end_, SortElementsInDim<Element, 1>());
            break;
        case 2:
            std::sort(begin_, end_, SortElementsInDim<Element, 2>());
            break;
        default:
            ;
            tgtAssert(false, "Invalid dimension");
    }
}

/// Node -----------------------------------------------------------------------
template<typename Element>
Node<Element>::Node(Element elm, size_t left_child, size_t right_child)
    : elm_(elm)
    , left_child_(left_child)
    , right_child_(right_child)
{
}

/// NodeStorageBuilder ---------------------------------------------------------
template<typename Element>
NodeStorageBuilder<Element>::NodeStorageBuilder(const std::string& filename)
    : file_(filename, std::ofstream::binary | std::ofstream::out | std::ofstream::trunc)
    , filename_(filename)
    , numNodes_(0)
{
}
template<typename Element>
NodeStorageBuilder<Element>::NodeStorageBuilder(NodeStorageBuilder&& other)
    : file_(std::move(other.file_))
    , filename_(std::move(other.filename_))
    , numNodes_(other.numNodes_)
{
}

template<typename Element>
size_t NodeStorageBuilder<Element>::push(const Node<Element>& node) {
    file_.write(reinterpret_cast<const char*>(&node), sizeof(Node<Element>));
    tgtAssert(file_.good(), "Write failed");
    return numNodes_++;
}

template<typename Element>
NodeStorage<Element> NodeStorageBuilder<Element>::finalize() && {
    auto tmp = std::move(*this);
    tmp.file_.flush();
    tmp.file_.close();
    return NodeStorage<Element>(tmp.filename_, tmp.numNodes_);
}

/// NodeStorage ----------------------------------------------------------------
template<typename Element>
NodeStorage<Element>::NodeStorage()
    : file_()
    , filename_()
    , numElements_(0)
{
}
template<typename Element>
NodeStorage<Element>::NodeStorage(const std::string& filename, size_t numElements)
    : file_(filename)
    , filename_(filename)
    , numElements_(numElements)
{
    tgtAssert(file_.is_open(), "File not open");
}

template<typename Element>
NodeStorage<Element>::NodeStorage(NodeStorage&& other)
    : file_(std::move(other.file_))
    , filename_(std::move(other.filename_))
    , numElements_(other.numElements_)
{
    //Somehow moving a file source only copies a reference to it, thus later closing our file in the destructor
    //REALLY make sure that the other does not happen here:
    other.file_ = boost::iostreams::mapped_file_source();
    tgtAssert(file_.is_open(), "File not open");
    tgtAssert(!other.file_.is_open(), "File not open");
}
template<typename Element>
NodeStorage<Element>& NodeStorage<Element>::operator=(NodeStorage&& other) {
    if(&other != this) {
        this->~NodeStorage();
        new(this) NodeStorage(std::move(other));
        tgtAssert(file_.is_open(), "File not open");
    }
    return *this;
}

template<typename Element>
NodeStorage<Element>::~NodeStorage() {
    file_.close();
    //tgt::FileSystem::deleteFile(filename_); //TODO
}
template<typename Element>
const Node<Element>& NodeStorage<Element>::operator[](size_t index) const {
    tgtAssert(index < numElements_, "Invalid index");
    const Node<Element>* data = reinterpret_cast<const Node<Element>*>(file_.data());
    return data[index];
}

/// Tree -----------------------------------------------------------------------
const size_t NO_NODE_ID = -1;

template<typename Element>
static size_t buildTree(NodeStorageBuilder<Element>& nodes, ElementArrayView<Element>& elements, int depth) {
    if(elements.size() == 0) {
        return NO_NODE_ID;
    } else {
        int dim = depth % 3;
        auto res = elements.split(dim);
        size_t left_node = buildTree(nodes, std::get<0>(res), depth+1);
        size_t right_node = buildTree(nodes, std::get<2>(res), depth+1);
        return nodes.push(Node<Element>(std::get<1>(res), left_node, right_node));
    }
}

template<typename Element>
Tree<Element>::Tree(const std::string& storagefilename, ElementArrayBuilder<Element>&& elements)
{
    auto elems = std::move(elements).finalize();
    auto elem_view = elems.view();

    NodeStorageBuilder<Element> nodeStorageBuilder(storagefilename);

    root_ = buildTree(nodeStorageBuilder, elem_view, 0);
    nodes_ = std::move(nodeStorageBuilder).finalize();
}

template<typename Element>
const Node<Element>& Tree<Element>::root() const {
    return nodes_[root_];
}


} //namespace KDTree
} //namespace voreen
