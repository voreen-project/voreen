#pragma once

#include "tgt/vector.h"
#include "tgt/filesystem.h"
#include <fstream>
#include <vector>
#include <boost/iostreams/device/mapped_file.hpp>

namespace voreen {
namespace static_kdtree {

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
    typedef tgt::Vector3<typename Element::CoordType> PosType;
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
    size_t size() const;
private:
    boost::iostreams::mapped_file_source file_;
    std::string filename_;
    size_t numNodes_;
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
struct SearchResult {
    typename Element::CoordType distSq_;
    const Element* element_;

    SearchResult(typename Element::CoordType distSq, const Element* element);
    static SearchResult none();

    bool found() const;
    typename Element::CoordType dist() const;

    void takeBetter(SearchResult&& other);
};

template<typename Element>
struct SearchResultSet {
    typename Element::CoordType distSq_;
    std::vector<const Element*> elements_;

    SearchResultSet(typename Element::CoordType distSq, const Element* element);
    static SearchResultSet none();

    bool found() const;
    typename Element::CoordType dist() const;

    void takeBetter(SearchResultSet&& other);
private:
    SearchResultSet();
};

template<typename Element>
class Tree {
public:
    static_assert(std::numeric_limits<typename Element::CoordType>::is_signed, "Element pos type must be signed");
    typedef typename Node<Element>::PosType PosType;
    Tree(const std::string& storagefilename, ElementArrayBuilder<Element>&& elements);

    SearchResult<Element> findNearest(const PosType& pos) const;
    SearchResultSet<Element> findAllNearest(const PosType& pos) const;

    const Node<Element>& root() const;
private:
    template<typename Result>
    void findNearestFrom(const PosType& pos, int depth, size_t node, Result& best_result) const;

    NodeStorage<Element> nodes_;
    size_t root_;
};

/// Implementation -------------------------------------------------------------

/// Impl: ElementArrayBuilder --------------------------------------------------------
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

/// Impl: ElementArray ---------------------------------------------------------------
template<typename Element>
ElementArray<Element>::ElementArray(const std::string& filename, size_t numElements)
    : file_()
    , filename_(filename)
    , numElements_(numElements)
{
    if(numElements_ > 0) {
        size_t fileSize = numElements_ * sizeof(Element);
        boost::iostreams::mapped_file_params openParams;
        openParams.path = filename_;
        openParams.mode = std::ios::in | std::ios::out;
        openParams.length = fileSize;

        file_.open(openParams);
        tgtAssert(file_.is_open(), "File not open");
    }
}
template<typename Element>
ElementArray<Element>::~ElementArray() {
    file_.close();
    tgt::FileSystem::deleteFile(filename_);
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

/// Impl: ElementArrayView -----------------------------------------------------------
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

/// Impl: Node -----------------------------------------------------------------------
template<typename Element>
Node<Element>::Node(Element elm, size_t left_child, size_t right_child)
    : elm_(elm)
    , left_child_(left_child)
    , right_child_(right_child)
{
}

/// Impl: NodeStorageBuilder ---------------------------------------------------------
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

/// Impl: NodeStorage ----------------------------------------------------------------
template<typename Element>
NodeStorage<Element>::NodeStorage()
    : file_()
    , filename_()
    , numNodes_(0)
{
}
template<typename Element>
NodeStorage<Element>::NodeStorage(const std::string& filename, size_t numNodes)
    : file_()
    , filename_(filename)
    , numNodes_(numNodes)
{
    if(numNodes_ > 0) {
        boost::iostreams::mapped_file_params openParams;
        openParams.path = filename_;

        file_.open(openParams);
        tgtAssert(file_.is_open(), "File not open");
    }
}

template<typename Element>
NodeStorage<Element>::NodeStorage(NodeStorage&& other)
    : file_(std::move(other.file_))
    , filename_(std::move(other.filename_))
    , numNodes_(other.numNodes_)
{
    //Somehow moving a file source only copies a reference to it, thus later closing our file in the destructor
    //REALLY make sure that the other does not happen here:
    other.file_ = boost::iostreams::mapped_file_source();
    tgtAssert(numNodes_ == 0 || file_.is_open(), "File not open");
    tgtAssert(!other.file_.is_open(), "File not open");
}
template<typename Element>
NodeStorage<Element>& NodeStorage<Element>::operator=(NodeStorage&& other) {
    if(&other != this) {
        this->~NodeStorage();
        new(this) NodeStorage(std::move(other));
        tgtAssert(numNodes_ == 0 || file_.is_open(), "File not open");
    }
    return *this;
}

template<typename Element>
NodeStorage<Element>::~NodeStorage() {
    file_.close();
    tgt::FileSystem::deleteFile(filename_);
}
template<typename Element>
const Node<Element>& NodeStorage<Element>::operator[](size_t index) const {
    tgtAssert(index < numNodes_, "Invalid index");
    const Node<Element>* data = reinterpret_cast<const Node<Element>*>(file_.data());
    return data[index];
}
template<typename Element>
size_t NodeStorage<Element>::size() const {
    return numNodes_;
}
/// Impl: SearchResult ---------------------------------------------------------------
template<typename Element>
SearchResult<Element> SearchResult<Element>::none() {
    return SearchResult (
        std::numeric_limits<typename Element::CoordType>::max(),
        nullptr
    );
}
template<typename Element>
SearchResult<Element>::SearchResult(typename Element::CoordType distSq, const Element* element)
    : distSq_(distSq)
    , element_(element)
{
}

template<typename Element>
bool SearchResult<Element>::found() const {
    bool res = element_ != nullptr;
    tgtAssert(res != (distSq_ == std::numeric_limits<typename Element::CoordType>::max()), "Valid distance for invalid element");
    return res;
}
template<typename Element>
typename Element::CoordType SearchResult<Element>::dist() const {
    return std::sqrt(distSq_);
}

template<typename Element>
void SearchResult<Element>::takeBetter(SearchResult&& other) {
    if(!found() || (other.found() && other.distSq_ < distSq_)) {
        distSq_ = other.distSq_;
        element_ = other.element_;
    }
}

/// Impl: SearchResultSet ------------------------------------------------------------

template<typename Element>
SearchResultSet<Element>::SearchResultSet()
    : distSq_(std::numeric_limits<typename Element::CoordType>::max())
    , elements_()
{
}
template<typename Element>
SearchResultSet<Element> SearchResultSet<Element>::none() {
    return SearchResultSet();
}
template<typename Element>
SearchResultSet<Element>::SearchResultSet(typename Element::CoordType distSq, const Element* element)
    : distSq_(distSq)
    , elements_()
{
    if(element) {
        elements_.push_back(element);
    }
}

template<typename Element>
bool SearchResultSet<Element>::found() const {
    bool res = !elements_.empty();
    tgtAssert(res != (distSq_ == std::numeric_limits<typename Element::CoordType>::max()), "Valid distance for invalid element");
    return res;
}
template<typename Element>
typename Element::CoordType SearchResultSet<Element>::dist() const {
    return std::sqrt(distSq_);
}

template<typename Element>
void SearchResultSet<Element>::takeBetter(SearchResultSet&& other) {
    if(other.distSq_ <= distSq_) {
        if(other.distSq_ == distSq_) {
            elements_.insert(elements_.end(), other.elements_.begin(), other.elements_.end());
        } else {
            distSq_ = other.distSq_;
            elements_ = std::move(other.elements_);
        }
    }
}

/// Impl: Tree -----------------------------------------------------------------------
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
    tgtAssert(nodes_.size() > 0, "No root in empty tree");
    return nodes_[root_];
}

template<typename Element>
SearchResult<Element> Tree<Element>::findNearest(const PosType& pos) const {
    if(nodes_.size() == 0) {
        return SearchResult<Element>::none();
    }
    SearchResult<Element> result = SearchResult<Element>::none();
    findNearestFrom(pos, 0, root_, result);
    return result;
}

template<typename Element>
SearchResultSet<Element> Tree<Element>::findAllNearest(const PosType& pos) const {
    if(nodes_.size() == 0) {
        return SearchResultSet<Element>::none();
    }
    SearchResultSet<Element> result = SearchResultSet<Element>::none();
    findNearestFrom(pos, 0, root_, result);
    return result;
}

template<typename Element>
template<typename Result>
void Tree<Element>::findNearestFrom(const PosType& pos, int depth, size_t node, Result& best_result) const {
    const Node<Element>& current = nodes_[node];
    size_t firstChild, secondChild;
    int dim = depth%3;
    typename Element::CoordType planeDist = pos.elem[dim] - current.elm_.getPos().elem[dim];
    if(planeDist < 0) {
        firstChild = current.left_child_;
        secondChild = current.right_child_;
    } else {
        firstChild = current.right_child_;
        secondChild = current.left_child_;
    }

    if(firstChild != NO_NODE_ID) {
        findNearestFrom(pos, depth+1, firstChild, best_result);
    }

    best_result.takeBetter(Result(tgt::distanceSq(pos,current.elm_.getPos()), &current.elm_));

    tgtAssert(best_result.found(), "Invalid result");
    if(planeDist*planeDist <= best_result.distSq_ && secondChild != NO_NODE_ID) {
        findNearestFrom(pos, depth+1, secondChild, best_result);
    }
}


} //namespace static_kdtree
} //namespace voreen
