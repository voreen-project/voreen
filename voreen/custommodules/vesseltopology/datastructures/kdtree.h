/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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
#pragma once

#include "tgt/vector.h"
#include "tgt/filesystem.h"
#include "diskarraystorage.h"
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
    void partitionAtMedian(int dimension);
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
class NodeStorage { //Memory mapped version of file created from NodeStorageBuilder
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
    size_t size() const;
    void flush();
private:
    std::ofstream file_;
    std::string filename_;
    size_t numNodes_;
};

template<typename Element>
using SharedNodeStorage = DiskArray<Node<Element>>;

template<typename Element, typename Storage = NodeStorage<Element>>
class Tree;

template<typename Element>
class SharedMemoryTreeBuilder {
public:
    SharedMemoryTreeBuilder(const std::string& storagefilename);

    Tree<Element, SharedNodeStorage<Element>> buildTree(ElementArrayBuilder<Element>&& elements);
private:
    DiskArrayStorage<Node<Element>> storage_;
};

template<typename Element>
struct SearchNearestResult {
    typename Element::CoordType distSq_;
    const Element* element_;

    static SearchNearestResult none();

    bool found() const;
    typename Element::CoordType dist() const;

    void tryInsert(typename Element::CoordType distSq, const Element* element);
};

template<typename Element>
struct SearchNearestResultSet {
    typename Element::CoordType distSq_;
    std::vector<const Element*> elements_;

    SearchNearestResultSet(typename Element::CoordType maxDistSq);

    bool found() const;
    typename Element::CoordType dist() const;

    void tryInsert(typename Element::CoordType distSq, const Element* element);
private:
    SearchNearestResultSet();
};

template<typename Element>
struct SearchAllWithinResultSet {
    typename Element::CoordType distSq_;
    std::vector<const Element*> elements_;

    SearchAllWithinResultSet(typename Element::CoordType maxDistSq);

    bool found() const;

    void tryInsert(typename Element::CoordType distSq, const Element* element);
};

template<typename Element, typename Storage>
class Tree {
public:
    static_assert(std::numeric_limits<typename Element::CoordType>::is_signed, "Element pos type must be signed");
    typedef typename Node<Element>::PosType PosType;
    Tree(const std::string& storagefilename, ElementArrayBuilder<Element>&& elements);

    SearchNearestResult<Element> findNearest(const PosType& pos) const;
    SearchNearestResultSet<Element> findAllNearest(const PosType& pos, typename Element::CoordType maxDistSq = std::numeric_limits<typename Element::CoordType>::max()) const;
    SearchAllWithinResultSet<Element> findAllWithin(const PosType& pos, typename Element::CoordType maxDistSq) const;

    const Node<Element>& root() const;
private:
    template<typename Result>
    void findNearestFrom(const PosType& pos, int depth, size_t node, Result& best_result) const;

    friend SharedMemoryTreeBuilder<Element>;
    Tree(size_t root, Storage&& elements);

    Storage nodes_;
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
    tgtAssert(tmp.file_.good(), "Flush failed");
    tgtAssert(tmp.file_.is_open(), "Flush failed");
    tmp.file_.flush();
    tgtAssert(tmp.file_.good(), "Flush failed");
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
    partitionAtMedian(dimension);
    Element* center = begin_ + dist/2;
    return std::tuple<ElementArrayView<Element>, Element, ElementArrayView<Element>>{
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
        static_assert(0 <= dim && dim <= 2, "invalid dim");
        return e1.getPos().elem[dim] < e2.getPos().elem[dim];
    }
};

template<typename Element>
void ElementArrayView<Element>::partitionAtMedian(int dimension)
{
    auto middle = begin_ + size()/2;
    switch(dimension) {
        case 0:
            std::nth_element(begin_, middle, end_, SortElementsInDim<Element, 0>());
            break;
        case 1:
            std::nth_element(begin_, middle, end_, SortElementsInDim<Element, 1>());
            break;
        case 2:
            std::nth_element(begin_, middle, end_, SortElementsInDim<Element, 2>());
            break;
        //default:
            //tgtAssert(false, "Invalid dimension");
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
    tgtAssert(tmp.file_.good(), "Flush failed");
    tgtAssert(tmp.file_.is_open(), "Flush failed");
    tmp.file_.flush();
    tgtAssert(tmp.file_.good(), "Flush failed");
    tmp.file_.close();
    return NodeStorage<Element>(tmp.filename_, tmp.numNodes_);
}
template<typename Element>
size_t NodeStorageBuilder<Element>::size() const {
    return numNodes_;
}
template<typename Element>
void NodeStorageBuilder<Element>::flush() {
    tgtAssert(file_.good(), "bad file");
    file_.flush();
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

/// Impl: SharedMemoryTreeBuilder ----------------------------------------------------
const size_t NO_NODE_ID = -1;

template<typename Element, typename Builder>
static size_t buildTreeRecursively(Builder& nodes, ElementArrayView<Element>& elements, int depth) {
    if(elements.size() == 0) {
        return NO_NODE_ID;
    } else {
        int dim = depth % 3;
        auto res = elements.split(dim);
        size_t left_node = buildTreeRecursively(nodes, std::get<0>(res), depth+1);
        size_t right_node = buildTreeRecursively(nodes, std::get<2>(res), depth+1);
        return nodes.push(Node<Element>(std::get<1>(res), left_node, right_node));
    }
}

template<typename Element>
SharedMemoryTreeBuilder<Element>::SharedMemoryTreeBuilder(const std::string& storagefilename)
    : storage_(storagefilename)
{
}

template<typename Element>
Tree<Element, SharedNodeStorage<Element>> SharedMemoryTreeBuilder<Element>::buildTree(ElementArrayBuilder<Element>&& elements) {
    auto elems = std::move(elements).finalize();
    auto elem_view = elems.view();

    DiskArrayBuilder<Node<Element>> builder = storage_.build();
    size_t root = buildTreeRecursively(builder, elem_view, 0);

    return Tree<Element, SharedNodeStorage<Element>>(root, SharedNodeStorage<Element>(std::move(builder).finalize()));
}

/// Impl: Tree -----------------------------------------------------------------------
template<typename Element, typename Storage>
Tree<Element, Storage>::Tree(const std::string& storagefilename, ElementArrayBuilder<Element>&& elements) {

    auto elems = std::move(elements).finalize();
    auto elem_view = elems.view();

    NodeStorageBuilder<Element> nodeStorageBuilder(storagefilename);

    root_ = buildTreeRecursively(nodeStorageBuilder, elem_view, 0);
    nodes_ = std::move(nodeStorageBuilder).finalize();
}

template<typename Element, typename Storage>
Tree<Element, Storage>::Tree(size_t root, Storage&& elements)
    : root_(root)
    , nodes_(std::move(elements))
{
}

template<typename Element, typename Storage>
const Node<Element>& Tree<Element, Storage>::root() const {
    tgtAssert(nodes_.size() > 0, "No root in empty tree");
    return nodes_[root_];
}

template<typename Element, typename Storage>
SearchNearestResult<Element> Tree<Element, Storage>::findNearest(const PosType& pos) const {
    SearchNearestResult<Element> result = SearchNearestResult<Element>::none();
    if(nodes_.size() == 0) {
        return result;
    }
    findNearestFrom(pos, 0, root_, result);
    return result;
}

template<typename Element, typename Storage>
SearchNearestResultSet<Element> Tree<Element, Storage>::findAllNearest(const PosType& pos, typename Element::CoordType maxDistSq) const {
    SearchNearestResultSet<Element> result(maxDistSq);
    if(nodes_.size() == 0) {
        return result;
    }
    findNearestFrom(pos, 0, root_, result);
    return result;
}

template<typename Element, typename Storage>
SearchAllWithinResultSet<Element> Tree<Element, Storage>::findAllWithin(const PosType& pos, typename Element::CoordType maxDistSq) const {
    SearchAllWithinResultSet<Element> result(maxDistSq);
    if(nodes_.size() == 0) {
        return result;
    }
    findNearestFrom(pos, 0, root_, result);
    return result;
}

template<typename Element, typename Storage>
template<typename Result>
void Tree<Element, Storage>::findNearestFrom(const PosType& pos, int depth, size_t node, Result& best_result) const {
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

    best_result.tryInsert(tgt::distanceSq(pos,current.elm_.getPos()), &current.elm_);

    if(planeDist*planeDist <= best_result.distSq_ && secondChild != NO_NODE_ID) {
        findNearestFrom(pos, depth+1, secondChild, best_result);
    }
}

/// Impl: SearchNearestResult ---------------------------------------------------------------
template<typename Element>
SearchNearestResult<Element> SearchNearestResult<Element>::none() {
    return SearchNearestResult {
        std::numeric_limits<typename Element::CoordType>::max(),
        nullptr
    };
}

template<typename Element>
bool SearchNearestResult<Element>::found() const {
    return element_ != nullptr;
}
template<typename Element>
typename Element::CoordType SearchNearestResult<Element>::dist() const {
    return std::sqrt(distSq_);
}

template<typename Element>
void SearchNearestResult<Element>::tryInsert(typename Element::CoordType distSq, const Element* element) {
    if(distSq < distSq_) {
        distSq_ = distSq;
        element_ = element;
    }
}

/// Impl: SearchNearestResultSet ------------------------------------------------------------

template<typename Element>
SearchNearestResultSet<Element>::SearchNearestResultSet()
    : distSq_(std::numeric_limits<typename Element::CoordType>::max())
    , elements_()
{
}

template<typename Element>
SearchNearestResultSet<Element>::SearchNearestResultSet(typename Element::CoordType maxDistSq)
    : distSq_(maxDistSq)
    , elements_()
{
}

template<typename Element>
bool SearchNearestResultSet<Element>::found() const {
    return !elements_.empty();
}
template<typename Element>
typename Element::CoordType SearchNearestResultSet<Element>::dist() const {
    return std::sqrt(distSq_);
}

template<typename Element>
void SearchNearestResultSet<Element>::tryInsert(typename Element::CoordType distSq, const Element* element) {
    if(distSq <= distSq_) {
        if(distSq != distSq_) {
            elements_.clear();
        }
        distSq_ = distSq;
        elements_.push_back(element);
    }
}

/// Impl: SearchAllWithinResultSet ----------------------------------------------------------
template<typename Element>
SearchAllWithinResultSet<Element>::SearchAllWithinResultSet(typename Element::CoordType maxDistSq)
    : distSq_(maxDistSq)
    , elements_()
{
}

template<typename Element>
bool SearchAllWithinResultSet<Element>::found() const {
    return !elements_.empty();
}

template<typename Element>
void SearchAllWithinResultSet<Element>::tryInsert(typename Element::CoordType distSq, const Element* element) {
    if(distSq <= distSq_) {
        elements_.push_back(element);
    }
}


} //namespace static_kdtree
} //namespace voreen
