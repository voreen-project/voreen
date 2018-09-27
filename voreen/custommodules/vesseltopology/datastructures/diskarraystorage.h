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
#include <string>
#include <vector>
#include <fstream>
#include <boost/iostreams/device/mapped_file.hpp>
#include "tgt/filesystem.h"

template<typename Element>
struct DiskArrayReverseConstIterator {
    DiskArrayReverseConstIterator(const Element* init)
        : current_(init)
    {
    }
    DiskArrayReverseConstIterator& operator++() {
        current_--;
        return *this;
    }
    bool operator==(const DiskArrayReverseConstIterator& other) {
        return current_ == other.current_;
    }
    bool operator!=(const DiskArrayReverseConstIterator& other) {
        return current_ != other.current_;
    }
    const Element& operator*() {
        return *current_;
    }
    const Element* operator->() {
        return current_;
    }
    const Element* current_;
};

template<typename Element>
class DiskArray {
    //typedef Element* iterator;
    typedef const Element* const_iterator;
    typedef DiskArrayReverseConstIterator<Element> const_reverse_iterator;
public:
    DiskArray(); //An empty disk array, not associated with any storage
    DiskArray(boost::iostreams::mapped_file* file_, size_t begin, size_t end);
    DiskArray(DiskArray&& other);
    DiskArray& operator=(DiskArray&& other);
    ~DiskArray() {}

    const Element& at(size_t index) const;
    Element& at(size_t index);
    const Element& operator[](size_t index) const;
    Element& operator[](size_t index);

    size_t size() const;
    bool empty() const;
    const Element& front() const;
    const Element& back() const;

    // May be invalidated by creating other diskarrays in the meantime!!
    const_iterator begin() const;
    const_iterator end() const;
    const_reverse_iterator rbegin() const;
    const_reverse_iterator rend() const;

    const void* storageIdentifier() const {
        return file_;
    }
private:
    boost::iostreams::mapped_file* file_;
    size_t begin_;
    size_t end_;
};

template<typename Element>
class DiskArrayBuilder;

template<typename Element>
class DiskArrayStorage {
public:
    DiskArrayStorage(const std::string& storagefilename);
    ~DiskArrayStorage();

    // Note: The store must live longer than the returned DiskArray!
    // Absolutely not threadsafe!
    DiskArray<Element> store(const std::vector<Element>& elements);
    DiskArray<Element> store(const DiskArray<Element>& elements);

    // Important! Do not use DiskArrayStorage until finalize is called on the builder!
    // (As a consequence) only one builder can be active at a time.
    DiskArrayBuilder<Element> build();

    // Only store a single Element and return the storage position with which it can be retrieved
    size_t storeElement(const Element& elm);
    size_t size() const;

    const void* identifier() const {
        return &file_;
    }
private:
    friend class DiskArrayBuilder<Element>;

    template<typename Arr>
    DiskArray<Element> store_internal(const Arr& array);
    void ensureFit(size_t numElements);

    boost::iostreams::mapped_file file_;
    size_t numElements_;
    const std::string storagefilename_;
    size_t physicalFileSize_;
};

template<typename Element>
class DiskArrayBuilder {
public:
    DiskArrayBuilder(DiskArrayStorage<Element>& storage);
    DiskArrayBuilder(DiskArrayBuilder&& other);
    size_t push(const Element&);
    DiskArray<Element> finalize() &&;

private:
    DiskArrayStorage<Element>& storage_;
    size_t begin_;
    size_t end_;
    size_t numElements_;
};


/// Impl: DiskArray ------------------------------------------------------------
template<typename Element>
DiskArray<Element>::DiskArray()
    : file_(nullptr)
    , begin_(0)
    , end_(0)
{
}

template<typename Element>
DiskArray<Element>::DiskArray(boost::iostreams::mapped_file* file, size_t begin, size_t end)
    : file_(file)
    , begin_(begin)
    , end_(end)
{
}

template<typename Element>
DiskArray<Element>::DiskArray(DiskArray&& other)
    : file_(other.file_)
    , begin_(other.begin_)
    , end_(other.end_)
{
}

template<typename Element>
DiskArray<Element>& DiskArray<Element>::operator=(DiskArray&& other)
{
    if(this != &other) {
        this->~DiskArray();
        new(this) DiskArray(std::move(other));
    }
    return *this;

}

template<typename Element>
const Element& DiskArray<Element>::at(size_t index) const {
    size_t fileIndex = index + begin_;
    tgtAssert(fileIndex < end_, "Invalid index");
    tgtAssert(file_, "No file");
    const Element* data = reinterpret_cast<const Element*>(file_->data());
    return data[fileIndex];
}

template<typename Element>
Element& DiskArray<Element>::at(size_t index) {
    size_t fileIndex = index + begin_;
    tgtAssert(fileIndex < end_, "Invalid index");
    tgtAssert(file_, "No file");
    Element* data = reinterpret_cast<Element*>(file_->data());
    return data[fileIndex];
}

template<typename Element>
const Element& DiskArray<Element>::operator[](size_t index) const {
    return at(index);
}

template<typename Element>
Element& DiskArray<Element>::operator[](size_t index) {
    return at(index);
}

template<typename Element>
size_t DiskArray<Element>::size() const {
    return end_ - begin_;
}

template<typename Element>
bool DiskArray<Element>::empty() const {
    return size() == 0;
}

template<typename Element>
const Element& DiskArray<Element>::front() const {
    tgtAssert(!empty(), "Empty array");
    return at(0);
}

template<typename Element>
const Element& DiskArray<Element>::back() const {
    tgtAssert(!empty(), "Empty array");
    return at(size() - 1);
}

template<typename Element>
typename DiskArray<Element>::const_iterator DiskArray<Element>::begin() const {
    if(empty()) {
        return nullptr;
    } else {
        return &operator[](0);
    }
}

template<typename Element>
typename DiskArray<Element>::const_iterator DiskArray<Element>::end() const {
    return begin() + size();
}

template<typename Element>
typename DiskArray<Element>::const_reverse_iterator DiskArray<Element>::rbegin() const {
    if(empty()) {
        return DiskArray<Element>::const_reverse_iterator(nullptr);
    } else {
        return DiskArray<Element>::const_reverse_iterator(&operator[](0) + size());
    }
}

template<typename Element>
typename DiskArray<Element>::const_reverse_iterator DiskArray<Element>::rend() const {
    if(empty()) {
        return DiskArray<Element>::const_reverse_iterator(nullptr);
    } else {
        return DiskArray<Element>::const_reverse_iterator(&operator[](0));
    }
}

/// Impl: DiskArrayStorage -----------------------------------------------------
template<typename Element>
DiskArrayStorage<Element>::DiskArrayStorage(const std::string& storagefilename)
    : file_()
    , numElements_(0)
    , storagefilename_(storagefilename)
    , physicalFileSize_(0)
{
}

template<typename Element>
DiskArrayStorage<Element>::~DiskArrayStorage() {
    file_.close();
    tgt::FileSystem::deleteFile(storagefilename_);
}

static void growFile(const std::string& fileName, size_t size) {
    {
        //Create file if it does not exist
        std::ofstream file(fileName, std::ios_base::binary | std::ios_base::in | std::ios_base::app);
    }
    // Now open the existing file and write to it
    std::ofstream file(fileName, std::ios_base::binary | std::ios_base::in /* Do not truncate file */);
    tgtAssert(file.good(), "Growing file failed");
    file.seekp(size);
    tgtAssert(file.good(), "Growing file failed");
    file << '\0';
    tgtAssert(file.good(), "Growing file failed");
}

static size_t nextPowerOfTwo(size_t input) {
    size_t out = 1;
    while(out < input) {
        out = out << 1;
    }
    return out;
}
template<typename Element>
void DiskArrayStorage<Element>::ensureFit(size_t numElements) {
    size_t requiredFileSize = numElements * sizeof(Element);
    if(requiredFileSize > physicalFileSize_) {
        file_.close();

        physicalFileSize_ = std::max(2*physicalFileSize_, nextPowerOfTwo(requiredFileSize));
        growFile(storagefilename_, physicalFileSize_);

        boost::iostreams::mapped_file_params openParams;
        openParams.path = storagefilename_;
        openParams.mode = std::ios::in | std::ios::out;

        file_.open(openParams);
        tgtAssert(file_.is_open(), "File not open");
    }
}

template<typename Element>
DiskArray<Element> DiskArrayStorage<Element>::store(const std::vector<Element>& elements) {
    return store_internal(elements);
}

template<typename Element>
DiskArray<Element> DiskArrayStorage<Element>::store(const DiskArray<Element>& elements) {
    return store_internal(elements);
}

template<typename Element>
template<typename Arr>
DiskArray<Element> DiskArrayStorage<Element>::store_internal(const Arr& array) {
    if(array.empty()) {
        return DiskArray<Element>(nullptr, 0, 0);
    } else {
        size_t oldNumElements = numElements_;
        numElements_ += array.size();

        ensureFit(numElements_);

        Element* data = reinterpret_cast<Element*>(file_.data());
        tgtAssert(file_.is_open(), "file not opened");
        tgtAssert(data, "Invalid data pointer");
        std::copy(&*array.begin(), &*array.end(), data + oldNumElements);

        return DiskArray<Element>(&file_, oldNumElements, numElements_);
    }
}

template<typename Element>
size_t DiskArrayStorage<Element>::storeElement(const Element& elm) {
    size_t oldNumElements = numElements_;
    numElements_ += 1;

    ensureFit(numElements_);

    Element* data = reinterpret_cast<Element*>(file_.data());
    tgtAssert(file_.is_open(), "file not opened");
    tgtAssert(data, "Invalid data pointer");
    data[oldNumElements] = elm;

    return oldNumElements;
}

template<typename Element>
size_t DiskArrayStorage<Element>::size() const {
    return numElements_;
}

template<typename Element>
DiskArrayBuilder<Element> DiskArrayStorage<Element>::build() {
    return DiskArrayBuilder<Element>(*this);
}

/// Impl: DiskArrayBuilder -----------------------------------------------------

template<typename Element>
DiskArrayBuilder<Element>::DiskArrayBuilder(DiskArrayStorage<Element>& storage)
    : storage_(storage)
    , begin_(storage.size())
    , end_(storage.size())
{
}

template<typename Element>
DiskArrayBuilder<Element>::DiskArrayBuilder(DiskArrayBuilder&& other)
    : storage_(other.storage_)
    , begin_(other.begin_)
    , end_(other.end_)
    , numElements_(other.numElements_)
{
}

template<typename Element>
size_t DiskArrayBuilder<Element>::push(const Element& elm) {
    size_t insertPos = storage_.storeElement(elm);
    (void)insertPos; //mark variable as used even in release builds...
    tgtAssert(insertPos == end_, "Invalid insert pos. Was storage used since builder was created?");
    ++end_;
    return insertPos - begin_;
}

template<typename Element>
DiskArray<Element> DiskArrayBuilder<Element>::finalize() && {
    auto tmp = std::move(*this);
    return DiskArray<Element>(&storage_.file_, begin_, end_);
}
