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
class DiskArray {
    //typedef Element* iterator;
    typedef const Element* const_iterator;
public:
    DiskArray(const boost::iostreams::mapped_file_source& file_, size_t begin, size_t end);
    DiskArray(const DiskArray& other);
    ~DiskArray() {}

    const Element& operator[](size_t index) const;
    size_t size() const;
    bool empty() const;
    const Element& front() const;
    const Element& back() const;

    // May be invalidated by creating other diskarrays in the meantime!!
    const_iterator begin() const;
    const_iterator end() const;

private:
    const boost::iostreams::mapped_file_source& file_;
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

private:
    friend class DiskArrayBuilder<Element>;
    size_t storeElement(const Element& elm);

    DiskArray<Element> store(const Element* elements, size_t len);
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
    void push(const Element&);
    DiskArray<Element> finalize() &&;

private:
    DiskArrayStorage<Element>& storage_;
    size_t begin_;
    size_t end_;
    size_t numElements_;
};


/// Impl: DiskArray ------------------------------------------------------------
template<typename Element>
DiskArray<Element>::DiskArray(const boost::iostreams::mapped_file_source& file, size_t begin, size_t end)
    : file_(file)
    , begin_(begin)
    , end_(end)
{
}

template<typename Element>
DiskArray<Element>::DiskArray(const DiskArray& other)
    : file_(other.file_)
    , begin_(other.begin_)
    , end_(other.end_)
{
}

template<typename Element>
const Element& DiskArray<Element>::operator[](size_t index) const {
    size_t fileIndex = index + begin_;
    tgtAssert(fileIndex < end_, "Invalid index");
    const Element* data = reinterpret_cast<const Element*>(file_.data());
    return data[fileIndex];
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
    return operator[](0);
}

template<typename Element>
const Element& DiskArray<Element>::back() const {
    tgtAssert(!empty(), "Empty array");
    return operator[](size() - 1);
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
    if(elements.empty()) {
        return store(nullptr, 0);
    } else {
        return store(&*elements.begin(), elements.size());
    }
}

template<typename Element>
DiskArray<Element> DiskArrayStorage<Element>::store(const DiskArray<Element>& elements) {
    if(elements.empty()) {
        return store(nullptr, 0);
    } else {
        return store(elements.begin(), elements.size());
    }
}

template<typename Element>
DiskArray<Element> DiskArrayStorage<Element>::store(const Element* elements, size_t len) {
    size_t oldNumElements = numElements_;
    numElements_ += len;

    ensureFit(numElements_);

    Element* data = reinterpret_cast<Element*>(file_.data());
    tgtAssert(file_.is_open(), "file not opened");
    tgtAssert(data, "Invalid data pointer");
    std::copy(elements, elements + len, data + oldNumElements);

    return DiskArray<Element>(file_, oldNumElements, numElements_);
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
DiskArrayBuilder<Element> DiskArrayStorage<Element>::build() {
    return DiskArrayBuilder<Element>(*this);
}

/// Impl: DiskArrayBuilder -----------------------------------------------------

template<typename Element>
DiskArrayBuilder<Element>::DiskArrayBuilder(DiskArrayStorage<Element>& storage)
    : storage_(storage)
    , begin_(storage.numElements_)
    , end_(storage.numElements_)
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
void DiskArrayBuilder<Element>::push(const Element& elm) {
    size_t insertPos = storage_.storeElement(elm);
    (void)insertPos; //mark variable as used even in release builds...
    tgtAssert(insertPos == end_, "Invalid insert pos. Was storage used since builder was created?");
    ++end_;
}

template<typename Element>
DiskArray<Element> DiskArrayBuilder<Element>::finalize() && {
    auto tmp = std::move(*this);
    return DiskArray<Element>(storage_.file_, begin_, end_);
}
