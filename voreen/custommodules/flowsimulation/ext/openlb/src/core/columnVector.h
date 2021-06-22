/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender
 *
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/


#ifndef COLUMN_VECTOR_H
#define COLUMN_VECTOR_H

#include <type_traits>
#include <memory>
#include <array>
#include <stdexcept>
#include <cstring>
#include <cstdlib>

#include "serializer.h"
#include "utilities/meta.h"
#include "utilities/aliases.h"

#include "genericVector.h"
#include "scalarVector.h"


namespace olb {


template<typename T>
struct MemoryView {
  T* const          data;
  const std::size_t size;
};

/// Array with componentwise arithmetic for usage in ColumnVector
/**
 * Storage for all field data encapsulated by FieldArrayD
 **/
template<typename T>
class Column : public Serializable {
private:
  std::size_t _count;
  T* _data;
  bool _owning;

public:
  Column(std::size_t count):
    _count(count),
    _data(new T[count]{}),
    _owning(true) { }

  Column(MemoryView<T>&& mem):
    _count(mem.size),
    _data(mem.data),
    _owning(false) { }

  Column(Column<T>&& rhs):
    _count(rhs._count),
    _data(rhs._data),
    _owning(rhs._owning) {
    if (_owning) {
      rhs._owning = false;
    }
  }

  ~Column() {
    if (_owning) {
      std::free(_data);
    }
  }

  void resize(std::size_t newCount) {
    if (!_owning) {
      throw std::logic_error("Can not resize a Column view");
    }
    if (void* data = std::realloc(_data, newCount*sizeof(T))) {
      _data = static_cast<T*>(data);
      _count = newCount;
    } else {
      throw std::bad_alloc();
    }
  }

  const T& operator[](std::size_t i) const {
    return _data[i];
  }

  T& operator[](std::size_t i) {
    return _data[i];
  }

  std::size_t size() const {
    return _count;
  }

  const T* data() const {
    return _data;
  }

  T* data() {
    return _data;
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template<typename T>
std::size_t Column<T>::getNblock() const
{
  return 1;
}

template<typename T>
std::size_t Column<T>::getSerializableSize() const
{
  return _count * sizeof(T);
}

template<typename T>
bool* Column<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *data(), _count);

  return dataPtr;
}


/// Base of all ColumnVector specializations
/**
 * Used in AnyFieldArray for anonymous reference to arbirary vector sizes.
 **/
struct ColumnVectorBase : public Serializable {
  virtual ~ColumnVectorBase() { };

  /// Store all columns in a single continguous buffer, no resizing
  struct compact_storage_tag { };
  /// Let columns maintain their own local buffers, support resizing
  struct default_storage_tag { };
};

/// Vector of columns
template<typename T, unsigned D>
class ColumnVector : public GenericVector<Column<T>,D,ColumnVector<T,D>>,
                     public ColumnVectorBase {
protected:
  /// Number of rows
  std::size_t _count;
  /// Column data buffer for ColumnVectorBase::compact_storage_tag
  std::unique_ptr<T[]> _data;
  std::array<Column<T>,D> _column;

  const Column<T>* getComponentPointer(unsigned iDim) const {
    return &_column[iDim];
  }
  Column<T>* getComponentPointer(unsigned iDim) {
    return &_column[iDim];
  }

  friend GenericVector<Column<T>,D,ColumnVector<T,D>>;

public:
  ColumnVector(std::size_t count, ColumnVectorBase::compact_storage_tag):
    _count(count),
    _data(new T[D * count]{}),
    _column(utilities::meta::make_array_f<Column<T>,D>([this,count](unsigned iDim) -> MemoryView<T> {
      return MemoryView<T>{_data.get() + iDim*count, count};
    })) { }

  ColumnVector(std::size_t count, ColumnVectorBase::default_storage_tag):
    _count(count),
    _data(),
    _column(utilities::meta::make_array_f<Column<T>,D>([count](unsigned iDim) -> std::size_t {
      return count;
    })) { }

  ColumnVector(std::size_t count):
    // Resizable storage by default, compact storage is a performance optimization for population storage
    ColumnVector(count, ColumnVectorBase::default_storage_tag())
  { }

  ColumnVector(ColumnVector&& rhs):
    _count(rhs._count),
    _data(rhs._data.release()),
    _column(std::move(rhs._column)) { }

  /// Resize columns, potentially invalidates any inbound pointers
  void resize(std::size_t newCount) {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _column[iDim].resize(newCount);
    }
    _count = newCount;
  }

  /// Return pointer to Dth column
  const T* data(unsigned iDim) const {
    return _column[iDim].data();
  }
  /// Return pointer to Dth column
  T* data(unsigned iDim) {
    return _column[iDim].data();
  }

  /// Swap contents of row i and row j
  /**
   * Required to realize sorted column stores in CellIndexListD
   **/
  void swap(std::size_t i, std::size_t j) {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      std::swap(_column[iDim][i], _column[iDim][j]);
    }
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template<typename T, unsigned D>
std::size_t ColumnVector<T,D>::getNblock() const
{
  return D;
}

template<typename T, unsigned D>
std::size_t ColumnVector<T,D>::getSerializableSize() const
{
  return D * _column[0].getSerializableSize();
}

template<typename T, unsigned D>
bool* ColumnVector<T,D>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (unsigned iDim=0; iDim < D; ++iDim) {
    this->registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _column[iDim], loadingMode);
  }

  return dataPtr;
}


/// Read-only proxy for accessing a column vector entry
template<typename T, unsigned D>
class ConstVectorPtr : public ScalarVector<const T,D,ConstVectorPtr<T,D>> {
private:
  const ColumnVector<T,D>& _data;
  std::size_t _index;

  friend typename ScalarVector<const T,D,ConstVectorPtr<T,D>>::type;

protected:
  const T* getComponentPointer(unsigned iDim) const {
    return _data.data(iDim) + _index;
  }

public:
  ConstVectorPtr(const ColumnVector<T,D>& columns, std::size_t index):
    _data(columns),
    _index(index) { }

  ConstVectorPtr(ConstVectorPtr&& rhs):
    _data(rhs._data),
    _index(rhs._index) { }

  std::size_t getIndex() const {
    return _index;
  }

  void setIndex(std::size_t index) {
    _index = index;
  }

};

/// Proxy for accessing a column vector entry
template<typename T, unsigned D>
class VectorPtr : public ScalarVector<T,D,VectorPtr<T,D>> {
private:
  ColumnVector<T,D>& _data;
  std::size_t _index;

  friend typename ScalarVector<T,D,VectorPtr<T,D>>::type;

protected:
  const T* getComponentPointer(unsigned iDim) const {
    return _data.data(iDim) + _index;
  }
  T* getComponentPointer(unsigned iDim) {
    return _data.data(iDim) + _index;
  }

public:
  VectorPtr(ColumnVector<T,D>& columns, std::size_t index):
    _data(columns),
    _index(index) { }

  VectorPtr(VectorPtr&& rhs):
    _data(rhs._data),
    _index(rhs._index) { }

  VectorPtr<T,D>& operator=(const ConstVectorPtr<T,D>& rhs) {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      this->operator[](iDim) = rhs[iDim];
    }
    return *this;
  }

  std::size_t getIndex() const {
    return _index;
  }

  void setIndex(std::size_t index) {
    _index = index;
  }

};


}

#endif
