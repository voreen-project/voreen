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
#include <vector>

#include "serializer.h"
#include "meta.h"
#include "communication/communicatable.h"

#include "genericVector.h"
#include "scalarVector.h"

namespace olb {


/// Base of all ColumnVector specializations
/**
 * Used in AnyFieldArray for anonymous reference to arbirary vector sizes.
 **/
struct ColumnVectorBase : public Serializable {
  virtual ~ColumnVectorBase() { };
};

/// Vector of columns
/**
 * Base storage of FieldArrayD and as such all field data including populations
 **/
template<typename COLUMN, unsigned D>
class ColumnVector : public ColumnVectorBase {
protected:
  /// Number of rows
  std::size_t _count;
  std::array<COLUMN,D> _column;

public:
  class const_ptr;
  class ptr;

  static constexpr unsigned d = D;

  ColumnVector(std::size_t count):
    _count(count),
    _column(meta::make_array_f<COLUMN,D>([count](unsigned iDim) -> std::size_t {
      return count;
    }))
  { }

  ColumnVector(ColumnVector&& rhs):
    _count(rhs._count),
    _column(std::move(rhs._column)) { }

  const COLUMN& operator[](unsigned iDim) const
  {
    return _column[iDim];
  }

  COLUMN& operator[](unsigned iDim)
  {
    return _column[iDim];
  }

  /// Return copy of data at index i
  auto getRow(std::size_t i) const
  {
    if constexpr (D == 1) {
      return operator[](0)[i];
    } else {
      return Vector<typename COLUMN::value_t,D>([this,i](unsigned iDim) {
        return operator[](iDim)[i];
      });
    }
	__builtin_unreachable();
  }

  auto setRow(std::size_t i, const Vector<typename COLUMN::value_t,D>& value)
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      operator[](iDim)[i] = value[iDim];
    }
  }

  ptr getRowPointer(std::size_t i)
  {
    return ptr(*this, i);
  }

  const_ptr getRowPointer(std::size_t i) const
  {
    return const_ptr(*this, i);
  }

  /// Resize columns, potentially invalidates any inbound pointers
  void resize(std::size_t newCount)
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _column[iDim].resize(newCount);
    }
    _count = newCount;
  }

  /// Swap contents of row i and row j
  /**
   * Required to realize sorted column stores in CellIndexListD
   **/
  void swap(std::size_t i, std::size_t j)
  {
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
  void postLoad() override;

};

template<typename COLUMN, unsigned D>
std::size_t ColumnVector<COLUMN,D>::getNblock() const
{
  return D * _column[0].getNblock();
}

template<typename COLUMN, unsigned D>
std::size_t ColumnVector<COLUMN,D>::getSerializableSize() const
{
  return D * _column[0].getSerializableSize();
}

template<typename COLUMN, unsigned D>
bool* ColumnVector<COLUMN,D>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (unsigned iDim=0; iDim < D; ++iDim) {
    this->registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _column[iDim], loadingMode);
  }

  return dataPtr;
}

template<typename COLUMN, unsigned D>
void ColumnVector<COLUMN,D>::postLoad()
{
  for (unsigned iD=0; iD < D; ++iD) {
    this->operator[](iD).postLoad();
  }
}


/// Read-only proxy for accessing a column vector entry
template<typename COLUMN, unsigned D>
class ColumnVector<COLUMN,D>::const_ptr : public ScalarVector<const typename COLUMN::value_t,D,const_ptr> {
private:
  const ColumnVector<COLUMN,D>& _data;
  std::size_t _index;

  friend typename ScalarVector<const typename COLUMN::value_t,D,const_ptr>::type;

protected:
  const typename COLUMN::value_t* getComponentPointer(unsigned iDim) const
  {
    return &_data[iDim][_index];
  }

public:
  const_ptr(const ColumnVector<COLUMN,D>& columns, std::size_t index):
    _data(columns),
    _index(index) { }

  const_ptr(const_ptr&& rhs):
    _data(rhs._data),
    _index(rhs._index) { }

  std::size_t getIndex() const
  {
    return _index;
  }

  void setIndex(std::size_t index)
  {
    _index = index;
  }

};

/// Proxy for accessing a column vector entry
template<typename COLUMN, unsigned D>
class ColumnVector<COLUMN,D>::ptr : public ScalarVector<typename COLUMN::value_t,D,ptr> {
private:
  ColumnVector<COLUMN,D>& _data;
  std::size_t _index;

  friend typename ScalarVector<typename COLUMN::value_t,D,ptr>::type;

protected:
  const typename COLUMN::value_t* getComponentPointer(unsigned iDim) const
  {
    return &_data[iDim][_index];
  }
  typename COLUMN::value_t* getComponentPointer(unsigned iDim)
  {
    return &_data[iDim][_index];
  }

public:
  ptr(ColumnVector<COLUMN,D>& columns, std::size_t index):
    _data(columns),
    _index(index) { }

  ptr(ptr&& rhs):
    _data(rhs._data),
    _index(rhs._index) { }

  template <typename U, typename IMPL>
  ptr& operator=(const GenericVector<U,D,IMPL>& rhs)
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      this->operator[](iDim) = rhs[iDim];
    }
    return *this;
  }

  std::size_t getIndex() const
  {
    return _index;
  }

  void setIndex(std::size_t index)
  {
    _index = index;
  }

  unsigned getSize() const
  {
    return D;
  }

};


template <typename COLUMN, unsigned D>
class ConcreteCommunicatable<ColumnVector<COLUMN,D>> final : public Communicatable {
private:
  ColumnVector<COLUMN,D>& _vector;

public:
  ConcreteCommunicatable(ColumnVector<COLUMN,D>& vector):
    _vector{vector} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const
  {
    std::size_t size = 0;
    for (unsigned iD=0; iD < D; ++iD) {
      size += ConcreteCommunicatable<COLUMN>(_vector[iD]).size(indices);
    }
    return size;
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const
  {
    std::size_t size = ConcreteCommunicatable<COLUMN>(_vector[0]).size(indices);
    std::uint8_t* curr = buffer;
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for schedule(static,1)
    #endif
    for (unsigned iD=0; iD < D; ++iD) {
      ConcreteCommunicatable<COLUMN>(_vector[iD]).serialize(indices, curr + iD*size);
    }
    return D * size;
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer)
  {
    std::size_t size = ConcreteCommunicatable<COLUMN>(_vector[0]).size(indices);
    const std::uint8_t* curr = buffer;
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for schedule(static,1)
    #endif
    for (unsigned iD=0; iD < D; ++iD) {
      ConcreteCommunicatable<COLUMN>(_vector[iD]).deserialize(indices, curr + iD*size);
    }
    return D * size;
  }

};


}

#endif
