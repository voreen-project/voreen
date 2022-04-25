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

#ifndef CPU_SISD_COLUMN_H
#define CPU_SISD_COLUMN_H

#include <memory>
#include <array>
#include <stdexcept>

#include "core/platform/platform.h"
#include "core/serializer.h"
#include "communication/communicatable.h"

namespace olb {

namespace cpu {

namespace sisd {

/// Plain column for SISD CPU targets (default)
template<typename T>
class Column final : public AbstractColumn<T>
                   , public Serializable {
private:
  std::size_t _count;
  std::unique_ptr<T[]> _data;

public:
  using value_t = T;

  Column(std::size_t count):
    _count(count),
    _data(new T[count] {})
  { }

  Column(Column<T>&& rhs):
    _count(rhs._count),
    _data(rhs._data.release())
  { }

  virtual ~Column() = default;

  void resize(std::size_t count)
  {
    std::unique_ptr<T[]> data = std::unique_ptr<T[]>(new T[count] { });
    std::copy(_data.get(), _data.get() + std::min(_count, count), data.get());
    _data.swap(data);
    _count = count;
  }

  const T& operator[](std::size_t i) const override
  {
    return _data[i];
  }

  T& operator[](std::size_t i) override
  {
    return _data[i];
  }

  std::size_t size() const
  {
    return _count;
  }

  const T* data() const
  {
    return _data.get();
  }

  T* data()
  {
    return _data.get();
  }

  void setProcessingContext(ProcessingContext) { };

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};


/// Cyclic column for usage in ColumnVector
/**
 * Column type used for propagatable population storage using
 * the branching PS pattern.
 **/
template<typename T>
class CyclicColumn final : public AbstractCyclicColumn<T>
                         , public Serializable {
private:
  const std::size_t    _count;
  std::unique_ptr<T[]> _data;

  std::ptrdiff_t   _shift;
  std::size_t      _remainder;
  std::array<T*,2> _start;

public:
  using value_t = T;

  CyclicColumn(std::size_t count):
    _count(count),
    _data(new T[count] {}),
    _shift(0),
    _remainder(count)
  {
    refresh();
  }

  CyclicColumn(CyclicColumn<T>&& rhs):
    _count(rhs._count),
    _data(rhs._data.release()),
    _shift(rhs._shift),
    _remainder(rhs._remainder)
  {
    refresh();
  }

  ~CyclicColumn() = default;

  const T& operator[](std::size_t i) const override
  {
    return (i > _remainder ? _start[1] : _start[0])[i];
  }

  T& operator[](std::size_t i) override
  {
    return (i > _remainder ? _start[1] : _start[0])[i];
  }

  std::size_t size() const
  {
    return _count;
  }

  void refresh()
  {
    const std::ptrdiff_t n = size();
    T* const base = _data.get();
    if (_shift >= 0) {
      _remainder = n - _shift - 1;
      _start[0] = base + _shift;
      _start[1] = base - (n - _shift);
    }
    else {
      _remainder = -_shift - 1;
      _start[0] = base + (n + _shift);
      _start[1] = base + _shift;
    }
  }

  void rotate(std::ptrdiff_t offset)
  {
    const std::ptrdiff_t n = size();
    _shift -= offset;
    if (_shift >= n) {
      _shift -= n;
    }
    else if (_shift <= -n) {
      _shift += n;
    }
    refresh();
  }

  void setProcessingContext(ProcessingContext) { };

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;
  void postLoad() override { refresh(); }

};

}

}

/// Declare cpu::sisd::Column as the AbstractColumn implementation for CPU SISD targets
template <typename T>
struct ImplementationOf<AbstractColumn<T>,Platform::CPU_SISD> {
  using type = cpu::sisd::Column<T>;
};

/// Declare cpu::sisd::CyclicColumn as the AbstractCyclicColumn implementation for CPU SISD targets
template <typename T>
struct ImplementationOf<AbstractCyclicColumn<T>,Platform::CPU_SISD> {
  using type = cpu::sisd::CyclicColumn<T>;
};

/// Use CPU SISD as default Column
template <typename T>
using Column = cpu::sisd::Column<T>;

/// Use CPU SISD as default CyclicColumn
template <typename T>
using CyclicColumn = cpu::sisd::CyclicColumn<T>;

}

#endif

#include "column.hh"
