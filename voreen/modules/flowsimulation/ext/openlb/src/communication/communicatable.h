/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef COMMUNICATABLE_H
#define COMMUNICATABLE_H

#include <vector>
#include <cstdint>
#include <typeindex>

#include "core/blockStructure.h"

namespace olb {

template <typename T>
class ConstSpan {
private:
  const T* const _base;
  const std::size_t _size;

public:
  ConstSpan(const T* base, std::size_t size):
    _base(base), _size(size) { }

  ConstSpan(const std::vector<T>& v):
    ConstSpan{v.data(), v.size()} { }

  const T* data() const {
    return _base;
  }

  T operator[](std::size_t i) const {
    return _base[i];
  }

  const T* begin() const {
    return _base;
  }
  const T* end() const {
    return _base + _size;
  }

  std::size_t size() const {
    return _size;
  }

};

struct Communicatable {
  virtual ~Communicatable() { };

  virtual std::size_t size(ConstSpan<CellID> indices) const = 0;
  /// Serialize data at locations `indices` to `buffer`
  /**
   * Used for preparing inter-block communication messages
   * \returns size of serialized data in bytes
   **/
  virtual std::size_t serialize(ConstSpan<CellID> indices,
                                std::uint8_t* buffer) const = 0;
  /// Deserialize data at locations `indices` to `buffer`
  /**
   * Used for applying inter-block communication messages
   * \returns size of serialized data in bytes
   **/
  virtual std::size_t deserialize(ConstSpan<CellID> indices,
                                  const std::uint8_t* buffer) = 0;
};

template <typename COMMUNICATEE>
class ConcreteCommunicatable final : public Communicatable {
private:
  COMMUNICATEE& _communicatee;

public:
  ConcreteCommunicatable(COMMUNICATEE& communicatee):
    _communicatee{communicatee} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const override
  {
    return indices.size() * sizeof(typename COMMUNICATEE::value_type);
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const override
  {
    if (meta::is_aligned<typename COMMUNICATEE::value_type>(buffer)) {
      auto* target = reinterpret_cast<typename COMMUNICATEE::value_type*>(buffer);
      for (CellID index : indices) {
        *(target++) = _communicatee[index];
      }
    } else {
      std::uint8_t* target = buffer;
      for (CellID index : indices) {
        std::memcpy(target,
                    reinterpret_cast<const void*>(&_communicatee[index]),
                    sizeof(typename COMMUNICATEE::value_type));
        target += sizeof(typename COMMUNICATEE::value_type);
      }
    }
    return indices.size() * sizeof(typename COMMUNICATEE::value_type);
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer) override
  {
    if (meta::is_aligned<typename COMMUNICATEE::value_type>(buffer)) {
      const auto* source = reinterpret_cast<const typename COMMUNICATEE::value_type*>(buffer);
      for (CellID index : indices) {
        _communicatee[index] = *(source++);
      }
    } else {
      const std::uint8_t* source = buffer;
      for (CellID index : indices) {
        std::memcpy(reinterpret_cast<void*>(&_communicatee[index]),
                    source,
                    sizeof(typename COMMUNICATEE::value_type));
        source += sizeof(typename COMMUNICATEE::value_type);
      }
    }
    return indices.size() * sizeof(typename COMMUNICATEE::value_type);
  }

};

template <typename COLUMN>
class ConcreteCommunicatable<std::vector<COLUMN>> final : public Communicatable {
private:
  std::vector<COLUMN>& _vector;

public:
  ConcreteCommunicatable(std::vector<COLUMN>& vector):
    _vector{vector} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const
  {
    std::size_t size = 0;
    for (unsigned iD=0; iD < _vector.size(); ++iD) {
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
    for (unsigned iD=0; iD < _vector.size(); ++iD) {
      ConcreteCommunicatable<COLUMN>(_vector[iD]).serialize(indices, curr + iD*size);
    }
    return _vector.size() * size;
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
    for (unsigned iD=0; iD < _vector.size(); ++iD) {
      ConcreteCommunicatable<COLUMN>(_vector[iD]).deserialize(indices, curr + iD*size);
    }
    return _vector.size() * size;
  }

};

template <typename COMMUNICATEE>
class MultiConcreteCommunicatable final : public Communicatable {
private:
  COMMUNICATEE& _communicatee;
  const std::vector<std::type_index>& _fields;

public:
  MultiConcreteCommunicatable(COMMUNICATEE& communicatee,
                              const std::vector<std::type_index>& fields):
    _communicatee{communicatee},
    _fields{fields} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const override
  {
    std::size_t size = 0;
    for (auto& field : _fields) {
      size += _communicatee.getCommunicatable(field).size(indices);
    }
    return size;
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const override
  {
    std::uint8_t* curr = buffer;
    for (auto& field : _fields) {
      curr += _communicatee.getCommunicatable(field).serialize(indices, curr);
    }
    return curr - buffer;
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer) override
  {
    const std::uint8_t* curr = buffer;
    for (auto& field : _fields) {
      curr += _communicatee.getCommunicatable(field).deserialize(indices, curr);
    }
    return curr - buffer;
  }
};



}

#endif
