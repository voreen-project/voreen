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

#ifndef CPU_SIMD_MASK_H_
#define CPU_SIMD_MASK_H_

#include "pack.h"
#include "column.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
class ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SIMD> final : public Serializable {
private:
  const std::size_t _size;
  std::size_t _weight;
  cpu::simd::Column<bool> _mask;

  const std::size_t _serializedSize;
  std::unique_ptr<typename cpu::simd::Mask<T>::storage_t[]> _serialized;

  bool _modified;

public:
  ConcreteBlockMask(std::size_t size):
    _size(size + cpu::simd::Pack<T>::size),
    _weight(0),
    _mask(_size),
    _serializedSize((_size + cpu::simd::Pack<T>::size - 1) / cpu::simd::Mask<T>::storage_size + 1),
    _serialized(new typename cpu::simd::Mask<T>::storage_t[_serializedSize] { }),
    _modified(true)
  { }

  bool operator[](std::size_t i) const {
    return _mask[i];
  }

  std::size_t weight() const {
    return _weight;
  }

  void set(std::size_t i, bool active) {
           if ( _mask[i] && !active) {
      _weight -= 1;
    } else if (!_mask[i] &&  active) {
      _weight += 1;
    }
    _mask[i] = active;
    _modified = true;
  }

  auto* raw()
  {
    if constexpr (cpu::simd::Mask<T>::storage_size == 1) {
      return _mask.data(); // Transparently convert bools to std::uint64_t for AVX2
    } else {
      return _serialized.get(); // Use serialized bit-mask for AVX-512
    }
  }

  void setProcessingContext(ProcessingContext) {
    if (_modified) {
      for (std::size_t i=0; i < _size; i += cpu::simd::Mask<T>::storage_size) {
        _serialized[i / cpu::simd::Mask<T>::storage_size] = cpu::simd::Mask<T>::encode(_mask.data() + i);
      }
      _modified = false;
    }
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template<typename T, typename DESCRIPTOR>
std::size_t ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SIMD>::getNblock() const
{
  return 1 + _mask.getNblock();
}

template<typename T, typename DESCRIPTOR>
std::size_t ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SIMD>::getSerializableSize() const
{
  return sizeof(_weight) + _mask.getSerializableSize();
}

template<typename T, typename DESCRIPTOR>
bool* ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SIMD>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _weight);
  registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _mask, loadingMode);

  return dataPtr;
}


}

#endif
