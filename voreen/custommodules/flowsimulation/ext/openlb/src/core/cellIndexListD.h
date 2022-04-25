/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender
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

#ifndef CELL_INDEX_LIST_D_H
#define CELL_INDEX_LIST_D_H

#include "fieldArrayD.h"
#include "utilities/aliases.h"
#include "serializer.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <type_traits>


namespace olb {


/// List of cell indices and associated field data
/**
 * To be used as a possible input for upcoming Operators.
 *
 * e.g. a velocity boundary operator accepts such a CellIndexListD
 *      containing the indices of boundary cells together with the
 *      velocity field to be applied.
 **/
template <typename T, typename DESCRIPTOR, typename... FIELDS>
class CellIndexListD : public Serializable {
private:
  std::size_t _count;
  std::size_t _capacity;
  MultiFieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,descriptors::CELL_ID,FIELDS...> _fields;

public:
  CellIndexListD(std::size_t capacity=64):
    _count(0),
    _capacity(capacity),
    _fields(capacity)
  { }

  CellIndexListD(std::vector<CellID>&& indices):
    _capacity(indices.size()),
    _fields(_capacity)
  {
    for (std::size_t i=0; i < indices.size(); ++i) {
      _fields.template get<descriptors::CELL_ID>()[0][i] = indices[i];
    }
  }

  std::size_t size() const {
    return _count;
  };

  /// Append cell index and allocate attached field data
  std::size_t append(std::size_t iCell);

  /// Ascending sort of cell indices
  /**
   * Attached data is moved accordingly.
   * Only useful for when processing chunks of sequential cell indices.
   **/
  void sort();

  /// Return pointer to FIELD at index
  template <typename FIELD>
  auto getFieldPointer(std::size_t index)
  {
    return _fields.template get<FIELD>().getRowPointer(index);
  }

  /// Return copy of FIELD data at index
  template <typename FIELD>
  auto getField(std::size_t index) const
  {
    return _fields.template getField<FIELD>(index);
  }
  /// Set FIELD data at index
  template <typename FIELD>
  void setField(std::size_t index, const FieldD<T,DESCRIPTOR,FIELD>& v)
  {
    return _fields.template setField<FIELD>(index, v);
  }

  template <typename FIELD>
  const auto& getFieldArray() const
  {
    return _fields.template get<FIELD>();
  }

  std::size_t getNblock() const override;
  std::size_t getSerializableSize() const override;
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template <typename T, typename DESCRIPTOR, typename... FIELDS>
std::size_t CellIndexListD<T,DESCRIPTOR,FIELDS...>::append(std::size_t iCell)
{
  if (_count == _capacity) {
    _capacity *= 2;
    _fields.resize(_capacity);
  }
  _fields.template get<descriptors::CELL_ID>()[0][_count] = iCell;
  return _count++;
}

template <typename T, typename DESCRIPTOR, typename... FIELDS>
void CellIndexListD<T,DESCRIPTOR,FIELDS...>::sort()
{
  auto& cell_id = _fields.template get<descriptors::CELL_ID>()[0];
  std::vector<std::size_t> permutation(_count);
  std::iota(permutation.begin(), permutation.end(), 0);
  std::sort(permutation.begin(), permutation.end(), [&cell_id](auto i, auto j) {
    return cell_id[i] < cell_id[j];
  });
  std::vector<bool> swapped(_count, false);
  for (std::size_t i=0; i < _count; ++i) {
    if (swapped[i]) {
      continue;
    }
    swapped[i] = true;
    std::size_t prev_j = i;
    std::size_t next_j = permutation[i];
    while (prev_j != next_j && !swapped[next_j]) {
      _fields.swap(prev_j, next_j);
      swapped[next_j] = true;
      prev_j = next_j;
      next_j = permutation[next_j];
    }
  }
}

template <typename T, typename DESCRIPTOR, typename... FIELDS>
std::size_t CellIndexListD<T,DESCRIPTOR,FIELDS...>::getNblock() const
{
  return 2
         + _fields.getNblock();
}

template <typename T, typename DESCRIPTOR, typename... FIELDS>
std::size_t CellIndexListD<T,DESCRIPTOR,FIELDS...>::getSerializableSize() const
{
  return 2 * sizeof(std::size_t)
         + _fields.getSerializableSize();
}

template <typename T, typename DESCRIPTOR, typename... FIELDS>
bool* CellIndexListD<T,DESCRIPTOR,FIELDS...>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _count);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _capacity);
  if (loadingMode && iBlock == 2) {
    _fields.resize(_capacity);
  }
  registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _fields, loadingMode);

  return dataPtr;
}

}

#endif
