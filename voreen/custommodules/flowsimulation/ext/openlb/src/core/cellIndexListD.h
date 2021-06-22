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
class CellIndexListD {
private:
  std::size_t _count;
  std::size_t _capacity;
  MultiFieldArrayD<T,DESCRIPTOR,descriptors::CELL_ID,FIELDS...> _fields;

public:
  CellIndexListD(std::size_t capacity=64):
    _count(0),
    _capacity(capacity),
    _fields(capacity)
  { }

  CellIndexListD(std::vector<std::size_t>&& indices):
    _capacity(indices.size()),
    _fields(_capacity)
  {
    for (std::size_t i=0; i < indices.size(); ++i) {
      _fields.template get<descriptors::CELL_ID>()[0][i] = indices[i];
    }
  }

  /// Append cell index and allocate attached field data
  void append(std::size_t iCell);

  template <typename FIELD>
  FieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer(std::size_t index) {
    return _fields.template get<FIELD>().getFieldPointer(index);
  }

  /// Ascending sort of cell indices
  /**
   * Attached data is moved accordingly.
   * Only useful for when processing chunks of sequential cell indices. 
   **/
  void sort();

  template <typename O>
  void forEach(BlockLattice<T,DESCRIPTOR>& block, O op);

};

template <typename T, typename DESCRIPTOR, typename... FIELDS>
void CellIndexListD<T,DESCRIPTOR,FIELDS...>::append(std::size_t iCell)
{
  if (_count == _capacity) {
    _capacity *= 2;
    _fields.resize(_capacity);
  }
  _fields.template get<descriptors::CELL_ID>()[0][_count] = iCell;
  _count += 1;
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
template <typename O>
void CellIndexListD<T,DESCRIPTOR,FIELDS...>::forEach(BlockLattice<T,DESCRIPTOR>& block, O op)
{
  auto cell = block.get(0ul);
  auto& indices = _fields.template get<descriptors::CELL_ID>()[0];
  for (std::size_t i=0; i < _count; ++i) {
    cell.setCellId(indices[i]);
    op(cell, getFieldPointer<FIELDS>(i)...);
  }
}
 
}

#endif
