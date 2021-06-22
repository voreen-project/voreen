/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef BLOCK_STATIC_POPULATION_D_2D_HH
#define BLOCK_STATIC_POPULATION_D_2D_HH

#include "blockStaticPopulationD2D.h"

#include "dynamics/latticeDescriptors.h"

#include <algorithm>

namespace olb {

template<typename T, typename DESCRIPTOR>
BlockStaticPopulationD2D<T,DESCRIPTOR>::BlockStaticPopulationD2D(int nX, int nY):
  BlockStructure2D(nX, nY),
  _size(this->getNcells()),
  _data(_size, ColumnVectorBase::compact_storage_tag())
{
  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    _remainder[iPop] = _size;
    _shift[iPop] = 0;

    _population[iPop][0] = _data.data(iPop);
    _population[iPop][1] = _data.data(iPop) - _size;
  }
}

template<typename T, typename DESCRIPTOR>
T* BlockStaticPopulationD2D<T,DESCRIPTOR>::getPopulationPointer(unsigned iPop, std::size_t iCell) {
  return (iCell > _remainder[iPop] ? _population[iPop][1] : _population[iPop][0]) + iCell;
}

template<typename T, typename DESCRIPTOR>
const T* BlockStaticPopulationD2D<T,DESCRIPTOR>::getPopulationPointer(unsigned iPop, std::size_t iCell) const {
  return (iCell > _remainder[iPop] ? _population[iPop][1] : _population[iPop][0]) + iCell;
}

template<typename T, typename DESCRIPTOR>
void BlockStaticPopulationD2D<T,DESCRIPTOR>::shift()
{
  const std::ptrdiff_t nCells = this->getNcells();

  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    _shift[iPop] -= this->getNeighborDistance(
      descriptors::c<DESCRIPTOR>(iPop,0),
      descriptors::c<DESCRIPTOR>(iPop,1)
    );

    if (_shift[iPop] >= nCells) {
      _shift[iPop] -= nCells;
    } else if (_shift[iPop] <= -nCells) {
      _shift[iPop] += nCells;
    }
  }

  refreshControlStructure();
}

template<typename T, typename DESCRIPTOR>
void BlockStaticPopulationD2D<T,DESCRIPTOR>::refreshControlStructure()
{
  const std::ptrdiff_t nCells = this->getNcells();

  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    T* const base = _data.data(iPop);

    if (_shift[iPop] >= 0) {
      _remainder[iPop] = nCells - _shift[iPop] - 1;
      _population[iPop][0] = base + _shift[iPop];
      _population[iPop][1] = base - (nCells - _shift[iPop]);
    } else {
      _remainder[iPop] = -_shift[iPop] - 1;
      _population[iPop][0] = base + (nCells + _shift[iPop]);
      _population[iPop][1] = base + _shift[iPop];
    }
  }
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockStaticPopulationD2D<T,DESCRIPTOR>::getNblock() const
{
  return 5 + _data.getNblock();
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockStaticPopulationD2D<T,DESCRIPTOR>::getSerializableSize() const
{
  return sizeof(decltype(_size))
       + sizeof(decltype(this->_nx))
       + sizeof(decltype(this->_ny))
       + DESCRIPTOR::q * sizeof(std::ptrdiff_t)
       + DESCRIPTOR::q * sizeof(std::size_t)
       + _data.getSerializableSize();
}

template<typename T, typename DESCRIPTOR>
bool* BlockStaticPopulationD2D<T,DESCRIPTOR>::getBlock(
  std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _size);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_nx);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_ny);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *_shift,     DESCRIPTOR::q);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *_remainder, DESCRIPTOR::q);
  registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _data, loadingMode);

  return dataPtr;
}


}

#endif
