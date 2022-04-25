/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef BLOCK_DATA_HH
#define BLOCK_DATA_HH

#include "blockData.h"

#include "geometry/cuboid2D.h"
#include "geometry/cuboid3D.h"

#include "functors/lattice/blockBaseF2D.h"
#include "functors/lattice/blockBaseF3D.h"

namespace olb {

template<unsigned D, typename T, typename U>
BlockData<D,T,U>::BlockData(Cuboid<T,D>& cuboid, int overlap, int size):
  BlockStructureD<D>(cuboid.getExtent(), overlap),
  _size(size),
  _communicatable(_data)
{
  for (unsigned iD=0; iD < _size; ++iD) {
    _data.emplace_back(this->getNcells());
  }
}

template<unsigned D, typename T, typename U>
BlockData<D,T,U>::BlockData(BlockStructureD<D>&& block, int size):
  BlockStructureD<D>(block),
  _size(size),
  _communicatable(_data)
{
  for (unsigned iD=0; iD < _size; ++iD) {
    _data.emplace_back(this->getNcells());
  }
}

template<unsigned D, typename T, typename U>
BlockData<D,T,U>::BlockData(BlockF<U,D>& blockF):
  BlockStructureD<D>(blockF.getBlockStructure()),
  _size(blockF.getTargetDim())
{
  for (unsigned iD=0; iD < _size; ++iD) {
    _data.emplace_back(this->getNcells());
  }
  int input[D];
  U output[_size];
  this->forCoreSpatialLocations([&](LatticeR<D> latticeR) {
    for (unsigned iD=0; iD < D; ++iD) {
      input[iD] = latticeR[iD];
    }
    blockF(output, input);
    for (unsigned iD=0; iD < _size; ++iD) {
      get(latticeR,iD) = output[iD];
    }
  });
}

template<unsigned D, typename T, typename U>
bool BlockData<D,T,U>::operator()(T output[], const int input[])
{
  const std::size_t iCell = this->getCellId(input);
  for (unsigned iD=0; iD < _size; ++iD) {
    output[iD] = _data[iD][iCell];
  }
  return true;
}

template<unsigned D, typename T, typename U>
U& BlockData<D,T,U>::get(std::size_t iCell, int iD)
{
  return _data[iD][iCell];
}

template<unsigned D, typename T, typename U>
U& BlockData<D,T,U>::get(LatticeR<D> latticeR, int iD)
{
  return get(this->getCellId(latticeR), iD);
}

template<unsigned D, typename T, typename U>
U BlockData<D,T,U>::get(LatticeR<D> latticeR, int iD) const
{
  return _data[iD][this->getCellId(latticeR)];
}

template<unsigned D, typename T, typename U>
Column<U>& BlockData<D,T,U>::getColumn(unsigned iD)
{
  return _data[iD];
}

template<unsigned D, typename T, typename U>
unsigned BlockData<D,T,U>::getSize() const
{
  return _size;
}

template<unsigned D, typename T, typename U>
std::size_t BlockData<D,T,U>::getNblock() const
{
  return _data.size() * _data[0].getNblock();
}

template<unsigned D, typename T, typename U>
std::size_t BlockData<D,T,U>::getSerializableSize() const
{
  return _data.size() * _data[0].getSerializableSize();
}

template<unsigned D, typename T, typename U>
bool* BlockData<D,T,U>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (unsigned iD=0; iD < _size; ++iD) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _data[iD], loadingMode);
  }

  return dataPtr;
}

}

#endif
