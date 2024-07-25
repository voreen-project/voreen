/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Benjamin Förster
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

#ifndef SUPER_DATA_HH
#define SUPER_DATA_HH

#include "superData.h"

#include "geometry/cuboidGeometry2D.h"
#include "geometry/cuboidGeometry3D.h"

#include "functors/lattice/superBaseF2D.h"
#include "functors/lattice/superBaseF3D.h"

namespace olb {


template<unsigned D, typename T, typename U>
SuperData<D,T,U>::SuperData(CuboidGeometry<T,D>& cuboidGeometry,
                            LoadBalancer<T>& loadBalancer,
                            int overlap, int size)
  : SuperStructure<T,D>(cuboidGeometry, loadBalancer, overlap),
    _size(size)
{
  auto& load = this->getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    _block.emplace_back(
      new BlockData<D,T,U>(cuboidGeometry.get(load.glob(iC)), overlap, size));
  }

  if (overlap >= 1) {
    _communicator = std::make_unique<SuperCommunicator<T,SuperData>>(*this);
    _communicator->template requestField<typename BlockData<D,T,U>::DUMMY_FIELD>();
    _communicator->requestOverlap(overlap);
    _communicator->exchangeRequests();
  }
}

template<unsigned D, typename T, typename U>
SuperData<D,T,U>::SuperData(SuperF<D,T,U>& rhs)
  : SuperData(rhs.getSuperStructure().getCuboidGeometry(),
              rhs.getSuperStructure().getLoadBalancer(),
              rhs.getSuperStructure().getOverlap(),
              rhs.getTargetDim())
{
  auto& load = this->getLoadBalancer();

  int input[D+1];
  U output[rhs.getTargetDim()];

  for (int iC=0; iC < load.size(); ++iC) {
    auto& block = getBlock(iC);
    input[0] = load.glob(iC);
    block.forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      for (unsigned iD=0; iD < D; ++iD) {
        input[1+iD] = latticeR[iD];
      }
      rhs(output, input);
      for (unsigned iD=0; iD < rhs.getTargetDim(); ++iD) {
        block.get(latticeR,iD) = output[iD];
      }
    });
  }
}

template<unsigned D, typename T, typename U>
const BlockData<D,T,U>& SuperData<D,T,U>::getBlock(int iC) const
{
  return *_block[iC];
}

template<unsigned D, typename T, typename U>
BlockData<D,T,U>& SuperData<D,T,U>::getBlock(int iC)
{
  return *_block[iC];
}

template<unsigned D, typename T, typename U>
template <typename BLOCK>
BLOCK& SuperData<D,T,U>::getBlock(int iC)
{
  return *_block[iC];
}

template<unsigned D, typename T, typename U>
template <typename BLOCK>
const BLOCK& SuperData<D,T,U>::getBlock(int iC) const
{
  return *_block[iC];
}

template<unsigned D, typename T, typename U>
void SuperData<D,T,U>::communicate()
{
  if (_communicator) {
    _communicator->communicate();
  }
}

template<unsigned D, typename T, typename U>
int SuperData<D,T,U>::getDataSize() const
{
  return _size;
}

template<unsigned D, typename T, typename U>
int SuperData<D,T,U>::getDataTypeSize() const
{
  return sizeof(U);
}

}

#endif
