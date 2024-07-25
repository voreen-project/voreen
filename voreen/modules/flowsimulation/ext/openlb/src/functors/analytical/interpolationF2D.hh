/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2018 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Adrian Kummerlaender
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

#ifndef INTERPOLATION_F_2D_HH
#define INTERPOLATION_F_2D_HH

#include "interpolationF2D.h"
#include "core/superLattice2D.h"
#include "dynamics/lbm.h"

namespace olb {


template <typename T, typename W>
AnalyticalFfromBlockF2D<T,W>::AnalyticalFfromBlockF2D(
  BlockF2D<W>& f, Cuboid2D<T>& cuboid)
  : AnalyticalF2D<T,W>(f.getTargetDim()),
    _f(f), _cuboid(cuboid)
{
  this->getName() = "fromBlockF";
}

template <typename T, typename W>
bool AnalyticalFfromBlockF2D<T,W>::operator()(W output[], const T physC[])
{
  int latticeC[2];
  int latticeR[2];
  _cuboid.getFloorLatticeR(latticeR, physC);

  auto& block = _f.getBlockStructure();
  auto padding = std::min(1, block.getPadding());

  if (LatticeR<2>(latticeR) >= -padding && LatticeR<2>(latticeR) < block.getExtent()+padding-1) {
    const int& locX = latticeR[0];
    const int& locY = latticeR[1];

    Vector<T,2> physRiC;
    Vector<T,2> physCv(physC);
    _cuboid.getPhysR(physRiC.data(), {locX, locY});

    // compute weights
    Vector<W,2> d = (physCv - physRiC) * (1. / _cuboid.getDeltaR());
    Vector<W,2> e = 1. - d;

    T output_tmp[_f.getTargetDim()];
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output_tmp[iD] = T();
    }

    latticeC[0] = locX;
    latticeC[1] = locY;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * e[0] * e[1];
      output_tmp[iD] = T();
    }

    latticeC[0] = locX;
    latticeC[1] = locY + 1;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * e[0] * d[1];
      output_tmp[iD] = T();
    }

    latticeC[0] = locX + 1;
    latticeC[1] = locY;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * d[0] * e[1];
      output_tmp[iD] = T();
    }

    latticeC[0] = locX + 1;
    latticeC[1] = locY + 1;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * d[0] * d[1];
      output_tmp[iD] = T();
    }

    return true;
  }
  else {
    return false;
  }
}

template <typename T, typename W>
AnalyticalFfromSuperF2D<T,W>::AnalyticalFfromSuperF2D(SuperF2D<T>& f,
    bool communicateToAll, bool communicateOverlap)
  : AnalyticalF2D<T,W>(f.getTargetDim()),
    _communicateToAll(communicateToAll),
    _communicateOverlap(communicateOverlap),
    _f(f),
    _cuboidGeometry(_f.getSuperStructure().getCuboidGeometry())
{
  this->getName() = "fromSuperF";

  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();
  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new AnalyticalFfromBlockF2D<T>(_f.getBlockF(iC),
                                     _cuboidGeometry.get(load.glob(iC)))
    );
  }
}

template <typename T, typename W>
bool AnalyticalFfromSuperF2D<T,W>::operator() (T output[], const T physC[])
{
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = W();
  }

  int latticeR[3];
  if (!_cuboidGeometry.getLatticeR(latticeR, physC)) {
    return false;
  }

  if (_communicateOverlap) {
    _f.getSuperStructure().communicate();
  }

  int dataSize = 0;
  int dataFound = 0;

  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();
  for (int iC = 0; iC < load.size(); ++iC) {
    if (_blockF[iC]->operator()(output, physC)) {
      dataSize += _f.getTargetDim();
      ++dataFound;
    }
  }

  if (_communicateToAll) {
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(dataFound, MPI_SUM);
    singleton::mpi().reduceAndBcast(dataSize, MPI_SUM);
#endif
    dataSize /= dataFound;
#ifdef PARALLEL_MODE_MPI
    for (int iD = 0; iD < dataSize; ++iD) {
      singleton::mpi().reduceAndBcast(output[iD], MPI_SUM);
    }
#endif
    for (int iD = 0; iD < dataSize; ++iD) {
      output[iD]/=dataFound;
    }
  }
  else {
    if (dataFound!=0) {
      dataSize /= dataFound;
      for (int iD = 0; iD < dataSize; ++iD) {
        output[iD]/=dataFound;
      }
    }
  }

  if (dataFound>0) {
    return true;
  }
  return false;
}

template <typename T, typename W>
int AnalyticalFfromSuperF2D<T,W>::getBlockFSize() const
{
  OLB_ASSERT(_blockF.size() < UINT32_MAX,
             "it is safe to cast std::size_t to int");
  return _blockF.size();
}

template <typename T, typename W>
AnalyticalFfromBlockF2D<T,W>& AnalyticalFfromSuperF2D<T,W>::getBlockF(int iCloc)
{
  OLB_ASSERT(size_t(iCloc) < _blockF.size() && iCloc >= 0,
             "block functor index within bounds");
  return *(_blockF[iCloc]);
}


} // end namespace olb

#endif
