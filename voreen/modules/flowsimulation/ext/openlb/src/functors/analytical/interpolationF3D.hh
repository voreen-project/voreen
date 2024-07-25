/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Fabian Klemens, Benjamin Förster, Marie-Luise Maier,
 *  Adrian Kummerlönder
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

#ifndef INTERPOLATION_F_3D_HH
#define INTERPOLATION_F_3D_HH

#include <algorithm>

#include "interpolationF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity

namespace olb {


/// trilinear interpolation for rectangular lattice with dimensions delta[i];
/// if the cuboid is a plane (e.g. nZ==1) this functor will convert to bilinear interpolation and to linear interpolation for a "line cuboid" (e.g. nY==nZ==1)
template <typename T, typename W>
SpecialAnalyticalFfromBlockF3D<T,W>::SpecialAnalyticalFfromBlockF3D(
  BlockF3D<W>& f, Cuboid3D<T>& cuboid,
  Vector<T,3> delta, T scale)
  : AnalyticalF3D<T,W>(f.getTargetDim()), _f(f), _cuboid(cuboid), _delta(delta), _scale(scale)
{
  this->getName() = "fromBlockF";
}


template <typename T, typename W>
bool SpecialAnalyticalFfromBlockF3D<T,W>::operator()(W output[], const T physC[])
{
  Vector<T,3> origin = _cuboid.getOrigin();

  // scale physC in all 3 dimensions
  Vector<T,3> physCv;
  for (int i=0; i<3; i++) {
    physCv[i] = origin[i] + (physC[i] - origin[i]) * ( _cuboid.getDeltaR() / _delta[i] );
  }

  int latticeR[3];
  for (int i=0; i<3; i++) {
    latticeR[i] = util::max((int)util::floor( (physCv[i] - origin[i])/
                            _cuboid.getDeltaR()), 0);
  }
  Vector<T,3> physRiC;
  Vector<W,3> d, e;
  W output_tmp[3];
  Vector<T,3> latticeRv;

  for (int i=0; i<3; i++) {
    latticeRv[i] = (T) latticeR[i];
  }
  physRiC = origin + latticeRv * _cuboid.getDeltaR();
  T dr = 1. / _cuboid.getDeltaR();

  // compute weights
  d = (physCv - physRiC) * dr;
  e = 1. - d;

  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = W();
    output_tmp[iD] = W();
  }

  //0=1=2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*e[1]*e[2];
  }

  if (_cuboid.getNy() != 1) {
    latticeR[1]++;
  }

  //0=1+2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*d[1]*e[2];
  }

  if (_cuboid.getNx() != 1) {
    latticeR[0]++;
  }
  if (_cuboid.getNy() != 1) {
    latticeR[1]--;
  }
  //0+1=2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*e[1]*e[2];
  }

  if (_cuboid.getNy() != 1) {
    latticeR[1]++;
  }
  //0+1+2=
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*d[1]*e[2];
  }

  if (_cuboid.getNx() != 1) {
    latticeR[0]--;
  }
  if (_cuboid.getNy() != 1) {
    latticeR[1]--;
  }
  if (_cuboid.getNz() != 1) {
    latticeR[2]++;
  }
  //0=1=2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*e[1]*d[2];
  }

  if (_cuboid.getNy() != 1) {
    latticeR[1]++;
  }
  //0=1+2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * e[0]*d[1]*d[2];
  }

  if (_cuboid.getNx() != 1) {
    latticeR[0]++;
  }
  if (_cuboid.getNy() != 1) {
    latticeR[1]--;
  }
  //0+1=2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*e[1]*d[2];
  }

  if (_cuboid.getNy() != 1) {
    latticeR[1]++;
  }
  //0+1+2+
  _f(output_tmp, latticeR);
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] += output_tmp[iD] * d[0]*d[1]*d[2];
  }

  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] *= _scale;
  }

  return true;
}


template <typename T, typename W>
AnalyticalFfromBlockF3D<T,W>::AnalyticalFfromBlockF3D(
  BlockF3D<W>& f, Cuboid3D<T>& cuboid)
  : AnalyticalF3D<T,W>(f.getTargetDim()),
    _f(f), _cuboid(cuboid)
{
  this->getName() = "fromBlockF";
}

/// trilinear interpolation on cubic lattice
template <typename T, typename W>
bool AnalyticalFfromBlockF3D<T,W>::operator()(W output[], const T physC[])
{
  int latticeC[3];
  int latticeR[3];
  _cuboid.getFloorLatticeR(latticeR, physC);

  auto& block = _f.getBlockStructure();
  auto padding = std::min(1, block.getPadding());

  if (LatticeR<3>(latticeR) >= -padding && LatticeR<3>(latticeR) < block.getExtent()+padding-1) {
    const int& locX = latticeR[0];
    const int& locY = latticeR[1];
    const int& locZ = latticeR[2];

    Vector<T,3> physRiC;
    Vector<T,3> physCv(physC);
    _cuboid.getPhysR(physRiC.data(), {locX, locY, locZ});

    // compute weights
    Vector<W,3> d = (physCv - physRiC) * (1. / _cuboid.getDeltaR());
    Vector<W,3> e = 1. - d;

    W output_tmp[_f.getTargetDim()];
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output_tmp[iD] = W();
    }

    latticeC[0] = locX;
    latticeC[1] = locY;
    latticeC[2] = locZ;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * e[0] * e[1] * e[2];
      output_tmp[iD] = W();
    }

    latticeC[0] = locX;
    latticeC[1] = locY + 1;
    latticeC[2] = locZ;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * e[0] * d[1] * e[2];
      output_tmp[iD] = W();
    }

    latticeC[0] = locX + 1;
    latticeC[1] = locY;
    latticeC[2] = locZ;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * d[0] * e[1] * e[2];
      output_tmp[iD] = W();
    }

    latticeC[0] = locX + 1;
    latticeC[1] = locY + 1;
    latticeC[2] = locZ;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * d[0] * d[1] * e[2];
      output_tmp[iD] = W();
    }

    latticeC[0] = locX;
    latticeC[1] = locY;
    latticeC[2] = locZ + 1;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * e[0] * e[1] * d[2];
      output_tmp[iD] = W();
    }

    latticeC[0] = locX;
    latticeC[1] = locY + 1;
    latticeC[2] = locZ + 1;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * e[0] * d[1] * d[2];
      output_tmp[iD] = W();
    }

    latticeC[0] = locX + 1;
    latticeC[1] = locY;
    latticeC[2] = locZ + 1;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * d[0] * e[1] * d[2];
      output_tmp[iD] = W();
    }

    latticeC[0] = locX + 1;
    latticeC[1] = locY + 1;
    latticeC[2] = locZ + 1;
    _f(output_tmp,latticeC);
    for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
      output[iD] += output_tmp[iD] * d[0] * d[1] * d[2];
      output_tmp[iD] = W();
    }

    return true;
  }
  else {
    return false;
  }
}


template <typename T, typename W>
AnalyticalFfromSuperF3D<T,W>::AnalyticalFfromSuperF3D(SuperF3D<T,W>& f,
    bool communicateToAll, bool communicateOverlap)
  : AnalyticalF3D<T,W>(f.getTargetDim()),
    _communicateToAll(communicateToAll),
    _communicateOverlap(communicateOverlap),
    _f(f),
    _cuboidGeometry(f.getSuperStructure().getCuboidGeometry())
{
  this->getName() = "fromSuperF";

  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();
  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new AnalyticalFfromBlockF3D<T>(_f.getBlockF(iC),
                                     _cuboidGeometry.get(load.glob(iC)))
    );
  }
}

template <typename T, typename W>
bool AnalyticalFfromSuperF3D<T,W>::operator()(W output[], const T physC[])
{
  for (int iD = 0; iD < _f.getTargetDim(); ++iD) {
    output[iD] = W();
  }

  int latticeR[4];
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
int AnalyticalFfromSuperF3D<T,W>::getBlockFSize() const
{
  OLB_ASSERT(_blockF.size() < INT32_MAX,
             "it is safe to cast std::size_t to int");
  return _blockF.size();
}

template <typename T, typename W>
AnalyticalFfromBlockF3D<T,W>& AnalyticalFfromSuperF3D<T,W>::getBlockF(int iCloc)
{
  OLB_ASSERT(iCloc < int(_blockF.size()) && iCloc >= 0,
             "block functor index within bounds");
  return *(_blockF[iCloc]);
}


}  // end namespace olb

#endif
