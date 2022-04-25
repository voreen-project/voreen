/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef EUKLID_NORM_3D_HH
#define EUKLID_NORM_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "euklidNorm3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperEuklidNorm3D<T, DESCRIPTOR>::SuperEuklidNorm3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), 1), _f(f)
{
  this->getName() = "EuklidNorm(" + _f.getName() + ")";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockEuklidNorm3D<T, DESCRIPTOR>(f.getBlockF(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
BlockEuklidNorm3D<T, DESCRIPTOR>::BlockEuklidNorm3D(BlockF3D<T>& f)
  : BlockF3D<T>(f.getBlockStructure(), 1),
    _f(f)
{
  this->getName() = "EuklidNorm(" + f.getName() + ")";
}

template<typename T, typename DESCRIPTOR>
bool BlockEuklidNorm3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();
  T data[_f.getTargetDim()];
  _f(data,input);
  for (int i = 0; i < _f.getTargetDim(); ++i) {
    output[0] += data[i] * data[i];
  }
  output[0] = util::sqrt(output[0]);
  return true;
}

}
#endif
