/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef EUKLID_NORM_2D_HH
#define EUKLID_NORM_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "euklidNorm2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T,typename DESCRIPTOR>
SuperEuklidNorm2D<T,DESCRIPTOR>::SuperEuklidNorm2D(
  SuperLatticeF2D<T,DESCRIPTOR>& f)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice(), 1), _f(f)
{
  this->getName() = "l2(" + _f.getName() + ")";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockEuklidNorm2D<T,DESCRIPTOR>(f.getBlockF(iC)));
  }
}

template <typename T, typename DESCRIPTOR>
BlockEuklidNorm2D<T,DESCRIPTOR>::BlockEuklidNorm2D(BlockF2D<T>& f)
  : BlockF2D<T>(f.getBlockStructure(),1), _f(f)
{
  this->getName() = "l2("+f.getName()+")";
}

template <typename T, typename DESCRIPTOR>
bool BlockEuklidNorm2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = T();  // flash output, since this methods adds values.
  T data[_f.getTargetDim()];
  _f(data,input);
  for ( int i = 0; i < _f.getTargetDim(); ++i) {
    output[0] += data[i]*data[i];
  }
  output[0] = util::sqrt(output[0]);
  return true;
}

}
#endif
