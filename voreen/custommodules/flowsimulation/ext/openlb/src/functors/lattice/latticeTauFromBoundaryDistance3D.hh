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

#ifndef LATTICE_TAU_FROM_BOUNDARY_DISTANCE_3D_HH
#define LATTICE_TAU_FROM_BOUNDARY_DISTANCE_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include "latticeTauFromBoundaryDistance3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T,typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysTauFromBoundaryDistance3D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticePhysTauFromBoundaryDistance3D(
  SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& sGeometry, XMLreader const& xmlReader, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter, const T p, const T T_avg, const T c_p, const T beta, const T lambda_0, const T sigma, const T p_0, const T n_0)
  : SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physTauFromBoundaryDistance";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysTauFromBoundaryDistance3D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlock(iC), sGeometry.getBlockGeometry(iC), xmlReader, this->_converter, p, T_avg, c_p, beta, lambda_0, sigma, p_0, n_0));
  }
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysTauFromBoundaryDistance3D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysTauFromBoundaryDistance3D(
  BlockLattice<T,DESCRIPTOR>& blockLattice, BlockGeometry<T,3>& blockGeometry, XMLreader const& xmlReader, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter, const T p, const T T_avg, const T c_p, const T beta, const T lambda_0, const T sigma, const T p_0, const T n_0)
  : BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,1), _blockGeometry(blockGeometry), _distanceFunctor(blockLattice, blockGeometry, xmlReader), _tmp1(lambda_0 * 287.058 * T_avg / p / c_p), _tmp2(2. * beta / (util::sqrt(2.) * M_PI * sigma * sigma * p * n_0 / p_0))
{
  this->getName() = "physTauFromBoundaryDistance";
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysTauFromBoundaryDistance3D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  T L[1] = {0.};
  _distanceFunctor(L, input);
  if ( L[0] < this->_converter.getPhysDeltaX() ) {
    L[0] = this->_converter.getPhysDeltaX();
  }

  const T alpha = _tmp1 / ( 1. + _tmp2 / L[0] );

  output[0] = alpha / this->_converter.getConversionFactorViscosity() * descriptors::invCs2<T,TDESCRIPTOR>() + 0.5;

  // std::cout << L[0] << " " << alpha << " " << output[0] << std::endl;

  return true;
}

}
#endif
