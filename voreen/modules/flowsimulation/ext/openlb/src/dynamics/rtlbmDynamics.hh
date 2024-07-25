/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017-2019 Albert Mink, Christopher McHardy
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

/** \file
 * A collection of radiative transport dynamics classes -- generic implementation.
 */

#ifndef RTLBM_DYNAMICS_HH
#define RTLBM_DYNAMICS_HH

#include "rtlbmDynamics.h"
#include "rtlbmDescriptors.h"
#include "lbm.h"

namespace olb {



//==================================================================//
//============= BGK Model for Advection diffusion anisotropic ===//
//==================================================================//

template<typename T, typename DESCRIPTOR, typename MOMENTA>
RTLBMdynamicsMcHardy<T, DESCRIPTOR, MOMENTA>::RTLBMdynamicsMcHardy
(T latticeAbsorption, T latticeScattering, std::array<std::array<T,DESCRIPTOR::q>, DESCRIPTOR::q>& anisoMatrix)
  : legacy::BasicDynamics<T, DESCRIPTOR, MOMENTA>(), _absorption(latticeAbsorption), _scattering(latticeScattering), _anisoMatrix(anisoMatrix)
{
  this->getName() = "RTLBMdynamicsMcHardy";
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T RTLBMdynamicsMcHardy<T, DESCRIPTOR, MOMENTA>::computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const
{
  return rho*descriptors::t<T,DESCRIPTOR>(iPop) - descriptors::t<T,DESCRIPTOR>(iPop);
}


template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> RTLBMdynamicsMcHardy<T, DESCRIPTOR, MOMENTA>::collide( Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics )
{
  std::array<double, DESCRIPTOR::q> feq = {};
  for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
    for ( int jPop = 0; jPop < DESCRIPTOR::q; ++jPop ) {
      feq[iPop] += (cell[jPop] + descriptors::t<T,DESCRIPTOR>(jPop)) * _anisoMatrix[jPop][iPop];
    }
    feq[iPop] *= descriptors::t<T,DESCRIPTOR>(iPop);
  }
  // execute collision
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = (cell[iPop]+descriptors::t<T,DESCRIPTOR>(iPop))
                 - descriptors::norm_c<T,DESCRIPTOR>(iPop)*(_absorption+_scattering) * ( (cell[iPop]+descriptors::t<T,DESCRIPTOR>(iPop))- feq[iPop] )
                 - _absorption*descriptors::norm_c<T,DESCRIPTOR>(iPop) *(cell[iPop]+descriptors::t<T,DESCRIPTOR>(iPop))
                 - descriptors::t<T,DESCRIPTOR>(iPop);
  }
  T temperature = lbm<DESCRIPTOR>::computeRho(cell);
  statistics.incrementStats( temperature, T() );
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T RTLBMdynamicsMcHardy<T, DESCRIPTOR, MOMENTA>::getOmega() const
{
  return -1;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void RTLBMdynamicsMcHardy<T, DESCRIPTOR, MOMENTA>::setOmega( T omega )
{
}

//==================================================================================//
template<typename T, typename DESCRIPTOR, typename MOMENTA>
RTLBMdynamicsMcHardyRK<T, DESCRIPTOR, MOMENTA>::RTLBMdynamicsMcHardyRK
(T latticeAbsorption, T latticeScattering, std::array<std::array<T,DESCRIPTOR::q>, DESCRIPTOR::q>& anisoMatrix)
  : legacy::BasicDynamics<T, DESCRIPTOR, MOMENTA>(), _absorption(latticeAbsorption), _scattering(latticeScattering), _anisoMatrix(anisoMatrix)
{
  this->getName() = "RTLBMdynamicsMcHardyRK";
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T RTLBMdynamicsMcHardyRK<T, DESCRIPTOR, MOMENTA>::computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const
{
  return rho*descriptors::t<T,DESCRIPTOR>(iPop) - descriptors::t<T,DESCRIPTOR>(iPop);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void RTLBMdynamicsMcHardyRK<T,DESCRIPTOR,MOMENTA>::computeEquilibriumAniso(Cell<T,DESCRIPTOR>& cell, std::array<T,DESCRIPTOR::q>& feq)
{
  feq.fill( T() );
  for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
    for ( int jPop = 0; jPop < DESCRIPTOR::q; ++jPop ) {
      feq[iPop] += cell[jPop] * _anisoMatrix[jPop][iPop];
    }
    feq[iPop] *= descriptors::t<T,DESCRIPTOR>(iPop);
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
std::array<T,DESCRIPTOR::q> RTLBMdynamicsMcHardyRK<T,DESCRIPTOR,MOMENTA>::doCollision(Cell<T,DESCRIPTOR>& cell, std::array<T,DESCRIPTOR::q>& feq)
{
  std::array<T,DESCRIPTOR::q> k;
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    k[iPop]  = - descriptors::norm_c<T,DESCRIPTOR>(iPop)*(_absorption+_scattering) * (cell[iPop])
               + descriptors::norm_c<T,DESCRIPTOR>(iPop)*_scattering * feq[iPop];
  }
  return k;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> RTLBMdynamicsMcHardyRK<T, DESCRIPTOR, MOMENTA>::collide( Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics )
{
  std::array<T,DESCRIPTOR::q> feq;
  std::array<T,DESCRIPTOR::q> f_pre_collision;
  // separate cell and precollision f_i
  for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
    f_pre_collision[iPop] = cell[iPop] + descriptors::t<T,DESCRIPTOR>(iPop);
  }

  // shift only first collision und equilibrium and then at the very end
  for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
    cell[iPop] += descriptors::t<T,DESCRIPTOR>(iPop);
  }
  computeEquilibriumAniso(cell,feq);
  std::array<T,DESCRIPTOR::q> k1 = doCollision(cell,feq);
  // update cell
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = f_pre_collision[iPop] + 0.5*k1[iPop];
  }

  computeEquilibriumAniso(cell,feq);
  std::array<T,DESCRIPTOR::q> k2 = doCollision(cell,feq);
  // update cell
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = f_pre_collision[iPop] + 0.5*k2[iPop];
  }

  computeEquilibriumAniso(cell,feq);
  std::array<T,DESCRIPTOR::q> k3 = doCollision(cell,feq);
  // update cell
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = f_pre_collision[iPop] + k3[iPop];
  }

  computeEquilibriumAniso(cell,feq);
  std::array<T,DESCRIPTOR::q> k4 = doCollision(cell,feq);
  // update cell
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = f_pre_collision[iPop] + 1/6.*(k1[iPop] + 2*k2[iPop] + 2*k3[iPop] + k4[iPop])
                 - descriptors::t<T,DESCRIPTOR>(iPop); // back shift for OpenLB
  }
  T temperature = lbm<DESCRIPTOR>::computeRho(cell);
  statistics.incrementStats( temperature, T() );
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T RTLBMdynamicsMcHardyRK<T, DESCRIPTOR, MOMENTA>::getOmega() const
{
  return -1;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void RTLBMdynamicsMcHardyRK<T, DESCRIPTOR, MOMENTA>::setOmega( T omega )
{
}


} // namespace olb


#endif

