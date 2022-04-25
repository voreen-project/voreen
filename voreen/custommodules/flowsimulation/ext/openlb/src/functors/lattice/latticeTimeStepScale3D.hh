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

#ifndef LATTICE_TIME_STEP_SCALE_3D_HH
#define LATTICE_TIME_STEP_SCALE_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticeTimeStepScale3D.h"
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
SuperLatticeTimeStepScale3D<T, DESCRIPTOR>::SuperLatticeTimeStepScale3D(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  T oldTau, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, DESCRIPTOR::q)
{
  this->getName() = "latticeTimeStepScale";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticeTimeStepScale3D<T, DESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        oldTau,
        converter)
    );
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticeTimeStepScale3D<T,DESCRIPTOR>::BlockLatticeTimeStepScale3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, T oldTau, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, DESCRIPTOR::q), _tau_old(oldTau), _converter(converter)
{
  this->getName() = "latticeTimeStepScale";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeTimeStepScale3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

  Cell<T,DESCRIPTOR> cell = this->_blockLattice.get(input[0], input[1], input[2]);
  Dynamics<T,DESCRIPTOR>* dynamics = this->_blockLattice.getDynamics(input[0], input[1], input[2]);
  T rho_old, rho_new, u_old[DESCRIPTOR::d], u_new[DESCRIPTOR::d];;
  cell.computeRhoU(rho_old,u_old);
  T fNeq[DESCRIPTOR::q];

  T uSqr_old = util::normSqr<T,DESCRIPTOR::d>(u_old);

  T tau_new = _converter.getLatticeRelaxationTime();

  T physDeltaT_new = _converter.getPhysDeltaT();
  T physDeltaT_old = physDeltaT_new * (_tau_old-0.5) / (tau_new-0.5);
  T neqScaling = ((physDeltaT_new*tau_new)/(physDeltaT_old*_tau_old));

  for (int iDim=0; iDim < DESCRIPTOR::d; iDim++) {
    u_new[iDim] = u_old[iDim] * physDeltaT_new / physDeltaT_old ;
  }

  rho_new = (rho_old-1.0) * (physDeltaT_new / physDeltaT_old) * (physDeltaT_new / physDeltaT_old) + 1.0;

  T uSqr_new = util::normSqr<T,DESCRIPTOR::d>(u_new);


  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    fNeq[iPop] = cell[iPop] - dynamics -> computeEquilibrium(iPop,rho_old,u_old,uSqr_old);
//    output[iPop] = cell[iPop];
    output[iPop] = (dynamics -> computeEquilibrium(iPop,rho_old,u_old,uSqr_old)+descriptors::t<T,DESCRIPTOR>(iPop)// * ((dynamics -> computeEquilibrium(iPop,rho_old,u_new,uSqr_new)+descriptors::t<T,DESCRIPTOR>(iPop))/(dynamics -> computeEquilibrium(iPop,rho_old,u_old,uSqr_old)+descriptors::t<T,DESCRIPTOR>(iPop)))
                    + fNeq[iPop]*neqScaling)
                   * ((dynamics -> computeEquilibrium(iPop,rho_new,u_new,uSqr_new)+descriptors::t<T,DESCRIPTOR>(iPop))/(dynamics -> computeEquilibrium(iPop,rho_old,u_old,uSqr_old)+descriptors::t<T,DESCRIPTOR>(iPop))) - descriptors::t<T,DESCRIPTOR>(iPop);
  }
  return true;
}

}
#endif
