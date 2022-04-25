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

#ifndef LATTICE_EFFEVTIVE_DISSIPATION_3D_HH
#define LATTICE_EFFEVTIVE_DISSIPATION_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticeEffevtiveDissipation3D.h"
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
SuperLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::SuperLatticeEffevtiveDissipation3D(
  SuperLattice<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst, LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1), _converter(converter)
{
  this->getName() = "EffevtiveDissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC),
                               this->_converter, smagoConst, LESdynamics));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::BlockLatticeEffevtiveDissipation3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
  LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter), _smagoConst(smagoConst), _LESdynamics(LESdynamics)
{
  this->getName() = "EffevtiveDissipation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  auto cell = this->_blockLattice.get(input[0], input[1], input[2]);
  T omegaEff = _LESdynamics.getEffectiveOmega(cell);
  T nuEff = ((1./omegaEff)-0.5)/descriptors::invCs2<T,DESCRIPTOR>();

  output[0] = PiNeqNormSqr * (nuEff)
              * util::pow(omegaEff * descriptors::invCs2<T,DESCRIPTOR>() / rho, 2)  / 2.;

  return true;
}

}
#endif
