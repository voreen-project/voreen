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

#ifndef LATTICE_PHYS_EFFEVTIVE_DISSIPATION_3D_HH
#define LATTICE_PHYS_EFFEVTIVE_DISSIPATION_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysEffevtiveDissipation3D.h"
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
SuperLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::SuperLatticePhysEffevtiveDissipation3D(
  SuperLattice<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
                                         std::function<T(Cell<T,DESCRIPTOR>&)> effectiveOmegaF)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physEffevtiveDissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC),
                               this->_converter, smagoConst, effectiveOmegaF));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::BlockLatticePhysEffevtiveDissipation3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
  std::function<T(Cell<T,DESCRIPTOR>&)> effectiveOmegaF)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter), _smagoConst(smagoConst), _effectiveOmegaF(effectiveOmegaF)
{
  this->getName() = "physEffevtiveDissipation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
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

  T dt = 1./ _converter.getConversionFactorTime();
  auto cell = this->_blockLattice.get(input[0], input[1], input[2]);
  T omegaEff = _effectiveOmegaF(cell);
  T nuEff = ((1./omegaEff)-0.5)/descriptors::invCs2<T,DESCRIPTOR>();  // BGK shear viscosity definition

  output[0] = PiNeqNormSqr * nuEff
              * util::pow(omegaEff * descriptors::invCs2<T,DESCRIPTOR>() / rho, 2) / 2.
              * _converter.getPhysViscosity() / _converter.getLatticeViscosity() / dt / dt;

  return true;
}

}
#endif
