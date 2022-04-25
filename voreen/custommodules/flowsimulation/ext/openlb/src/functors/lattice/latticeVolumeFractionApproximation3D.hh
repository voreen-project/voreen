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

#ifndef LATTICE_VOLUME_FRACTION_APPROXIMATION_3D_HH
#define LATTICE_VOLUME_FRACTION_APPROXIMATION_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticeVolumeFractionApproximation3D.h"
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
SuperLatticeVolumeFractionApproximation3D<T, DESCRIPTOR>::SuperLatticeVolumeFractionApproximation3D(
  SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry,
  IndicatorF3D<T>& indicator, int refinementLevel,
  const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "volumeFractionApproximation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVolumeFractionApproximation3D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC),
                               superGeometry.getBlockGeometry(iC),
                               indicator, refinementLevel,
                               converter, insideOut));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeVolumeFractionApproximation3D<T, DESCRIPTOR>::BlockLatticeVolumeFractionApproximation3D(BlockLattice<T,DESCRIPTOR>& blockLattice,
    BlockGeometry<T,3>& blockGeometry,
    IndicatorF3D<T>& indicator,
    int refinementLevel,
    const UnitConverter<T,DESCRIPTOR>& converter,
    bool insideOut)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _indicator(indicator), _refinementLevel(refinementLevel), _converter(converter), _insideOut(insideOut),
    _physSubGridMinPhysRshift((1./ T(_refinementLevel * 2.) - 0.5) * _converter.getPhysDeltaX()),
    _physSubGridDeltaX(1. / T(_refinementLevel) * _converter.getPhysDeltaX()),
    _latticeSubGridVolume(1. / (_refinementLevel * _refinementLevel * _refinementLevel))
{
  this->getName() = "volumeFractionApproximation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeVolumeFractionApproximation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  T physR[3];
  bool inside[1];
  _blockGeometry.getPhysR(physR, {input[0], input[1], input[2]});

  T subGridMinPhysR[3];
  subGridMinPhysR[0] = physR[0] + _physSubGridMinPhysRshift;
  subGridMinPhysR[1] = physR[1] + _physSubGridMinPhysRshift;
  subGridMinPhysR[2] = physR[2] + _physSubGridMinPhysRshift;

  for (int x = 0; x < _refinementLevel; x++) {
    for (int y = 0; y < _refinementLevel; y++) {
      for (int z = 0; z < _refinementLevel; z++) {
        physR[0] = subGridMinPhysR[0] + x * _physSubGridDeltaX;
        physR[1] = subGridMinPhysR[1] + y * _physSubGridDeltaX;
        physR[2] = subGridMinPhysR[2] + z * _physSubGridDeltaX;
        _indicator(inside, physR);
        if (!_insideOut) {
          if (inside[0]) {
            output[0] += _latticeSubGridVolume;
          }
        }
        else {
          if (!inside[0]) {
            output[0] += _latticeSubGridVolume;
          }
        }
      }
    }
  }
  return true;
}

}
#endif
