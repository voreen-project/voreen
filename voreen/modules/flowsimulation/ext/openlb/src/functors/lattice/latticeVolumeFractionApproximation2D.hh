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

#ifndef LATTICE_VOLUME_FRACTION_APPROXIMATION_2D_HH
#define LATTICE_VOLUME_FRACTION_APPROXIMATION_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticeVolumeFractionApproximation2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticeVolumeFractionApproximation2D<T, DESCRIPTOR>::SuperLatticeVolumeFractionApproximation2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
  IndicatorF2D<T>& indicator, int refinementLevel,
  const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "volumeFractionApproximation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVolumeFractionApproximation2D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC),
                               superGeometry.getBlockGeometry(iC),
                               indicator, refinementLevel,
                               converter, insideOut));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeVolumeFractionApproximation2D<T, DESCRIPTOR>::BlockLatticeVolumeFractionApproximation2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
    BlockGeometry<T,2>& blockGeometry,
    IndicatorF2D<T>& indicator,
    int refinementLevel,
    const UnitConverter<T,DESCRIPTOR>& converter,
    bool insideOut)
  : BlockLatticeF2D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _indicator(indicator), _refinementLevel(refinementLevel), _converter(converter), _insideOut(insideOut),
    _physSubGridMinPhysRshift((1./ T(_refinementLevel * 2.) - 0.5) * _converter.getPhysDeltaX()),
    _physSubGridDeltaX(1. / T(_refinementLevel) * _converter.getPhysDeltaX()),
    _latticeSubGridVolume(1. / (_refinementLevel * _refinementLevel))
{
  this->getName() = "volumeFractionApproximation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeVolumeFractionApproximation2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  T physR[2];
  bool inside[1];
  _blockGeometry.getPhysR(physR, {input[0], input[1]});

  T subGridMinPhysR[2];
  subGridMinPhysR[0] = physR[0] + _physSubGridMinPhysRshift;
  subGridMinPhysR[1] = physR[1] + _physSubGridMinPhysRshift;

  for (int x = 0; x < _refinementLevel; x++) {
    for (int y = 0; y < _refinementLevel; y++) {
      physR[0] = subGridMinPhysR[0] + x * _physSubGridDeltaX;
      physR[1] = subGridMinPhysR[1] + y * _physSubGridDeltaX;
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
  return true;
}

}
#endif
