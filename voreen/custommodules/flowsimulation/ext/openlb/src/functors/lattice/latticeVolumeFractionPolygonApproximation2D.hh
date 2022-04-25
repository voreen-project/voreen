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

#ifndef LATTICE_VOLUME_FRACTION_POLYGON_APPROXIMATION_2D_HH
#define LATTICE_VOLUME_FRACTION_POLYGON_APPROXIMATION_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticeVolumeFractionPolygonApproximation2D.h"
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
SuperLatticeVolumeFractionPolygonApproximation2D<T, DESCRIPTOR>::SuperLatticeVolumeFractionPolygonApproximation2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
  IndicatorF2D<T>& indicator, const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "volumeFractionPolygonApproximation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVolumeFractionPolygonApproximation2D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC),
                               superGeometry.getBlockGeometry(iC),
                               indicator, converter, insideOut));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeVolumeFractionPolygonApproximation2D<T, DESCRIPTOR>::BlockLatticeVolumeFractionPolygonApproximation2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
    BlockGeometry<T,2>& blockGeometry,
    IndicatorF2D<T>& indicator,
    const UnitConverter<T,DESCRIPTOR>& converter,
    bool insideOut)
  : BlockLatticeF2D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _indicator(indicator), _converter(converter), _insideOut(insideOut)
{
  this->getName() = "volumeFractionPolygonApproximation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeVolumeFractionPolygonApproximation2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  T physR[2];

  bool cornerXMYM_inside[1];
  bool cornerXMYP_inside[1];
  bool cornerXPYM_inside[1];
  bool cornerXPYP_inside[1];
  _blockGeometry.getPhysR(physR, {input[0], input[1]});

  T cornerXMYM[2];
  T cornerXMYP[2];
  T cornerXPYM[2];
  T cornerXPYP[2];

  cornerXMYM[0] = physR[0] - 0.5*_converter.getPhysDeltaX();
  cornerXMYM[1] = physR[1] - 0.5*_converter.getPhysDeltaX();

  cornerXMYP[0] = physR[0] - 0.5*_converter.getPhysDeltaX();
  cornerXMYP[1] = physR[1] + 0.5*_converter.getPhysDeltaX();

  cornerXPYM[0] = physR[0] + 0.5*_converter.getPhysDeltaX();
  cornerXPYM[1] = physR[1] - 0.5*_converter.getPhysDeltaX();

  cornerXPYP[0] = physR[0] + 0.5*_converter.getPhysDeltaX();
  cornerXPYP[1] = physR[1] + 0.5*_converter.getPhysDeltaX();

  _indicator(cornerXMYM_inside, cornerXMYM);
  _indicator(cornerXMYP_inside, cornerXMYP);
  _indicator(cornerXPYM_inside, cornerXPYM);
  _indicator(cornerXPYP_inside, cornerXPYP);

  if ((cornerXMYM_inside[0] && cornerXMYP_inside[0] && cornerXPYM_inside[0] && cornerXPYP_inside[0]) ||
      (!cornerXMYM_inside[0] && !cornerXMYP_inside[0] && !cornerXPYM_inside[0] && !cornerXPYP_inside[0]) ) {
    if (!_insideOut) {
      if (cornerXMYM_inside[0]) {
        output[0] = 1.0;
      }
    }
    else {
      if (!cornerXMYM_inside[0]) {
        output[0] = 1.0;
      }
    }
  }
  else {
    Vector<T,2> cornerXMYM_vec(physR[0] - 0.5*_converter.getPhysDeltaX(), physR[1] - 0.5*_converter.getPhysDeltaX());
    Vector<T,2> cornerXPYP_vec(physR[0] + 0.5*_converter.getPhysDeltaX(), physR[1] + 0.5*_converter.getPhysDeltaX());

    Vector<T,2> directionXP(_converter.getPhysDeltaX(), 0);
    Vector<T,2> directionXM(-_converter.getPhysDeltaX(), 0);
    Vector<T,2> directionYP(0, _converter.getPhysDeltaX());
    Vector<T,2> directionYM(0, -_converter.getPhysDeltaX());

    T distanceXP = 1.01 * _converter.getPhysDeltaX();
    T distanceXM = 1.01 * _converter.getPhysDeltaX();
    T distanceYP = 1.01 * _converter.getPhysDeltaX();
    T distanceYM = 1.01 * _converter.getPhysDeltaX();

    if ((cornerXMYM_inside[0] && !cornerXMYP_inside[0]) ||
        (!cornerXMYM_inside[0] && cornerXMYP_inside[0]) ) {
      _indicator.distance(distanceYP, cornerXMYM, directionYP);
    }

    if ((cornerXMYM_inside[0] && !cornerXPYM_inside[0]) ||
        (!cornerXMYM_inside[0] && cornerXPYM_inside[0]) ) {
      _indicator.distance(distanceXP, cornerXMYM, directionXP);
    }

    if ((cornerXPYP_inside[0] && !cornerXMYP_inside[0]) ||
        (!cornerXPYP_inside[0] && cornerXMYP_inside[0]) ) {
      _indicator.distance(distanceXM, cornerXPYP, directionXM);
    }

    if ((cornerXPYP_inside[0] && !cornerXPYM_inside[0]) ||
        (!cornerXPYP_inside[0] && cornerXPYM_inside[0]) ) {
      _indicator.distance(distanceYM, cornerXPYP, directionYM);
    }

    T volumeFraction = 0.0;

    if (distanceXP < _converter.getPhysDeltaX() && distanceXM < _converter.getPhysDeltaX()) {
      volumeFraction = distanceXP * _converter.getPhysDeltaX();
      volumeFraction += 0.5 * (_converter.getPhysDeltaX() - distanceXM - distanceXP) * _converter.getPhysDeltaX();
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXMYM_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceYP < _converter.getPhysDeltaX() && distanceYM < _converter.getPhysDeltaX()) {
      volumeFraction = distanceYP * _converter.getPhysDeltaX();
      volumeFraction += 0.5 * (_converter.getPhysDeltaX() - distanceYM - distanceYP) * _converter.getPhysDeltaX();
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXMYM_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceXP < _converter.getPhysDeltaX() && distanceYP < _converter.getPhysDeltaX()) {
      volumeFraction = 0.5 * distanceXP * distanceYP;
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXMYM_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceXM < _converter.getPhysDeltaX() && distanceYM < _converter.getPhysDeltaX()) {
      volumeFraction = 0.5 * distanceXM * distanceYM;
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXPYP_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceXM < _converter.getPhysDeltaX() && distanceYP < _converter.getPhysDeltaX()) {
      volumeFraction = 0.5 * (_converter.getPhysDeltaX() - distanceXM) * (_converter.getPhysDeltaX() - distanceYP);
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXMYP_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceYM < _converter.getPhysDeltaX() && distanceXP < _converter.getPhysDeltaX()) {
      volumeFraction = 0.5 * (_converter.getPhysDeltaX() - distanceXP) * (_converter.getPhysDeltaX() - distanceYM);
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXPYM_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (!_insideOut) {
      output[0] = volumeFraction;
    }
    else {
      output[0] = 1.0 - volumeFraction;
    }

  }
  return true;
}

}
#endif
