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

#ifndef LATTICE_PHYS_HEAT_FLUX_BOUNDARY_3D_HH
#define LATTICE_PHYS_HEAT_FLUX_BOUNDARY_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysHeatFluxBoundary3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysHeatFluxBoundary3D<T, DESCRIPTOR, TDESCRIPTOR>::SuperLatticePhysHeatFluxBoundary3D(
  SuperLattice<T, TDESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry,
  const int material, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter,
  IndicatorF3D<T>& indicator)
  : SuperLatticeThermalPhysF3D<T, DESCRIPTOR, TDESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physHeatFluxBoundary";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysHeatFluxBoundary3D<T, DESCRIPTOR, TDESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        _superGeometry.getBlockGeometry(iC),
        _material,
        this->_converter,
        indicator)
    );
  }
}

// constructor calculates the normal of the boundary on every lattice for one specific material number
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysHeatFluxBoundary3D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysHeatFluxBoundary3D(
  BlockLattice<T,TDESCRIPTOR>& blockLattice,
  BlockGeometry<T,3>& blockGeometry,
  int material,
  const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter,
  IndicatorF3D<T>& indicator)
  : BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,1),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "physHeatFluxBoundary";
  const T scaling = this->_converter.getConversionFactorLength() * 0.1;
  std::vector<int> discreteNormalOutwards(4, 0);

  for (int iX = 1 ; iX < _blockGeometry.getNx() - 1; iX++) {
    _discreteNormal.resize(_blockGeometry.getNx() - 2);
    _normal.resize(_blockGeometry.getNx() - 2);

    for (int iY = 1; iY < _blockGeometry.getNy() - 1; iY++) {
      _discreteNormal[iX-1].resize(_blockGeometry.getNy() - 2);
      _normal[iX-1].resize(_blockGeometry.getNy() - 2);

      for (int iZ = 1; iZ < _blockGeometry.getNz() - 1; iZ++) {
        _discreteNormal[iX-1][iY-1].resize(_blockGeometry.getNz() - 2);
        _normal[iX-1][iY-1].resize(_blockGeometry.getNz() - 2);

        if (_blockGeometry.get({iX, iY, iZ}) == _material) {
          discreteNormalOutwards = _blockGeometry.getStatistics().getType(iX, iY, iZ);
          _discreteNormal[iX - 1][iY - 1][iZ - 1].resize(3);
          _normal[iX - 1][iY - 1][iZ - 1].resize(3);

          _discreteNormal[iX - 1][iY- 1][iZ- 1][0] = -discreteNormalOutwards[1];
          _discreteNormal[iX- 1][iY- 1][iZ- 1][1] = -discreteNormalOutwards[2];
          _discreteNormal[iX- 1][iY- 1][iZ- 1][2] = -discreteNormalOutwards[3];

          T physR[3];
          _blockGeometry.getPhysR(physR,{iX, iY, iZ});
          Vector<T,3> origin(physR[0],physR[1],physR[2]);
          Vector<T,3> direction(-_discreteNormal[iX- 1][iY- 1][iZ- 1][0] * scaling,
                                -_discreteNormal[iX- 1][iY- 1][iZ- 1][1] * scaling,
                                -_discreteNormal[iX- 1][iY- 1][iZ- 1][2] * scaling);
          Vector<T,3> normal(0.,0.,0.);
          origin[0] = physR[0];
          origin[1] = physR[1];
          origin[2] = physR[2];

          indicator.normal(normal, origin, direction);
          normalize(normal);

          _normal[iX- 1][iY- 1][iZ- 1][0] = normal[0];
          _normal[iX- 1][iY- 1][iZ- 1][1] = normal[1];
          _normal[iX- 1][iY- 1][iZ- 1][2] = normal[2];
        }
      }
    }
  }
}

// functor calculates specific heat flux for one lattice, which is selected with the input
template<typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysHeatFluxBoundary3D<T, DESCRIPTOR, TDESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();

  if (this->_blockLattice.getNeighborhoodRadius(input) < 1) {
#ifdef OLB_DEBUG
    std::cout << "Input address not mapped by _discreteNormal, overlap too small" << std::endl;
#endif
    return true;
  }

  if (_blockGeometry.get(input) == _material) {

    // lattice temperature next to the boundary in the direction of the normal
    T temp1 = this->_blockLattice.get(
                input[0] + _discreteNormal[input[0]-1][input[1]-1][input[2]-1][0],
                input[1] + _discreteNormal[input[0]-1][input[1]-1][input[2]-1][1],
                input[2] + _discreteNormal[input[0]-1][input[1]-1][input[2]-1][2]).computeRho();

    // second lattice temperature in the direction of the normal
    T temp2 = this->_blockLattice.get(
                input[0] + 2*_discreteNormal[input[0]-1][input[1]-1][input[2]-1][0],
                input[1] + 2*_discreteNormal[input[0]-1][input[1]-1][input[2]-1][1],
                input[2] + 2*_discreteNormal[input[0]-1][input[1]-1][input[2]-1][2]).computeRho();

    // calculation of deltaX with the normal, yields higher accuracy for deltaX
    T deltaX = util::sqrt( _normal[input[0]-1][input[1]-1][input[2]-1][0] *
                           _normal[input[0]-1][input[1]-1][input[2]-1][0] +
                           _normal[input[0]-1][input[1]-1][input[2]-1][1] *
                           _normal[input[0]-1][input[1]-1][input[2]-1][1] +
                           _normal[input[0]-1][input[1]-1][input[2]-1][2] *
                           _normal[input[0]-1][input[1]-1][input[2]-1][2] );

    // lattice temperature on the boundary
    // ATTENTION: here the temperature is hardcoded as the characteristic high temperature
    T temp0 = 1.5;

    //return of specific heat flux in direction of the normal
    output[0] = this->_converter.getThermalConductivity() *
                this->_converter.getCharPhysTemperatureDifference() *
                (1.5 * temp0 - 2.0 * temp1 + 0.5 * temp2) /
                this->_converter.getPhysLength(deltaX);

    return true;
  }
  else {
    return true;
  }
}

}
#endif
