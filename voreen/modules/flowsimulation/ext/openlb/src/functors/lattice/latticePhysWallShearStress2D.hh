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

#ifndef LATTICE_PHYS_WALL_SHEAR_STRESS_2D_HH
#define LATTICE_PHYS_WALL_SHEAR_STRESS_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePhysWallShearStress2D.h"
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
SuperLatticePhysWallShearStress2D<T, DESCRIPTOR>::SuperLatticePhysWallShearStress2D(
  SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter,
  IndicatorF2D<T>& indicator)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physWallShearStress";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysWallShearStress2D<T, DESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        _superGeometry.getBlockGeometry(iC),
        _material,
        this->_converter,
        indicator)
    );
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysWallShearStress2D<T,DESCRIPTOR>::BlockLatticePhysWallShearStress2D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockGeometry<T,2>& blockGeometry,
  int material,
  const UnitConverter<T,DESCRIPTOR>& converter,
  IndicatorF2D<T>& indicator)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physWallShearStress";
  const T scaling = this->_converter.getConversionFactorLength() * 0.1;
  const T omega = 1. / this->_converter.getLatticeRelaxationTime();
  const T dt = this->_converter.getConversionFactorTime();
  _physFactor = -omega * descriptors::invCs2<T,DESCRIPTOR>() / dt * this->_converter.getPhysDensity() * this->_converter.getPhysViscosity();
  std::vector<int> discreteNormalOutwards(3, 0);

  for (int iX = 1; iX < _blockGeometry.getNx() - 1; iX++) {
    _discreteNormal.resize(_blockGeometry.getNx() - 2);
    _normal.resize(_blockGeometry.getNx() - 2);

    for (int iY = 1; iY < _blockGeometry.getNy() - 1; iY++) {
      _discreteNormal[iX-1].resize(_blockGeometry.getNy() - 2);
      _normal[iX-1].resize(_blockGeometry.getNy() - 2);

      if (_blockGeometry.get({iX, iY}) == _material) {
        discreteNormalOutwards = _blockGeometry.getStatistics().getType(iX, iY);
        _discreteNormal[iX-1][iY-1].resize(2);
        _normal[iX-1][iY-1].resize(2);

        _discreteNormal[iX-1][iY-1][0] = -discreteNormalOutwards[1];
        _discreteNormal[iX-1][iY-1][1] = -discreteNormalOutwards[2];

        T physR[2];
        _blockGeometry.getPhysR(physR, {iX, iY});
        Vector<T,2> origin(physR[0],physR[1]);
        Vector<T,2> direction(-_discreteNormal[iX-1][iY-1][0] * scaling,
                              -_discreteNormal[iX-1][iY-1][1] * scaling);
        Vector<T,2> normal(0.,0.);
        origin[0] = physR[0];
        origin[1] = physR[1];

        indicator.normal(normal, origin, direction);
        normal = normalize(normal);

        _normal[iX-1][iY-1][0] = normal[0];
        _normal[iX-1][iY-1][1] = normal[0];
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysWallShearStress2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = T();

  if (_blockGeometry.getNeighborhoodRadius(input) < 1) {
#ifdef OLB_DEBUG
    std::cout << "Input address not mapped by _discreteNormal, overlap too small" << std::endl;
#endif
    return true;
  }

  if (this->_blockGeometry.get(input)==_material) {

    T traction[2];
    T stress[3];
    T rho = this->_blockLattice.get(input[0] + _discreteNormal[input[0]-1][input[1]-1][0],
                                    input[1] + _discreteNormal[input[0]-1][input[1]-1][1]).computeRho();

    this->_blockLattice.get(input[0] +   _discreteNormal[input[0]-1][input[1]-1][0],
                            input[1] +   _discreteNormal[input[0]-1][input[1]-1][1]).computeStress(stress);

    traction[0] = stress[0]*_physFactor/rho*_normal[input[0]-1][input[1]-1][0] +
                  stress[1]*_physFactor/rho*_normal[input[0]-1][input[1]-1][1];
    traction[1] = stress[1]*_physFactor/rho*_normal[input[0]-1][input[1]-1][0] +
                  stress[2]*_physFactor/rho*_normal[input[0]-1][input[1]-1][1];

    T traction_normal_SP;
    T tractionNormalComponent[2];
    // scalar product of traction and normal vector
    traction_normal_SP = traction[0] * _normal[input[0]-1][input[1]-1][0] +
                         traction[1] * _normal[input[0]-1][input[1]-1][1];
    tractionNormalComponent[0] = traction_normal_SP * _normal[input[0]-1][input[1]-1][0];
    tractionNormalComponent[1] = traction_normal_SP * _normal[input[0]-1][input[1]-1][1];

    T WSS[2];
    WSS[0] = traction[0] - tractionNormalComponent[0];
    WSS[1] = traction[1] - tractionNormalComponent[1];
    // magnitude of the wall shear stress vector
    output[0] = util::sqrt(WSS[0]*WSS[0] + WSS[1]*WSS[1]);
    return true;
  }
  else {
    return true;
  }
}

}
#endif
