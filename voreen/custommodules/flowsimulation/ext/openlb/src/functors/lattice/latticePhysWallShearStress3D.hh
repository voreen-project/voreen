/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Jonathan Jeppener-Haltenhoff
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

#ifndef LATTICE_PHYS_WALL_SHEAR_STRESS_3D_HH
#define LATTICE_PHYS_WALL_SHEAR_STRESS_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysWallShearStress3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry3D.h"
#include "blockBaseF3D.h"
#include "core/blockLatticeStructure3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticePhysWallShearStress3D<T, DESCRIPTOR>::SuperLatticePhysWallShearStress3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter,
  IndicatorF3D<T>& indicator)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physWallShearStress";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysWallShearStress3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        _superGeometry.getExtendedBlockGeometry(iC),
        this->_sLattice.getOverlap(),
        _material,
        this->_converter,
        indicator)
    );
  }
}

//SL
template <typename T, typename DESCRIPTOR>
size_t BlockLatticePhysWallShearStress3D<T,DESCRIPTOR>::index(size_t x, size_t y, size_t z) const {
    return (x-1)*(_blockGeometry.getNz()-2)*(_blockGeometry.getNy()-2)*3 + (y-1) * (_blockGeometry.getNz()-2)*3 + (z-1)*3;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysWallShearStress3D<T,DESCRIPTOR>::BlockLatticePhysWallShearStress3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
  BlockGeometryStructure3D<T>& blockGeometry,
  int overlap,
  int material,
  const UnitConverter<T,DESCRIPTOR>& converter,
  IndicatorF3D<T>& indicator)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,1),
    _blockGeometry(blockGeometry),
    _overlap(overlap),
    _material(material)
{
  this->getName() = "physWallShearStress";
  const T scaling = this->_converter.getConversionFactorLength() * 0.1;
  const T omega = 1. / this->_converter.getLatticeRelaxationTime();
  const T dt = this->_converter.getConversionFactorTime();
  _physFactor = -omega * descriptors::invCs2<T,DESCRIPTOR>() / dt * this->_converter.getPhysDensity() * this->_converter.getPhysViscosity();
  std::vector<int> discreteNormalOutwards(4, 0);

  _discreteNormal.resize((_blockGeometry.getNx()-2)*(_blockGeometry.getNy()-2)*(_blockGeometry.getNz()-2)*3); //SL
  _normal.resize((_blockGeometry.getNx()-2)*(_blockGeometry.getNy()-2)*(_blockGeometry.getNz()-2)*3);         //SL

  _blockGeometry.getStatistics().update(true); //SL

#pragma omp parallel for private(discreteNormalOutwards) //SL
  for (int iX = 1 ; iX < _blockGeometry.getNx() - 1; iX++) {
    //_discreteNormal.resize(_blockGeometry.getNx() - 2); //SL
    //_normal.resize(_blockGeometry.getNx() - 2); //SL

    for (int iY = 1; iY < _blockGeometry.getNy() - 1; iY++) {
      //_discreteNormal[iX-1].resize(_blockGeometry.getNy() - 2); //SL
      //_normal[iX-1].resize(_blockGeometry.getNy() - 2); //SL

      for (int iZ = 1; iZ < _blockGeometry.getNz() - 1; iZ++) {
        //_discreteNormal[iX-1][iY-1].resize(_blockGeometry.getNz() - 2); //SL
        //_normal[iX-1][iY-1].resize(_blockGeometry.getNz() - 2); //SL

        if (_blockGeometry.get(iX, iY, iZ) == _material) {
          const size_t idx = index(iX, iY, iZ); //SL

          discreteNormalOutwards = const_cast<const BlockGeometryStructure3D<T>&>(_blockGeometry).getStatistics().getType(iX, iY, iZ); //SL

          _discreteNormal[idx+0] = -discreteNormalOutwards[1];
          _discreteNormal[idx+1] = -discreteNormalOutwards[2];
          _discreteNormal[idx+2] = -discreteNormalOutwards[3];

          T physR[3];
          _blockGeometry.getPhysR(physR,iX, iY, iZ);
          Vector<T,3> origin(physR[0],physR[1],physR[2]);
          Vector<T,3> direction(-_discreteNormal[idx+0] * scaling,
                                -_discreteNormal[idx+1] * scaling,
                                -_discreteNormal[idx+2] * scaling);
          Vector<T,3> normal(0.,0.,0.);
          origin[0] = physR[0];
          origin[1] = physR[1];
          origin[2] = physR[2];

          indicator.normal(normal, origin, direction);
          normal = normalize(normal);

          _normal[idx+0] = normal[0];
          _normal[idx+1] = normal[1];
          _normal[idx+2] = normal[2];
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysWallShearStress3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();

  if (input[0] + _overlap < 1 ||
      input[1] + _overlap < 1 ||
      input[2] + _overlap < 1 ||
      input[0] + _overlap >= _blockGeometry.getNx()-1 ||
      input[1] + _overlap >= _blockGeometry.getNy()-1 ||
      input[2] + _overlap >= _blockGeometry.getNz()-1 ) {

#ifdef OLB_DEBUG
    std::cout << "Input address not mapped by _discreteNormal, overlap too small" << std::endl;
#endif
    return true;
  }

  if (_blockGeometry.get(input[0]+_overlap,input[1]+_overlap,input[2]+_overlap) == _material) {

    const size_t idx = index(input[0]+_overlap-1, input[1]+_overlap-1, input[2]+_overlap-1); //SL

    T traction[3];
    T stress[6];
    T rho = this->_blockLattice.get(input[0] + _overlap + _discreteNormal[idx+0],
                                    input[1] + _overlap + _discreteNormal[idx+1],
                                    input[2] + _overlap + _discreteNormal[idx+2]).computeRho();
    this->_blockLattice.get(input[0] + _overlap +   _discreteNormal[idx+0],
                            input[1] + _overlap +   _discreteNormal[idx+1],
                            input[2] + _overlap +   _discreteNormal[idx+2]).computeStress(stress);

    traction[0] = stress[0]*_physFactor/rho*_normal[idx+0] +
                  stress[1]*_physFactor/rho*_normal[idx+1] +
                  stress[2]*_physFactor/rho*_normal[idx+2];
    traction[1] = stress[1]*_physFactor/rho*_normal[idx+0] +
                  stress[3]*_physFactor/rho*_normal[idx+1] +
                  stress[4]*_physFactor/rho*_normal[idx+2];
    traction[2] = stress[2]*_physFactor/rho*_normal[idx+0] +
                  stress[4]*_physFactor/rho*_normal[idx+1] +
                  stress[5]*_physFactor/rho*_normal[idx+2];

    T traction_normal_SP;
    T tractionNormalComponent[3];
    // scalar product of traction and normal vector
    traction_normal_SP = traction[0] * _normal[idx+0] +
                         traction[1] * _normal[idx+1] +
                         traction[2] * _normal[idx+2];
    tractionNormalComponent[0] = traction_normal_SP * _normal[idx+0];
    tractionNormalComponent[1] = traction_normal_SP * _normal[idx+1];
    tractionNormalComponent[2] = traction_normal_SP * _normal[idx+2];

    T WSS[3];
    WSS[0] = traction[0] - tractionNormalComponent[0];
    WSS[1] = traction[1] - tractionNormalComponent[1];
    WSS[2] = traction[2] - tractionNormalComponent[2];
    // magnitude of the wall shear stress vector
    output[0] = sqrt(WSS[0]*WSS[0] + WSS[1]*WSS[1] + WSS[2]*WSS[2]);

    return true;
  }
  else {
    return true;
  }
}

}
#endif
