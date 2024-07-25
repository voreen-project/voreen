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

#ifndef LATTICE_PHYS_VELOCITY_3D_HH
#define LATTICE_PHYS_VELOCITY_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysVelocity3D.h"
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
SuperLatticePhysVelocity3D<T, DESCRIPTOR>::SuperLatticePhysVelocity3D(
  SuperLattice<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter, bool print)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3), _print(print)
{
  this->getName() = "physVelocity";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysVelocity3D<T, DESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        this->_converter,
        _print)
    );
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysVelocity3D<T, DESCRIPTOR>::BlockLatticePhysVelocity3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  bool print)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _print(print)
{
  this->getName() = "physVelocity";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (_print) {
    std::cout << input[0] << " " << input[1] << " " << input[2] << " | "
              << singleton::mpi().getRank() << std::endl;
  }

  this->_blockLattice.get(input).computeU(output);
  output[0] = this->_converter.getPhysVelocity(output[0]);
  output[1] = this->_converter.getPhysVelocity(output[1]);
  output[2] = this->_converter.getPhysVelocity(output[2]);

  return true;
}

}
#endif
