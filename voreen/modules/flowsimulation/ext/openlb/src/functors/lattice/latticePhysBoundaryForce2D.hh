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

#ifndef LATTICE_PHYS_BOUNDARY_FORCE_2D_HH
#define LATTICE_PHYS_BOUNDARY_FORCE_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePhysBoundaryForce2D.h"
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
SuperLatticePhysBoundaryForce2D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce2D(
  SuperLattice<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 2),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "physBoundaryForce";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysBoundaryForce2D<T, DESCRIPTOR>(
        this->_sLattice.getBlock(iC),
        _indicatorF->getBlockIndicatorF(iC),
        this->_converter));
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysBoundaryForce2D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce2D(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  SuperGeometry<T,2>& superGeometry, const int material,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysBoundaryForce2D(sLattice,
                                    superGeometry.getMaterialIndicator(material),
                                    converter)
{ }

template <typename T, typename DESCRIPTOR>
BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR>::BlockLatticePhysBoundaryForce2D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockIndicatorF2D<T>&                  indicatorF,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2),
    _indicatorF(indicatorF),
    _blockGeometry(indicatorF.getBlockGeometry())
{
  this->getName() = "physBoundaryForce";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T();
  }

  if (_indicatorF(input)) {
    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
      // Get direction
      // Get next cell located in the current direction
      // Check if the next cell is a fluid node
      if (_blockGeometry.get({input[0] + descriptors::c<DESCRIPTOR >(iPop,0), input[1] + descriptors::c<DESCRIPTOR >(iPop,1)}) == 1) {
        // Get f_q of next fluid cell where l = opposite(q)
        T f = this->_blockLattice.get(input[0] + descriptors::c<DESCRIPTOR >(iPop,0), input[1] + descriptors::c<DESCRIPTOR >(iPop,1))[iPop];
        // Get f_l of the boundary cell
        // Add f_q and f_opp
        f += this->_blockLattice.get(input[0], input[1])[util::opposite<DESCRIPTOR >(iPop)];
        // Update force
        for (int i = 0; i < this->getTargetDim(); ++i) {
          output[i] -= descriptors::c<DESCRIPTOR >(iPop,i) * f;
        }
      }
    }
    for (int i = 0; i < this->getTargetDim(); ++i) {
      output[i] = this->_converter.getPhysForce(output[i]);
    }
  }
  return true;
}

}
#endif
