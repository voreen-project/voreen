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

#ifndef LATTICE_PHYS_CORR_BOUNDARY_FORCE_2D_HH
#define LATTICE_PHYS_CORR_BOUNDARY_FORCE_2D_HH

#include <vector>
#include <cmath>
#include <limits>

#include "latticePhysCorrBoundaryForce2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry2D.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "core/blockLattice2D.h"
#include "communication/mpiManager.h"
#include "core/blockLatticeStructure2D.h"


namespace olb {

template<typename T,typename DESCRIPTOR>
SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::SuperLatticePhysCorrBoundaryForce2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physCorrBoundaryForce";
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  std::vector<T> force(3, T());
  //  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
  //    int globX = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosX() + locix;
  //    int globY = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosY() + lociy;
  //    int globZ = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosZ() + lociz;
  //    if(superGeometry.getMaterial(globX, globY, globZ) == material) {
  //      for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
  //        // Get direction
  //        const int* c = DESCRIPTOR::c[iPop];
  //        // Get next cell located in the current direction
  //        // Check if the next cell is a fluid node
  //        if (superGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
  //          int overlap = this->sLattice.getOverlap();
  //          // Get f_q of next fluid cell where l = opposite(q)
  //          T f = this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
  //          // Get f_l of the boundary cell
  //          // Add f_q and f_opp
  //          f += this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR >(iPop)];
  //          // Update force
  //          force[0] -= c[0]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //          force[1] -= c[1]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //          force[2] -= c[2]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //        }
  //      }
  //      force[0]=this->converter.physForce(force[0]);
  //      force[1]=this->converter.physForce(force[1]);
  //      force[2]=this->converter.physForce(force[2]);
  //      return force;
  //    }
  //    else {
  //      return force;
  //    }
  //  }
  //  else {
  //    return force;
  //  }
  return false;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::BlockLatticePhysCorrBoundaryForce2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,BlockGeometry2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physCorrBoundaryForce";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  std::vector<T> force(3, T());
  //  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
  //    int globX = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosX() + locix;
  //    int globY = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosY() + lociy;
  //    int globZ = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosZ() + lociz;
  //    if(BlockGeometry.getMaterial(globX, globY, globZ) == material) {
  //      for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
  //        // Get direction
  //        const int* c = DESCRIPTOR::c[iPop];
  //        // Get next cell located in the current direction
  //        // Check if the next cell is a fluid node
  //        if (BlockGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
  //          int overlap = this->_blockLattice.getOverlap();
  //          // Get f_q of next fluid cell where l = opposite(q)
  //          T f = this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
  //          // Get f_l of the boundary cell
  //          // Add f_q and f_opp
  //          f += this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR >(iPop)];
  //          // Update force
  //          force[0] -= c[0]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //          force[1] -= c[1]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //          force[2] -= c[2]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //        }
  //      }
  //      force[0]=this->_converter.physForce(force[0]);
  //      force[1]=this->_converter.physForce(force[1]);
  //      force[2]=this->_converter.physForce(force[2]);
  //      return force;
  //    }
  //    else {
  //      return force;
  //    }
  //  }
  //  else {
  //    return force;
  //  }
  return false;
}

}
#endif
