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

#ifndef LATTICE_INDICATOR_SMOOTH_INDICATOR_INTERSECTION_2D_HH
#define LATTICE_INDICATOR_SMOOTH_INDICATOR_INTERSECTION_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticeIndicatorSmoothIndicatorIntersection2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template <typename T, typename DESCRIPTOR, bool HLBM>
SuperLatticeIndicatorSmoothIndicatorIntersection2D<T,DESCRIPTOR,HLBM>::SuperLatticeIndicatorSmoothIndicatorIntersection2D (
  SuperLattice<T,DESCRIPTOR>& sLattice,
  SuperGeometry<T,2>& superGeometry,
  IndicatorF2D<T>& normalInd, SmoothIndicatorF2D<T,T,HLBM>& smoothInd )
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "Indicator-SmoothIndicator Intersection";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeIndicatorSmoothIndicatorIntersection2D<T,DESCRIPTOR,HLBM>(this->_sLattice.getBlock(iC), superGeometry.getBlockGeometry(iC), normalInd, smoothInd));
  }
}

template <typename T, typename DESCRIPTOR, bool HLBM>
bool SuperLatticeIndicatorSmoothIndicatorIntersection2D<T,DESCRIPTOR,HLBM>::operator() (T output[], const int input[])
{
  output[0] = 0.;
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    int globiC = this->_sLattice.getLoadBalancer().glob(iC);
    if ( this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank() ) {
      this->getBlockF(iC)(output,&input[1]);
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(output[0], MPI_MAX);
#endif
  return true;

}

template<typename T, typename DESCRIPTOR, bool HLBM>
BlockLatticeIndicatorSmoothIndicatorIntersection2D<T,DESCRIPTOR,HLBM>::BlockLatticeIndicatorSmoothIndicatorIntersection2D (
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  BlockGeometry<T,2>& blockGeometry,
  IndicatorF2D<T>& normalInd, SmoothIndicatorF2D<T,T,HLBM>& smoothInd )
  : BlockLatticeF2D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _normalInd(normalInd), _smoothInd(smoothInd)
{
  this->getName() = "Indicator-SmoothIndicator Intersection";
}

template<typename T, typename DESCRIPTOR, bool HLBM>
bool BlockLatticeIndicatorSmoothIndicatorIntersection2D<T, DESCRIPTOR,HLBM>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  int start[2] = {0};
  int end[2] = {0};
  // check for intersection of cuboid and smoothIndicator
  Cuboid2D<T> tmpCuboid(_blockGeometry.getOrigin()[0], _blockGeometry.getOrigin()[1], _blockGeometry.getDeltaR(), _blockGeometry.getNx(), _blockGeometry.getNy());
  T posXmin = _smoothInd.getPos()[0] - _smoothInd.getCircumRadius();
  T posXmax = _smoothInd.getPos()[0] + _smoothInd.getCircumRadius();
  T posYmin = _smoothInd.getPos()[1] - _smoothInd.getCircumRadius();
  T posYmax = _smoothInd.getPos()[1] + _smoothInd.getCircumRadius();
  if (tmpCuboid.checkInters(posXmin, posXmax, posYmin, posYmax, start[0], end[0], start[1], end[1])) {

    for (int k=0; k<2; k++) {
      start[k] -= 1;
      if (start[k] < 0) {
        start[k] = 0;
      }
      end[k] += 2;
      if (end[k] > _blockGeometry.getExtent()[k]) {
        end[k] = _blockGeometry.getExtent()[k];
      }
    }

    // iterate over cells in the constructed intersection box
    for (int iX = start[0]; iX < end[0]; iX++) {
      for (int iY = start[1]; iY < end[1]; iY++) {

        // check if cell belongs to particle
        T insideT[1] = {0.};
        T posIn[2] = {0.};
        _blockGeometry.getPhysR(posIn, {iX, iY});
        _smoothInd( insideT, posIn);
        if ( !util::nearZero(insideT[0]) && this->_blockGeometry.get({iX,iY})==1) {
          // Return 1 if at least one cell is found to be inside both A and B
          bool insideBool[1] = {false};
          _normalInd(insideBool, posIn);
          if (insideBool[0]) {
            output[0] = 1.;
            return true;
          }
        }
      }
    }
  }

  return true;
}

}
#endif
