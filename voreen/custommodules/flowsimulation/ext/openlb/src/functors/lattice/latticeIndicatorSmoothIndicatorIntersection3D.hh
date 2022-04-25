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

#ifndef LATTICE_INDICATOR_SMOOTH_INDICATOR_INTERSECTION_3D_HH
#define LATTICE_INDICATOR_SMOOTH_INDICATOR_INTERSECTION_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticeIndicatorSmoothIndicatorIntersection3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template <typename T, typename DESCRIPTOR, bool HLBM>
SuperLatticeIndicatorSmoothIndicatorIntersection3D<T,DESCRIPTOR,HLBM>::SuperLatticeIndicatorSmoothIndicatorIntersection3D (
  SuperLattice<T,DESCRIPTOR>& sLattice,
  SuperGeometry<T,3>& superGeometry,
  IndicatorF3D<T>& normalInd, SmoothIndicatorF3D<T,T,HLBM>& smoothInd )
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "Indicator-SmoothIndicator Intersection";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeIndicatorSmoothIndicatorIntersection3D<T,DESCRIPTOR,HLBM>(this->_sLattice.getBlock(iC), superGeometry.getBlockGeometry(iC), normalInd, smoothInd));
  }
}

template <typename T, typename DESCRIPTOR, bool HLBM>
bool SuperLatticeIndicatorSmoothIndicatorIntersection3D<T,DESCRIPTOR,HLBM>::operator() (T output[], const int input[])
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
BlockLatticeIndicatorSmoothIndicatorIntersection3D<T, DESCRIPTOR, HLBM>::BlockLatticeIndicatorSmoothIndicatorIntersection3D (
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  BlockGeometry<T,3>& blockGeometry,
  IndicatorF3D<T>& normalInd,
  SmoothIndicatorF3D<T,T,HLBM>& smoothInd )
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _normalInd(normalInd), _smoothInd(smoothInd)
{
  this->getName() = "Indicator-SmoothIndicator Intersection";
}

template<typename T, typename DESCRIPTOR, bool HLBM>
bool BlockLatticeIndicatorSmoothIndicatorIntersection3D<T, DESCRIPTOR, HLBM>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  int start[3] = {0};
  int end[3] = {0};
  // check for intersection of cuboid and smoothIndicator
  Cuboid3D<T> tmpCuboid(_blockGeometry.getOrigin()[0], _blockGeometry.getOrigin()[1], _blockGeometry.getOrigin()[2], _blockGeometry.getDeltaR(), _blockGeometry.getNx(), _blockGeometry.getNy(), _blockGeometry.getNz());
  T posXmin = _smoothInd.getPos()[0] - _smoothInd.getCircumRadius();
  T posXmax = _smoothInd.getPos()[0] + _smoothInd.getCircumRadius();
  T posYmin = _smoothInd.getPos()[1] - _smoothInd.getCircumRadius();
  T posYmax = _smoothInd.getPos()[1] + _smoothInd.getCircumRadius();
  T posZmin = _smoothInd.getPos()[2] - _smoothInd.getCircumRadius();
  T posZmax = _smoothInd.getPos()[2] + _smoothInd.getCircumRadius();
  if (tmpCuboid.checkInters(posXmin, posXmax, posYmin, posYmax, posZmin, posZmax, start[0],
                            end[0], start[1], end[1], start[2], end[2])) {

    for (int k=0; k<3; k++) {
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
        for (int iZ = start[2]; iZ < end[2]; iZ++) {

          // check if cell belongs to particle
          T insideT[1] = {0.};
          T posIn[3] = {0.};
          _blockGeometry.getPhysR(posIn, {iX, iY, iZ});
          _smoothInd( insideT, posIn);
          if ( !util::nearZero(insideT[0]) && this->_blockGeometry.get({iX,iY,iZ})==1) {
            // Return true if at least one cell is found to be inside both A and B
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
  }

  return true;
}

}
#endif
