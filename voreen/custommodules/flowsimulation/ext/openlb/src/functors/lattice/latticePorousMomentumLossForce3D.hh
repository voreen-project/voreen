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

#ifndef LATTICE_POROUS_MOMENTUM_LOSS_FORCE_3D_HH
#define LATTICE_POROUS_MOMENTUM_LOSS_FORCE_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePorousMomentumLossForce3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {
/*
template <typename T, typename DESCRIPTOR>
SuperLatticePorousMomentumLossForce3D<T,DESCRIPTOR>::SuperLatticePorousMomentumLossForce3D
(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry,
 std::vector<SmoothIndicatorF3D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,7*indicator.size())
{
  this->getName() = "physPorousMomentumLossForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePorousMomentumLossForce3D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC), superGeometry.getBlockGeometry(iC), indicator, converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePorousMomentumLossForce3D<T,DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  for (int i=0; i<this->getTargetDim(); i++) {
    output[i] = 0.;
  }
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    int globiC = this->_sLattice.getLoadBalancer().glob(iC);
    if ( this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank() ) {
      this->getBlockF(iC)(output,&input[1]);
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim(); ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
#endif
  return true;

}

template<typename T, typename DESCRIPTOR>
BlockLatticePorousMomentumLossForce3D<T, DESCRIPTOR>::BlockLatticePorousMomentumLossForce3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, BlockGeometry<T,3>& blockGeometry,
  std::vector<SmoothIndicatorF3D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 7*indicator.size()), _blockGeometry(blockGeometry), _vectorOfIndicator(indicator)
{
  this->getName() = "physPorousMomentumLossForce";
}


template<typename T, typename DESCRIPTOR>
bool BlockLatticePorousMomentumLossForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  // iterate over all particles in _indicator
  for (size_t iInd=0; iInd!=_vectorOfIndicator.size(); iInd++) {

    int numVoxels = 0;
    std::vector<int> start{0,0,0};
    std::vector<int> end{0,0,0};
    T invDeltaX = 1./this->_converter.getPhysDeltaX();

    // check for intersection of cuboid and indicator
    if (getRangeBlockGeometrySmoothIndicatorIntersection3D(_blockGeometry, *(_vectorOfIndicator[iInd]), invDeltaX, start, end)) {

      // iterate over cells in the constructed intersection box
      for (int iX = start[0]; iX < end[0]; iX++) {
        for (int iY = start[1]; iY < end[1]; iY++) {
          for (int iZ = start[2]; iZ < end[2]; iZ++) {

            // check if cell belongs to particle
            T inside[1] = {0.};
            T posIn[3] = {0.};
            _blockGeometry.getPhysR(posIn, {iX, iY, iZ});
            (*(_vectorOfIndicator[iInd]))( inside, posIn);
            if ( !util::nearZero(inside[0]) && this->_blockGeometry.get({iX,iY,iZ})==1) {
              // compute momentum exchange force on particle
              T tmpForce[3] = {0.,0.,0.};
              tmpForce[0] += this->_blockLattice.get(iX, iY, iZ).template getFieldComponent<descriptors::VELOCITY_NUMERATOR>(0);
              tmpForce[1] += this->_blockLattice.get(iX, iY, iZ).template getFieldComponent<descriptors::VELOCITY_NUMERATOR>(1);
              tmpForce[2] += this->_blockLattice.get(iX, iY, iZ).template getFieldComponent<descriptors::VELOCITY_NUMERATOR>(2);
              // reset external field for next timestep
              T reset_to_zero[3] = {0.,0.,0.};
              this->_blockLattice.get(iX, iY, iZ).template setField<descriptors::VELOCITY_NUMERATOR>(reset_to_zero);
              // convert force to SI units and compute torque
              numVoxels++;
              // division bei length of lattice cell necessary due to converter handling of force
              tmpForce[0] = this->_converter.getPhysForce(tmpForce[0]);
              tmpForce[1] = this->_converter.getPhysForce(tmpForce[1]);
              tmpForce[2] = this->_converter.getPhysForce(tmpForce[2]);
              output[0+iInd*7] += tmpForce[0];
              output[1+iInd*7] += tmpForce[1];
              output[2+iInd*7] += tmpForce[2];
              output[3+iInd*7] += (posIn[1]-_vectorOfIndicator[iInd]->getPos()[1])*tmpForce[2]
                                  - (posIn[2]-_vectorOfIndicator[iInd]->getPos()[2])*tmpForce[1];
              output[4+iInd*7] += (posIn[2]-_vectorOfIndicator[iInd]->getPos()[2])*tmpForce[0]
                                  - (posIn[0]-_vectorOfIndicator[iInd]->getPos()[0])*tmpForce[2];
              output[5+iInd*7] += (posIn[0]-_vectorOfIndicator[iInd]->getPos()[0])*tmpForce[1]
                                  - (posIn[1]-_vectorOfIndicator[iInd]->getPos()[1])*tmpForce[0];
            }
          }
        }
      }

    }
    output[6+iInd*7] = numVoxels;

  }
  return true;
}
*/
}
#endif
