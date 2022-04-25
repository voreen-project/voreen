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

#ifndef LATTICE_POROUS_MOMENTUM_LOSS_FORCE_2D_HH
#define LATTICE_POROUS_MOMENTUM_LOSS_FORCE_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePorousMomentumLossForce2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {
/*
template <typename T, typename DESCRIPTOR>
SuperLatticePorousMomentumLossForce2D<T,DESCRIPTOR>::SuperLatticePorousMomentumLossForce2D
(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
 std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,4*indicator.size())
{
  this->getName() = "physPorousMomentumLossForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePorousMomentumLossForce2D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC), superGeometry.getBlockGeometry(iC), indicator, converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePorousMomentumLossForce2D<T,DESCRIPTOR>::operator() (T output[],
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
BlockLatticePorousMomentumLossForce2D<T, DESCRIPTOR>::BlockLatticePorousMomentumLossForce2D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, BlockGeometry<T,2>& blockGeometry,
  std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T, DESCRIPTOR>(blockLattice, converter, 4*indicator.size()), _blockGeometry(blockGeometry), _vectorOfIndicator(indicator)
{
  this->getName() = "physPorousMomentumLossForce";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePorousMomentumLossForce2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  // iterate over all particles in _indicator
  for (size_t iInd=0; iInd!=_vectorOfIndicator.size(); iInd++) {

    int numVoxels = 0;
    std::vector<int> start{0,0};
    std::vector<int> end{0,0};
    T invDeltaX = 1./this->_converter.getPhysDeltaX();

    if (getRangeBlockGeometrySmoothIndicatorIntersection2D(_blockGeometry, *(_vectorOfIndicator[iInd]), invDeltaX, start, end)) {

      // iterate over cells in the constructed intersection box
      for (int iX = start[0]; iX < end[0]; iX++) {
        for (int iY = start[1]; iY < end[1]; iY++) {
          // check if cell belongs to particle
          T inside[1] = {0.};
          T posIn[2] = {0.};
          _blockGeometry.getPhysR(posIn, {iX, iY});
          (*(_vectorOfIndicator[iInd]))( inside, posIn);
          if ( !util::nearZero(inside[0]) && this->_blockGeometry.get({iX,iY})==1) {
            // compute momentum exchange force on particle
            T tmpForce[2] = {0.,0.};
            tmpForce[0] += this->_blockLattice.get(iX, iY).template getFieldComponent<descriptors::VELOCITY_NUMERATOR>(0);
            tmpForce[1] += this->_blockLattice.get(iX, iY).template getFieldComponent<descriptors::VELOCITY_NUMERATOR>(1);
            // reset external field for next timestep
            T reset_to_zero[2] = {0.,0.};
            this->_blockLattice.get(iX, iY).template setField<descriptors::VELOCITY_NUMERATOR>(reset_to_zero);
            // convert force to SI units and compute torque
            numVoxels++;
            // division bei length of lattice cell necessary due to converter handling of force
            tmpForce[0] = this->_converter.getPhysForce(tmpForce[0])*invDeltaX;
            tmpForce[1] = this->_converter.getPhysForce(tmpForce[1])*invDeltaX;
            output[0+iInd*4] += tmpForce[0];
            output[1+iInd*4] += tmpForce[1];
            output[2+iInd*4] += (posIn[0]-_vectorOfIndicator[iInd]->getPos()[0])*tmpForce[1] - (posIn[1]-_vectorOfIndicator[iInd]->getPos()[1])*tmpForce[0];
          }
        }
      }

    }
    output[3+iInd*4] = numVoxels;

  }
  return true;
}
*/
}
#endif
