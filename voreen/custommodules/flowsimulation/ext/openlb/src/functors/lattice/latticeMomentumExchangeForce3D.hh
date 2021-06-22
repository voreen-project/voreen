/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Robin Trunk
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

#ifndef LATTICE_MOMENTUM_EXCHANGE_FORCE_3D_HH
#define LATTICE_MOMENTUM_EXCHANGE_FORCE_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePorousMomentumLossForce3D.h"
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

template <typename T, typename DESCRIPTOR>
SuperLatticeMomentumExchangeForce3D<T,DESCRIPTOR>::SuperLatticeMomentumExchangeForce3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
 std::vector<SmoothIndicatorF3D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter,
 Vector<bool,3> periodicity)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,7*indicator.size())
{
  this->getName() = "physMomentumExchangeForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  Vector<T,3> min = superGeometry.getStatistics().getMinPhysR( 1 );
  Vector<T,3> max = superGeometry.getStatistics().getMaxPhysR( 1 );
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeMomentumExchangeForce3D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), superGeometry.getBlockGeometry(iC), indicator, converter, min, max, periodicity));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeMomentumExchangeForce3D<T,DESCRIPTOR>::operator() (T output[],
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
BlockLatticeMomentumExchangeForce3D<T, DESCRIPTOR>::BlockLatticeMomentumExchangeForce3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>& blockGeometry,
  std::vector<SmoothIndicatorF3D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter,
  Vector<T,3>cellMin, Vector<T,3> cellMax,
  Vector<bool,3> periodic)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 7*indicator.size()), _blockGeometry(blockGeometry), _vectorOfIndicator(indicator), _cellMin(cellMin), _cellMax(cellMax), _periodic(periodic)
{
  this->getName() = "physMomentumExchangeForce";
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeMomentumExchangeForce3D<T, DESCRIPTOR>::evaluate(T output[], SmoothIndicatorF3D<T,T,true>& sIndicator, int iInd)
{
  int numVoxels = 0;
  std::vector<int> start{0,0,0};
  std::vector<int> end{0,0,0};
  T invDeltaX = 1./this->_converter.getPhysDeltaX();

  // check for intersection of cuboid and indicator
  if (getRangeBlockGeometrySmoothIndicatorIntersection3D(_blockGeometry, sIndicator, invDeltaX, start, end)) {

    // iterate over cells in the constructed intersection box
    for (int iX = start[0]; iX < end[0]; iX++) {
      for (int iY = start[1]; iY < end[1]; iY++) {
        for (int iZ = start[2]; iZ < end[2]; iZ++) {

          // check if cell belongs to particle
          T inside[1] = {0.};
          T posIn[3] = {0.};
          _blockGeometry.getPhysR(posIn, iX, iY, iZ);
          sIndicator(inside, posIn);
          T tmpForce[3] = {0.,0.,0.};
          if ( !util::nearZero(inside[0]) && this->_blockGeometry.get(iX,iY,iZ)==1) {

            for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
              const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
              // Get next cell located in the current direction
              T inside2[1] = {0.};
              T posOut[3] = {0.};
              _blockGeometry.getPhysR(posOut, iX+c[0], iY+c[1], iZ+c[2]);
              sIndicator(inside2, posOut);

              // if not both cells are in the full solid domain calculate force
              if ( !(inside[0]==1 && inside2[0]==1) ) {
                T f1 = this->_blockLattice.get( iX+c[0], iY+c[1], iZ+c[2])[iPop];
                T f2 = this->_blockLattice.get( iX, iY, iZ)[util::opposite<DESCRIPTOR>(iPop)];
                // Update force
#ifdef FEATURE_HLBM_MEA_LADD
                tmpForce[0] -= this->_converter.getPhysForce( (f1+f2) * c[0] );
                tmpForce[1] -= this->_converter.getPhysForce( (f1+f2) * c[1] );
                tmpForce[2] -= this->_converter.getPhysForce( (f1+f2) * c[2] );
#else
                T pVel[3] = {0., 0., 0.};
                for (int iDim=0; iDim<3; iDim++) {
                  pVel[iDim] = this->_blockLattice.get(iX, iY, iZ).template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[iDim];
                }
                tmpForce[0] -= this->_converter.getPhysForce( f1*(c[0]-pVel[0]) + f2*(c[0]+pVel[0]) );
                tmpForce[1] -= this->_converter.getPhysForce( f1*(c[1]-pVel[1]) + f2*(c[1]+pVel[1]) );
                tmpForce[2] -= this->_converter.getPhysForce( f1*(c[2]-pVel[2]) + f2*(c[2]+pVel[2]) );
#endif
              }
            }

            // reset external field for next timestep
            T reset_to_zero[3] = {0.,0.,0.};
            this->_blockLattice.get(iX, iY, iZ).template setField<descriptors::VELOCITY_NUMERATOR>(reset_to_zero);
            // count cells of considered object
            numVoxels++;
            // output force and torque
            output[0+iInd*7] += tmpForce[0];
            output[1+iInd*7] += tmpForce[1];
            output[2+iInd*7] += tmpForce[2];
            output[3+iInd*7] += (posIn[1]-sIndicator.getPos()[1])*tmpForce[2]
                                - (posIn[2]-sIndicator.getPos()[2])*tmpForce[1];
            output[4+iInd*7] += (posIn[2]-sIndicator.getPos()[2])*tmpForce[0]
                                - (posIn[0]-sIndicator.getPos()[0])*tmpForce[2];
            output[5+iInd*7] += (posIn[0]-sIndicator.getPos()[0])*tmpForce[1]
                                - (posIn[1]-sIndicator.getPos()[1])*tmpForce[0];

          }

        }
      }
    }

  }

  output[6+iInd*7] = numVoxels;

}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeMomentumExchangeForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  // iterate over all particles in _indicator
  for (long unsigned int iInd=0; iInd!=_vectorOfIndicator.size(); iInd++) {

    if (_periodic[0]||_periodic[1]||_periodic[2]) {
      bool outOfGeometry = false;
      Vector<T,3> ghostPos = Vector<T,3> (0.,0.,0.);
      checkSmoothIndicatorOutOfGeometry(outOfGeometry, ghostPos, *(_vectorOfIndicator[iInd]), _cellMin, _cellMax, _periodic);

      Vector<T,3>particlePosition = _vectorOfIndicator[iInd]->getPos();
      // Sets the particle to the other domainside if it leaves the domain
      if ( particlePosition[0]<_cellMin[0] || particlePosition[0]>_cellMax[0] || particlePosition[1]<_cellMin[1] || particlePosition[1]>_cellMax[1]
           || particlePosition[2]<_cellMin[2] || particlePosition[2]>_cellMax[2] ) {
        _vectorOfIndicator[iInd]->setPos(ghostPos);
        ghostPos = particlePosition;
        particlePosition = _vectorOfIndicator[iInd]->getPos();
      }

      if (!outOfGeometry) {
        evaluate(output, *(_vectorOfIndicator[iInd]), iInd);
      } else {
        // Calculate force for the ghost particle
        _vectorOfIndicator[iInd]->setPos(ghostPos);
        evaluate(output, *(_vectorOfIndicator[iInd]), iInd);
        // Calculate force of actual particle
        _vectorOfIndicator[iInd]->setPos(particlePosition);
        evaluate(output, *(_vectorOfIndicator[iInd]), iInd);
      }
    } else {
      evaluate(output, *(_vectorOfIndicator[iInd]), iInd);
    }
  }
  return true;
}

}
#endif
