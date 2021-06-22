/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

///This file contains the Bouzidi Zero Velocity Boundary
///This is an offLattice Boundary
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_BOUZIDI_ZERO_VELOCITY_BOUNDARY_2D_HH
#define SET_BOUZIDI_ZERO_VELOCITY_BOUNDARY_2D_HH

#include "setBouzidiZeroVelocityBoundary2D.h"


namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the BouzidiZeroVelocityBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiZeroVelocityBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry, int material,
                                    IndicatorF2D<T>& geometryIndicator,
                                    std::vector<int> bulkMaterials)
{
  setBouzidiZeroVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, superGeometry.getMaterialIndicator(material),
      superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
      geometryIndicator);


}

///Initialising the BouzidiZeroVelocityBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiZeroVelocityBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                                    FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                                    IndicatorF2D<T>&                   geometryIndicator)
{
  //out of superOffBoundary2D
  T _epsFraction = 0.0001;
  OstreamManager clout(std::cout, "setBouzidiZeroVelocityBoundary");
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */
  int _overlap = 1;

  clout << "epsFraction=" << _epsFraction << std::endl;
  clout.setMultiOutput(true);
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    clout << "Cuboid globiC " << sLattice.getLoadBalancer().glob(iCloc)
          << " starts to read distances for ZeroVelocity Boundary..." << std::endl;
    setBouzidiZeroVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getExtendedBlockLattice(iCloc),
        boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
        bulkIndicator->getExtendedBlockIndicatorF(iCloc),
        geometryIndicator);
    clout << "Cuboid globiC " << sLattice.getLoadBalancer().glob(iCloc)
          << " finished reading distances for ZeroVelocity Boundary." << std::endl;
  }
  clout.setMultiOutput(false);
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(boundaryIndicator)>(boundaryIndicator), _overlap);
}


////////// BlockLattice Domain  /////////////////////////////////////////
//the functions below set the boundary on indicated cells inside the block domain

template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiZeroVelocityBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator)
{
  if ( !boundaryIndicator.isEmpty() ) {
    const Vector<int,2> min = boundaryIndicator.getMin();
    const Vector<int,2> max = boundaryIndicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        if (boundaryIndicator(iX, iY)) {
          setBouzidiZeroVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(block, bulkIndicator.getBlockGeometryStructure(), iX, iY,
              bulkIndicator, geometryIndicator);
        }
      }
    }
  }
}


template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiZeroVelocityBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator)
{
  OstreamManager clout(std::cout, "setBouzidiZeroVelocityBoundary");
  T _epsFraction = 0.0001;
  T distances[DESCRIPTOR::q];
  std::fill(std::begin(distances), std::end(distances), -1);

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    const Vector<int,2> c = descriptors::c<DESCRIPTOR>(iPop);
    const int iXn = iX + c[0];
    const int iYn = iY + c[1];
    if (blockGeometryStructure.isInside(iXn,iYn)) {
      if (bulkIndicator(iXn,iYn)) {
        T dist = -1;
        T physR[2];
        blockGeometryStructure.getPhysR(physR,iXn,iYn);
        T voxelSize=blockGeometryStructure.getDeltaR();
        Vector<T,2> physC(physR);
        
        Vector<T,2> direction(-voxelSize*descriptors::c<DESCRIPTOR >(iPop,0),-voxelSize*descriptors::c<DESCRIPTOR >(iPop,1));
        T cPhysNorm = voxelSize*sqrt(descriptors::c<DESCRIPTOR >(iPop,0)*descriptors::c<DESCRIPTOR >(iPop,0)+descriptors::c<DESCRIPTOR >(iPop,1)*descriptors::c<DESCRIPTOR >(iPop,1));
        
        if (!geometryIndicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
          T epsX = voxelSize*descriptors::c<DESCRIPTOR >(iPop,0)*_epsFraction;
          T epsY = voxelSize*descriptors::c<DESCRIPTOR >(iPop,1)*_epsFraction;
        
          Vector<T,2> physC2(physC);
          physC2[0] += epsX;
          physC2[1] += epsY;
          Vector<T,2> direction2(direction);
          direction2[0] -= 2.*epsX;
          direction2[1] -= 2.*epsY;
        
          if ( !geometryIndicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
            clout << "ERROR: no boundary found at (" << iXn << "," << iYn <<") ~ ("
                  << physR[0] << "," << physR[1] << "), "
                  << "in direction " << util::opposite<DESCRIPTOR >(iPop)
                  << std::endl;
          }
          T distNew = (dist - sqrt(epsX*epsX+epsY*epsY))/cPhysNorm;
          if (distNew < 0.5) {
            dist = 0;
          }
          else {
            dist = 0.5 * cPhysNorm;
            clout << "WARNING: distance at (" << iXn << "," << iYn <<") ~ ("
                  << physR[0] << "," << physR[1] <<"), "
                  << "in direction " << util::opposite<DESCRIPTOR >(iPop) << ": "
                  << distNew
                  << " rounded to "
                  << dist/cPhysNorm
                  << std::endl;
          }
        }
        distances[util::opposite<DESCRIPTOR >(iPop)] = dist/cPhysNorm;
      } else {
        if (blockGeometryStructure.getMaterial(iXn,iYn) != 0) {
          auto postProcessor = new ZeroVelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>(iXn, iYn, util::opposite<DESCRIPTOR>(iPop), 0);
          block.addPostProcessor(*postProcessor);
        }
      }
    } // bulk indicator
  } // iPop loop
  setBouzidiZeroVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(block, blockGeometryStructure, iX, iY, distances);
}


template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiZeroVelocityBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[DESCRIPTOR::q])
{
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      const Vector<int,2> c = descriptors::c<DESCRIPTOR>(iPop);
      setBouzidiZeroVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(block, blockGeometryStructure, x-c[0], y-c[1], iPop, distances[iPop]);
    }
  }
}

//set postProcessor on indicated cells
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBouzidiZeroVelocityBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist)
{
  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = nullptr;
  const Vector<int,2> c = descriptors::c<DESCRIPTOR>(iPop);

  if (blockGeometryStructure.getMaterial(x-c[0], y-c[1]) != 1) {
    postProcessor = new ZeroVelocityBounceBackPostProcessorGenerator2D<T,DESCRIPTOR>(x, y, iPop, dist);
  }
  else {
    postProcessor = new ZeroVelocityBouzidiLinearPostProcessorGenerator2D<T,DESCRIPTOR>(x, y, iPop, dist);
  }

  if (postProcessor) {
    block.addPostProcessor(*postProcessor);
  }
}

}//namespace olb

#endif
