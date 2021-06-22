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
#ifndef SET_BOUNCE_BACK_VELOCITY_BOUNDARY_2D_HH
#define SET_BOUNCE_BACK_VELOCITY_BOUNDARY_2D_HH

#include "setBounceBackVelocityBoundary2D.h"


namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the BounceBackVelocityBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry, int material, IndicatorF2D<T>& geometryIndicator,
                                   std::vector<int> bulkMaterials)
{

  setBounceBackVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, superGeometry.getMaterialIndicator(material),
      superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
      geometryIndicator);
}

///Initialising the BounceBackVelocityBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                                   FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                                   IndicatorF2D<T>& geometryIndicator)
{
  //out of superOffBoundary2D
  T _epsFraction = 0.0001;
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */
  int _overlap = 1;

  OstreamManager clout(std::cout, "setBounceBackVelocityBoundary");
  clout << "epsFraction=" << _epsFraction << std::endl;
  clout.setMultiOutput(true);
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    clout << "Cuboid globiC " << sLattice.getLoadBalancer().glob(iCloc)
          << " starts to read distances for Velocity Boundary..." << std::endl;
    setBounceBackVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getExtendedBlockLattice(iCloc),boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
        bulkIndicator->getExtendedBlockIndicatorF(iCloc),
        geometryIndicator);
    clout << "Cuboid globiC " << sLattice.getLoadBalancer().glob(iCloc)
          << " finished reading distances for Velocity Boundary." << std::endl;
  }
  clout.setMultiOutput(false);
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(boundaryIndicator)>(boundaryIndicator), _overlap);
}


////////// BlockLattice Domain  /////////////////////////////////////////


template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator)
{
  if ( !boundaryIndicator.isEmpty() ) {
    /*
     *Range of cells to be traversed
     **/
    const Vector<int,2> min = boundaryIndicator.getMin();
    const Vector<int,2> max = boundaryIndicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        if (boundaryIndicator(iX, iY)) {
          setBounceBackVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(block, bulkIndicator.getBlockGeometryStructure(), iX, iY,
              bulkIndicator, geometryIndicator);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator)
{
  T _epsFraction = 0.0001;
  OstreamManager clout(std::cout, "setBounceBackVelocityBoundary");
  T distances[DESCRIPTOR::q];
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    const int iXn = iX + descriptors::c<DESCRIPTOR >(iPop,0);
    const int iYn = iY + descriptors::c<DESCRIPTOR >(iPop,1);
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
    } // bulk indicator
  } // iPop loop
  setBounceBackVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(block, blockGeometryStructure, iX, iY, distances);

}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[DESCRIPTOR::q])
{
  typedef DESCRIPTOR L;
  T location[DESCRIPTOR::d];
  blockGeometryStructure.getPhysR(location, x,y);

  T distancesCopy[L::q];
  T spacing = blockGeometryStructure.getDeltaR();
  for (int iPop = 1; iPop < L::q ; ++iPop) {
    distancesCopy[iPop] = spacing*(1.-distances[iPop]);
    if ( !util::nearZero(distances[iPop]+1) ) {
      distancesCopy[iPop] = -1;
    }
  }
  //set dynamics on indicated cells
  //function used for bouzidi and bounceBackVelocityBoundary
  setOffDynamics<T,DESCRIPTOR,MixinDynamics>(block, x, y, location, distancesCopy);

  for (int iPop = 1; iPop < L::q ; ++iPop) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      setBounceBackVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(block, blockGeometryStructure, x-descriptors::c<L>(iPop,0), y-descriptors::c<L>(iPop,1), iPop, distances[iPop]);
    }
  }

}



template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist)
{
  if (blockGeometryStructure.getMaterial(x-descriptors::c<DESCRIPTOR >(iPop,0), y-descriptors::c<DESCRIPTOR >(iPop,1)) != 1) {
    /// Using Bouzidi BC OnePoint. This corresponds to Bounce Back and TwoPoint to linear interpolation
    PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = new VelocityBounceBackPostProcessorGenerator2D
    <T, DESCRIPTOR>(x, y, iPop, dist);
    if (postProcessor) {
      block.addPostProcessor(*postProcessor);
    }

  }
  else {
    ///Using Bouzidi BC TwoPoint. This corresponds to Linear Interpolation
    PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =new VelocityBouzidiLinearPostProcessorGenerator2D
    <T, DESCRIPTOR>(x, y, iPop, dist);
    if (postProcessor) {
      block.addPostProcessor(*postProcessor);
    }
  }
}



}//namespace olb

#endif
