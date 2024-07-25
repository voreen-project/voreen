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

//This file contains the Interpolated Velocity Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_INTERPOLATED_VELOCITY_BOUNDARY_2D_HH
#define SET_INTERPOLATED_VELOCITY_BOUNDARY_2D_HH

#include "setInterpolatedVelocityBoundary2D.h"

namespace olb {
///Initialising the setInterpolatedVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setInterpolatedVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry<T,2>& superGeometry, int material)
{
  setInterpolatedVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));
}

///Initialising the setInterpolatedVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setInterpolatedVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  bool includeOuterCells = false;
  int _overlap = 1;
  OstreamManager clout(std::cout, "setInterpolatedVelocityBoundary");
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setInterpolatedVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iCloc), omega,
        indicator->getBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  //the addPoints2CommBC function is initialised inside setLocalVelocityBoundary2D.h/hh
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}



/// Set Interpolated velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setInterpolatedVelocityBoundary(BlockLattice<T,DESCRIPTOR>& block, T omega,BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  using namespace boundaryhelper;
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(3, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY}) >= margin
        && indicator(iX, iY)) {
      Dynamics<T, DESCRIPTOR>* dynamics = nullptr;
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY);
      if (discreteNormal[0] == 0) {
        dynamics = block.getDynamics(MixinDynamicsExchangeDirectionOrientationMomenta<T,DESCRIPTOR,
          MixinDynamics,momenta::BasicDirichletVelocityBoundaryTuple
        >::construct(Vector<int,2>(discreteNormal.data() + 1)));
        block.addPostProcessor(
          typeid(PostStream), {iX, iY},
          promisePostProcessorForDirectionOrientation<T,DESCRIPTOR,StraightFdBoundaryProcessor2D>(
            Vector<int,2>(discreteNormal.data() + 1)));
      }
      else if (discreteNormal[0] == 1) {
        dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
          momenta::FixedVelocityBoundaryTuple
        >>();
        block.addPostProcessor(
          typeid(PostStream), {iX, iY},
          promisePostProcessorForNormal<T,DESCRIPTOR,OuterVelocityCornerProcessor2D>(
            Vector<int,2>(discreteNormal.data() + 1)));
      }
      else if (discreteNormal[0] == 2) {
        dynamics = block.getDynamics(PlainMixinDynamicsForNormalMomenta<T,DESCRIPTOR,
          CombinedRLBdynamics,MixinDynamics,momenta::InnerCornerVelocityTuple2D
        >::construct(Vector<int,2>(discreteNormal.data() + 1)));
      }
      if (dynamics) {
        dynamics->getParameters(block).template set<descriptors::OMEGA>(omega);
      }
      setBoundary(block, iX, iY, dynamics);
    }
  });
}


}//namespace olb

#endif
