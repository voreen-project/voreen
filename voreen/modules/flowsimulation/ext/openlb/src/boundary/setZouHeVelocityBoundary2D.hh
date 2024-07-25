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

//This file contains the ZouHe Velocity Boundary
//This is a new version of the Boundary, which only contains free floating functions

#ifndef SET_ZOUHE_VELOCITY_BOUNDARY_2D_HH
#define SET_ZOUHE_VELOCITY_BOUNDARY_2D_HH

#include "setZouHeVelocityBoundary2D.h"

namespace olb {
///Initialising the setZouHeVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setZouHeVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,T omega, SuperGeometry<T,2>& superGeometry, int material)
{
  setZouHeVelocityBoundary<T,DESCRIPTOR, MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));
}
///Initialising the setZouHeVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setZouHeVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setZouHeVelocityBoundary2D");
  int _overlap = 1;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setZouHeVelocityBoundary<T,DESCRIPTOR, MixinDynamics>(sLattice.getBlock(iCloc), omega,
        indicator->getBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}



/// Set ZouHe velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setZouHeVelocityBoundary(BlockLattice<T,DESCRIPTOR>& block, T omega, BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  using namespace boundaryhelper;
  OstreamManager clout(std::cout, "setZouHeVelocityBoundary2D");
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  /*
   *x0,x1,y0,y1 Range of cells to be traversed
   **/
  int x0 = margin;
  int y0 = margin;
  int x1 = blockGeometryStructure.getNx()-1 -margin;
  int y1 = blockGeometryStructure.getNy()-1 -margin;
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      Dynamics<T, DESCRIPTOR>* dynamics = nullptr;
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {
          //set the momenta,dynamics and post processor on the indicated ZouHe velocity boundary cells
          dynamics = block.getDynamics(DirectionOrientationMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
            ZouHeDynamics,MixinDynamics,momenta::BasicDirichletVelocityBoundaryTuple
          >::construct(Vector<int,2>(discreteNormal.data() + 1)));
        }
        else if (discreteNormal[0] == 1) {
          //sets the momenta, dynamics and post processors on indicated ZouHeVelocityCornerBoundary cells
          dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::FixedVelocityBoundaryTuple
          >>();
          block.addPostProcessor(
            typeid(PostStream), {iX, iY},
            promisePostProcessorForNormal<T,DESCRIPTOR,OuterVelocityCornerProcessor2D>(
              Vector<int,2>(discreteNormal.data() + 1)));
        }
        //sets momenta, dynamics and postProcessors on ZouHe InnerVelocityCornerBoundary Cells
        else if (discreteNormal[0] == 2) {
          dynamics = block.getDynamics(PlainMixinDynamicsForNormalMomenta<T,DESCRIPTOR,
            CombinedRLBdynamics,MixinDynamics,momenta::InnerCornerVelocityTuple2D
          >::construct(Vector<int,2>(discreteNormal.data() + 1)));
        }
        //sets the boundary on any indicated cell
        setBoundary(block, iX,iY, dynamics);
      }
    }
  }
}


}//namespace olb
#endif
