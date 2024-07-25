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

//This file contains the Regularized Temperature Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_REGULARIZED_TEMPERATURE_BOUNDARY_2D_HH
#define SET_REGULARIZED_TEMPERATURE_BOUNDARY_2D_HH

#include "setRegularizedTemperatureBoundary2D.h"

namespace olb {
///Initialising the RegularizedTemperatureBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setRegularizedTemperatureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,T omega,
                                       SuperGeometry<T,2>& superGeometry, int material)
{
  setRegularizedTemperatureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));
}

///Initialising the RegularizedTemperatureBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setRegularizedTemperatureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega,
                                       FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  int _overlap = 1;
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setRegularizedTemperatureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iCloc),omega,
        indicator->getBlockIndicatorF(iCloc));
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T, DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}


/// Set RegularizedTemperatureBoundary for indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setRegularizedTemperatureBoundary(BlockLattice<T,DESCRIPTOR>& block, T omega, BlockIndicatorF2D<T>& indicator)
{
  using namespace boundaryhelper;
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  /*
   *x0,x1,y0,y1 Range of cells to be traversed
   **/
  int x0 = 0;
  int y0 = 0;
  int x1 = blockGeometryStructure.getNx()-1;
  int y1 = blockGeometryStructure.getNy()-1;
  std::vector<int> discreteNormal(3,0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      Dynamics<T,DESCRIPTOR>* dynamics = nullptr;
      PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor = nullptr;
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX,iY);
        if (discreteNormal[0] == 0) {   ////set momenta, dynamics and postProcessors on indicated RegularizedTemperatureBoundary cells
          dynamics = block.getDynamics(PlainMixinDynamicsForDirectionOrientationMomenta<T,DESCRIPTOR,
            CombinedAdvectionDiffusionRLBdynamics,
            MixinDynamics,momenta::RegularizedTemperatureBoundaryTuple
          >::construct(Vector<int,2>(discreteNormal.data() + 1)));
          dynamics->getParameters(block).template set<descriptors::OMEGA>(omega);
          setBoundary(block, iX,iY, dynamics, postProcessor);
        }
        else if (discreteNormal[0] == 1) {////set momenta, dynamics and postProcessors on indicated RegularizedTemperatureBoundary corner cells
          dynamics = block.getDynamics(NormalMixinDynamicsForPlainMomenta<T,DESCRIPTOR,
            AdvectionDiffusionCornerDynamics2D,
            MixinDynamics,momenta::EquilibriumBoundaryTuple
          >::construct(Vector<int,2>(discreteNormal.data() + 1)));
          dynamics->getParameters(block).template set<descriptors::OMEGA>(omega);
          setBoundary(block, iX,iY, dynamics, postProcessor);
        }
      }
    }
  }
}

}//namespace olb


#endif
