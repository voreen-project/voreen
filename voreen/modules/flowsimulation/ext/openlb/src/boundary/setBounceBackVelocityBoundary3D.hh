/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Max Gaedtke
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

#ifndef SET_BOUNCE_BACK_VELOCITY_BOUNDARY_3D_HH
#define SET_BOUNCE_BACK_VELOCITY_BOUNDARY_3D_HH

#include "setBounceBackVelocityBoundary3D.h"

namespace olb {


///Initialising the setLocalVelocityBoundary function on the superLattice domain
template<typename T,typename DESCRIPTOR>
void setBounceBackVelocityBoundary(SuperGeometry<T,3>& superGeometry, int material, T omega, SuperLattice<T, DESCRIPTOR>& sLattice)
{
  setBounceBackVelocityBoundary<T,DESCRIPTOR>(superGeometry.getMaterialIndicator(material),omega,sLattice);
}

template<typename T, typename DESCRIPTOR>
void setBounceBackVelocityBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega,SuperLattice<T, DESCRIPTOR>& sLattice)
{
  OstreamManager clout(std::cout, "BounceBackVelocityBoundary");
  
  int _overlap = 1;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  // addPoints2CommBC(sLattice,std::forward<decltype(indicator)>(indicator), _overlap);
  clout << sLattice.getLoadBalancer().size() <<"sLattice.getLoadBalancer.size()" << std::endl;
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setBounceBackVelocityBoundary<T,DESCRIPTOR>(indicator->getBlockIndicatorF(iC),omega,includeOuterCells,sLattice.getBlock(iC));
  }
}



template<typename T, typename DESCRIPTOR>
void setBounceBackVelocityBoundary(BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells,BlockLattice<T,DESCRIPTOR>& _block)
{
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4,0);
  T default_rho = 1.0;
  T default_u[] = {0.0, 0.0, 0.0};
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      Dynamics<T,DESCRIPTOR>* dynamics = nullptr;
      PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor = nullptr;
      dynamics = new BounceBackVelocity<T,DESCRIPTOR>(default_rho, default_u);

      setBoundary(_block, iX,iY,iZ, dynamics, postProcessor);
    }
  });
}

}//namespace olb
#endif
