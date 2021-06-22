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

////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the setLocalVelocityBoundary function on the superLattice domain
template<typename T,typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega, SuperLattice3D<T, DESCRIPTOR>& sLattice)
{
  setBounceBackVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(superGeometry.getMaterialIndicator(material),omega,sLattice);
}
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega,SuperLattice3D<T, DESCRIPTOR>& sLattice)
{
  OstreamManager clout(std::cout, "BounceBackVelocityBoundary");
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */
  int _overlap = 1;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  // addPoints2CommBC(sLattice,std::forward<decltype(indicator)>(indicator), _overlap);
  clout << sLattice.getLoadBalancer().size() <<"sLattice.getLoadBalancer.size()" << std::endl;
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setBounceBackVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(indicator->getExtendedBlockIndicatorF(iC),omega,includeOuterCells,sLattice.getExtendedBlockLattice(iC));
  }
}

////////// BlockLattice Domain  /////////////////////////////////////////


template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells,BlockLatticeStructure3D<T,DESCRIPTOR>& _block)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  /*
   *x0,x1,y0,y1,z0,z1 Range of cells to be traversed
   **/
  int x0 = margin;
  int y0 = margin;
  int z0 = margin;
  int x1 = blockGeometryStructure.getNx()-1 -margin;
  int y1 = blockGeometryStructure.getNy()-1 -margin;
  int z1 = blockGeometryStructure.getNz()-1 -margin;
  std::vector<int> discreteNormal(4,0);
  T default_rho = 1.0;
  T default_u[] = {0.0, 0.0, 0.0};
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {
        Momenta<T,DESCRIPTOR>* momenta = nullptr;
        Dynamics<T,DESCRIPTOR>* dynamics = nullptr;
        PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor = nullptr;
        if (indicator(iX, iY, iZ)) {
          momenta = nullptr;
          dynamics = new BounceBackVelocity<T,DESCRIPTOR>(default_rho, default_u);
          postProcessor = nullptr;
          //Defined in setLocalVelocityBoundary3d
          setBoundary<T, DESCRIPTOR, MixinDynamics>(_block, iX,iY,iZ, omega, momenta, dynamics, postProcessor);

        }
      }
    }
  }

}

}
#endif
