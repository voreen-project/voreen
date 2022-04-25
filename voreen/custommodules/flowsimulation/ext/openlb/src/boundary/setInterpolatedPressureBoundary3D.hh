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

//This file contains the Interpolated Pressure Boundary
//This is an onLattice Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_INTERPOLATED_PRESSURE_BOUNDARY_HH
#define SET_INTERPOLATED_PRESSURE_BOUNDARY_HH

#include "setInterpolatedPressureBoundary3D.h"

namespace olb {

///Initialising the setInterpolatedPressureBoundary function on the superLattice domain
template<typename T,typename DESCRIPTOR, class MixinDynamics>
void setInterpolatedPressureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry<T,3>& superGeometry,
                                     int material)
{
  setInterpolatedPressureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));
}

///Initialising the setInterpolatedPressureBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setInterpolatedPressureBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega,FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setInterpolatedPressureBoundary");
  int _overlap = 1;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  //clout << sLattice.getLoadBalancer().size() <<"sLattice.getLoadBalancer.size()" << std::endl;
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setInterpolatedPressureBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iC), omega, indicator->getBlockIndicatorF(iC),
        includeOuterCells);
  }
  //defined in setLocalVelocityBoundary
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);

}


/// Add interpolated pressure boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setInterpolatedPressureBoundary(BlockLattice<T,DESCRIPTOR>& _block, T omega, BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      Dynamics<T,DESCRIPTOR>* dynamics = nullptr;
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);

      if (discreteNormal[0] == 0) {
        if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
          dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletPressureBoundaryTuple<0,-1>
          >>();
          _block.addPostProcessor(
            typeid(PostStream), {iX, iY, iZ},
            meta::id<PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,0,-1>>());
        }
        else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
          dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletPressureBoundaryTuple<0,1>
          >>();
          _block.addPostProcessor(
            typeid(PostStream), {iX, iY, iZ},
            meta::id<PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,0,1>>());
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
          dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletPressureBoundaryTuple<1,-1>
          >>();
          _block.addPostProcessor(
            typeid(PostStream), {iX, iY, iZ},
            meta::id<PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,1,-1>>());
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
          dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletPressureBoundaryTuple<1,1>
          >>();
          _block.addPostProcessor(
            typeid(PostStream), {iX, iY, iZ},
            meta::id<PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,1,1>>());
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
          dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletPressureBoundaryTuple<2,-1>
          >>();
          _block.addPostProcessor(
            typeid(PostStream), {iX, iY, iZ},
            meta::id<PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,2,-1>>());
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
          dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletPressureBoundaryTuple<2,1>
          >>();
          _block.addPostProcessor(
            typeid(PostStream), {iX, iY, iZ},
            meta::id<PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,2,1>>());
        }
        else{
          std::cout << "Discrete Normal is not valid for interpolated pressure boundary!" << std::endl;
          singleton::exit(1);
        }
      }
      else {
        // When this error occurs, check your geometry for stray cells or try slightly variating your resolution
        std::cout << "Discrete Normal is not valid for interpolated pressure boundary!" << std::endl;
        singleton::exit(1);
      }
      if (dynamics) {
        dynamics->getParameters(_block).template set<descriptors::OMEGA>(omega);
      }
      setBoundary(_block, iX, iY, iZ, dynamics);
    }
  });
}


}

#endif
