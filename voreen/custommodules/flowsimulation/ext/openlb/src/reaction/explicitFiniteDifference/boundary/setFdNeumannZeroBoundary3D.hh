/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz, Davide Dapelo
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
#ifndef SET_NO_PENETRATION_BOUNDARY_HH
#define SET_NO_PENETRATION_BOUNDARY_HH

#include "setFdNeumannZeroBoundary3D.h"

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the setFdNeumannZeroBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
void setFdNeumannZeroBoundary(SuperLattice<T,DESCRIPTOR>& sLattice, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model,
                                SuperGeometry<T,3>& superGeometry, int material)
{
  setFdNeumannZeroBoundary<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(sLattice, iT, model, superGeometry.getMaterialIndicator(material));
}

///Initialising the setFdNeumannZeroBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
void setFdNeumannZeroBoundary(SuperLattice<T,DESCRIPTOR>& sLattice, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model,
                                FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setFdNeumannZeroBoundary");
  int _overlap = indicator->getSuperGeometry().getOverlap();
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setFdNeumannZeroBoundary<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(sLattice.getBlock(iC), iT, model, indicator->getBlockIndicatorF(iC),includeOuterCells);
    /// Adds needed Cells to the Communicator _commBC in SuperLattice
    //the addPoints2CommBC function is initialised inside setLocalVelocityBoundary3D.h/hh
    addPoints2CommBC(sLattice,std::forward<decltype(indicator)>(indicator), _overlap);
  }
}

////////// BlockLattice Domain  /////////////////////////////////////////

//set FdNeumannZeroBoundary on indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
void setFdNeumannZeroBoundary(BlockLattice<T,DESCRIPTOR>& block, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model,
                                BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      std::vector<int> discreteNormal(4,0);
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY, iZ);
      PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor;
      postProcessor = new FdNeumannZeroBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(
        iX,iX,iY,iY,iZ,iZ,iT,-discreteNormal[1],-discreteNormal[2],-discreteNormal[3],model);
      block.addPostProcessor(*postProcessor);
    }
  });
}

}
#endif
