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
#ifndef SET_NO_PENETRATION_BOUNDARY_2D_HH
#define SET_NO_PENETRATION_BOUNDARY_2D_HH

#include "setFdNeumannZeroBoundary2D.h"

namespace olb {
///Initialising the setFdNeumannZeroBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
void setFdNeumannZeroBoundary(SuperLattice<T,DESCRIPTOR>& sLattice, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model,
                                SuperGeometry<T,2>& superGeometry, int material)
{
  setFdNeumannZeroBoundary<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(sLattice, iT, model, superGeometry.getMaterialIndicator(material));
}

///Initialising the setFdNeumannZeroBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
void setFdNeumannZeroBoundary(SuperLattice<T,DESCRIPTOR>& sLattice, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model,
                                FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  bool includeOuterCells = false;
  int _overlap = indicator->getSuperGeometry().getOverlap();
  OstreamManager clout(std::cout, "setOnBCInterpolatedBoundary");
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setFdNeumannZeroBoundary<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(sLattice.getBlock(iCloc), iT, model,
        indicator->getBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  //the addPoints2CommBC function is initialised inside setLocalVelocityBoundary2D.h/hh
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);


}
////////// BlockLattice Domain  /////////////////////////////////////////

/// Set Interpolated velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
void setFdNeumannZeroBoundary(BlockLattice<T,DESCRIPTOR>& block, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model,
                                BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setFdNeumannZeroBoundary");
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(3, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY}) >= margin
        && indicator(iX, iY)) {
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY);
      PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor;
      postProcessor = new FdNeumannZeroBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(
        iX,iX,iY,iY,iT,-discreteNormal[1],-discreteNormal[2],model);
      block.addPostProcessor(*postProcessor);
    }
  });
}


}//namespace olb

#endif
