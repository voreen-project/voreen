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

//This file contains the ExternalFieldBoundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_EXT_FIELD_BOUNDARY_3D_HH
#define SET_EXT_FIELD_BOUNDARY_3D_HH

#include "setExtFieldBoundary3D.h"
namespace olb {

///Initialising the ExternalFieldBoundary on the superLattice domain
//advection diffusion boundary. Therefore mostly --> MixinDynamics = AdvectionDiffusionRLBdynamics<T,DESCRIPTOR>>
template <typename T, typename DESCRIPTOR, typename FIELD_A, typename FIELD_B>
void setExtFieldBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,SuperGeometry<T,3>& superGeometry, int material)
{
  setExtFieldBoundary<T,DESCRIPTOR,FIELD_A,FIELD_B>(sLattice,
      superGeometry.getMaterialIndicator(material));
}

///Initialising the ExternalFieldBoundary on the superLattice domain
template <typename T, typename DESCRIPTOR, typename FIELD_A, typename FIELD_B>
void setExtFieldBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  int _overlap = 1;
  OstreamManager clout(std::cout, "setExtFieldBoundary");
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setExtFieldBoundary<T,DESCRIPTOR,FIELD_A,FIELD_B>(
      sLattice.getBlock(iCloc),
      indicator->getBlockIndicatorF(iCloc),
      includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}


/// Set externalFieldBoundary for any indicated cells inside the block domain
template <typename T, typename DESCRIPTOR, typename FIELD_A, typename FIELD_B>
void setExtFieldBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator,
                         bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setExtFieldBoundary");
  const auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
      if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
        auto postProcessor = std::unique_ptr<PostProcessorGenerator3D<T, DESCRIPTOR>>{
        new ExtFieldBoundaryProcessorGenerator3D<T,DESCRIPTOR,FIELD_A,FIELD_B>(iX, iX, iY, iY, iZ, iZ, -discreteNormal[1], -discreteNormal[2], -discreteNormal[3]) };
        if (postProcessor) {
          _block.addPostProcessor(*postProcessor);
        }
      }
      else {
        clout << "Warning: Could not setExternalFieldBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] << "," << discreteNormal[3] << ")" << std::endl;
      }
    }
  });
}


}//namespace olb


#endif

