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

//This file contains the Interpolated Convection Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_INTERPOLATED_CONVECTION_BOUNDARY_2D_HH
#define SET_INTERPOLATED_CONVECTION_BOUNDARY_2D_HH

#include "setInterpolatedConvectionBoundary2D.h"

namespace olb {

///Initialising the InterpolatedConvectionBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setInterpolatedConvectionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,T omega, SuperGeometry<T,2>& superGeometry, int material, T* uAv)
{
  setInterpolatedConvectionBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry.getMaterialIndicator(material), uAv);
}
///Initialising the InterpolatedConvectionBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setInterpolatedConvectionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega,FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T* uAv)
{
  OstreamManager clout(std::cout, "setInterpolatedConvectionBoundary");
  int _overlap = 1;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setInterpolatedConvectionBoundary<T,DESCRIPTOR>(sLattice.getBlock(iCloc), omega,
        indicator->getBlockIndicatorF(iCloc), uAv, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T, DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

///Set InterpolatedConvectionBoundary for indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setInterpolatedConvectionBoundary(BlockLattice<T,DESCRIPTOR>& block, T omega,BlockIndicatorF2D<T>& indicator,
                                       T* uAv, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  OstreamManager clout(std::cout, "setInterpolatedConvectionBoundary");
  bool _output = false;
  std::vector<int> discreteNormal(3, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY}) >= margin
        && indicator(iX, iY)) {
      PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = nullptr;
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY);
      if (discreteNormal[0] == 0) {
        if (discreteNormal[1] == -1) {
          if (_output) {
            clout << "setInterpolatedConvectionBoundary<" << 0 << ","<< -1 << ">("  << iX << ", "<< iY << ", " << omega << " )" << std::endl;
          }
          postProcessor = new StraightConvectionBoundaryProcessorGenerator2D
          < T,DESCRIPTOR,0,-1 >  (iX, iX, iY, iY, uAv);
        }
        else if (discreteNormal[1] == 1) {
          if (_output) {
            clout << "setInterpolatedConvectionBoundary<" << 0 << "," << 1 << ">("  << iX << ", "<< iY << ", " << omega << " )" << std::endl;
          }
          postProcessor = new StraightConvectionBoundaryProcessorGenerator2D
          < T,DESCRIPTOR,0,1 >  (iX, iX, iY, iY, uAv);
        }
        else if (discreteNormal[2] == -1) {
          if (_output) {
            clout << "setInterpolatedConvectionBoundary<" << 1 << "," << -1 << ">("  << iX << ", "<< iY << ", " << omega << " )" << std::endl;
          }
          postProcessor = new StraightConvectionBoundaryProcessorGenerator2D
          < T,DESCRIPTOR,1,-1 >  (iX, iX, iY, iY, uAv);
        }
        else if (discreteNormal[2] == 1) {
          if (_output) {
            clout << "setInterpolatedConvectionBoundary<" << 1 << "," << 1 << ">("  << iX << ", "<< iY << ", " << omega << " )" << std::endl;
          }
          postProcessor = new StraightConvectionBoundaryProcessorGenerator2D
          < T,DESCRIPTOR,1,1 >  (iX, iX, iY, iY, uAv);
        }
        if (postProcessor) {
          block.addPostProcessor(*postProcessor);
        }
      }
    }
  });
}




}//namespace olb


#endif
