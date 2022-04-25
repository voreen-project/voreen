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

//This file contains the Slip Boundary
//This is an onLattice boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_SLIP_BOUNDARY_HH
#define SET_SLIP_BOUNDARY_HH

#include "setSlipBoundary3D.h"


namespace olb {


///Initialising the setslipBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material)
{
  setSlipBoundary<T, DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material));
}

///Initialising the setslipBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setSlipBoundary( SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  bool includeOuterCells = false;
  int _overlap = 1;
  OstreamManager clout(std::cout, "setslipBoundary");
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setSlipBoundary(sLattice.getBlock(iCloc), indicator->getBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

/// Set slip boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setSlipBoundary(BlockLattice<T,DESCRIPTOR>& block,BlockIndicatorF3D<T>& indicator, bool includeOuterCells )
{
  OstreamManager clout(std::cout, "setslipBoundary");
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
      if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {//set postProcessors for indicated cells
        OstreamManager clout(std::cout, "setslipBoundary");
        bool _output = false;
        if (_output) {
          clout << "setSlipBoundary<" << discreteNormal[1] << ","<< discreteNormal[2] << ","<< discreteNormal[3] << ">("  << iX << ", "<< iX << ", " << iY << ", " << iY << ", " << iZ << ", " << iZ << " )" << std::endl;
        }
        PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new SlipBoundaryProcessorGenerator3D<T, DESCRIPTOR>(iX, iX, iY, iY, iZ, iZ, discreteNormal[1], discreteNormal[2], discreteNormal[3]);
        if (postProcessor) {
          block.addPostProcessor(*postProcessor);
        }
      }
      else {//define dynamics for indicated cells
        clout << "Warning: Could not setSlipBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<","<< discreteNormal[3] <<"), set to bounceBack" << std::endl;
        block.template defineDynamics<BounceBack>({iX, iY, iZ});
      }
    }
  });
}
}//namespace olb

#endif
