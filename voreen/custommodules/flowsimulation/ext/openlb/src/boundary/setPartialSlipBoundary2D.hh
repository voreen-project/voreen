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

///This file contains the Partial Slip Boundary
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_PARTIAL_SLIP_BOUNDARY_2D_HH
#define SET_PARTIAL_SLIP_BOUNDARY_2D_HH


#include "setPartialSlipBoundary2D.h"


namespace olb {

///Initialising the Partial slip boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setPartialSlipBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, T tuner, SuperGeometry2D<T>& superGeometry, int material)
{

  setPartialSlipBoundary<T,DESCRIPTOR>(sLattice, tuner, superGeometry.getMaterialIndicator(material));

}

///Initialising the Partial slip boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setPartialSlipBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, T tuner, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setPartialSlipBoundary");
  bool includeOuterCells = false;
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */
  int _overlap = 1;

  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setPartialSlipBoundary<T,DESCRIPTOR>(sLattice.getExtendedBlockLattice(iCloc),
                                         tuner, indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  //the addPoints2CommBC function is initialised inside setLocalVelocityBoundary2D.h/hh
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

/// Set Partial Slip boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setPartialSlipBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, T tuner, BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setPartialSlipBoundary");
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
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
      if (indicator(iX, iY)) {
        if (tuner < 0. || tuner > 1.) {
          clout << "Warning: Could not setPartialSlipBoundary (" << iX << ", " << iY << "), tuner must be between 0.1 and instead is=" << tuner <<", set to bounceBack" << std::endl;
          block.defineDynamics(iX, iY, &instances::getBounceBack<T, DESCRIPTOR>() );
        }
        else {
          discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0) {
            //set partial slip boundary on indicated cells
            bool _output = false;
            if (_output) {
              clout << "setPartialSlipBoundary<" << discreteNormal[1] << ","<< discreteNormal[2] << ">("  << iX << ", "<< iX << ", " << iY << ", " << iY << " )" << std::endl;
            }
            PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = new PartialSlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(tuner, iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);
            if (postProcessor) {
              block.addPostProcessor(*postProcessor);
            }
          }
          else {
            clout << "Warning: Could not setPartialSlipBoundary (" << iX << ", " << iY << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<"), set to bounceBack" << std::endl;
            block.defineDynamics(iX, iY, &instances::getBounceBack<T, DESCRIPTOR>() );
          }
        }
      }
    }
  }
}
}//namespace olb


#endif
