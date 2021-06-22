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

///This file contains the Slip Boundary
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_SLIP_BOUNDARY_2D_HH
#define SET_SLIP_BOUNDARY_2D_HH

#include "setSlipBoundary2D.h"



namespace olb {
///Initialising the SlipBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setSlipBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry, int material)
{
  setSlipBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material));

}
///Initialising the SlipBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setSlipBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setSlipBoundary");
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
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setSlipBoundary<T, DESCRIPTOR>(sLattice.getExtendedBlockLattice(iCloc),
                                   indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  //defined inside setLocalVelocityBoundary2D.h/hh
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T, DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}
///Set slip boundary on blockLattice domain
template<typename T, typename DESCRIPTOR>
void setSlipBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setSlipBoundary");
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
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
        if (discreteNormal[1]!=0 || discreteNormal[2]!=0) {
          //set slip boundary on indicated cells
          bool _output = false;
          if (_output) {
            clout << "setSlipBoundary<" << discreteNormal[1] << ","<< discreteNormal[2] << ">("  << iX << ", "<< iX << ", " << iY << ", " << iY << " )" << std::endl;
          }
          PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);
          if (postProcessor) {
            block.addPostProcessor(*postProcessor);
          }
        }
        else {
          clout << "Warning: Could not addSlipBoundary (" << iX << ", " << iY << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<"), set to bounceBack" << std::endl;
          block.defineDynamics(iX, iY, &instances::getBounceBack<T, DESCRIPTOR>() );
        }
      }
    }
  }
}


}//namespace olb

#endif
