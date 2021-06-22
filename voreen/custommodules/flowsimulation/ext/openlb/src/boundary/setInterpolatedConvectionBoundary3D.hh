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
#ifndef SET_INTERPOLATED_CONVECTION_BOUNDARY_3D_HH
#define SET_INTERPOLATED_CONVECTION_BOUNDARY_3D_HH

#include "setInterpolatedConvectionBoundary3D.h"

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the setInterpolatedConvectionBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setInterpolatedConvectionBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice,T omega, SuperGeometry3D<T>& superGeometry, int material,
                                       T* uAv)
{
  setInterpolatedConvectionBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material), uAv);

}

///Initialising the setInterpolatedConvectionBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setInterpolatedConvectionBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                                       T* uAv)
{
  OstreamManager clout(std::cout, "setInterpolatedConvectionBoundary");
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
    setInterpolatedConvectionBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getExtendedBlockLattice(iCloc),omega, indicator->getExtendedBlockIndicatorF(iCloc),
        uAv, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

////////// BlockLattice Domain  /////////////////////////////////////////


/// set InterpolatedConvection boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setInterpolatedConvectionBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, T omega, BlockIndicatorF3D<T>& indicator, T* uAv, bool includeOuterCells)
{

  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  /*
   *x0,x1,y0,y1, z0, z1 Range of cells to be traversed
   **/
  int x0 = margin;
  int y0 = margin;
  int z0 = margin;
  int x1 = blockGeometryStructure.getNx()-1 -margin;
  int y1 = blockGeometryStructure.getNy()-1 -margin;
  int z1 = blockGeometryStructure.getNz()-1 -margin;
  std::vector<int> discreteNormal(4,0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor = nullptr;
        if (indicator(iX, iY, iZ)) {//set postProcessors for indicated boundary cells
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);

          if (discreteNormal[0] == 0) {
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
              postProcessor = new StraightConvectionBoundaryProcessorGenerator3D
              <T,DESCRIPTOR,0,-1>(x0, x1, y0, y1, z0, z1, uAv);
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              postProcessor = new StraightConvectionBoundaryProcessorGenerator3D
              <T,DESCRIPTOR,0,1>(x0, x1, y0, y1, z0, z1, uAv);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              postProcessor = new StraightConvectionBoundaryProcessorGenerator3D
              <T,DESCRIPTOR,1,-1>(x0, x1, y0, y1, z0, z1, uAv);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              postProcessor = new StraightConvectionBoundaryProcessorGenerator3D
              <T,DESCRIPTOR,1,1>(x0, x1, y0, y1, z0, z1, uAv);
            }

            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              postProcessor = new StraightConvectionBoundaryProcessorGenerator3D
              <T,DESCRIPTOR,2,-1>(x0, x1, y0, y1, z0, z1, uAv);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              postProcessor = new StraightConvectionBoundaryProcessorGenerator3D
              <T,DESCRIPTOR,2,1>(x0, x1, y0, y1, z0, z1, uAv);

            }
            if (postProcessor) {
              _block.addPostProcessor(*postProcessor);
            }
          }
        }
      }
    }
  }
}


}//namespace olb
#endif

