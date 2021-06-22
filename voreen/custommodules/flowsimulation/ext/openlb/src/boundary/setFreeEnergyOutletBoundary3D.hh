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

///This file contains the FreeEnergyOutletBoundary
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_FREE_ENERGY_OUTLET_BOUNDARY_3D_HH
#define SET_FREE_ENERGY_OUTLET_BOUNDARY_3D_HH

#include "setFreeEnergyOutletBoundary3D.h"



namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////

/// Implementation of a outlet boundary condition for the partner lattices of the binary or the ternary free energy model.
///Initialising the FreeEnergyOutletBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setFreeEnergyOutletBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry3D<T>& superGeometry, int material, std::string type, int latticeNumber)
{
  setFreeEnergyOutletBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material), type, latticeNumber);
}


/// Implementation of a outlet boundary condition for the partner lattices of the binary or the ternary free energy model.
///Initialising the FreeEnergyOutletBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setFreeEnergyOutletBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice,T omega,FunctorPtr<SuperIndicatorF3D<T>>&& indicator, std::string type, int latticeNumber)
{
  OstreamManager clout(std::cout, "setFreeEnergyOutletBoundary");
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
    setFreeEnergyOutletBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getExtendedBlockLattice(iCloc),
        indicator->getExtendedBlockIndicatorF(iCloc),  omega, type, latticeNumber, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}


////////// BlockLattice Domain  /////////////////////////////////////////

/// Set FreeEnergyOutletBoundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setFreeEnergyOutletBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, T omega, std::string type,
                                 int latticeNumber, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setFreeEnergyOutletBoundary");
  bool _output = false;
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  //setFreeEnergyInletBoundary on indicated cells inside the blockLattice domain
  setFreeEnergyInletBoundary<T,DESCRIPTOR,MixinDynamics>(_block, omega, indicator, type, latticeNumber,includeOuterCells);
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
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX,iY,iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ, true);
          if (discreteNormal[0] == 0) {//set postProcessor on indicated cells

            PostProcessorGenerator3D<T, DESCRIPTOR>* convectivePostProcessor =
              new FreeEnergyConvectiveProcessorGenerator3D<T, DESCRIPTOR> (
              iX, iX, iY, iY, iZ, iZ, discreteNormal[1], discreteNormal[2], discreteNormal[3] );
            if (convectivePostProcessor) {
              _block.addPostProcessor(*convectivePostProcessor);
            }

            if (_output) {
              clout << "setFreeEnergyOutletBoundary<" << "," << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << ")" << std::endl;
            }
          }
        }
      }
    }
  }


}


}//namespace olb
#endif
