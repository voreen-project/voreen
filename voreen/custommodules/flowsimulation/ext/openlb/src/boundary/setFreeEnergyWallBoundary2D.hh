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

///This file contains the Free Energy Wall Boundary
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_FREE_ENERGY_WALL_BOUNDARY_2D_HH
#define SET_FREE_ENERGY_WALL_BOUNDARY_2D_HH

#include "setFreeEnergyWallBoundary2D.h"


namespace olb {
/// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
/// dynamics and an FreeEnergyWall PostProcessor.
/// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
/// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
/// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]

////////// SuperLattice Domain  /////////////////////////////////////////
///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
                               int material, T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber)
{
  setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material),
                                          alpha, kappa1, kappa2, h1, h2, latticeNumber);
}
///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                               T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber)
{
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */
  int _overlap = 1;
  OstreamManager clout(std::cout, "setFreeEnergyWallBoundary");
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  T addend = 0;
  if (latticeNumber==1) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) );
  }
  else if (latticeNumber==2) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (-h2/kappa2) );
  }
  else if (latticeNumber==3) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) );
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice.getExtendedBlockLattice(iCloc),
                                            indicator->getExtendedBlockIndicatorF(iCloc), addend, latticeNumber, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}




/// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
/// dynamics and an FreeEnergyWall PostProcessor.
/// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
/// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa3_ - Parameter related to surface tension. [lattice units]
/// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h3_ - Parameter related to resulting contact angle of the boundary. [lattice units]

////////// SuperLattice Domain  /////////////////////////////////////////
///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
                               int material, T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber)
{
  setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material),
                                          alpha, kappa1, kappa2, kappa3, h1, h2, h3, latticeNumber);
}
///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                               T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber)
{

  int _overlap = indicator->getSuperGeometry().getOverlap();
  OstreamManager clout(std::cout, "setFreeEnergyWallBoundary");
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  T addend = 0;
  if (latticeNumber==1) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) + (h3/kappa3) );
  }
  else if (latticeNumber==2) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (-h2/kappa2) );
  }
  else if (latticeNumber==3) {
    addend = 1./(alpha*alpha) * ( (h3/kappa3) );
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice.getExtendedBlockLattice(iCloc),
                                            indicator->getExtendedBlockIndicatorF(iCloc), addend, latticeNumber, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

////////// BlockLattice Domain  /////////////////////////////////////////
//set FreeEnergyWallBoundary on block domain.
//This function works for the setFreeEnergyWallBoundary with h1,h2,h3 Parameters and h1,h2 Parameters
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator,
                               T addend, int latticeNumber, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setFreeEnergyWallBoundary");
  bool _output = false;
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
      ///set dynamics and post Processors on indicated cells
      if (indicator(iX,iY)) {
        discreteNormal = blockGeometryStructure.getStatistics().getType(iX,iY);
        if (discreteNormal[1]!=0 || discreteNormal[2]!=0) {

          Dynamics<T, DESCRIPTOR>* dynamics = NULL;
          if (latticeNumber == 1) {
            dynamics = &instances::getBounceBack<T, DESCRIPTOR>();
          }
          else {
            dynamics = new FreeEnergyWallDynamics<T, DESCRIPTOR>;
            block.dynamicsVector.push_back(dynamics);
          }
          block.get(iX,iY).defineDynamics(dynamics);
          block.dynamicsVector.push_back(dynamics);

          PostProcessorGenerator2D<T, DESCRIPTOR>* wettingPostProcessor =
            new FreeEnergyWallProcessorGenerator2D<T, DESCRIPTOR> (
            iX, iX, iY, iY, discreteNormal[1], discreteNormal[2], addend );
          PostProcessorGenerator2D<T, DESCRIPTOR>* chemPotPostProcessor =
            new FreeEnergyChemPotBoundaryProcessorGenerator2D<T, DESCRIPTOR> (
            iX, iX, iY, iY, discreteNormal[1], discreteNormal[2], latticeNumber );
          if (wettingPostProcessor) {
            block.addPostProcessor(*wettingPostProcessor);
          }
          if (chemPotPostProcessor) {
            block.addPostProcessor(*chemPotPostProcessor);
          }
        }
        if (_output) {
          clout << "setFreeEnergyWallBoundary<" << "," << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << " )" << std::endl;
        }
      }
    }
  }

}

}//namespace olb
#endif
