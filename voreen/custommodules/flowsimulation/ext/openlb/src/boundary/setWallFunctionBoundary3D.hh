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

///This file contains the WallFunctionBoundary
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_WALL_FUNCTION_BOUNDARY_3D_HH
#define SET_WALL_FUNCTION_BOUNDARY_3D_HH

#include "setWallFunctionBoundary3D.h"

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the WallFunctionBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setWallFunctionBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry, int material,
                             UnitConverter<T, DESCRIPTOR> const& converter,
                             wallFunctionParam<T> const& wallFunctionParam,
                             IndicatorF3D<T>* geoIndicator)
{

  setWallFunctionBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, superGeometry.getMaterialIndicator(material),
      converter,
      wallFunctionParam,
      geoIndicator);
}

///Initialising the WallFunctionBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setWallFunctionBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                             UnitConverter<T, DESCRIPTOR> const& converter,
                             wallFunctionParam<T> const& wallFunctionParam,
                             IndicatorF3D<T>* geoIndicator)
{
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */
  int _overlap = 1;
  bool includeOuterCells = false;
  OstreamManager clout(std::cout, "setWallFunctionBoundary");
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setWallFunctionBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getExtendedBlockLattice(iCloc),
        indicator->getExtendedBlockIndicatorF(iCloc),
        converter, wallFunctionParam, geoIndicator, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);

}

////////// BlockLattice Domain  /////////////////////////////////////////

/// Set WallFunction boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setWallFunctionBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator,
                             UnitConverter<T, DESCRIPTOR> const& converter,
                             wallFunctionParam<T> const& wallFunctionParam,
                             IndicatorF3D<T>*            geoIndicator,
                             bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setWallFunctionBoundary");
  bool _output = false;
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
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY, iZ);
          std::vector<int> missingIndices;
          for (int x = -1 ; x < 2; ++x) {
            for (int y = -1 ; y < 2; ++y) {
              for (int z = -1 ; z < 2; ++z) {
                if (blockGeometryStructure.getMaterial(iX + x, iY + y, iZ + z) == 0) {
                  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
                    if (descriptors::c<DESCRIPTOR>(iPop,0) == x &&
                        descriptors::c<DESCRIPTOR>(iPop,1) == y &&
                        descriptors::c<DESCRIPTOR>(iPop,2) == z) {
                      missingIndices.push_back(descriptors::opposite<DESCRIPTOR>(iPop));
                    }
                  }
                }
              }
            }
          }
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
            discreteNormal.erase(discreteNormal.begin());
            //OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
            if (_output) {
              clout << "setWallFunctionBoundary<" << discreteNormal[0] << ","<< discreteNormal[1] << ","<< discreteNormal[2] << ">("  << iX << ", "<< iX << ", " << iY << ", " << iY << ", " << iZ << ", " << iZ << " )" << std::endl;
            }
            PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new WallFunctionBoundaryProcessorGenerator3D<T, DESCRIPTOR>(iX, iX, iY, iY, iZ, iZ, blockGeometryStructure, discreteNormal, missingIndices,
                converter, wallFunctionParam, geoIndicator);
            if (postProcessor) {
              _block.addPostProcessor(*postProcessor);
            }
          }
          else {
            clout << "Warning: Could not setWallFunctionBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<","<< discreteNormal[3] <<"), set to bounceBack" << std::endl;
            _block.defineDynamics(iX, iY, iZ, &instances::getBounceBack<T, DESCRIPTOR>());
          }
        }
      }
    }
  }
}


} //namespace olb


#endif

