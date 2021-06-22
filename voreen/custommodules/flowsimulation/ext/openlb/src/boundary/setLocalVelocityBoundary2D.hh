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

///This file contains the Local Velocity Boundary
///This is a new version of the Boundary, which only contains free floating functions

#ifndef SET_LOCAL_VELOCITY_BOUNDARY_2D_HH
#define SET_LOCAL_VELOCITY_BOUNDARY_2D_HH

#include "setLocalVelocityBoundary2D.h"

namespace olb {
///Initialising the setLocalVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setLocalVelocityBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice,T omega, SuperGeometry2D<T>& superGeometry, int material)
{
  setLocalVelocityBoundary<T,DESCRIPTOR, MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));


}
///Initialising the setLocalVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setLocalVelocityBoundary(SuperLattice2D<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setLocalVelocityBoundary2D");
  /*  local boundaries: _overlap = 0;
   *  interp boundaries: _overlap = 1;
   *  bouzidi boundaries: _overlap = 1;
   *  extField boundaries: _overlap = 1;
   *  advectionDiffusion boundaries: _overlap = 1;
   */
  int _overlap = 0;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setLocalVelocityBoundary<T,DESCRIPTOR, MixinDynamics>(sLattice.getExtendedBlockLattice(iCloc), omega,
        indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

/// Adds needed Cells to the Communicator _commBC in SuperLattice
template<typename T, typename DESCRIPTOR>
void addPoints2CommBC(SuperLattice2D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,int _overlap)
{
  if (_overlap == 0) {
    return;
  }

  SuperGeometry2D<T>& superGeometry = indicator->getSuperGeometry();
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    const int nX = superGeometry.getBlockGeometry(iCloc).getNx();
    const int nY = superGeometry.getBlockGeometry(iCloc).getNy();

    for (int iX = -_overlap; iX < nX+_overlap; ++iX) {
      for (int iY = -_overlap; iY < nY+_overlap; ++iY) {
        if (iX < 0 || iX > nX - 1 ||
            iY < 0 || iY > nY - 1 ) { // if within overlap
          if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY) != 0) {
            bool found = false;
            for (int iXo = -_overlap; iXo <= _overlap && !found; ++iXo) {
              for (int iYo = -_overlap; iYo <= _overlap && !found; ++iYo) {
                const int nextX = iXo + iX;
                const int nextY = iYo + iY;
                if (indicator->getBlockIndicatorF(iCloc)(nextX, nextY)) {
                  sLattice.get_commBC().add_cell(iCloc, iX, iY);
                  found = true;
                }
              }
            }
          }
        }
      }
    }
  }

}


/// Set local velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setLocalVelocityBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block, T omega, BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setLocalVelocityBoundary2D");
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
      //momenta vector provisionally inside src/core/blockLatticeStructure3D.h
      Momenta<T, DESCRIPTOR>* momenta = nullptr;
      //dynamics vector provisionally inside src/core/blockLatticeStructure3D.h
      Dynamics<T, DESCRIPTOR>* dynamics = nullptr;
      PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = nullptr;
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {
          //set the momenta,dynamics and post processor on the indicated local velocity boundary cells
          if (discreteNormal[1] == 1) {
            momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 0,1>;
            dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            postProcessor = nullptr;
          }
          else if (discreteNormal[1] == -1) {
            momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 0,-1>;
            dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            postProcessor = nullptr;
          }
          else if (discreteNormal[2] == 1) {
            momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 1,1>;
            dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            postProcessor = nullptr;
          }
          else if (discreteNormal[2] == -1) {
            momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 1,-1>;
            dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            postProcessor = nullptr;
          }
          else {
            clout << "Could not setLocalVelocityBoundary2D (" << iX
                  << ", " << iY << ")" << std::endl;
          }
        }
        else if (discreteNormal[0] == 1) {
          //sets the momenta, dynamics and post processors on indicated localVelocityCornerBoundary cells
          if (discreteNormal[1] == 1) {
            if (discreteNormal[2] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator2D<T,DESCRIPTOR, 1,1> (iX,iY);
            }
            else if (discreteNormal[2] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator2D<T,DESCRIPTOR, 1,-1> (iX,iY);
            }
            else {
              clout << "Could not setLocalVelocityBoundary2D (" << iX
                    << ", " << iY << ")" << std::endl;
            }
          }
          else if (discreteNormal[1] == -1) {
            if (discreteNormal[2] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator2D<T,DESCRIPTOR, -1,1> (iX,iY);
            }
            else if (discreteNormal[2] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator2D<T,DESCRIPTOR, -1,-1> (iX,iY);
            }
            else {
              clout << "Could not setLocalVelocityBoundary2D (" << iX
                    << ", " << iY << ")" << std::endl;
            }
          }
        }
        //sets momenta, dynamics and postProcessors on local InnerVelocityCornerBoundary Cells
        else if (discreteNormal[0] == 2) {
          if (discreteNormal[1] == 1) {
            if (discreteNormal[2] == 1) {
              momenta = new InnerCornerVelBM2D<T,DESCRIPTOR, 1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[2] == -1) {
              momenta = new InnerCornerVelBM2D<T,DESCRIPTOR, 1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else {
              clout << "Could not setLocalVelocityBoundary2D (" << iX
                    << ", " << iY << ")" << std::endl;
            }
          }
          else if (discreteNormal[1] == -1) {
            if (discreteNormal[2] == 1) {
              momenta = new InnerCornerVelBM2D<T,DESCRIPTOR, -1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[2] == -1) {
              momenta = new InnerCornerVelBM2D<T,DESCRIPTOR, -1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else {
              clout << "Could not setLocalVelocityBoundary2D (" << iX
                    << ", " << iY << ")" << std::endl;
            }
          }
        }
        //sets the boundary on any indicated cell
        setBoundary<T, DESCRIPTOR, MixinDynamics>(block, omega, iX,iY, momenta, dynamics, postProcessor);

      }
    }
  }

}

//sets the boundary on any indicated cell
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBoundary(BlockLatticeStructure2D<T,DESCRIPTOR>& block,T omega, int iX,int iY,
                 Momenta<T,DESCRIPTOR>* momenta,Dynamics<T,DESCRIPTOR>* dynamics,PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor)
{
  OstreamManager clout(std::cout, "setBoundary");
  block.defineDynamics(iX,iX,iY,iY, dynamics);
  //momenta vector provisionally inside src/core/blockLatticeStructure3D.h
  block.momentaVector.push_back(momenta);
  //dynamics vector provisionally inside src/core/blockLatticeStructure3D.h
  block.dynamicsVector.push_back(dynamics);
  if (postProcessor) {
    block.addPostProcessor(*postProcessor);
  }
}


}//namespace olb
#endif
