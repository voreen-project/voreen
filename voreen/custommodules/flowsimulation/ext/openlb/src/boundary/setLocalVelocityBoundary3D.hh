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
#ifndef SET_LOCAL_VELOCITY_BOUNDARY_HH
#define SET_LOCAL_VELOCITY_BOUNDARY_HH

#include "setLocalVelocityBoundary3D.h"

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////


///Initialising the setLocalVelocityBoundary function on the superLattice domain
template<typename T,typename DESCRIPTOR, class MixinDynamics>
void setLocalVelocityBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry3D<T>& superGeometry, int material)
{
  setLocalVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));
}
///Initialising the setLocalVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setLocalVelocityBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setLocalVelocityBoundary");
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

  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setLocalVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getExtendedBlockLattice(iC), omega, indicator->getExtendedBlockIndicatorF(iC),includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}

/// Adds needed Cells to the Communicator _commBC in SuperLattice
template<typename T, typename DESCRIPTOR>
void addPoints2CommBC(SuperLattice3D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator, int _overlap)
{
  if (_overlap == 0) {
    return;
  }

  SuperGeometry3D<T>& superGeometry = indicator->getSuperGeometry();
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    const int nX = superGeometry.getBlockGeometry(iCloc).getNx();
    const int nY = superGeometry.getBlockGeometry(iCloc).getNy();
    const int nZ = superGeometry.getBlockGeometry(iCloc).getNz();

    for (int iX = -_overlap; iX < nX+_overlap; ++iX) {
      for (int iY = -_overlap; iY < nY+_overlap; ++iY) {
        for (int iZ = -_overlap; iZ < nZ+_overlap; ++iZ) {
          if (iX < 0 || iX > nX - 1 ||
              iY < 0 || iY > nY - 1 ||
              iZ < 0 || iZ > nZ - 1 ) { // if within overlap
            if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY,iZ) != 0) {
              bool found = false;
              for (int iXo = -_overlap; iXo <= _overlap && !found; ++iXo) {
                for (int iYo = -_overlap; iYo <= _overlap && !found; ++iYo) {
                  for (int iZo = -_overlap; iZo <= _overlap && !found; ++iZo) {
                    const int nextX = iXo + iX;
                    const int nextY = iYo + iY;
                    const int nextZ = iZo + iZ;
                    if (indicator->getBlockIndicatorF(iCloc)(nextX, nextY, nextZ)
                        && nextX >= -_overlap && nextX < nX+_overlap
                        && nextY >= -_overlap && nextY < nY+_overlap
                        && nextZ >= -_overlap && nextZ < nZ+_overlap) {
                      sLattice.get_commBC().add_cell(iCloc, iX, iY, iZ);
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
  }
}

////////// BlockLattice Domain  /////////////////////////////////////////


/// Set local velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setLocalVelocityBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block,T omega, BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
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
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {
        Momenta<T,DESCRIPTOR>* momenta = nullptr;
        Dynamics<T,DESCRIPTOR>* dynamics = nullptr;
        PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor = nullptr;
        if (indicator(iX, iY, iZ)) {
          discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY, iZ);
          if (discreteNormal[0] == 0) {
            //setVelocityBoundary
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {//set momenta, dynamics and postProcessors on indicated cells
              momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 0,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 0,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 2,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              momenta = new RegularizedVelocityBM<T,DESCRIPTOR, 2,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
          }
          //ExternalVelocityCorner
          else if (discreteNormal[0] == 1) {
            if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {//set momenta, dynamics and postProcessors on indicated external corner cells
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, 1,1,1> (iX,iY,iZ);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, 1,-1,1> (iX,iY,iZ);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, 1,1,-1> (iX,iY,iZ);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, 1,-1,-1> (iX,iY,iZ);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, -1,1,1> (iX,iY,iZ);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, -1,-1,1> (iX,iY,iZ);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, -1,1,-1> (iX,iY,iZ);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, -1,-1,-1> (iX,iY,iZ);
            }
          }
          //InternalVelocityCorner
          else if (discreteNormal[0] == 2) {//set momenta, dynamics and postProcessors on indicated internal corner cells

            if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, 1,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, 1,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, 1,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, 1,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, -1,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, -1,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, -1,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, -1,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            /// addInternalVelocityCorner<discreteNormal[1],discreteNormal[2],discreteNormal[3]>(iX,iY,iZ, omega);

          }
          //ExternalVelocityEdge
          else if (discreteNormal[0] == 3) {//set momenta, dynamics and postProcessors on indicated external edge cells
            if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 0,1,1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 0,-1,1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 0,1,-1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 0,-1,-1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 1,1,1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 1,1,-1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 1,-1,1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 1,-1,-1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 2,1,1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 2,-1,1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 2,1,-1>(iX,iX,iY,iY,iZ,iZ);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, 2,-1,-1>(iX,iX,iY,iY,iZ,iZ);
            }

          }
          //InternalVelocityEdge
          else if (discreteNormal[0] == 4) {//set momenta, dynamics and postProcessors on indicated internal edge cells
            if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 0,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 0,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 0,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 0,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 1,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 1,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 1,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 1,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 2,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 2,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 2,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 2,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
              postProcessor = nullptr;
            }
          }
          //setBoundary on indicated cells
          setBoundary<T, DESCRIPTOR, MixinDynamics>(_block, omega, iX,iY,iZ, momenta, dynamics, postProcessor);

        }
      }
    }
  }

}


///sets boundary on indicated cells. This is a function, which can be used on many boundaries. It is not specific to the local velocity boundary
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block,T omega, int iX,int iY,int iZ,
                 Momenta<T,DESCRIPTOR>* momenta,Dynamics<T,DESCRIPTOR>* dynamics,PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor)
{
  _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
  //momentaVector and dynamicsVector are two vectors, which were previously located in the BoundaryInstantiator3D
  _block.momentaVector.push_back(momenta);
  _block.dynamicsVector.push_back(dynamics);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }

}

}//namespace olb
#endif

