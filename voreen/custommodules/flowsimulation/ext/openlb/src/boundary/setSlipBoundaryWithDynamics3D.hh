/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz, Davide Dapelo
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

///This file contains the slip boundary with dynamics
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_SLIP_BOUNDARY_WITH_DYNAMICS_HH
#define SET_SLIP_BOUNDARY_WITH_DYNAMICS_HH

#include "setSlipBoundaryWithDynamics2D.h"

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the setSlipBoundaryWithDynamics function on the superLattice domain
template<typename T,typename DESCRIPTOR, class MixinDynamics>
void setSlipBoundaryWithDynamics(SuperLattice3D<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry3D<T>& superGeometry, int material)
{
  setSlipBoundaryWithDynamics<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));
}

///Initialising the setSlipBoundaryWithDynamics function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setSlipBoundaryWithDynamics(SuperLattice3D<T, DESCRIPTOR>& sLattice, T omega,FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setSlipBoundaryWithDynamics");
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
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setSlipBoundaryWithDynamics<T,DESCRIPTOR,MixinDynamics>(sLattice.getExtendedBlockLattice(iC),omega, indicator->getExtendedBlockIndicatorF(iC),includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  //the addPoints2CommBC function is initialised inside setLocalVelocityBoundary3D.h/hh
  addPoints2CommBC(sLattice,std::forward<decltype(indicator)>(indicator), _overlap);
}

////////// BlockLattice Domain  /////////////////////////////////////////

//set SlipBoundaryWithDynamics on indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setSlipBoundaryWithDynamics(BlockLatticeStructure3D<T,DESCRIPTOR>& block, T omega,BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setslipBoundaryWithDynamics");
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  /*
   *x0,x1,y0,y1,z0z1 Range of cells to be traversed
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
        if (indicator(iX, iY, iZ)) {
          discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY, iZ);
          // Setting postProcessor as in setSlipBoundary3D
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {//set postProcessors for indicated cells
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
            block.defineDynamics(iX, iY, iZ, &instances::getBounceBack<T, DESCRIPTOR>());
          }
          // Setting dynamics and momenta as in interpolatedVelocityBoundary3D
          if (discreteNormal[0] == 0) {
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {//set momenta, dynamics and postProcessors on indicated velocityBoundaryCells
              //momenta vector provisionally inside src/core/blockLatticeStructure3D.h
              momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM, 0,-1>;
              //dynamics vector provisionally inside src/core/blockLatticeStructure3D.h
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM, 0,1>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM, 1,-1>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM, 1,1>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM, 2,-1>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM, 2,1>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
          }

          else if (discreteNormal[0] == 1) {//set momenta,dynamics and postProcessors on indicated velocityBoundary External Corner cells
            //ExternalVelocityCorner
            if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            /// addExternalVelocityCorner<discreteNormal[1],discreteNormal[2],discreteNormal[3]>(iX,iY,iZ, omega);
          }

          else if (discreteNormal[0] == 2) {//
            //Internalvelocitycorner
            if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, 1,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, 1,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, 1,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, 1,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, -1,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, 1,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, -1,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new InnerCornerVelBM3D<T,DESCRIPTOR, -1,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }

          }
          //ExternalVelocityEdge
          else if (discreteNormal[0] == 3) {//set momenta,dynamics and postProcessors on indicated velocityBoundary External Edge cells
            if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
            }

          }
          //InternalVelocityEdge
          else if (discreteNormal[0] == 4) {//set momenta,dynamics and postProcessors on indicated velocityBoundary Inner Edge cells
            if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 0,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 0,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 0,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 0,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 1,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 1,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 1,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 1,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 2,1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 2,-1,1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 2,1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              momenta = new InnerEdgeVelBM3D<T,DESCRIPTOR, 2,-1,-1>;
              dynamics = new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, *momenta);
            }
          }
        }
      }
    }
  }

}

}
#endif
