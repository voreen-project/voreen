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
#ifndef SET_SLIP_BOUNDARY_WITH_DYNAMICS_2D_HH
#define SET_SLIP_BOUNDARY_WITH_DYNAMICS_2D_HH

#include "setSlipBoundaryWithDynamics2D.h"

namespace olb {
///Initialising the setSlipBoundaryWithDynamics function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setSlipBoundaryWithDynamics(SuperLattice2D<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry2D<T>& superGeometry, int material)
{

  setSlipBoundaryWithDynamics<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));
}

///Initialising the setSlipBoundaryWithDynamics function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setSlipBoundaryWithDynamics(SuperLattice2D<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  bool includeOuterCells = false;
  int _overlap = indicator->getSuperGeometry().getOverlap();
  OstreamManager clout(std::cout, "setOnBCInterpolatedBoundary");
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setSlipBoundaryWithDynamics<T,DESCRIPTOR,MixinDynamics>(sLattice.getExtendedBlockLattice(iCloc), omega,
        indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  //the addPoints2CommBC function is initialised inside setLocalVelocityBoundary2D.h/hh
  addPoints2CommBC<T,DESCRIPTOR>(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);


}
////////// BlockLattice Domain  /////////////////////////////////////////

/// Set Interpolated velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setSlipBoundaryWithDynamics(BlockLatticeStructure2D<T,DESCRIPTOR>& block, T omega,BlockIndicatorF2D<T>& indicator, bool includeOuterCells)
{
  OstreamManager clout(std::cout, "setSlipBoundaryWithDynamics");
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
      Momenta<T, DESCRIPTOR>* momenta = nullptr;
      Dynamics<T, DESCRIPTOR>* dynamics = nullptr;
      PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = nullptr;
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {
          ///sets momenta, dynamics and post processor on velocity boundary cells
          if (discreteNormal[1] == 1) {
            //momenta vector provisionally inside src/core/blockLatticeStructure2D.h
            momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM,0,1>;
            //dynamics vector provisionally inside src/core/blockLatticeStructure2D.h
            dynamics = new MixinDynamics(omega, *momenta);
            postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);

          }
          else if (discreteNormal[1] == -1) {
            momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM,0,-1>;
            dynamics = new MixinDynamics(omega, *momenta);
            postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);


          }
          else if (discreteNormal[2] == 1) {
            momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM,1,1>;
            dynamics = new MixinDynamics(omega, *momenta);
            postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);

          }
          else if (discreteNormal[2] == -1) {
            momenta = new BasicDirichletBM<T,DESCRIPTOR,VelocityBM,1,-1>;
            dynamics = new MixinDynamics(omega, *momenta);
            postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);
          }
          else {
            clout << "Could not setSlipBoundaryWithDynamics (" << iX
                  << ", " << iY << ")" << std::endl;
          }
        }

        else if (discreteNormal[0] == 1) {
          ///sets momenta, dynamics and post processors on externalVelocityCorner Boundary cells
          if (discreteNormal[1] == 1) {
            if (discreteNormal[2] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);
            }
            else if (discreteNormal[2] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);

            }
            else {
              clout << "Could not setSlipBoundaryWithDynamics ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          }
          else if (discreteNormal[1] == -1) {
            if (discreteNormal[2] == 1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);

            }
            else if (discreteNormal[2] == -1) {
              momenta = new FixedVelocityBM<T,DESCRIPTOR>;
              dynamics = new MixinDynamics(omega, *momenta);
              postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);


            }
            else {
              clout << "Could not setSlipBoundaryWithDynamics ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          }
        }
        else if (discreteNormal[0] == 2) {
          ///sets momenta, dynamics and post processors on internalVelocityCorner Boundary cells
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
              clout << "Could not setSlipBoundaryWithDynamics ("
                    << iX << ", " << iY << ")" << std::endl;
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
              clout << "Could not setSlipBoundaryWithDynamics ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          }
        }
        //this function is initialised inside setLocalVelocityBoundary2D.h/hh
        //setting of boundary on cells
        setBoundary<T, DESCRIPTOR, MixinDynamics>(block, omega, iX,iY, momenta, dynamics, postProcessor);
      }
    }
  }

}


}//namespace olb

#endif
