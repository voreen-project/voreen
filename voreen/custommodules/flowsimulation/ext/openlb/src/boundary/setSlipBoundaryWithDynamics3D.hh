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

//This file contains the slip boundary with dynamics
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_SLIP_BOUNDARY_WITH_DYNAMICS_HH
#define SET_SLIP_BOUNDARY_WITH_DYNAMICS_HH

#include "setSlipBoundaryWithDynamics3D.h"

namespace olb {

///Initialising the setSlipBoundaryWithDynamics function on the superLattice domain
template<typename T,typename DESCRIPTOR, class MixinDynamics>
void setSlipBoundaryWithDynamics(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry<T,3>& superGeometry, int material)
{
  setSlipBoundaryWithDynamics<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material));
}

///Initialising the setSlipBoundaryWithDynamics function on the superLattice domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setSlipBoundaryWithDynamics(SuperLattice<T, DESCRIPTOR>& sLattice, T omega,FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  OstreamManager clout(std::cout, "setSlipBoundaryWithDynamics");
  int _overlap = 1;
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setSlipBoundaryWithDynamics<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iC),omega, indicator->getBlockIndicatorF(iC),includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  //the addPoints2CommBC function is initialised inside setLocalVelocityBoundary3D.h/hh
  addPoints2CommBC(sLattice,std::forward<decltype(indicator)>(indicator), _overlap);
}


//set SlipBoundaryWithDynamics on indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setSlipBoundaryWithDynamics(BlockLattice<T,DESCRIPTOR>& block, T omega,BlockIndicatorF3D<T>& indicator, bool includeOuterCells)
{
  using namespace boundaryhelper;
  OstreamManager clout(std::cout, "setslipBoundaryWithDynamics");
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  std::vector<int> discreteNormal(4,0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      Dynamics<T,DESCRIPTOR>* dynamics = nullptr;
      discreteNormal = indicator.getBlockGeometry().getStatistics().getType(iX, iY, iZ);
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
        block.template defineDynamics<BounceBack>({iX, iY, iZ});
      }
      // Setting dynamics and momenta as in interpolatedVelocityBoundary3D
      if (discreteNormal[0] == 0) {
        if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {//set momenta, dynamics and postProcessors on indicated velocityBoundaryCells
          dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletVelocityBoundaryTuple<0,-1>
          >>();
        }
        else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
          dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletVelocityBoundaryTuple<0,1>
          >>();
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
          dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletVelocityBoundaryTuple<1,-1>
          >>();
        }
        else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
          dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletVelocityBoundaryTuple<1,1>
          >>();
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
          dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletVelocityBoundaryTuple<2,-1>
          >>();
        }
        else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
          dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
            momenta::BasicDirichletVelocityBoundaryTuple<2,1>
          >>();
        }
      }
      else if (discreteNormal[0] == 1) {//set momenta,dynamics and postProcessors on indicated velocityBoundary External Corner cells
        dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
          momenta::FixedVelocityBoundaryTuple
        >>();
      }
      else if (discreteNormal[0] == 2) {//Internalvelocitycorner
        dynamics = block.getDynamics(PlainMixinDynamicsForNormalMomenta<T,DESCRIPTOR,
          CombinedRLBdynamics,MixinDynamics,momenta::InnerCornerVelocityTuple3D
        >::construct(Vector<int,3>(discreteNormal.data() + 1)));
      }
      //ExternalVelocityEdge
      else if (discreteNormal[0] == 3) {//set momenta,dynamics and postProcessors on indicated velocityBoundary External Edge cells
        dynamics = block.template getDynamics<typename MixinDynamics::template exchange_momenta<
          momenta::FixedVelocityBoundaryTuple
        >>();
      }
      //InternalVelocityEdge
      else if (discreteNormal[0] == 4) {//set momenta,dynamics and postProcessors on indicated velocityBoundary Inner Edge cells
            dynamics = block.getDynamics(PlainMixinDynamicsForNormalSpecialMomenta<T,DESCRIPTOR,
              CombinedRLBdynamics,MixinDynamics,momenta::InnerEdgeVelocityTuple3D
            >::construct(Vector<int,3>(discreteNormal.data() + 1)));
      }
      dynamics->getParameters(block).template set<descriptors::OMEGA>(omega);
      setBoundary(block, iX,iY,iZ, dynamics);
    }
  });
}

}
#endif
