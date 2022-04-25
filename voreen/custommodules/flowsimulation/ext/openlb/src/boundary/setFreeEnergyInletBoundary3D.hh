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

//This file contains the FreeEnergyInlet Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_FREE_ENERGY_INLET_BOUNDARY_3D_HH
#define SET_FREE_ENERGY_INLET_BOUNDARY_3D_HH

#include "setFreeEnergyInletBoundary3D.h"

namespace olb {

///Initialising the FreeEnergyInletBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setFreeEnergyInletBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,T omega, SuperGeometry<T,3>& superGeometry, int material, std::string type, int latticeNumber)
{
  setFreeEnergyInletBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice, omega, superGeometry.getMaterialIndicator(material), type, latticeNumber);
}

///Initialising the FreeEnergyInletBoundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setFreeEnergyInletBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,T omega,FunctorPtr<SuperIndicatorF3D<T>>&& indicator,std::string type, int latticeNumber)
{
  bool includeOuterCells = false;
  int _overlap = 1;
  OstreamManager clout(std::cout, "setFreeEnergyInletBoundary");
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    setFreeEnergyInletBoundary<T,DESCRIPTOR,MixinDynamics>(sLattice.getBlock(iCloc),omega,
        indicator->getBlockIndicatorF(iCloc), type, latticeNumber, includeOuterCells);
  }
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  addPoints2CommBC(sLattice, std::forward<decltype(indicator)>(indicator), _overlap);
}


/// Set FreeEnergyInletBoundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setFreeEnergyInletBoundary(BlockLattice<T,DESCRIPTOR>& _block, T omega, BlockIndicatorF3D<T>& indicator, std::string type,
                                int latticeNumber, bool includeOuterCells)
{
  auto& blockGeometryStructure = indicator.getBlockGeometry();
  const int margin = includeOuterCells ? 0 : 1;
  OstreamManager clout(std::cout, "setFreeEnergyInletBoundary");
  bool _output = false;
  std::vector<int> discreteNormal(4, 0);
  blockGeometryStructure.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (blockGeometryStructure.getNeighborhoodRadius({iX, iY, iZ}) >= margin
        && indicator(iX, iY, iZ)) {
      discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ, true);

      if (discreteNormal[0] == 0) {
        Dynamics<T, DESCRIPTOR>* dynamics = nullptr;
        if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
          if (latticeNumber == 1) {//set PressureBoundary momenta and dynamics on indicated cells
            if (type == "density") {
              dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
                MixinDynamics,momenta::RegularizedPressureBoundaryTuple<0,-1>>>();
            }
            else {  //set VelocityBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
                momenta::BasicDirichletVelocityBoundaryTuple<0,-1>
              >>();
            }
          }
          else {
            dynamics = _block.template getDynamics<FreeEnergyInletOutletDynamics<T,DESCRIPTOR,0,-1>>();
          }
        }

        else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
          if (latticeNumber == 1) {
            if (type == "density") {//set PressureBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
                MixinDynamics,momenta::RegularizedPressureBoundaryTuple<0,1>>>();
            }
            else {  //set VelocityBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
                momenta::BasicDirichletVelocityBoundaryTuple<0,1>
              >>();
            }
          }
          else {
            dynamics = _block.template getDynamics<FreeEnergyInletOutletDynamics<T,DESCRIPTOR,0,1>>();
          }
        }

        else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
          if (latticeNumber == 1) {
            if (type == "density") {//set PressureBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
                MixinDynamics,momenta::RegularizedPressureBoundaryTuple<1,-1>>>();
            }
            else {  //set VelocityBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
                momenta::BasicDirichletVelocityBoundaryTuple<1,-1>
              >>();
            }
          }
          else {
            dynamics = _block.template getDynamics<FreeEnergyInletOutletDynamics<T,DESCRIPTOR,1,-1>>();
          }
        }

        else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
          if (latticeNumber == 1) {
            if (type == "density") {//set PressureBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
                MixinDynamics,momenta::RegularizedPressureBoundaryTuple<1,1>>>();
            }
            else {  //set VelocityBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
                momenta::BasicDirichletVelocityBoundaryTuple<1,1>
              >>();
            }
          }
          else {
            dynamics = _block.template getDynamics<FreeEnergyInletOutletDynamics<T,DESCRIPTOR,1,1>>();
          }
        }

        else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
          if (latticeNumber == 1) {
            if (type == "density") {//set PressureBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
                MixinDynamics,momenta::RegularizedPressureBoundaryTuple<2,-1>>>();
            }
            else {  //set VelocityBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
                momenta::BasicDirichletVelocityBoundaryTuple<2,-1>
              >>();
            }
          }
          else {
            dynamics = _block.template getDynamics<FreeEnergyInletOutletDynamics<T,DESCRIPTOR,2,-1>>();
          }
        }

        else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
          if (latticeNumber == 1) {
            if (type == "density") {//set PressureBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<CombinedRLBdynamics<T,DESCRIPTOR,
                MixinDynamics,momenta::RegularizedPressureBoundaryTuple<2,1>>>();
            }
            else {  //set VelocityBoundary momenta and dynamics on indicated cells
              dynamics = _block.template getDynamics<typename MixinDynamics::template exchange_momenta<
                momenta::BasicDirichletVelocityBoundaryTuple<2,1>
              >>();
            }
          }
          else {
            dynamics = _block.template getDynamics<FreeEnergyInletOutletDynamics<T,DESCRIPTOR,2,1>>();
          }
        }
        dynamics->getParameters(_block).template set<descriptors::OMEGA>(omega);
        if (latticeNumber != 1) {
          _block.defineDynamics({iX, iY, iZ}, dynamics);
          auto cell = _block.get(iX,iY,iZ);
          dynamics->initialize(cell);
        }
      }

      PostProcessorGenerator3D<T, DESCRIPTOR>* chemPotPostProcessor =
        new FreeEnergyChemPotBoundaryProcessorGenerator3D<T, DESCRIPTOR> ( iX, iX, iY, iY, iZ, iZ,
            discreteNormal[1], discreteNormal[2], discreteNormal[3], latticeNumber );
      if (chemPotPostProcessor) {
        _block.addPostProcessor(*chemPotPostProcessor);
      }

      if (_output) {
        clout << "setFreeEnergyInletBoundary<" << "," << ">("  << iX << ", "<< iY << ", " << iZ << ")" << std::endl;
      }
    }
  });
}



}//namespace olb
#endif
