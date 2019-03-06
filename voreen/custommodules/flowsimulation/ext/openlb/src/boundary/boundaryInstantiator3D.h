/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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

/** \file A helper for initialising 3D boundaries -- header file.  */
#ifndef BOUNDARY_INSTANTIATOR_3D_H
#define BOUNDARY_INSTANTIATOR_3D_H

#include "boundaryCondition3D.h"
#include "boundaryPostProcessors3D.h"
#include "io/ostreamManager.h"

namespace olb {


template<typename T, template<typename U> class Lattice, class BoundaryManager>
class BoundaryConditionInstantiator3D : public OnLatticeBoundaryCondition3D<T,Lattice> {
public:
  BoundaryConditionInstantiator3D( BlockLatticeStructure3D<T,Lattice>& block_ );
  ~BoundaryConditionInstantiator3D() override;

  void addVelocityBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;

  void addSlipBoundary(int x0, int x1, int y0, int y1, int z0, int z1, int discreteNormalX, int discreteNormalY, int discreteNormalZ);

  void addPressureBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;

  void addConvectionBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;

  void addExternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;

  void addInternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;

  void addExternalVelocityCornerNNN(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerNNP(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerNPN(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerNPP(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerPNN(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerPNP(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerPPN(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerPPP(int x, int y, int z, T omega) override;

  void addInternalVelocityCornerNNN(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerNNP(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerNPN(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerNPP(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerPNN(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerPNP(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerPPN(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerPPP(int x, int y, int z, T omega) override;

  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1,
                           T omega) override;
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega) override;

  void addSlipBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1) override;
  void addSlipBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material) override;

  void addPressureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1,
                           T omega) override;
  void addPressureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega) override;

  void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1,
                             T omega, T* uAv=NULL) override;
  void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega, T* uAv=NULL) override;

  void outputOn() override;
  void outputOff() override;

private:
  template<int direction, int orientation>
  void addVelocityBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int direction, int orientation>
  void addPressureBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int direction, int orientation>
  void addConvectionBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL);
  template<int plane, int normal1, int normal2>
  void addExternalVelocityEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int plane, int normal1, int normal2>
  void addInternalVelocityEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int normalX, int normalY, int normalZ>
  void addExternalVelocityCorner(int x, int y, int z, T omega);
  template<int normalX, int normalY, int normalZ>
  void addInternalVelocityCorner(int x, int y, int z, T omega);
private:
  BlockLatticeStructure3D<T,Lattice>& _block;
  std::vector<Momenta<T,Lattice>*>  momentaVector;
  std::vector<Dynamics<T,Lattice>*> dynamicsVector;
  bool _output;
  mutable OstreamManager clout;
};


///////// class BoundaryConditionInstantiator3D ////////////////////////

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::BoundaryConditionInstantiator3D (
  BlockLatticeStructure3D<T,Lattice>& block)
  : _block(block), _output(false), clout(std::cout,"BoundaryConditionInstantiator3D")
{ }

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::~BoundaryConditionInstantiator3D()
{
  for (auto &iDynamics : dynamicsVector) {
    delete iDynamics;
  }
  for (auto &iMomenta : momentaVector) {
    delete iMomenta;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addVelocityBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  OLB_PRECONDITION( x0==x1 || y0==y1 || z0==z1 );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,Lattice>* momenta
          = BoundaryManager::template getVelocityBoundaryMomenta<direction,orientation>();
        Dynamics<T,Lattice>* dynamics
          = BoundaryManager::template getVelocityBoundaryDynamics<direction,orientation>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
        if (_output) {
          clout << "addVelocityBoundary<" << direction << ", " << orientation << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,Lattice>* postProcessor
    = BoundaryManager::template getVelocityBoundaryProcessor<direction,orientation>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

// Slip BC

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addSlipBoundary(
  int x0, int x1, int y0, int y1, int z0, int z1, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (_output) {
          clout << "addSlipBoundary<" << discreteNormalX << ","<< discreteNormalY << ","<< discreteNormalZ << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T, Lattice>* postProcessor = new SlipBoundaryProcessorGenerator3D<T, Lattice>(x0, x1, y0, y1, z0, z1, discreteNormalX, discreteNormalY, discreteNormalZ);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

// Pressure BC

template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addPressureBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  OLB_PRECONDITION( x0==x1 || y0==y1 || z0==z1 );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,Lattice>* momenta
          = BoundaryManager::template getPressureBoundaryMomenta<direction,orientation>();
        Dynamics<T,Lattice>* dynamics
          = BoundaryManager::template getPressureBoundaryDynamics<direction,orientation>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
        if (_output) {
          clout << "addPressureBoundary<" << direction << ", " << orientation << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,Lattice>* postProcessor
    = BoundaryManager::template getPressureBoundaryProcessor<direction,orientation>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

// Convection BC

template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addConvectionBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  OLB_PRECONDITION( x0==x1 || y0==y1 || z0==z1 );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        if (_output) {
          clout << "addConvectionBoundary<" << direction << ", " << orientation << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,Lattice>* postProcessor
    = BoundaryManager::template getConvectionBoundaryProcessor<direction,orientation>(x0,x1, y0,y1, z0,z1, uAv);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int plane, int normal1, int normal2>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  OLB_PRECONDITION(
    ( x0==x1 && y0==y1 ) ||
    ( x0==x1 && z0==z1 ) ||
    ( y0==y1 && z0==z1 ) );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,Lattice>* momenta
          = BoundaryManager::template getExternalVelocityEdgeMomenta<plane,normal1,normal2>();
        Dynamics<T,Lattice>* dynamics
          = BoundaryManager::template getExternalVelocityEdgeDynamics<plane,normal1,normal2>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
        if (_output) {
          clout << "addExternalVelocityEdge<" << plane << ", " << normal1 << ", " << normal2 << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,Lattice>* postProcessor
    = BoundaryManager::template getExternalVelocityEdgeProcessor<plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int plane, int normal1, int normal2>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  if (!(( x0==x1 && y0==y1 ) ||
        ( x0==x1 && z0==z1 ) ||
        ( y0==y1 && z0==z1 ) )) {
    clout << x0 <<" "<< x1 <<" "<< y0 <<" "<< y1 <<" "<< z0 <<" "<< z1 << std::endl;
  }

  OLB_PRECONDITION(
    ( x0==x1 && y0==y1 ) ||
    ( x0==x1 && z0==z1 ) ||
    ( y0==y1 && z0==z1 ) );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,Lattice>* momenta
          = BoundaryManager::template getInternalVelocityEdgeMomenta<plane,normal1,normal2>();
        Dynamics<T,Lattice>* dynamics
          = BoundaryManager::template getInternalVelocityEdgeDynamics<plane,normal1,normal2>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
        if (_output) {
          clout << "addInternalVelocityEdge<" << plane << ", " << normal1 << ", " << normal2 << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,Lattice>* postProcessor
    = BoundaryManager::template getInternalVelocityEdgeProcessor<plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int xNormal, int yNormal, int zNormal>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityCorner(int x, int y, int z, T omega)
{
  Momenta<T,Lattice>* momenta
    = BoundaryManager::template getExternalVelocityCornerMomenta<xNormal,yNormal,zNormal>();
  Dynamics<T,Lattice>* dynamics
    = BoundaryManager::template getExternalVelocityCornerDynamics<xNormal,yNormal,zNormal>(omega, *momenta);

  _block.defineDynamics(x,x,y,y,z,z, dynamics);

  momentaVector.push_back(momenta);
  dynamicsVector.push_back(dynamics);

  PostProcessorGenerator3D<T,Lattice>* postProcessor
    = BoundaryManager::template getExternalVelocityCornerProcessor<xNormal,yNormal,zNormal>(x, y, z);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
  if (_output) {
    clout << "addExternalVelocityCorner<" << xNormal << ", " << yNormal << ", " << zNormal << ">(" << x << ", " << y << ", "<< z << omega << " )" << std::endl;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int xNormal, int yNormal, int zNormal>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityCorner(int x, int y, int z, T omega)
{
  Momenta<T,Lattice>* momenta
    = BoundaryManager::template getInternalVelocityCornerMomenta<xNormal,yNormal,zNormal>();
  Dynamics<T,Lattice>* dynamics
    = BoundaryManager::template getInternalVelocityCornerDynamics<xNormal,yNormal,zNormal>(omega, *momenta);

  _block.defineDynamics(x,x,y,y,z,z, dynamics);

  momentaVector.push_back(momenta);
  dynamicsVector.push_back(dynamics);

  PostProcessorGenerator3D<T,Lattice>* postProcessor
    = BoundaryManager::template getInternalVelocityCornerProcessor<xNormal,yNormal,zNormal>(x, y, z);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
  if (_output) {
    clout << "addInternalVelocityCorner<" << xNormal << ", " << yNormal << ", " << zNormal << ">(" << x << ", " << y << ", "<< z << omega << " )" << std::endl;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  std::vector<int> discreteNormal(4,0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {

        if (blockGeometryStructure.getMaterial(iX, iY, iZ)==material) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX,iY,iZ);
          if (discreteNormal[0] == 0) {
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
              addVelocityBoundary<0,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              addVelocityBoundary<0,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              addVelocityBoundary<1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              addVelocityBoundary<1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              addVelocityBoundary<2,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              addVelocityBoundary<2,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }

          else if (discreteNormal[0] == 1) {
            if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addExternalVelocityCorner<1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addExternalVelocityCorner<1,-1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<-1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<-1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addExternalVelocityCorner<-1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addExternalVelocityCorner<-1,-1,-1>(iX,iY,iZ, omega);
            }
            ///                     addExternalVelocityCorner<discreteNormal[1],discreteNormal[2],discreteNormal[3]>(iX,iY,iZ, omega);
          }

          else if (discreteNormal[0] == 2) {
            if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addInternalVelocityCorner<1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addInternalVelocityCorner<1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addInternalVelocityCorner<1,-1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addInternalVelocityCorner<-1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addInternalVelocityCorner<-1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addInternalVelocityCorner<-1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addInternalVelocityCorner<-1,-1,-1>(iX,iY,iZ, omega);
            }
            ///                     addInternalVelocityCorner<discreteNormal[1],discreteNormal[2],discreteNormal[3]>(iX,iY,iZ, omega);
          }

          else if (discreteNormal[0] == 3) {
            if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addExternalVelocityEdge<0,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addExternalVelocityEdge<0,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addExternalVelocityEdge<0,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addExternalVelocityEdge<0,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addExternalVelocityEdge<1,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addExternalVelocityEdge<1,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addExternalVelocityEdge<1,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addExternalVelocityEdge<1,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addExternalVelocityEdge<2,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addExternalVelocityEdge<2,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addExternalVelocityEdge<2,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addExternalVelocityEdge<2,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }

          else if (discreteNormal[0] == 4) {
            if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addInternalVelocityEdge<0,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addInternalVelocityEdge<0,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addInternalVelocityEdge<0,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addInternalVelocityEdge<0,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addInternalVelocityEdge<1,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addInternalVelocityEdge<1,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addInternalVelocityEdge<1,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addInternalVelocityEdge<1,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addInternalVelocityEdge<2,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addInternalVelocityEdge<2,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addInternalVelocityEdge<2,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addInternalVelocityEdge<2,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega)
{
  addVelocityBoundary(blockGeometryStructure, material, 0, blockGeometryStructure.getNx()-1, 0, blockGeometryStructure.getNy()-1, 0, blockGeometryStructure.getNz()-1, omega);
}

// Slip BC

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addSlipBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1)
{
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {
        const BlockGeometryStructure3D<T>& bgs = blockGeometryStructure;
        if (bgs.get(iX, iY, iZ) == material) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
            addSlipBoundary(iX, iX, iY, iY, iZ, iZ, discreteNormal[1], discreteNormal[2], discreteNormal[3]);
          } else {
            clout << "Warning: Could not addSlipBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<","<< discreteNormal[3] <<"), set to bounceBack" << std::endl;
            _block.defineDynamics(iX, iY, iZ, &instances::getBounceBack<T, Lattice>() );
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addSlipBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material)
{
  addSlipBoundary(blockGeometryStructure, material, 0,
                  blockGeometryStructure.getNx()-1, 0,
                  blockGeometryStructure.getNy()-1, 0,
                  blockGeometryStructure.getNz()-1);
}


// Pressure BC

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addPressureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  std::vector<int> discreteNormal(4,0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {

        if (blockGeometryStructure.getMaterial(iX, iY, iZ)==material) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX,iY,iZ);

          if (discreteNormal[0] == 0) {
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
              addPressureBoundary<0,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              addPressureBoundary<0,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              addPressureBoundary<1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              addPressureBoundary<1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              addPressureBoundary<2,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              addPressureBoundary<2,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addPressureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega)
{

  addPressureBoundary(blockGeometryStructure, material, 0, blockGeometryStructure.getNx()-1, 0, blockGeometryStructure.getNy()-1, 0, blockGeometryStructure.getNz()-1, omega);
}

// Convection BC

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  std::vector<int> discreteNormal(4,0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {

        if (blockGeometryStructure.getMaterial(iX, iY, iZ)==material) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX,iY,iZ);

          if (discreteNormal[0] == 0) {
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
              addConvectionBoundary<0,-1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              addConvectionBoundary<0,1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              addConvectionBoundary<1,-1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              addConvectionBoundary<1,1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }

            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              addConvectionBoundary<2,-1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              addConvectionBoundary<2,1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
          }
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega, T* uAv)
{

  addConvectionBoundary(blockGeometryStructure, material, 0, blockGeometryStructure.getNx()-1, 0, blockGeometryStructure.getNy()-1, 0, blockGeometryStructure.getNz()-1, omega, uAv);
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addVelocityBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<0,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addVelocityBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<0,1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addVelocityBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addVelocityBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addVelocityBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<2,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addVelocityBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<2, 1>(x0,x1,y0,y1,z0,z1, omega);
}

// Pressure BC

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addPressureBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<0,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addPressureBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<0,1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addPressureBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addPressureBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addPressureBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<2,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addPressureBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<2, 1>(x0,x1,y0,y1,z0,z1, omega);
}

// Convection BC

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addConvectionBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<0,-1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addConvectionBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<0,1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addConvectionBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<1,-1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addConvectionBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<1, 1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addConvectionBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<2,-1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addConvectionBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<2, 1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}


// Velocity BC

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<0,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<0,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<0, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<0, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<1,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<1,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<1, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<1, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<2,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<2,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<2, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<2, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}



template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<0,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<0,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<0, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<0, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<1,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<1,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<1, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<1, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<2,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<2,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<2, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<2, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityCornerNNN(int x, int y, int z, T omega)
{
  addExternalVelocityCorner<-1,-1,-1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityCornerNNP(int x, int y, int z, T omega)
{
  addExternalVelocityCorner<-1,-1, 1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityCornerNPN(int x, int y, int z, T omega)
{
  addExternalVelocityCorner<-1, 1,-1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityCornerNPP(int x, int y, int z, T omega)
{
  addExternalVelocityCorner<-1, 1, 1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityCornerPNN(int x, int y, int z, T omega)
{
  addExternalVelocityCorner< 1,-1,-1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityCornerPNP(int x, int y, int z, T omega)
{
  addExternalVelocityCorner< 1,-1, 1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityCornerPPN(int x, int y, int z, T omega)
{
  addExternalVelocityCorner< 1, 1,-1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addExternalVelocityCornerPPP(int x, int y, int z, T omega)
{
  addExternalVelocityCorner< 1, 1, 1>(x,y,z, omega);
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityCornerNNN(int x, int y, int z, T omega)
{
  addInternalVelocityCorner<-1,-1,-1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityCornerNNP(int x, int y, int z, T omega)
{
  addInternalVelocityCorner<-1,-1, 1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityCornerNPN(int x, int y, int z, T omega)
{
  addInternalVelocityCorner<-1, 1,-1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityCornerNPP(int x, int y, int z, T omega)
{
  addInternalVelocityCorner<-1, 1, 1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityCornerPNN(int x, int y, int z, T omega)
{
  addInternalVelocityCorner< 1,-1,-1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityCornerPNP(int x, int y, int z, T omega)
{
  addInternalVelocityCorner< 1,-1, 1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityCornerPPN(int x, int y, int z, T omega)
{
  addInternalVelocityCorner< 1, 1,-1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::
addInternalVelocityCornerPPP(int x, int y, int z, T omega)
{
  addInternalVelocityCorner< 1, 1, 1>(x,y,z, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::outputOn()
{
  _output = true;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,Lattice,BoundaryManager>::outputOff()
{
  _output = false;
}


}  // namespace olb

#endif
