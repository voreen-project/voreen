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

/** \file A helper for initialising 2D boundaries -- header file.  */
#ifndef BOUNDARY_INSTANTIATOR_2D_H
#define BOUNDARY_INSTANTIATOR_2D_H

#include "boundaryCondition2D.h"
#include "boundaryPostProcessors2D.h"
#include "geometry/blockGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "io/ostreamManager.h"

namespace olb {

template<typename T, template<typename U> class Lattice, class BoundaryManager>
class BoundaryConditionInstantiator2D: public OnLatticeBoundaryCondition2D<T,
  Lattice> {
public:
  BoundaryConditionInstantiator2D(BlockLatticeStructure2D<T, Lattice>& block_);
  ~BoundaryConditionInstantiator2D() override;

  void addVelocityBoundary0N(int x0, int x1, int y0, int y1, T omega) override;
  void addVelocityBoundary0P(int x0, int x1, int y0, int y1, T omega) override;
  void addVelocityBoundary1N(int x0, int x1, int y0, int y1, T omega) override;
  void addVelocityBoundary1P(int x0, int x1, int y0, int y1, T omega) override;

  void addSlipBoundary(int x0, int x1, int y0, int y1, int discreteNormalX, int discreteNormalY);

  void addPressureBoundary0N(int x0, int x1, int y0, int y1, T omega) override;
  void addPressureBoundary0P(int x0, int x1, int y0, int y1, T omega) override;
  void addPressureBoundary1N(int x0, int x1, int y0, int y1, T omega) override;
  void addPressureBoundary1P(int x0, int x1, int y0, int y1, T omega) override;

  void addConvectionBoundary0N(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary0P(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary1N(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary1P(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) override;

  void addExternalVelocityCornerNN(int x, int y, T omega) override;
  void addExternalVelocityCornerNP(int x, int y, T omega) override;
  void addExternalVelocityCornerPN(int x, int y, T omega) override;
  void addExternalVelocityCornerPP(int x, int y, T omega) override;

  void addInternalVelocityCornerNN(int x, int y, T omega) override;
  void addInternalVelocityCornerNP(int x, int y, T omega) override;
  void addInternalVelocityCornerPN(int x, int y, T omega) override;
  void addInternalVelocityCornerPP(int x, int y, T omega) override;

  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega) override;
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega) override;
  void addSlipBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1) override;
  void addSlipBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material) override;
  void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega) override;
  void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega) override;
  void addConvectionBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, T* uAv=NULL) override;

  void outputOn() override;
  void outputOff() override;

  BlockLatticeStructure2D<T, Lattice>& getBlock() override;
  BlockLatticeStructure2D<T, Lattice> const& getBlock() const override;
private:
  template<int direction, int orientation>
  void addVelocityBoundary(int x0, int x1, int y0, int y1, T omega);
  template<int direction, int orientation>
  void addPressureBoundary(int x0, int x1, int y0, int y1, T omega);
  template<int direction, int orientation>
  void addConvectionBoundary(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL);
  template<int normalX, int normalY>
  void addExternalVelocityCorner(int x, int y, T omega);
  template<int normalX, int normalY>
  void addInternalVelocityCorner(int x, int y, T omega);
private:
  BlockLatticeStructure2D<T, Lattice>& block;
  std::vector<Momenta<T, Lattice>*> momentaVector;
  std::vector<Dynamics<T, Lattice>*> dynamicsVector;
  bool _output;
  mutable OstreamManager clout;
};

///////// class BoundaryConditionInstantiator2D ////////////////////////

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::BoundaryConditionInstantiator2D(
  BlockLatticeStructure2D<T, Lattice>& block_) :
  block(block_), _output(false), clout(std::cout,"BoundaryConditionInstantiator2D")
{
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::~BoundaryConditionInstantiator2D()
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
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary(
  int x0, int x1, int y0, int y1, T omega)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      Momenta<T, Lattice>* momenta =
        BoundaryManager::template getVelocityBoundaryMomenta<
          direction, orientation>();
      Dynamics<T, Lattice>* dynamics =
        BoundaryManager::template getVelocityBoundaryDynamics<
          direction, orientation>(omega, *momenta);
      this->getBlock().defineDynamics(iX, iX, iY, iY, dynamics);
      momentaVector.push_back(momenta);
      dynamicsVector.push_back(dynamics);
      if (_output) {
        clout << "addVelocityBoundary<" << direction << ","<< orientation << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << omega << " )" << std::endl;
      }
    }
  }

  PostProcessorGenerator2D<T, Lattice>* postProcessor =
    BoundaryManager::template getVelocityBoundaryProcessor<direction,
        orientation>(x0, x1, y0, y1);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addSlipBoundary(
  int x0, int x1, int y0, int y1, int discreteNormalX, int discreteNormalY)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (_output) {
        clout << "addSlipBoundary<" << discreteNormalX << ","<< discreteNormalY << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << " )" << std::endl;
      }
    }
  }

  PostProcessorGenerator2D<T, Lattice>* postProcessor = new SlipBoundaryProcessorGenerator2D<T, Lattice>(x0, x1, y0, y1, discreteNormalX, discreteNormalY);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary(
  int x0, int x1, int y0, int y1, T omega)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      Momenta<T, Lattice>* momenta =
        BoundaryManager::template getPressureBoundaryMomenta<
          direction, orientation>();
      Dynamics<T, Lattice>* dynamics =
        BoundaryManager::template getPressureBoundaryDynamics<
          direction, orientation>(omega, *momenta);
      this->getBlock().defineDynamics(iX, iX, iY, iY, dynamics);
      momentaVector.push_back(momenta);
      dynamicsVector.push_back(dynamics);
      if (_output) {
        clout << "addPressureBoundary<" << direction << ","<< orientation << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << omega << " )" << std::endl;
      }
    }
  }

  PostProcessorGenerator2D<T, Lattice>* postProcessor =
    BoundaryManager::template getPressureBoundaryProcessor<direction,
        orientation>(x0, x1, y0, y1);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addConvectionBoundary(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (_output) {
        clout << "addConvectionBoundary<" << direction << ","<< orientation << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << omega << " )" << std::endl;
      }
    }
  }

  PostProcessorGenerator2D<T, Lattice>* postProcessor =
    BoundaryManager::template getConvectionBoundaryProcessor<direction,orientation>(x0, x1, y0, y1, uAv);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addExternalVelocityCorner(
  int x, int y, T omega)
{
  Momenta<T, Lattice>* momenta =
    BoundaryManager::template getExternalVelocityCornerMomenta<xNormal,
        yNormal>();
  Dynamics<T, Lattice>* dynamics =
    BoundaryManager::template getExternalVelocityCornerDynamics<
      xNormal, yNormal>(omega, *momenta);

  this->getBlock().defineDynamics(x, x, y, y, dynamics);
  momentaVector.push_back(momenta);
  dynamicsVector.push_back(dynamics);

  PostProcessorGenerator2D<T, Lattice>* postProcessor =
    BoundaryManager::template getExternalVelocityCornerProcessor<
      xNormal, yNormal>(x, y);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
  if (_output) {
    clout << "addExternalVelocityCorner<" << xNormal << ","<< yNormal << ">("  << x << ", "<< y << omega << " )" << std::endl;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
template<int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addInternalVelocityCorner(
  int x, int y, T omega)
{
  Momenta<T, Lattice>* momenta =
    BoundaryManager::template getInternalVelocityCornerMomenta<xNormal,
        yNormal>();
  Dynamics<T, Lattice>* dynamics =
    BoundaryManager::template getInternalVelocityCornerDynamics<
      xNormal, yNormal>(omega, *momenta);

  this->getBlock().defineDynamics(x, x, y, y, dynamics);
  momentaVector.push_back(momenta);
  dynamicsVector.push_back(dynamics);

  PostProcessorGenerator2D<T, Lattice>* postProcessor =
    BoundaryManager::template getInternalVelocityCornerProcessor<
      xNormal, yNormal>(x, y);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
  if (_output) {
    clout << "addInternalVelocityCorner<" << xNormal << ","<< yNormal << ">("  << x << ", "<< y << omega << " )" << std::endl;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega)
{
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {

      if (blockGeometryStructure.get(iX, iY)
          == material) {
        discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {

          if (discreteNormal[1] == 1) {
            addVelocityBoundary<0, 1> (iX, iX, iY, iY, omega);
          } else if (discreteNormal[1] == -1) {
            addVelocityBoundary<0, -1> (iX, iX, iY, iY, omega);
          } else if (discreteNormal[2] == 1) {
            addVelocityBoundary<1, 1> (iX, iX, iY, iY, omega);
          } else if (discreteNormal[2] == -1) {
            addVelocityBoundary<1, -1> (iX, iX, iY, iY, omega);
          } else {
            clout << "Could not addVelocityBoundary (" << iX
                  << ", " << iY << ")" << std::endl;
          }
        } else if (discreteNormal[0] == 1) {
          if (discreteNormal[1] == 1) {
            if (discreteNormal[2] == 1) {
              addExternalVelocityCorner<1, 1> (iX, iY, omega);
            } else if (discreteNormal[2] == -1) {
              addExternalVelocityCorner<1, -1> (iX, iY, omega);
            } else {
              clout << "Could not addVelocityBoundary ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          } else if (discreteNormal[1] == -1) {
            if (discreteNormal[2] == 1) {
              addExternalVelocityCorner<-1, 1> (iX, iY, omega);
            } else if (discreteNormal[2] == -1) {
              addExternalVelocityCorner<-1, -1> (iX, iY, omega);
            } else {
              clout << "Could not addVelocityBoundary ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          }
        } else if (discreteNormal[0] == 2) {
          if (discreteNormal[1] == 1) {
            if (discreteNormal[2] == 1) {
              addInternalVelocityCorner<1, 1> (iX, iY, omega);
            } else if (discreteNormal[2] == -1) {
              addInternalVelocityCorner<1, -1> (iX, iY, omega);
            } else {
              clout << "Could not addVelocityBoundary ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          } else if (discreteNormal[1] == -1) {
            if (discreteNormal[2] == 1) {
              addInternalVelocityCorner<-1, 1> (iX, iY, omega);
            } else if (discreteNormal[2] == -1) {
              addInternalVelocityCorner<-1, -1> (iX, iY, omega);
            } else {
              clout << "Could not addVelocityBoundary ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addSlipBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1)
{
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {

      if (blockGeometryStructure.get(iX, iY) == material) {
        discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY);
        if (discreteNormal[1]!=0 || discreteNormal[2]!=0) {
          addSlipBoundary(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);
        } else {
          clout << "Warning: Could not addSlipBoundary (" << iX << ", " << iY << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<"), set to bounceBack" << std::endl;
          this->getBlock().defineDynamics(iX, iY, &instances::getBounceBack<T, Lattice>() );
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega)
{
  addVelocityBoundary(blockGeometryStructure, material, 0,
                      blockGeometryStructure.getNx()-1, 0,
                      blockGeometryStructure.getNy()-1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addSlipBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material)
{
  addSlipBoundary(blockGeometryStructure, material, 0,
                  blockGeometryStructure.getNx()-1, 0,
                  blockGeometryStructure.getNy()-1);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega)
{
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      if (blockGeometryStructure.get(iX, iY)
          == material) {
        discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {
          if (discreteNormal[1] == -1) {
            addPressureBoundary0N(iX, iX, iY, iY, omega);
          } else if (discreteNormal[1] == 1) {
            addPressureBoundary0P(iX, iX, iY, iY, omega);
          } else if (discreteNormal[2] == -1) {
            addPressureBoundary1N(iX, iX, iY, iY, omega);
          } else if (discreteNormal[2] == 1) {
            addPressureBoundary1P(iX, iX, iY, iY, omega);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega)
{
  addPressureBoundary(blockGeometryStructure, material, 0,
                      blockGeometryStructure.getNx()-1, 0,
                      blockGeometryStructure.getNy()-1, omega);
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addConvectionBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1,
  T omega, T* uAv)
{
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      if (blockGeometryStructure.get(iX, iY)
          == material) {
        discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {
          if (discreteNormal[1] == -1) {
            addConvectionBoundary0N(iX, iX, iY, iY, omega, uAv);
          } else if (discreteNormal[1] == 1) {
            addConvectionBoundary0P(iX, iX, iY, iY, omega, uAv);
          } else if (discreteNormal[2] == -1) {
            addConvectionBoundary1N(iX, iX, iY, iY, omega, uAv);
          } else if (discreteNormal[2] == 1) {
            addConvectionBoundary1P(iX, iX, iY, iY, omega, uAv);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addConvectionBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, T* uAv)
{
  addConvectionBoundary(blockGeometryStructure, material, 0,
                        blockGeometryStructure.getNx()-1, 0,
                        blockGeometryStructure.getNy()-1, omega, uAv);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary0N(
  int x0, int x1, int y0, int y1, T omega)
{
  addVelocityBoundary<0, -1> (x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary0P(
  int x0, int x1, int y0, int y1, T omega)
{
  addVelocityBoundary<0, 1> (x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary1N(
  int x0, int x1, int y0, int y1, T omega)
{
  addVelocityBoundary<1, -1> (x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary1P(
  int x0, int x1, int y0, int y1, T omega)
{
  addVelocityBoundary<1, 1> (x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary0N(
  int x0, int x1, int y0, int y1, T omega)
{
  addPressureBoundary<0, -1> (x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary0P(
  int x0, int x1, int y0, int y1, T omega)
{
  addPressureBoundary<0, 1> (x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary1N(
  int x0, int x1, int y0, int y1, T omega)
{
  addPressureBoundary<1, -1> (x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary1P(
  int x0, int x1, int y0, int y1, T omega)
{
  addPressureBoundary<1, 1> (x0, x1, y0, y1, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addConvectionBoundary0N(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  addConvectionBoundary<0, -1> (x0, x1, y0, y1, omega, uAv);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addConvectionBoundary0P(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  addConvectionBoundary<0, 1> (x0, x1, y0, y1, omega, uAv);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addConvectionBoundary1N(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  addConvectionBoundary<1, -1> (x0, x1, y0, y1, omega, uAv);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addConvectionBoundary1P(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  addConvectionBoundary<1, 1> (x0, x1, y0, y1, omega, uAv);
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addExternalVelocityCornerNN(
  int x, int y, T omega)
{
  addExternalVelocityCorner<-1, -1> (x, y, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addExternalVelocityCornerNP(
  int x, int y, T omega)
{
  addExternalVelocityCorner<-1, 1> (x, y, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addExternalVelocityCornerPN(
  int x, int y, T omega)
{
  addExternalVelocityCorner<1, -1> (x, y, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addExternalVelocityCornerPP(
  int x, int y, T omega)
{
  addExternalVelocityCorner<1, 1> (x, y, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addInternalVelocityCornerNN(
  int x, int y, T omega)
{
  addInternalVelocityCorner<-1, -1> (x, y, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addInternalVelocityCornerNP(
  int x, int y, T omega)
{
  addInternalVelocityCorner<-1, 1> (x, y, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addInternalVelocityCornerPN(
  int x, int y, T omega)
{
  addInternalVelocityCorner<1, -1> (x, y, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addInternalVelocityCornerPP(
  int x, int y, T omega)
{
  addInternalVelocityCorner<1, 1> (x, y, omega);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BlockLatticeStructure2D<T, Lattice>& BoundaryConditionInstantiator2D<T, Lattice,
                        BoundaryManager>::getBlock()
{
  return block;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BlockLatticeStructure2D<T, Lattice> const& BoundaryConditionInstantiator2D<T, Lattice,
                        BoundaryManager>::getBlock() const
{
  return block;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::outputOn()
{
  _output = true;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::outputOff()
{
  _output = false;
}

}

#endif
