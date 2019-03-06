/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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
#ifndef OFF_BOUNDARY_INSTANTIATOR_2D_H
#define OFF_BOUNDARY_INSTANTIATOR_2D_H

#include "offBoundaryCondition2D.h"
#include "offBoundaryPostProcessors2D.h"
#include "geometry/blockGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "dynamics/dynamics.h"
#include "core/cell.h"
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/stlReader.h"

namespace olb {

/**
* This class gets the needed processors from BoundaryManager and adds them
* to the Processor Vector of the Lattice
*/

template<typename T, template<typename U> class Lattice, class BoundaryManager>
class OffBoundaryConditionInstantiator2D: public OffLatticeBoundaryCondition2D<T,
  Lattice> {
public:
  OffBoundaryConditionInstantiator2D(BlockLatticeStructure2D<T, Lattice>& block_, T epsFraction_ = 0.0001);
  ~OffBoundaryConditionInstantiator2D() override;

  void addOnePointZeroVelocityBoundary(int x, int y, int iPop, T dist) override;
  void addTwoPointZeroVelocityBoundary(int x, int y, int iPop, T dist) override;

  void addOnePointVelocityBoundary(int x, int y, int iPop, T dist) override;
  void addTwoPointVelocityBoundary(int x, int y, int iPop, T dist) override;

  void addOffDynamics(int x, int y, T location[Lattice<T>::d]) override;
  void addOffDynamics(int x, int y, T location[Lattice<T>::d], T distances[Lattice<T>::q]) override;
  void addOffDynamics(BlockGeometryStructure2D<T>& blockGeometryStructure, int material) override;

  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist) override;
  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[Lattice<T>::q]);
  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1));
  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) override;
  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::list<int> bulkMaterials = std::list<int>(1,1)) override;

  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist);
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[Lattice<T>::q]);
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1));
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) override;
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::list<int> bulkMaterials = std::list<int>(1,1)) override;

  void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) override;
  void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::list<int> bulkMaterials = std::list<int>(1,1)) override;

  void setBoundaryIntersection(int iX, int iY, int iPop, T distance);
  bool getBoundaryIntersection(int iX, int iY, int iPop, T point[Lattice<T>::d]);

  void defineU(int iX, int iY, int iPop, const T u[Lattice<T>::d]) override;
  void defineU(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, AnalyticalF2D<T,T>& u, std::list<int> bulkMaterials = std::list<int>(1,1) ) override;

  void defineRho(int iX, int iY, int iPop, const T rho) override;
  void defineRho(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, AnalyticalF2D<T,T>& rho, std::list<int> bulkMaterials = std::list<int>(1,1) ) override;

  void outputOn() override;
  void outputOff() override;

  BlockLatticeStructure2D<T, Lattice>& getBlock() override;
  BlockLatticeStructure2D<T, Lattice> const& getBlock() const override;
private:
  BlockLatticeStructure2D<T, Lattice>& block;
  //std::vector<Momenta<T, Lattice>*> momentaVector;
  std::vector<Dynamics<T, Lattice>*> dynamicsVector;
  bool _output;
  mutable OstreamManager clout;
};

///////// class OffBoundaryConditionInstantiator2D ////////////////////////

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BlockLatticeStructure2D<T, Lattice>& OffBoundaryConditionInstantiator2D<T, Lattice,
                        BoundaryManager>::getBlock()
{
  return block;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BlockLatticeStructure2D<T, Lattice> const& OffBoundaryConditionInstantiator2D<T, Lattice,
                        BoundaryManager>::getBlock() const
{
  return block;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::OffBoundaryConditionInstantiator2D(
  BlockLatticeStructure2D<T, Lattice>& block_, T epsFraction_) :
  block(block_), _output(false), clout(std::cout,"BoundaryConditionInstantiator2D")
{
  this->_epsFraction = epsFraction_;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::~OffBoundaryConditionInstantiator2D()
{
  for (unsigned iDynamics = 0; iDynamics < dynamicsVector.size(); ++iDynamics) {
    delete dynamicsVector[iDynamics];
  }
  /*
  for (unsigned iMomenta = 0; iMomenta < dynamicsVector.size(); ++iMomenta) {
    delete momentaVector[iMomenta];
  }*/
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addOnePointZeroVelocityBoundary(
  int x, int y, int iPop, T dist)
{
  PostProcessorGenerator2D<T, Lattice>* postProcessor =
    BoundaryManager::getOnePointZeroVelocityBoundaryProcessor
    (x, y, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addTwoPointZeroVelocityBoundary(
  int x, int y, int iPop, T dist)
{
  PostProcessorGenerator2D<T, Lattice>* postProcessor =
    BoundaryManager::getTwoPointZeroVelocityBoundaryProcessor
    (x, y, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addOnePointVelocityBoundary(
  int x, int y, int iPop, T dist)
{
  PostProcessorGenerator2D<T, Lattice>* postProcessor =
    BoundaryManager::getOnePointVelocityBoundaryProcessor
    (x, y, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addTwoPointVelocityBoundary(
  int x, int y, int iPop, T dist)
{
  PostProcessorGenerator2D<T, Lattice>* postProcessor =
    BoundaryManager::getTwoPointVelocityBoundaryProcessor
    (x, y, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addOffDynamics(
  int x, int y, T location[Lattice<T>::d])
{
  Dynamics<T,Lattice>* dynamics
    = BoundaryManager::getOffDynamics(location);
  this->getBlock().defineDynamics(x,x,y,y, dynamics);
  dynamicsVector.push_back(dynamics);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addOffDynamics(
  int x, int y, T location[Lattice<T>::d], T distances[Lattice<T>::q])
{
  Dynamics<T,Lattice>* dynamics
    = BoundaryManager::getOffDynamics(location, distances);
  this->getBlock().defineDynamics(x,x,y,y, dynamics);
  dynamicsVector.push_back(dynamics);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addOffDynamics(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material)
{
  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];
    for (int ix = x0; ix <= x1; ix++)
      for (int iy = y0; iy <= y1; iy++)
        if (blockGeometryStructure.getMaterial(ix,iy) == material ) {
          T location[Lattice<T>::d];
          blockGeometryStructure.getPhysR(location, ix,iy);
          addOffDynamics(ix, iy, location);
        }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist)
{
  const int* c = Lattice<T>::c[iPop];
  if (blockGeometryStructure.getMaterial(x-c[0], y-c[1]) != 1) {
    addOnePointZeroVelocityBoundary(x, y, iPop, dist);
  } else {
    addTwoPointZeroVelocityBoundary(x, y, iPop, dist);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::
addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[Lattice<T>::q])
{
  typedef Lattice<T> L;
  //T location[Lattice<T>::d];
  //location[0] = blockGeometryStructure.physCoordX(x);
  //location[1] = blockGeometryStructure.physCoordY(y);
  //location[2] = blockGeometryStructure.physCoordZ(z);
  //T distancesCopy[L::q];
  //T spacing = blockGeometryStructure.getDeltaR();
  //for (int iPop = 1; iPop < L::q ; ++iPop) {
  //  distancesCopy[iPop] = spacing*(1.-distances[iPop]);
  //  if (distances[iPop] == -1)
  //    distancesCopy[iPop] = -1;
  //}
  //addOffDynamics(x, y, z, location, distancesCopy);

  for (int iPop = 1; iPop < L::q ; ++iPop) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      const int* c = L::c[iPop];
      addZeroVelocityBoundary(blockGeometryStructure, x-c[0], y-c[1], iPop, distances[iPop]);
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials)
{

  T distances[Lattice<T>::q];
  for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
    const int* c = Lattice<T>::c[iPop];
    const int iXn = iX + c[0];
    const int iYn = iY + c[1];
    std::list<int>::iterator mat;
    for (mat=bulkMaterials.begin(); mat!=bulkMaterials.end(); ++mat) {
      if (blockGeometryStructure.getMaterial(iXn,iYn) == *mat ) {
        T dist = -1;

        T physR[2];
        blockGeometryStructure.getPhysR(physR,iXn,iYn);
        T voxelSize=blockGeometryStructure.getDeltaR();

        Vector<T,2> physC(physR);

        Vector<T,2> direction(-voxelSize*c[0],-voxelSize*c[1]) ;
        T cPhysNorm = voxelSize*sqrt(c[0]*c[0]+c[1]*c[1]);

        if (!indicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
          T epsX = voxelSize*c[0]*this->_epsFraction;
          T epsY = voxelSize*c[1]*this->_epsFraction;

          Vector<T,2> physC2(physC);
          physC2[0] += epsX;
          physC2[1] += epsY;
          Vector<T,2> direction2(direction);
          direction2[0] -= 2.*epsX;
          direction2[1] -= 2.*epsY;

          if ( !indicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
            clout << "ERROR: no boundary found at (" << iXn << "," << iYn <<") ~ ("
                  << physR[0] << "," << physR[1] << "), "
                  << "in direction " << util::opposite<Lattice<T> >(iPop)
                  << std::endl;
            //blockGeometryStructure.printNode(iXn,iYn);
            //  exit(-1);
          }
          T distNew = (dist - sqrt(epsX*epsX+epsY*epsY) )/cPhysNorm;
          if (distNew < 0.5) {
            dist = 0;
          } else {
            dist = 0.5 * cPhysNorm;
            //dist = cPhysNorm;
            clout << "WARNING: distance at (" << iXn << "," << iYn << ") ~ ("
                  << physR[0] << "," << physR[1] << "), "
                  << "in direction " << util::opposite<Lattice<T> >(iPop) << ": "
                  << distNew
                  << " rounded to "
                  << dist/cPhysNorm
                  << std::endl;
          } // else if
        } // if
        distances[util::opposite<Lattice<T> >(iPop)] = dist/cPhysNorm;
      } // bulkMaterials if
    } // bulkMaterials loop
  } // iPop
  addZeroVelocityBoundary(blockGeometryStructure, iX, iY, distances);
} //method

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials)
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        if (blockGeometryStructure.getMaterial(iX,iY) == material ) {
          addZeroVelocityBoundary(blockGeometryStructure, iX, iY, indicator, bulkMaterials);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::list<int> bulkMaterials)
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        if (blockGeometryStructure.getMaterial(iX,iY) == material ) {
          bool streamDirections[Lattice<T>::q];
          // check if (ix,iY) is really a boundary node and compute the stream direction
          if (blockGeometryStructure.template findStreamDirections<T,Lattice>(iX, iY, material, bulkMaterials, streamDirections) ) {
            T distances[Lattice<T>::q];
            for (int iPop=0; iPop<Lattice<T>::q; iPop++) {
              if (streamDirections[iPop]) {
                distances[Lattice<T>::opposite[iPop]] = .5;
              } else {
                distances[Lattice<T>::opposite[iPop]] = -1;
              }
            }
            addZeroVelocityBoundary(blockGeometryStructure, iX, iY, distances);
          }
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist)
{
  const int* c = Lattice<T>::c[iPop];
  if (blockGeometryStructure.getMaterial(x-c[0], y-c[1]) != 1) {
    addOnePointVelocityBoundary(x, y, iPop, dist);
  } else {
    addTwoPointVelocityBoundary(x, y, iPop, dist);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::
addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[Lattice<T>::q])
{
  typedef Lattice<T> L;
  T location[Lattice<T>::d];
  blockGeometryStructure.getPhysR(location, x,y);

  T distancesCopy[L::q];
  T spacing = blockGeometryStructure.getDeltaR();
  for (int iPop = 1; iPop < L::q ; ++iPop) {
    distancesCopy[iPop] = spacing*(1.-distances[iPop]);
    if ( !util::nearZero(distances[iPop]+1) ) {
      distancesCopy[iPop] = -1;
    }
  }
  addOffDynamics(x, y, location, distancesCopy);

  for (int iPop = 1; iPop < L::q ; ++iPop) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      const int* c = L::c[iPop];
      addVelocityBoundary(blockGeometryStructure, x-c[0], y-c[1], iPop, distances[iPop]);
    }
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials)
{

  T distances[Lattice<T>::q];
  for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
    const int* c = Lattice<T>::c[iPop];
    const int iXn = iX + c[0];
    const int iYn = iY + c[1];
    std::list<int>::iterator mat;
    for (mat=bulkMaterials.begin(); mat!=bulkMaterials.end(); ++mat) {
      if (blockGeometryStructure.getMaterial(iXn,iYn) == *mat ) {
        T dist = -1;
        T physR[2];
        blockGeometryStructure.getPhysR(physR,iXn,iYn);
        T voxelSize=blockGeometryStructure.getDeltaR();
        Vector<T,2> physC(physR);

        Vector<T,2> direction(-voxelSize*c[0],-voxelSize*c[1]);
        T cPhysNorm = voxelSize*sqrt(c[0]*c[0]+c[1]*c[1]);

        if (!indicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
          T epsX = voxelSize*c[0]*this->_epsFraction;
          T epsY = voxelSize*c[1]*this->_epsFraction;

          Vector<T,2> physC2(physC);
          physC2[0] += epsX;
          physC2[1] += epsY;
          Vector<T,2> direction2(direction);
          direction2[0] -= 2.*epsX;
          direction2[1] -= 2.*epsY;

          if ( !indicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
            clout << "ERROR: no boundary found at (" << iXn << "," << iYn <<") ~ ("
                  << physR[0] << "," << physR[1] << "), "
                  << "in direction " << util::opposite<Lattice<T> >(iPop)
                  << std::endl;
            //blockGeometryStructure.printNode(iXn,iYn);
            //  exit(-1);
          }
          T distNew = (dist - sqrt(epsX*epsX+epsY*epsY))/cPhysNorm;
          if (distNew < 0.5) {
            dist = 0;
          } else {
            dist = 0.5 * cPhysNorm;
            //dist = cPhysNorm;
            clout << "WARNING: distance at (" << iXn << "," << iYn <<") ~ ("
                  << physR[0] << "," << physR[1] <<"), "
                  << "in direction " << util::opposite<Lattice<T> >(iPop) << ": "
                  << distNew
                  << " rounded to "
                  << dist/cPhysNorm
                  << std::endl;
          } // else if
        } // if
        distances[util::opposite<Lattice<T> >(iPop)] = dist/cPhysNorm;
      } // bulkMaterials if
    } // bulkMaterials loop
  } // iPop
  addVelocityBoundary(blockGeometryStructure, iX, iY, distances);
} //method

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials)
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        if (blockGeometryStructure.getMaterial(iX,iY) == material ) {
          addVelocityBoundary(blockGeometryStructure, iX, iY, indicator, bulkMaterials);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::list<int> bulkMaterials)
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        if (blockGeometryStructure.getMaterial(iX,iY) == material ) {
          bool streamDirections[Lattice<T>::q];
          // check if (ix,iY) is really a boundary node and compute the stream direction
          if (blockGeometryStructure.template findStreamDirections<T,Lattice>(iX, iY, material, bulkMaterials, streamDirections) ) {
            T distances[Lattice<T>::q];
            for (int iPop=0; iPop<Lattice<T>::q; iPop++) {
              if (streamDirections[iPop]) {
                distances[Lattice<T>::opposite[iPop]] = .5;
              } else {
                distances[Lattice<T>::opposite[iPop]] = -1;
              }
            }
            addVelocityBoundary(blockGeometryStructure, iX, iY, distances);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials)
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        bool streamDirections[Lattice<T>::q];
        // check if (ix,iY) is really a boundary node and compute the stream direction
        if (blockGeometryStructure.template findStreamDirections<T,Lattice>(iX, iY, material, bulkMaterials, streamDirections) ) {

          T location[Lattice<T>::d];
          blockGeometryStructure.getPhysR(location, iX, iX);

          block.defineDynamics(iX,iX,iY,iY,&instances::getBounceBackAnti<T, Lattice>(1.));
          //addOffDynamics(iX, iX, location, iPopStream);
          //addStreamBoundary(iX, iY, iPopStream);
        }
      }
    }
  }
  clout << "ERROR: OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary NOT YET IMPLEMENTED" << std::endl;
  exit(-1);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::addPressureBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::list<int> bulkMaterials)
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        bool streamDirections[Lattice<T>::q];
        // check if (ix,iY) is really a boundary node and compute the stream direction
        if (blockGeometryStructure.template findStreamDirections<T,Lattice>(iX, iY, material, bulkMaterials, streamDirections) ) {

          T location[Lattice<T>::d];
          blockGeometryStructure.getPhysR(location, iX, iX);

          block.defineDynamics(iX,iX,iY,iY,&instances::getBounceBackAnti<T, Lattice>(1.));

          //addOffDynamics(iX, iX, location, iPopStream);

          for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
            if (streamDirections[iPop]) {
              //std::cout << Lattice<T>::c[iPop][1] << std::endl;
              PostProcessorGenerator2D<T, Lattice>* postProcessor = new AntiBounceBackPostProcessorGenerator2D<T, Lattice>(iX+Lattice<T>::c[iPop][0], iY+Lattice<T>::c[iPop][1], Lattice<T>::opposite[iPop]);

              this->getBlock().addPostProcessor(*postProcessor);
            }
          }
          //addStreamBoundary(iX, iY, streamDirections);
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::
setBoundaryIntersection(int iX, int iY, int iPop, T distance)
{

  this->getBlock().getDynamics(iX, iY)->setBoundaryIntersection(iPop, distance);
  if (_output) {
    clout << "setBoundaryIntersection(" << iX << ", " << iY << " )" << std::endl;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
bool OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::
getBoundaryIntersection(int iX, int iY, int iPop, T point[Lattice<T>::d])
{

  return this->getBlock().getDynamics(iX, iY)->getBoundaryIntersection(iPop, point);
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::defineU(int iX, int iY, int iPop, const T u[Lattice<T>::d])
{

  this->getBlock().getDynamics(iX, iY)->defineU(iPop, u);
  if (_output) {
    clout << "defineU(" << iX << ", " << iY << " )" << std::endl;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::defineRho(int iX, int iY, int iPop, const T rho)
{

  this->getBlock().getDynamics(iX, iY)->defineRho(iPop, rho);
  if (_output) {
    clout << "defineRho(" << iX << ", " << iY << " )" << std::endl;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::defineU(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, AnalyticalF2D<T,T>& u, std::list<int> bulkMaterials )
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {

        if (blockGeometryStructure.getMaterial(iX, iY) == material) {
          for (int q = 1; q < Lattice<T>::q ; ++q) {
            // Get direction
            const int* c = Lattice<T>::c[q];
            int iXn = iX + c[0];
            int iYn = iY + c[1];
            std::list<int>::iterator i;
            for (i=bulkMaterials.begin(); i!=bulkMaterials.end(); ++i) {
              if (blockGeometryStructure.getMaterial(iXn, iYn) == *i) {
                T intersection[] = { T(), T() }; // coord. of intersection
                int opp = util::opposite<Lattice<T> >(q);
                if (this->getBoundaryIntersection(iX, iY, opp, intersection) ) {
                  //                  std::vector<double> intersection2;
                  //                  intersection2.push_back(intersection[0]);
                  //                  intersection2.push_back(intersection[1]);
                  //                  T vel[] = {u(intersection2)[0], u(intersection2)[1] };
                  T vel[]= {T(),T()};
                  u(vel,intersection);
                  this->defineU(iX, iY, opp, vel);
                }
              }
            }
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::defineRho(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, AnalyticalF2D<T,T>& rho, std::list<int> bulkMaterials )
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {

        if (blockGeometryStructure.getMaterial(iX, iY) == material) {
          for (int q = 1; q < Lattice<T>::q ; ++q) {
            // Get direction
            const int* c = Lattice<T>::c[q];
            int iXn = iX + c[0];
            int iYn = iY + c[1];
            std::list<int>::iterator i;
            for (i=bulkMaterials.begin(); i!=bulkMaterials.end(); ++i) {
              if (blockGeometryStructure.getMaterial(iXn, iYn) == *i) {
                T intersection[] = { T(), T() }; // coord. of intersection
                int opp = util::opposite<Lattice<T> >(q);
                if (this->getBoundaryIntersection(iX, iY, opp, intersection) ) {
                  //                  std::vector<double> intersection2;
                  //                  intersection2.push_back(intersection[0]);
                  //                  intersection2.push_back(intersection[1]);
                  //                  T vel[] = {u(intersection2)[0], u(intersection2)[1] };
                  T rhoLocal[]= {T(1)};
                  rho(rhoLocal,intersection);
                  this->defineRho(iX, iY, opp, rhoLocal[0]);
                }
              }
            }
          }
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::outputOn()
{
  _output = true;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, Lattice, BoundaryManager>::outputOff()
{
  _output = false;
}

}

#endif
