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

/** \file A helper for initialising 3D boundaries -- header file.  */
#ifndef OFF_BOUNDARY_INSTANTIATOR_3D_H
#define OFF_BOUNDARY_INSTANTIATOR_3D_H

#include "offBoundaryCondition3D.h"
#include "geometry/blockGeometry3D.h"
#include "geometry/blockGeometryStatistics3D.h"
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
class OffBoundaryConditionInstantiator3D: public OffLatticeBoundaryCondition3D<T,
  Lattice> {
public:
  OffBoundaryConditionInstantiator3D(BlockLatticeStructure3D<T, Lattice>& block_, T epsFraction_ = 0.0001);
  ~OffBoundaryConditionInstantiator3D() override;

  void addOnePointZeroVelocityBoundary(int x, int y, int z, int iPop, T dist) override;
  void addTwoPointZeroVelocityBoundary(int x, int y, int z, int iPop, T dist) override;

  void addOnePointVelocityBoundary(int x, int y, int z, int iPop, T dist) override;
  void addTwoPointVelocityBoundary(int x, int y, int z, int iPop, T dist) override;

  void addOffDynamics(int x, int y, int z, T location[Lattice<T>::d]) override;
  void addOffDynamics(int x, int y, int z, T location[Lattice<T>::d], T distances[Lattice<T>::q]) override;
  void addOffDynamics(BlockGeometryStructure3D<T>& blockGeometryStructure, int material) override;

  void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist) override;
  void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, T distances[Lattice<T>::q]);
  void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1));
  void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) override;

  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist);
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, T distances[Lattice<T>::q]);
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1));
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) override;

  void setBoundaryIntersection(int iX, int iY, int iZ, int iPop, T distance);
  bool getBoundaryIntersection(int iX, int iY, int iZ, int iPop, T point[Lattice<T>::d]);

  void defineU(int iX, int iY, int iZ, int iPop, const T u[Lattice<T>::d]) override;
  void defineU(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, AnalyticalF3D<T,T>& u, std::list<int> bulkMaterials = std::list<int>(1,1) ) override;

  void outputOn() override;
  void outputOff() override;

  BlockLatticeStructure3D<T, Lattice>& getBlock() override;
  BlockLatticeStructure3D<T, Lattice> const& getBlock() const override;
private:
  BlockLatticeStructure3D<T, Lattice>& block;
  //std::vector<Momenta<T, Lattice>*> momentaVector;
  std::vector<Dynamics<T, Lattice>*> dynamicsVector;
  bool _output;
  mutable OstreamManager clout;
};

///////// class OffBoundaryConditionInstantiator3D ////////////////////////

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BlockLatticeStructure3D<T, Lattice>& OffBoundaryConditionInstantiator3D<T, Lattice,
                        BoundaryManager>::getBlock()
{
  return block;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
BlockLatticeStructure3D<T, Lattice> const& OffBoundaryConditionInstantiator3D<T, Lattice,
                        BoundaryManager>::getBlock() const
{
  return block;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::OffBoundaryConditionInstantiator3D(
  BlockLatticeStructure3D<T, Lattice>& block_, T epsFraction_) :
  block(block_), _output(false), clout(std::cout,"BoundaryConditionInstantiator3D")
{
  this->_epsFraction = epsFraction_;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::~OffBoundaryConditionInstantiator3D()
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
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addOnePointZeroVelocityBoundary(
  int x, int y, int z, int iPop, T dist)
{
  PostProcessorGenerator3D<T, Lattice>* postProcessor =
    BoundaryManager::getOnePointZeroVelocityBoundaryProcessor
    (x, y, z, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addTwoPointZeroVelocityBoundary(
  int x, int y, int z, int iPop, T dist)
{
  PostProcessorGenerator3D<T, Lattice>* postProcessor =
    BoundaryManager::getTwoPointZeroVelocityBoundaryProcessor
    (x, y, z, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addOnePointVelocityBoundary(
  int x, int y, int z, int iPop, T dist)
{
  PostProcessorGenerator3D<T, Lattice>* postProcessor =
    BoundaryManager::getOnePointVelocityBoundaryProcessor
    (x, y, z, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addTwoPointVelocityBoundary(
  int x, int y, int z, int iPop, T dist)
{
  PostProcessorGenerator3D<T, Lattice>* postProcessor =
    BoundaryManager::getTwoPointVelocityBoundaryProcessor
    (x, y, z, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addOffDynamics(
  int x, int y, int z, T location[Lattice<T>::d])
{
  Dynamics<T,Lattice>* dynamics
    = BoundaryManager::getOffDynamics(location);
  this->getBlock().defineDynamics(x,x,y,y,z,z, dynamics);
  dynamicsVector.push_back(dynamics);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addOffDynamics(
  int x, int y, int z, T location[Lattice<T>::d], T distances[Lattice<T>::q])
{
  Dynamics<T,Lattice>* dynamics
    = BoundaryManager::getOffDynamics(location, distances);
  this->getBlock().defineDynamics(x,x,y,y,z,z, dynamics);
  dynamicsVector.push_back(dynamics);
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addOffDynamics(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material)
{
  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int z0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[2];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];
    const int z1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[2];
    for (int ix = x0; ix <= x1; ix++)
      for (int iy = y0; iy <= y1; iy++)
        for (int iz = z0; iz <= z1; iz++)
          if (blockGeometryStructure.getMaterial(ix,iy,iz) == material ) {
            T location[Lattice<T>::d];
            blockGeometryStructure.getPhysR(location, ix,iy,iz);
            addOffDynamics(ix, iy, iz, location);
          }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist)
{
  const int* c = Lattice<T>::c[iPop];
  if (blockGeometryStructure.getMaterial(x-c[0], y-c[1], z-c[2]) != 1) {
    addOnePointZeroVelocityBoundary(x, y, z, iPop, dist);
  } else {
    addTwoPointZeroVelocityBoundary(x, y, z, iPop, dist);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::
addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, T distances[Lattice<T>::q])
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
      addZeroVelocityBoundary(blockGeometryStructure, x-c[0], y-c[1], z-c[2], iPop, distances[iPop]);
    }
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials)
{

  T distances[Lattice<T>::q];
  for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
    const int* c = Lattice<T>::c[iPop];
    const int iXn = iX + c[0];
    const int iYn = iY + c[1];
    const int iZn = iZ + c[2];
    std::list<int>::iterator mat;
    for (mat=bulkMaterials.begin(); mat!=bulkMaterials.end(); ++mat) {
      if (blockGeometryStructure.getMaterial(iXn,iYn,iZn) == *mat ) {
        T dist = -1;
        T physR[3];
        blockGeometryStructure.getPhysR(physR,iXn,iYn,iZn);
        T voxelSize=blockGeometryStructure.getDeltaR();

        Vector<T,3> physC(physR);

        Vector<T,3> direction(-voxelSize*c[0],-voxelSize*c[1],-voxelSize*c[2]);
        T cPhysNorm = voxelSize*sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

        if (!indicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
          T epsX = voxelSize*c[0]*this->_epsFraction;
          T epsY = voxelSize*c[1]*this->_epsFraction;
          T epsZ = voxelSize*c[2]*this->_epsFraction;

          Vector<T,3> physC2(physC);
          physC2[0] += epsX;
          physC2[1] += epsY;
          physC2[2] += epsZ;
          Vector<T,3> direction2(direction);
          direction2[0] -= 2.*epsX;
          direction2[1] -= 2.*epsY;
          direction2[2] -= 2.*epsZ;

          if ( !indicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
            clout << "ERROR: no boundary found at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                  << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                  << "in direction " << util::opposite<Lattice<T> >(iPop)
                  << std::endl;
            //blockGeometryStructure.printNode(iXn,iYn,iZn);
            //  exit(-1);
          }
          T distNew = (dist - sqrt(epsX*epsX+epsY*epsY+epsZ*epsZ))/cPhysNorm;
          if (distNew < 0.5) {
            dist = 0;
          } else {
            dist = 0.5 * cPhysNorm;
            //dist = cPhysNorm;
            clout << "WARNING: distance at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                  << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
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
  addZeroVelocityBoundary(blockGeometryStructure, iX, iY, iZ, distances);
} //method

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials)
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int z0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[2];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];
    const int z1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[2];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        for (int iZ = z0; iZ <= z1; iZ++) {
          if (blockGeometryStructure.getMaterial(iX,iY,iZ) == material ) {
            addZeroVelocityBoundary(blockGeometryStructure, iX, iY, iZ, indicator, bulkMaterials);
          }
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist)
{
  const int* c = Lattice<T>::c[iPop];
  if (blockGeometryStructure.getMaterial(x-c[0], y-c[1], z-c[2]) != 1) {
    addOnePointVelocityBoundary(x, y, z, iPop, dist);
  } else {
    addTwoPointVelocityBoundary(x, y, z, iPop, dist);
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::
addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, T distances[Lattice<T>::q])
{
  typedef Lattice<T> L;
  T location[Lattice<T>::d];
  blockGeometryStructure.getPhysR(location, x,y,z);
  T distancesCopy[L::q];
  T spacing = blockGeometryStructure.getDeltaR();
  for (int iPop = 1; iPop < L::q ; ++iPop) {
    distancesCopy[iPop] = spacing*(1.-distances[iPop]);
    if ( util::nearZero(distances[iPop]+1) ) {
      distancesCopy[iPop] = -1;
    }
  }
  addOffDynamics(x, y, z, location, distancesCopy);

  for (int iPop = 1; iPop < L::q ; ++iPop) {
    if (!util::nearZero(distances[iPop]+1)) {
      const int* c = L::c[iPop];
      addVelocityBoundary(blockGeometryStructure, x-c[0], y-c[1], z-c[2], iPop, distances[iPop]);
    }
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials)
{

  T distances[Lattice<T>::q];
  for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
    const int* c = Lattice<T>::c[iPop];
    const int iXn = iX + c[0];
    const int iYn = iY + c[1];
    const int iZn = iZ + c[2];
    std::list<int>::iterator mat;
    for (mat=bulkMaterials.begin(); mat!=bulkMaterials.end(); ++mat) {
      if (blockGeometryStructure.getMaterial(iXn,iYn,iZn) == *mat ) {
        T dist = -1;
        T physR[3];
        blockGeometryStructure.getPhysR(physR,iXn,iYn,iZn);
        T voxelSize=blockGeometryStructure.getDeltaR();

        Vector<T,3> physC(physR);
        Vector<T,3> direction(-voxelSize*c[0],-voxelSize*c[1],-voxelSize*c[2]);
        T cPhysNorm = voxelSize*sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

        if (!indicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
          T epsX = voxelSize*c[0]*this->_epsFraction;
          T epsY = voxelSize*c[1]*this->_epsFraction;
          T epsZ = voxelSize*c[2]*this->_epsFraction;

          Vector<T,3> physC2(physC);
          physC2[0] += epsX;
          physC2[1] += epsY;
          physC2[2] += epsZ;
          Vector<T,3> direction2(direction);
          direction2[0] -= 2.*epsX;
          direction2[1] -= 2.*epsY;
          direction2[2] -= 2.*epsZ;

          if ( !indicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
            clout << "ERROR: no boundary found at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                  << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
                  << "in direction " << util::opposite<Lattice<T> >(iPop)
                  << std::endl;
            //blockGeometryStructure.printNode(iXn,iYn,iZn);
            //  exit(-1);
          }
          T distNew = (dist - sqrt(epsX*epsX+epsY*epsY+epsZ*epsZ))/cPhysNorm;
          if (distNew < 0.5) {
            dist = 0;
          } else {
            dist = 0.5 * cPhysNorm;
            //dist = cPhysNorm;
            clout << "WARNING: distance at (" << iXn << "," << iYn << "," << iZn <<") ~ ("
                  << physR[0] << "," << physR[1] << "," << physR[2] <<"), "
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
  addVelocityBoundary(blockGeometryStructure, iX, iY, iZ, distances);
} //method

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure3D<T>& blockGeometryStructure, int material, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials)
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int z0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[2];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];
    const int z1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[2];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        for (int iZ = z0; iZ <= z1; iZ++) {
          if (blockGeometryStructure.getMaterial(iX,iY,iZ) == material ) {
            addVelocityBoundary(blockGeometryStructure, iX, iY, iZ, indicator, bulkMaterials);
          }
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::
setBoundaryIntersection(int iX, int iY, int iZ, int iPop, T distance)
{

  this->getBlock().getDynamics(iX, iY, iZ)->setBoundaryIntersection(iPop, distance);
  if (_output) {
    clout << "setBoundaryIntersection(" << iX << ", " << iY << ", " << iZ << " )" << std::endl;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
bool OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::
getBoundaryIntersection(int iX, int iY, int iZ, int iPop, T point[Lattice<T>::d])
{

  return this->getBlock().getDynamics(iX, iY, iZ)->getBoundaryIntersection(iPop, point);
}


template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::defineU(int iX, int iY, int iZ, int iPop, const T u[Lattice<T>::d])
{

  this->getBlock().getDynamics(iX, iY, iZ)->defineU(iPop, u);
  if (_output) {
    clout << "defineU(" << iX << ", " << iY << ", " << iZ << " )" << std::endl;
  }
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::defineU(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, AnalyticalF3D<T,T>& u, std::list<int> bulkMaterials )
{

  if (blockGeometryStructure.getStatistics().getNvoxel(material)!=0) {
    const int x0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[0];
    const int y0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[1];
    const int z0 = blockGeometryStructure.getStatistics().getMinLatticeR(material)[2];
    const int x1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[0];
    const int y1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[1];
    const int z1 = blockGeometryStructure.getStatistics().getMaxLatticeR(material)[2];

    for (int iX = x0; iX <= x1; iX++) {
      for (int iY = y0; iY <= y1; iY++) {
        for (int iZ = z0; iZ <= z1; iZ++) {

          if (blockGeometryStructure.getMaterial(iX, iY, iZ) == material) {
            for (int q = 1; q < Lattice<T>::q ; ++q) {
              // Get direction
              const int* c = Lattice<T>::c[q];
              int iXn = iX + c[0];
              int iYn = iY + c[1];
              int iZn = iZ + c[2];
              std::list<int>::iterator i;
              for (i=bulkMaterials.begin(); i!=bulkMaterials.end(); ++i) {
                if (blockGeometryStructure.getMaterial(iXn, iYn, iZn) == *i) {
                  T intersection[3] = { T(), T() ,T() }; // coord. of intersection
                  int opp = util::opposite<Lattice<T> >(q);
                  if (this->getBoundaryIntersection(iX, iY, iZ, opp, intersection) ) {
                    //                    std::vector<double> intersection2;
                    //                    intersection2.push_back(intersection[0]);
                    //                    intersection2.push_back(intersection[1]);
                    //                    intersection2.push_back(intersection[2]);
                    //                    T vel[] = {u(intersection2)[0], u(intersection2)[1], u(intersection2)[2]};
                    T vel[3]= {T(),T(),T()};
                    u(vel,intersection);
                    this->defineU(iX, iY, iZ, opp, vel);
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

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::outputOn()
{
  _output = true;
}

template<typename T, template<typename U> class Lattice, class BoundaryManager>
void OffBoundaryConditionInstantiator3D<T, Lattice, BoundaryManager>::outputOff()
{
  _output = false;
}

}

#endif
