/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2016 Jonas Kratzke, Mathias J. Krause
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

#ifndef OFF_BOUNDARY_CONDITION_2D_H
#define OFF_BOUNDARY_CONDITION_2D_H

#include <list>
#include "core/blockLatticeStructure2D.h"
#include "core/blockLatticeStructure2D.h"
#include "offBoundaryCondition2D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "io/stlReader.h"
namespace olb {

/**
* This class provides a general off lattice boundary condition
*/

template<typename T, template<typename U> class Lattice>
class OffLatticeBoundaryCondition2D {

protected:
  T _epsFraction;

public:
  virtual ~OffLatticeBoundaryCondition2D() { }

  /// Using Bouzidi BC OnePoint corresponds to Bounce Back and TwoPoint to linear interpolation
  virtual void addOnePointZeroVelocityBoundary(int iX, int iY, int iPop, T dist) =0;
  virtual void addTwoPointZeroVelocityBoundary(int iX, int iY, int iPop, T dist) =0;
  virtual void addOnePointVelocityBoundary(int iX, int iY, int iPop, T dist) =0;
  virtual void addTwoPointVelocityBoundary(int iX, int iY, int iPop, T dist) =0;


  virtual void addOffDynamics(int iX, int iY, T location[Lattice<T>::d]) =0;
  virtual void addOffDynamics(int iX, int iY, T location[Lattice<T>::d], T distances[Lattice<T>::q]) =0;
  virtual void addOffDynamics(BlockGeometryStructure2D<T>& blockGeometryStructure, int material) =0;


  virtual void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, int iPop, T dist) =0;
  virtual void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) =0;

  //virtual void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator) =0;

  //virtual void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist) =0;
  //virtual void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int z, T distances[Lattice<T>::q]) =0;
  virtual void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) =0;

  virtual void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) =0;


  virtual void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::list<int> bulkMaterials = std::list<int>(1,1)) =0;

  virtual void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::list<int> bulkMaterials = std::list<int>(1,1)) =0;

  virtual void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, std::list<int> bulkMaterials = std::list<int>(1,1)) =0;


  virtual void defineU(int iX, int iY, int iPop, const T u[Lattice<T>::d]) =0;
  virtual void defineU(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, AnalyticalF2D<T,T>& u, std::list<int> bulkMaterials = std::list<int>(1,1) ) =0;
  virtual void defineRho(int iX, int iY, int iPop, const T rho) =0;
  virtual void defineRho(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, AnalyticalF2D<T,T>& rho, std::list<int> bulkMaterials = std::list<int>(1,1) ) =0;

  virtual BlockLatticeStructure2D<T,Lattice>& getBlock() =0;
  virtual BlockLatticeStructure2D<T,Lattice> const& getBlock() const =0;

  virtual void outputOn() =0;
  virtual void outputOff() =0;
};

////////// Factory functions //////////////////////////////////////////////////

/**
* Create specific off lattice boundary conditions
*/

template<typename T, template<typename U> class Lattice, typename MixinDynamics=BGKdynamics<T,Lattice> >
OffLatticeBoundaryCondition2D<T,Lattice>*
createBouzidiBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block);


template<typename T, template<typename U> class Lattice>
OffLatticeBoundaryCondition2D<T,Lattice>*
createBounceBackBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block)
{
  return createBounceBackBoundaryCondition2D<T,Lattice>(block);
}

}

#endif
