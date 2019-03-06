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

#ifndef OFF_BOUNDARY_CONDITION_3D_H
#define OFF_BOUNDARY_CONDITION_3D_H

#include <list>
#include "core/blockLatticeStructure3D.h"
#include "core/blockLatticeStructure3D.h"
#include "offBoundaryCondition3D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometryStatistics3D.h"
#include "io/stlReader.h"
namespace olb {

/**
* This class provides a general off lattice boundary condition
*/

template<typename T, template<typename U> class Lattice>
class OffLatticeBoundaryCondition3D {

protected:
  T _epsFraction;

public:
  virtual ~OffLatticeBoundaryCondition3D() { }

  /// Using Bouzidi BC OnePoint corresponds to Bounce Back and TwoPoint to linear interpolation
  virtual void addOnePointZeroVelocityBoundary(int x, int y, int z, int iPop, T dist) =0;
  virtual void addTwoPointZeroVelocityBoundary(int x, int y, int z, int iPop, T dist) =0;
  virtual void addOnePointVelocityBoundary(int x, int y, int z, int iPop, T dist) =0;
  virtual void addTwoPointVelocityBoundary(int x, int y, int z, int iPop, T dist) =0;


  virtual void addOffDynamics(int x, int y, int z, T location[Lattice<T>::d]) =0;
  virtual void addOffDynamics(int x, int y, int z, T location[Lattice<T>::d], T distances[Lattice<T>::q]) =0;
  virtual void addOffDynamics(BlockGeometryStructure3D<T>& blockGeometryStructure, int material) =0;


  virtual void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist) =0;
  virtual void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) =0;

  //virtual void addZeroVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, IndicatorF3D<T>& indicator) =0;

  //virtual void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, int iPop, T dist) =0;
  //virtual void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int x, int y, int z, T distances[Lattice<T>::q]) =0;
  virtual void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1)) =0;

  virtual void defineU(int iX, int iY, int iZ, int iPop, const T u[Lattice<T>::d]) =0;
  virtual void defineU(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, AnalyticalF3D<T,T>& u, std::list<int> bulkMaterials = std::list<int>(1,1) ) =0;

  virtual BlockLatticeStructure3D<T,Lattice>& getBlock() =0;
  virtual BlockLatticeStructure3D<T,Lattice> const& getBlock() const =0;

  virtual void outputOn() =0;
  virtual void outputOff() =0;
};

////////// Factory functions //////////////////////////////////////////////////

/**
* Create specific off lattice boundary conditions
*/

template<typename T, template<typename U> class Lattice, typename MixinDynamics=BGKdynamics<T,Lattice> >
OffLatticeBoundaryCondition3D<T,Lattice>*
createBouzidiBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block);

}

#endif
