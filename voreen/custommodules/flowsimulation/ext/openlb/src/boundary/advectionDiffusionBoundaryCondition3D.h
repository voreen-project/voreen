/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_H
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_H


#include "dynamics/advectionDiffusionDynamics.h"

namespace olb {

template<typename T, template<typename U> class Lattice>
class sOnLatticeBoundaryCondition3D;

/* Interface class whose childs are boundary instantiators.
*  Instantiators have a third template parameter BoundaryManager that finaly implements the interfaces from here.
*/
template<typename T, template<typename U> class Lattice>
class OnLatticeAdvectionDiffusionBoundaryCondition3D {
public:
  virtual ~OnLatticeAdvectionDiffusionBoundaryCondition3D() { }

  /// adds a temperature boundary for one material and a range (x0-x1, y0-y1, z0-z1) or the whole geometry
  virtual void addTemperatureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addTemperatureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega) =0;
  /// adds a convection boundary for one material and a range (x0-x1, y0-y1, z0-z1) or the whole geometry
  virtual void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1) =0;
  virtual void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material) =0;
  /// adds a boundary on the external field for one material and a range (x0-x1, y0-y1, z0-z1) or the whole geometry
  virtual void addExtFieldBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int offset, int x0, int x1, int y0, int y1, int z0, int z1) =0;
  virtual void addExtFieldBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int offset) =0;
  /// adds a boundary that initializes zero distributions and computes the density that entered the boundary for one material and a range (x0-x1, y0-y1, z0-z1) or the whole geometry
  virtual void addZeroDistributionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, int z0, int z1) =0;
  virtual void addZeroDistributionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material) =0;
};

//////  Factory function for Regularized Thermal BC

/// blockLattice creator
template<typename T, template<typename U> class Lattice, typename MixinDynamics=AdvectionDiffusionRLBdynamics<T,Lattice> >
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Lattice>*
createAdvectionDiffusionBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block);

/// superLattice creator, calls createAdvectionDiffusionBoundaryCondidtion3D from above.
template<typename T, template<typename U> class Lattice, typename MixinDynamics=AdvectionDiffusionRLBdynamics<T,Lattice> >
void createAdvectionDiffusionBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, Lattice>& sBC);


} //namespace olb


#endif
