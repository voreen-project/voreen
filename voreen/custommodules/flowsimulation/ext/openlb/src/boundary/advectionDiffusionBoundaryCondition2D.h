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

/** \file A helper for initialising 2D boundaries -- header file.  */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_H
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_H

#include "momentaOnBoundaries2D.h"
#include "dynamics/dynamics.h"
#include "dynamics/advectionDiffusionDynamics.h"
#include "superBoundaryCondition2D.h"

#include <vector>
#include <list>

namespace olb {

template<typename T, template<typename U> class Lattice>
class sOnLatticeBoundaryCondition2D;

template<typename T, template<typename U> class Lattice>
class OnLatticeAdvectionDiffusionBoundaryCondition2D {
public:
  virtual ~OnLatticeAdvectionDiffusionBoundaryCondition2D() { }

  virtual void addTemperatureBoundary0N(int x0, int x1, int y0, int y1,T omega) =0;
  virtual void addTemperatureBoundary0P(int x0, int x1, int y0, int y1,T omega) =0;
  virtual void addTemperatureBoundary1N(int x0, int x1, int y0, int y1,T omega) =0;
  virtual void addTemperatureBoundary1P(int x0, int x1, int y0, int y1,T omega) =0;

  virtual void addTemperatureCornerNN(int x, int y, T omega) =0;
  virtual void addTemperatureCornerNP(int x, int y, T omega) =0;
  virtual void addTemperatureCornerPN(int x, int y, T omega) =0;
  virtual void addTemperatureCornerPP(int x, int y, T omega) =0;

  BlockLatticeStructure2D<T,Lattice>& getBlock();
  BlockLatticeStructure2D<T,Lattice> const& getBlock() const;

  /// adds a temperature boundary for one material or a range (x0-x1, y0-y1, z0-z1)
  virtual void addTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega) =0;
};

//////  Factory function for Regularized Thermal BC

/// blockLattice creator
template<typename T, template<typename U> class Lattice, typename MixinDynamics=AdvectionDiffusionRLBdynamics<T,Lattice> >
OnLatticeAdvectionDiffusionBoundaryCondition2D<T,Lattice>*
createAdvectionDiffusionBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block);

/// superLattice creator, calls createAdvectionDiffusionBoundaryCondidtion3D from above.
template<typename T, template<typename U> class Lattice, typename MixinDynamics=AdvectionDiffusionRLBdynamics<T,Lattice> >
void createAdvectionDiffusionBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T, Lattice>& sBC);


} //namespace olb


#endif
