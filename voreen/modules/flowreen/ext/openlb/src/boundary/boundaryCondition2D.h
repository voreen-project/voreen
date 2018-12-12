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

#ifndef BOUNDARY_CONDITION_2D_H
#define BOUNDARY_CONDITION_2D_H

#include "core/blockLatticeStructure2D.h"
#include "momentaOnBoundaries2D.h"
#include "boundaryPostProcessors2D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometry2D.h"
namespace olb {

template<typename T, template<typename U> class Lattice>
class OnLatticeBoundaryCondition2D {
public:
  virtual ~OnLatticeBoundaryCondition2D() { }

  virtual void addVelocityBoundary0N(int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addVelocityBoundary0P(int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addVelocityBoundary1N(int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addVelocityBoundary1P(int x0, int x1, int y0, int y1, T omega) =0;

  virtual void addPressureBoundary0N(int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addPressureBoundary0P(int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addPressureBoundary1N(int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addPressureBoundary1P(int x0, int x1, int y0, int y1, T omega) =0;

  virtual void addConvectionBoundary0N(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary0P(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary1N(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary1P(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) =0;

  virtual void addExternalVelocityCornerNN(int x, int y, T omega) =0;
  virtual void addExternalVelocityCornerNP(int x, int y, T omega) =0;
  virtual void addExternalVelocityCornerPN(int x, int y, T omega) =0;
  virtual void addExternalVelocityCornerPP(int x, int y, T omega) =0;

  virtual void addInternalVelocityCornerNN(int x, int y, T omega) =0;
  virtual void addInternalVelocityCornerNP(int x, int y, T omega) =0;
  virtual void addInternalVelocityCornerPN(int x, int y, T omega) =0;
  virtual void addInternalVelocityCornerPP(int x, int y, T omega) =0;

  /// adds a pressure or velocity boundary for one material and a range (x0-x1, y0-y1, z0-z1) or the whole geometry
  virtual void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega) =0;
  virtual void addSlipBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1) =0;
  virtual void addSlipBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material) =0;
  virtual void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega) =0;
  virtual void addConvectionBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, T* uAv=NULL) =0;

  virtual BlockLatticeStructure2D<T,Lattice>& getBlock() =0;
  virtual BlockLatticeStructure2D<T,Lattice> const& getBlock() const =0;

  virtual void outputOn() =0;
  virtual void outputOff() =0;
};

////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics=RLBdynamics<T,Lattice> >
OnLatticeBoundaryCondition2D<T,Lattice>*
createLocalBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block);

template<typename T, template<typename U> class Lattice, typename MixinDynamics=BGKdynamics<T,Lattice> >
OnLatticeBoundaryCondition2D<T,Lattice>*
createInterpBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block);


}

#endif
