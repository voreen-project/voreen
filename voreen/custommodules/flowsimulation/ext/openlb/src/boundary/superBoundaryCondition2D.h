/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause
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

/** \file
 * A helper for initialising 2D boundaries -- header file.
 */


#ifndef SUPER_BOUNDARY_CONDITION_2D_H
#define SUPER_BOUNDARY_CONDITION_2D_H

#include <vector>
#include "boundaryCondition2D.h"
#include "advectionDiffusionBoundaryCondition2D.h"     // -> for AdvectionDiffusion
#include "geometry/blockGeometryStatistics2D.h"
#include "core/superLattice2D.h"
#include "io/ostreamManager.h"
#include "geometry/superGeometry2D.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

/// A helper for initialising 2D boundaries for super lattices.
/** Here we have methods that initializes for a given global
 * point or global range the local postprocessors and the
 * communicator (_commBC in SuperLattice) for boundary conditions.
 *
 * This class is not intended to be derived from.
 */



template<typename T, template<typename U> class Lattice>
class sOnLatticeBoundaryCondition2D {

public:
  /// Constructor
  sOnLatticeBoundaryCondition2D(SuperLattice2D<T,Lattice>& sLattice);
  /// Copy construction
  sOnLatticeBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T,Lattice> const& rhs);
  /// Copy assignment
  sOnLatticeBoundaryCondition2D operator=(sOnLatticeBoundaryCondition2D<T,Lattice> rhs);
  /// Destructor
  ~sOnLatticeBoundaryCondition2D();

  void addVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega);
  void addSlipBoundary(SuperGeometry2D<T>& superGeometry, int material);
  void addPressureBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega);
  void addConvectionBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega, T* uAv=NULL);
  void addTemperatureBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega);

  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(SuperGeometry2D<T>& superGeometry, int material);

  SuperLattice2D<T,Lattice>& getSuperLattice()
  {
    return _sLattice;
  };
  std::vector<OnLatticeBoundaryCondition2D<T,Lattice>* >& getBlockBCs()
  {
    return _blockBCs;
  };
  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Lattice>*>& getADblockBCs()
  {
    return _ADblockBCs;
  };
  int getOverlap()
  {
    return _overlap;
  };
  void setOverlap(int overlap)
  {
    _overlap = overlap;
  };

  void outputOn();
  void outputOff();

private:
  mutable OstreamManager clout;
  SuperLattice2D<T,Lattice>& _sLattice;
  std::vector<OnLatticeBoundaryCondition2D<T,Lattice>* > _blockBCs;
  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Lattice>*> _ADblockBCs;       // -> for AdvectionDiffusion
  int _overlap;
  bool _output;
};


template<typename T, template<typename U> class Lattice, typename MixinDynamics=RLBdynamics<T,Lattice> >
void createLocalBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T,Lattice>& sBC);

template<typename T, template<typename U> class Lattice, typename MixinDynamics=BGKdynamics<T,Lattice> >
void createInterpBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T,Lattice>& sBC);

template<typename T, template<typename U> class Lattice, typename MixinDynamics=BGKdynamics<T,Lattice> >
void createExtFdBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T, Lattice>& sBC);

}  // namespace olb

#endif
