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
 * A helper for initialising 3D boundaries -- header file.
 */

#ifndef SUPER_BOUNDARY_CONDITION_3D_H
#define SUPER_BOUNDARY_CONDITION_3D_H

#include <vector>
#include "io/ostreamManager.h"
#include "extendedFiniteDifferenceBoundary3D.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T, template<typename U> class Lattice> class OnLatticeAdvectionDiffusionBoundaryCondition3D;
template<typename T, template<typename U> class Lattice> class OnLatticeBoundaryCondition3D;
template<typename T, template<typename U> class Lattice> class SuperLattice3D;
template<typename T> class SuperGeometry3D;

/// A helper for initialising 3D boundaries for super lattices.
/** Here we have methods that initializes the local postprocessors and the
 * communicator (_commBC in SuperLattice) for boundary conditions
 * for a given global point or global range.
 *
 * This class is not intended to be derived from.
 */
template<typename T, template<typename U> class Lattice>
class sOnLatticeBoundaryCondition3D {
public:
  sOnLatticeBoundaryCondition3D(SuperLattice3D<T, Lattice>& sLattice);
  sOnLatticeBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, Lattice> const& rhs);
  sOnLatticeBoundaryCondition3D operator=(sOnLatticeBoundaryCondition3D<T,Lattice> rhs);
  ~sOnLatticeBoundaryCondition3D();

  void addVelocityBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega);
  void addSlipBoundary(SuperGeometry3D<T>& superGeometry, int material);
  void addPressureBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega);
  void addConvectionBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega, T* uAv=NULL);
  void addTemperatureBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega);
  void addConvectionBoundary(SuperGeometry3D<T>& superGeometry, int material);
  void addExtFieldBoundary(SuperGeometry3D<T>& superGeometry, int material, int offset);
  void addZeroDistributionBoundary(SuperGeometry3D<T>& superGeometry, int material);

  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(SuperGeometry3D<T>& superGeometry, int material);

  SuperLattice3D<T, Lattice>& getSuperLattice();
  std::vector<OnLatticeBoundaryCondition3D<T, Lattice>*>& getBlockBCs();
  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Lattice>*>& getADblockBCs();
  int getOverlap();
  void setOverlap(int overlap);

  void outputOn();
  void outputOff();

private:
  mutable OstreamManager clout;
  SuperLattice3D<T, Lattice>& _sLattice;
  std::vector<OnLatticeBoundaryCondition3D<T, Lattice>*> _blockBCs;
  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Lattice>*> _ADblockBCs;
  int _overlap;
  bool _output;
};


////////////////// Factory functions //////////////////////////////////

template<typename T, template<typename U> class Lattice>
void createLocalBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, Lattice>& sBC);

template<typename T, template<typename U> class Lattice, typename MixinDynamics=BGKdynamics<T,Lattice> >
void createInterpBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, Lattice>& sBC);

template<typename T, template<typename U> class Lattice, typename MixinDynamics=BGKdynamics<T,Lattice> >
void createExtFdBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, Lattice>& sBC);


} // namespace olb

#endif
