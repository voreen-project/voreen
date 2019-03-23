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

/** \file
 * A helper for initialising 2D boundaries -- header file.
 */


#ifndef SUPER_OFF_BOUNDARY_CONDITION_2D_H
#define SUPER_OFF_BOUNDARY_CONDITION_2D_H

#include <vector>
#include <list>

#include "offBoundaryCondition2D.h"
#include "geometry/superGeometry2D.h"
#include "core/superLattice2D.h"
#include "io/ostreamManager.h"
#include "functors/analytical/analyticalF.h"

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
class sOffLatticeBoundaryCondition2D {

public:
  /// Constructor
  sOffLatticeBoundaryCondition2D(SuperLattice2D<T,Lattice>& sLattice, T epsFraction_ = 0.0001);
  /// Copy construction
  sOffLatticeBoundaryCondition2D(sOffLatticeBoundaryCondition2D<T,Lattice> const& rhs);
  /// Copy assignment
  sOffLatticeBoundaryCondition2D operator=(sOffLatticeBoundaryCondition2D<T,Lattice> rhs);
  /// Destructor
  ~sOffLatticeBoundaryCondition2D();

  /// Automatically sets offDynamics with boundary links and post processors using an indicator function and material number
  /// Add offDynamics with initialisation of boundary links and the corresponding post processors
  /// Note: Uses information of the second neighbours of the cell (x,y,z)
  /// Add post processors. Ensure that offDynamics are defined!
  void addVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1));
  void addZeroVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1));
  void addPressureBoundary(SuperGeometry2D<T>& superGeometry, int material, IndicatorF2D<T>& indicator, std::list<int> bulkMaterials = std::list<int>(1,1));

  void addVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material, std::list<int> bulkMaterials = std::list<int>(1,1));
  void addZeroVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material,  std::list<int> bulkMaterials = std::list<int>(1,1));
  void addPressureBoundary(SuperGeometry2D<T>& superGeometry, int material,  std::list<int> bulkMaterials = std::list<int>(1,1));

  void defineRho(SuperGeometry2D<T>& superGeometry, int material, AnalyticalF2D<T,T>& rho, std::list<int> bulkMaterials = std::list<int>(1,1) );
  void defineU(SuperGeometry2D<T>& superGeometry, int material, AnalyticalF2D<T,T>& u, std::list<int> bulkMaterials = std::list<int>(1,1) );

  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(SuperGeometry2D<T>& superGeometry, int material);

  SuperLattice2D<T,Lattice>& getSuperLattice()
  {
    return _sLattice;
  };
  std::vector<OffLatticeBoundaryCondition2D<T,Lattice>* >& getBlockBCs()
  {
    return _blockBCs;
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
  std::vector<OffLatticeBoundaryCondition2D<T,Lattice>* > _blockBCs;
  T _epsFraction;
  int _overlap;
  bool _output;
};

template<typename T, template<typename U> class Lattice, typename MixinDynamics=BGKdynamics<T,Lattice> >
void createBouzidiBoundaryCondition2D(sOffLatticeBoundaryCondition2D<T,Lattice>& sBC);

template<typename T, template<typename U> class Lattice>
void createBounceBackBoundaryCondition2D(sOffLatticeBoundaryCondition2D<T,Lattice>& sBC);

}  // namespace olb

#endif
