/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_PHYS_BOUNDARY_FORCE_2D_H
#define LATTICE_PHYS_BOUNDARY_FORCE_2D_H

#include <vector>

#include "superBaseF2D.h"
#include "core/superLattice2D.h"
#include "indicator/superIndicatorBaseF2D.h"
#include "utilities/functorPtr.h"
#include "blockBaseF2D.h"
#include "geometry/blockGeometry.h"
#include "indicator/blockIndicatorF2D.h"
#include "dynamics/porousBGKdynamics.h"

namespace olb {

/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysBoundaryForce2D : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperIndicatorF2D<T>> _indicatorF;
public:
  SuperLatticePhysBoundaryForce2D(SuperLattice<T,DESCRIPTOR>&      sLattice,
                                  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticePhysBoundaryForce2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                                  SuperGeometry<T,2>& superGeometry, const int material,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
};

/// BlockLatticePhysBoundaryForce2D returns pointwise phys force acting on a boundary
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysBoundaryForce2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockIndicatorF2D<T>&        _indicatorF;
  BlockGeometry<T,2>& _blockGeometry;
public:
  BlockLatticePhysBoundaryForce2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                  BlockIndicatorF2D<T>& indicatorF,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
