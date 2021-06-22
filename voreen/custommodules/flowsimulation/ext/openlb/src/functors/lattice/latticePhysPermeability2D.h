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

#ifndef LATTICE_PHYS_PERMEABILITY_2D_H
#define LATTICE_PHYS_PERMEABILITY_2D_H

#include <vector>

#include "superBaseF2D.h"
#include "core/superLattice2D.h"
#include "indicator/superIndicatorBaseF2D.h"
#include "utilities/functorPtr.h"
#include "blockBaseF2D.h"
#include "geometry/blockGeometry2D.h"
#include "core/blockLattice2D.h"
#include "core/blockLatticeStructure2D.h"
#include "indicator/blockIndicatorF2D.h"
#include "dynamics/porousBGKdynamics.h"

namespace olb {

/// functor to get pointwise mesh-independent permeability values in (0,inf)  in combination with (Extended)PorousBGKdynamics
/// note: result is cropped to 999999
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysPermeability2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysPermeability2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                 SuperGeometry2D<T>& superGeometry,
                                 const int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/**
 *  BlockLatticePhysPermeability2D returns pointwise mesh-independent permeability
 *  values in (0,inf) in combination with (Extended)PorousBGKdynamics
 *  note: result is cropped to 999999.
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysPermeability2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysPermeability2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                 BlockGeometryStructure2D<T>& blockGeometry,
                                 int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
