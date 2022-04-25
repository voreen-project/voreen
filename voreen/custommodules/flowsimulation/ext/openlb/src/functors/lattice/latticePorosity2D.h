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

#ifndef LATTICE_POROSITY_2D_H
#define LATTICE_POROSITY_2D_H

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

/// functor to get pointwise, lattice-dependent porosity values in [0,1]
/// in combination with (Extended)PorousBGKdynamics: 0->solid, 1->fluid
template <typename T, typename DESCRIPTOR>
class SuperLatticePorosity2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry<T,2>& _superGeometry;
  const int _material;
public:
  SuperLatticePorosity2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                         SuperGeometry<T,2>& superGeometry, const int material,
                         const UnitConverter<T,DESCRIPTOR>& converter);
};

/**
 *  BlockLatticePorosity2D returns pointwise, lattice-dependent porosity values in [0,1]
 *  in combination with (Extended)PorousBGKdynamics: 0->solid, 1->fluid.
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticePorosity2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometry<T,2>& _blockGeometry;
  int _material;
public:
  BlockLatticePorosity2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                         BlockGeometry<T,2>& blockGeometry, int material,
                         const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
