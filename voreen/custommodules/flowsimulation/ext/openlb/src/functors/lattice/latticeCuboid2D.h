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

#ifndef LATTICE_CUBOID_2D_H
#define LATTICE_CUBOID_2D_H

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

/// functor to get pointwise the cuboid no. + 1 on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeCuboid2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeCuboid2D(SuperLattice<T,DESCRIPTOR>& sLattice);
};

/// BlockLatticeCuboid2D returns pointwise the cuboid no. + 1 on local lattice.
template <typename T, typename DESCRIPTOR>
class BlockLatticeCuboid2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
public:
  BlockLatticeCuboid2D(BlockLattice<T,DESCRIPTOR>& blockLattice, const int iC);
  bool operator() (T output[], const int input[]) override;
private:
  // holds cuboid nmb of current block
  const int _iC;
};

}
#endif
