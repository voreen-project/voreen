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

#ifndef LATTICE_GEOMETRY_2D_H
#define LATTICE_GEOMETRY_2D_H

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

template<typename T, unsigned D> class SuperGeometry;

/// functor to get pointwise the material no. presenting the geometry on local lattice
template <typename T, typename DESCRIPTOR=void>
class SuperLatticeGeometry2D final : public SuperF2D<T> {
private:
  SuperGeometry<T,2>& _superGeometry;
  const int _material;
public:
  SuperLatticeGeometry2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                         SuperGeometry<T,2>& superGeometry, const int material = -1);
  SuperLatticeGeometry2D(SuperGeometry<T,2>& superGeometry, const int material = -1);
};

/// BlockLatticeGeometry2D returns pointwise the material no. presenting the geometry on local lattice.
template <typename T, typename DESCRIPTOR=void>
class BlockLatticeGeometry2D final : public BlockF2D<T> {
private:
  BlockGeometry<T,2>& _blockGeometry;
  const int _material;
public:
  BlockLatticeGeometry2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                         BlockGeometry<T,2>& blockGeometry, int material = -1);
  BlockLatticeGeometry2D(BlockGeometry<T,2>& blockGeometry, int material = -1);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
