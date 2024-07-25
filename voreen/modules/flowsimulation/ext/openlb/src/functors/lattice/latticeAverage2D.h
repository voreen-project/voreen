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

#ifndef LATTICE_AVERAGE_2D_H
#define LATTICE_AVERAGE_2D_H

#include <vector>

#include "blockBaseF2D.h"
#include "geometry/blockGeometry.h"
#include "indicator/blockIndicatorF2D.h"
#include "dynamics/porousBGKdynamics.h"

namespace olb {

/**
 *  BlockLatticeAverage2D returns pointwise local average of a passed functor with
 *  a given material and radius on local lattice.
 *  the output data must be of the same size and dimension like f.
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticeAverage2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
  BlockGeometry<T,2>& _blockGeometry;
  int _material;
  T _radius;
public:
  BlockLatticeAverage2D(BlockLatticeF2D<T,DESCRIPTOR>& f,
                        BlockGeometry<T,2>& blockGeometry, int material, T radius);
  bool operator() (T output[], const int input[]) override;
};

}
#endif


