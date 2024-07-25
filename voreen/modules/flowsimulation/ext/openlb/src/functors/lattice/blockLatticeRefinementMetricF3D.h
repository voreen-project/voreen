/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef BLOCK_LATTICE_REFINEMENT_METRIC_F_3D_H
#define BLOCK_LATTICE_REFINEMENT_METRIC_F_3D_H

#include "blockBaseF3D.h"

namespace olb {


template <typename T, typename DESCRIPTOR>
class BlockLatticeKnudsen3D : public BlockLatticeF3D<T, DESCRIPTOR> {
public:
  BlockLatticeKnudsen3D(BlockLattice<T, DESCRIPTOR>& blockLattice);

  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeRefinementMetricKnudsen3D final : public BlockLatticeKnudsen3D<T, DESCRIPTOR> {
private:
  const T _knudsen;
public:
  BlockLatticeRefinementMetricKnudsen3D(
    BlockLattice<T, DESCRIPTOR>& blockLattice,
    const UnitConverter<T, DESCRIPTOR>&     converter);

  bool operator() (T output[]);
  bool operator() (T output[], const int input[]) override;
};


template <typename T, typename DESCRIPTOR>
class BlockLatticeHighOrderKnudsen3D : public BlockLatticeF3D<T, DESCRIPTOR> {
public:
  BlockLatticeHighOrderKnudsen3D(BlockLattice<T, DESCRIPTOR>& blockLattice);

  bool operator() (T output[], const int input[]) override;
};


}

#endif
