/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron, Davide Dapelo
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

#ifndef LATTICE_EXTERNAL_2D_H
#define LATTICE_EXTERNAL_2D_H

#include <vector>

namespace olb {

/// functor to get pointwise density rho on local lattices
template <typename T, typename DESCRIPTOR, typename FIELD>
class SuperLatticeExternal2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeExternal2D(SuperLattice<T,DESCRIPTOR>& sLattice, size_t& iT);
};

/// BlockLatticeExternal2D returns pointwise density rho on local lattices.
template <typename T, typename DESCRIPTOR, typename FIELD>
class BlockLatticeExternal2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
public:
  BlockLatticeExternal2D(BlockLattice<T,DESCRIPTOR>& blockLattice, size_t& iT);
  bool operator() (T output[], const int input[]) override;
private:
  size_t& _iT;
};

}
#endif
