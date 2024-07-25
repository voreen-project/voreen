/*  This file is part of the porous model described in
 *  Guo and Zhao (2002) and implemented for OpenLB
 *
 *  Copyright (C) 2016-2017 Davide Dapelo, Mathias J. Krause
 *  E-mail contact: dapelod@bham.ac.uk
 *
 *  OpenLB e-mail contact: info@openlb.net
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
 * Class to provide access to the external field for
 * Lagrangian particle density, converted to Eulerian -- header file
 */

#ifndef EUL2LAGR_DENSITY_H
#define EUL2LAGR_DENSITY_H

#include "functors/lattice/functors3D.h"

namespace olb {

/// functor returns pointwise external field for Lagrangian particle density, converted to Eulerian
template <typename T, typename DESCRIPTOR>
class BlockLatticeEul2LagrDensity3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeEul2LagrDensity3D(BlockLattice<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise external field for Lagrangian particle density, converted to Eulerian
template <typename T, typename DESCRIPTOR>
class SuperLatticeEul2LagrDensity3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeEul2LagrDensity3D(SuperLattice<T,DESCRIPTOR>& sLattice);
};


}

#endif
