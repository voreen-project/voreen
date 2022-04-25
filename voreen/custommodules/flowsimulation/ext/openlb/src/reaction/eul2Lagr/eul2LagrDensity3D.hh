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
 * Lagrangian particle density, converted to Eulerian -- generic file
 */

#ifndef EUL2LAGR_DENSITY_HH
#define EUL2LAGR_DENSITY_HH

namespace olb {

template<typename T, typename DESCRIPTOR>
BlockLatticeEul2LagrDensity3D<T, DESCRIPTOR>::BlockLatticeEul2LagrDensity3D(BlockLattice<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "eul2LagrDensity";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeEul2LagrDensity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  auto rhoP = this->_blockLattice.get(input[0], input[1], input[2]).template getFieldPointer<descriptors::EUL2LAGR>();
  output[0] = rhoP[0];
  return true;
}

template<typename T, typename DESCRIPTOR>
SuperLatticeEul2LagrDensity3D<T, DESCRIPTOR>::SuperLatticeEul2LagrDensity3D(
  SuperLattice<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "eul2LagrDensity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeEul2LagrDensity3D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC)));
  }
}


}

#endif
