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

#ifndef LATTICE_EXTERNAL_3D_HH
#define LATTICE_EXTERNAL_3D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

namespace olb {

template<typename T,typename DESCRIPTOR, typename FIELD>
SuperLatticeExternal3D<T,DESCRIPTOR,FIELD>::SuperLatticeExternal3D(
  SuperLattice<T,DESCRIPTOR>& sLattice, size_t& iT)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "external";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeExternal3D<T,DESCRIPTOR,FIELD>(this->_sLattice.getBlock(iC), iT) );
  }
}

template <typename T, typename DESCRIPTOR, typename FIELD>
BlockLatticeExternal3D<T,DESCRIPTOR,FIELD>::BlockLatticeExternal3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, size_t& iT)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice,1),
    _iT(iT)
{
  this->getName() = "external";
}


template <typename T, typename DESCRIPTOR, typename FIELD>
bool BlockLatticeExternal3D<T,DESCRIPTOR,FIELD>::operator() (T output[], const int input[])
{
  T* value0 = fd::accessNew<T,DESCRIPTOR,FIELD>( this->_blockLattice.get( input[0], input[1], input[2] ), this->_iT );
  //T* value0 = fd::accessOld<T,DESCRIPTOR,FIELD>( this->_blockLattice.get( input[0], input[1], input[2] ), this->_iT );
  output[0] = *value0;
  return true;
}

}
#endif
