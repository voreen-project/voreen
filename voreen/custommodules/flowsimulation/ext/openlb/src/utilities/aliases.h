/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender
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

#ifndef ALIASES_H
#define ALIASES_H

#include <type_traits>

namespace olb {

template <typename T, typename DESCRIPTOR> class BlockStaticPopulationD2D;
template <typename T, typename DESCRIPTOR> class BlockStaticPopulationD3D;

template <typename T, typename DESCRIPTOR>
using BlockStaticPopulationD = std::conditional_t<
  DESCRIPTOR::d == 2,
  BlockStaticPopulationD2D<T,DESCRIPTOR>,
  BlockStaticPopulationD3D<T,DESCRIPTOR>
>;

template <typename T, typename DESCRIPTOR> class BlockStaticFieldsD2D;
template <typename T, typename DESCRIPTOR> class BlockStaticFieldsD3D;

template <typename T, typename DESCRIPTOR>
using BlockStaticFieldsD = std::conditional_t<
  DESCRIPTOR::d == 2,
  BlockStaticFieldsD2D<T,DESCRIPTOR>,
  BlockStaticFieldsD3D<T,DESCRIPTOR>
>;

template <typename T, typename DESCRIPTOR> class BlockDynamicFieldsD2D;
template <typename T, typename DESCRIPTOR> class BlockDynamicFieldsD3D;

template <typename T, typename DESCRIPTOR>
using BlockDynamicFieldsD = std::conditional_t<
  DESCRIPTOR::d == 2,
  BlockDynamicFieldsD2D<T,DESCRIPTOR>,
  BlockDynamicFieldsD3D<T,DESCRIPTOR>
>;

template <typename T, typename DESCRIPTOR> class BlockLattice2D;
template <typename T, typename DESCRIPTOR> class BlockLattice3D;

template <typename T, typename DESCRIPTOR>
using BlockLattice = std::conditional_t<
  DESCRIPTOR::d == 2,
  BlockLattice2D<T,DESCRIPTOR>,
  BlockLattice3D<T,DESCRIPTOR>
>;

template <typename T, typename DESCRIPTOR> class SuperLattice2D;
template <typename T, typename DESCRIPTOR> class SuperLattice3D;

template <typename T, typename DESCRIPTOR>
using SuperLattice = std::conditional_t<
  DESCRIPTOR::d == 2,
  SuperLattice2D<T,DESCRIPTOR>,
  SuperLattice3D<T,DESCRIPTOR>
>;


}

#endif
