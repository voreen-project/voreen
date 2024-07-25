/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender, Lennart Neukamm
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

#ifndef LATTICE_PHYS_FIELD_2D_H
#define LATTICE_PHYS_FIELD_2D_H

#include "superBaseF2D.h"
#include "core/superLattice2D.h"
#include "utilities/functorDsl2D.h"

namespace olb {

template <typename T, typename DESCRIPTOR, typename FIELD>
struct SuperLatticePhysField2D final : public SuperIdentity2D<T> {
  SuperLatticePhysField2D(SuperLattice<T,DESCRIPTOR>& sLattice,
                          T convFactorToPhysUnits):
    SuperIdentity2D<T>([&](){
      return functor_dsl::field<T,DESCRIPTOR,FIELD>(sLattice) * convFactorToPhysUnits;
    }())
  {
    this->getName() = "physField";
  }
};

}
#endif
