/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
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

#pragma once

#include "dynamics/descriptorField.h"

#include <cstdint>
#include <array>

namespace olb {

namespace FreeSurface {

struct Stage0 { };
struct Stage1 { };
struct Stage2 { };
struct Stage3 { };
struct Stage4 { };

enum class Type : std::uint8_t {
  Gas = 0,
  Interface = 1,
  Fluid = 2,
  Solid = 4
};

enum class Flags : std::uint8_t {
  ToGas = 1,
  ToFluid = 2,
  NewInterface = 4
};

struct CELL_TYPE : public descriptors::TYPED_FIELD_BASE<Type,1> { };
struct CELL_FLAGS : public descriptors::TYPED_FIELD_BASE<Flags,1> { };
struct MASS : public descriptors::FIELD_BASE<1> { };
struct EPSILON : public descriptors::FIELD_BASE<1> { };
struct PREVIOUS_VELOCITY : public descriptors::FIELD_BASE<0,1,0>{};

struct TEMP_MASS_EXCHANGE : public descriptors::FIELD_BASE<0,0,1>{};

template<typename T, size_t S>
static std::array<T,S> solvePivotedLU(std::array<std::array<T,S>,S>& matrix, const std::array<T,S>& b, size_t N = S) any_platform;

template<typename T, typename DESCRIPTOR>
void initialize(SuperLattice<T,DESCRIPTOR>& lattice);

}

any_platform FreeSurface::Flags operator&(FreeSurface::Flags lhs, FreeSurface::Flags rhs) {
  return static_cast<FreeSurface::Flags>(static_cast<std::uint8_t>(lhs) & static_cast<std::uint8_t>(rhs));
}
any_platform FreeSurface::Flags operator|(FreeSurface::Flags lhs, FreeSurface::Flags rhs) {
  return static_cast<FreeSurface::Flags>(static_cast<std::uint8_t>(lhs) | static_cast<std::uint8_t>(rhs));
}


}

#include "freeSurfaceHelpers.hh"
