/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

/** \file
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef LATTICE_DESCRIPTORS_H
#define LATTICE_DESCRIPTORS_H

#include "descriptorBase.h"
#include "descriptorField.h"
#include "descriptorFunction.h"

namespace olb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors {


/// D2Q9 lattice
template <typename... FIELDS>
struct D2Q9 : public LATTICE_DESCRIPTOR<2,9,POPULATION,FIELDS...> {
  D2Q9() = delete;
};

namespace data {

template <>
platform_constant_definition int vicinity<2,9> = 1;

template <>
platform_constant_definition int c<2,9>[9][2] = {
  { 0, 0},
  {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
  { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
};

template <>
platform_constant_definition int opposite<2,9>[9] = {
  0, 5, 6, 7, 8, 1, 2, 3, 4
};

template <>
platform_constant_definition Fraction t<2,9>[9] = {
  {4, 9}, {1, 36}, {1, 9}, {1, 36}, {1, 9},
  {1, 36}, {1, 9}, {1, 36}, {1, 9}
};

template <>
platform_constant_definition Fraction cs2<2,9> = {1, 3};

}


/// D2Q5 lattice
template <typename... FIELDS>
struct D2Q5 : public LATTICE_DESCRIPTOR<2,5,POPULATION,FIELDS...> {
  D2Q5() = delete;
};

namespace data {

template <>
platform_constant_definition int vicinity<2,5> = 1;

template <>
platform_constant_definition int c<2,5>[5][2] = {
  { 0, 0},
  {-1, 0}, {0, -1}, {1,0}, { 0,1}
};

template <>
platform_constant_definition int opposite<2,5>[5] = {
  0, 3, 4, 1, 2
};

template <>
platform_constant_definition Fraction t<2,5>[5] = {
  {1, 3},
  {1, 6}, {1, 6},
  {1, 6}, {1, 6}
};

template <>
platform_constant_definition Fraction cs2<2,5> = {1, 3};

}


/// D3Q19 lattice
template <typename... FIELDS>
struct D3Q19 : public LATTICE_DESCRIPTOR<3,19,POPULATION,FIELDS...> {
  D3Q19() = delete;
};

namespace data {

template <>
platform_constant_definition int vicinity<3,19> = 1;

template <>
platform_constant_definition int c<3,19>[19][3] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
};

template <>
platform_constant_definition int opposite<3,19>[19] = {
  0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9
};

template <>
platform_constant_definition Fraction t<3,19>[19] = {
  {1, 3},

  {1, 18}, {1, 18}, {1, 18},
  {1, 36}, {1, 36}, {1, 36},
  {1, 36}, {1, 36}, {1, 36},

  {1, 18}, {1, 18}, {1, 18},
  {1, 36}, {1, 36}, {1, 36},
  {1, 36}, {1, 36}, {1, 36}
};

template <>
platform_constant_definition Fraction cs2<3,19> = {1, 3};

}


/// D3Q7 lattice
template <typename... FIELDS>
struct D3Q7 : public LATTICE_DESCRIPTOR<3,7,POPULATION,FIELDS...> {
  D3Q7() = delete;
};

namespace data {

template <>
platform_constant_definition int vicinity<3,7> = 1;

template <>
platform_constant_definition int c<3,7>[7][3] = {
  { 0, 0, 0},

  {-1, 0, 0}, {0,-1, 0},
  { 0, 0,-1}, {1, 0, 0},
  { 0, 1, 0}, {0, 0, 1},
};

template <>
platform_constant_definition int opposite<3,7>[7] = {
  0, 4, 5, 6, 1, 2, 3
};

template <>
platform_constant_definition Fraction cs2<3,7> = {1, 4};

template <>
platform_constant_definition Fraction t<3,7>[7] = {
  {1, 4},

  {1, 8}, {1, 8}, {1, 8},
  {1, 8}, {1, 8}, {1, 8}
};

}


/// D3Q13 lattice
template <typename... FIELDS>
struct D3Q13 : public LATTICE_DESCRIPTOR<3,13,POPULATION,FIELDS...> {
  D3Q13() = delete;
};

namespace data {

template <>
platform_constant_definition int vicinity<3,13> = 1;

template <>
platform_constant_definition int c<3,13>[13][3] = {
  { 0, 0, 0},

  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
};

template <>
platform_constant_definition int opposite<3,13>[13] = {
  0, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6
};

template <>
platform_constant_definition Fraction cs2<3,13> = {1, 3};

template <>
platform_constant_definition Fraction t<3,13>[13] = {
  {1, 2},

  {1, 24}, {1, 24}, {1, 24},
  {1, 24}, {1, 24}, {1, 24},

  {1, 24}, {1, 24}, {1, 24},
  {1, 24}, {1, 24}, {1, 24}
};

template <>
platform_constant_definition Fraction lambda_e<3,13> = {3, 2};

template <>
platform_constant_definition Fraction lambda_h<3,13> = {9, 5};

}


/// D3Q15 lattice
template <typename... FIELDS>
struct D3Q15 : public LATTICE_DESCRIPTOR<3,15,POPULATION,FIELDS...> {
  D3Q15() = delete;
};

namespace data {

template <>
platform_constant_definition int vicinity<3,15> = 1;

template <>
platform_constant_definition int c<3,15>[15][3] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1}
};

template <>
platform_constant_definition int opposite<3,15>[15] = {
  0, 8, 9, 10, 11, 12, 13, 14, 1, 2, 3, 4, 5, 6, 7
};

template <>
platform_constant_definition Fraction cs2<3,15> = {1, 3};

template <>
platform_constant_definition Fraction t<3,15>[15] = {
  {2, 9},

  {1, 9}, {1, 9}, {1, 9},
  {1, 72}, {1, 72}, {1, 72}, {1, 72},

  {1, 9}, {1, 9}, {1, 9},
  {1, 72}, {1, 72}, {1, 72}, {1, 72}
};

}


/// D3Q27 lattice
template <typename... FIELDS>
struct D3Q27 : public LATTICE_DESCRIPTOR<3,27,POPULATION,FIELDS...> {
  D3Q27() = delete;
};

namespace data {

template <>
platform_constant_definition int vicinity<3,27> = 1;

template <>
platform_constant_definition int c<3,27>[27][3] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},
  {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1},
  { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1}
};

template <>
platform_constant_definition int opposite<3,27>[27] = {
  0, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
};

template <>
platform_constant_definition Fraction cs2<3,27> = {1, 3};

template <>
platform_constant_definition Fraction t<3,27>[27] = {
  {8, 27},

  {2, 27},  {2, 27},  {2, 27},
  {1, 54},  {1, 54},  {1, 54},
  {1, 54},  {1, 54},  {1, 54},
  {1, 216}, {1, 216}, {1, 216}, {1, 216},

  {2, 27},  {2, 27},  {2, 27},
  {1, 54},  {1, 54},  {1, 54},
  {1, 54},  {1, 54},  {1, 54},
  {1, 216}, {1, 216}, {1, 216}, {1, 216}
};

}


}  // namespace descriptors

}  // namespace olb

#endif
