/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#ifndef BASETYPE_H
#define BASETYPE_H

namespace olb {

namespace util {

template <typename T>
struct BaseTypeHelper {
  using type = typename T::base_t;
};

template <>
struct BaseTypeHelper<long double> {
  using type = long double;
};

template <>
struct BaseTypeHelper<double> {
  using type = double;
};

template <>
struct BaseTypeHelper<float> {
  using type = float;
};

template <>
struct BaseTypeHelper<int> {
  using type = int;
};

} // namespace util

template <typename T>
using BaseType = typename util::BaseTypeHelper<T>::type;

} // namespace olb

#endif