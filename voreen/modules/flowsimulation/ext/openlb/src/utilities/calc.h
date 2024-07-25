/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
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
 * Some arithmetic helper functions.
 */

#ifndef CALC_H
#define CALC_H

#include "omath.h"

namespace olb {

namespace util {

/// Variant of fmod (floating point modulo) that always returns positive values
template<typename T, typename S>
inline auto fmod_pos(T a, S b) {
  const auto res = util::fmod(a, b);
  if (res < 0) {
    return (res + b);
  } else {
    return res;
  }
}


/** Accurate summation of floating point numbers with the Kahan algorithm
 * Reduces round-off effects which arise if the total sum is significantly
 * larger than the single summands. Works very well for adding many numbers of
 * the same sign. Does not really help if sign changes, e.g. if the summands
 * are randomly distributed around 0.
 * Cf. https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
template<typename T>
class KahanSummator
{
private:
  T sum {0};
  T c {0};   // tracks the cancellation errors

public:
  void add(T summand) {
    const T y {summand - c};
    const T t {sum + y};
    c = (t - sum) - y;
    sum = t;
  }

  T getSum() const {
    return sum;
  }

  void initialize(T value = 0.) {
    sum = value;
    c = T(0);
  }
};

} // namespace util

} // namespace olb

#endif
