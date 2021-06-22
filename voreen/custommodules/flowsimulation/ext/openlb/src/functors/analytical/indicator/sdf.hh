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

#ifndef SDF_HH
#define SDF_HH

#include "sdf.h"

#include <algorithm>

namespace olb {

template <typename T>
IndicatorSDF2D<T>::IndicatorSDF2D(std::function<T(Vector<T,2>)> f):
  _f(f) { }

template <typename T>
bool IndicatorSDF2D<T>::operator()(bool output[], const T input[]) {
  output[0] = _f(input) < 0.0;
  return true;
}

template <typename T>
IndicatorSDF3D<T>::IndicatorSDF3D(std::function<T(Vector<T,3>)> f):
  _f(f) { }

template <typename T>
bool IndicatorSDF3D<T>::operator()(bool output[], const T input[]) {
  output[0] = _f(input) < 0.0;
  return true;
}

namespace sdf {

template <typename T>
T mix(T a, T b, T h) {
  return b*(1.0-h) + a*h;
}

template <typename T>
T clamp(T x, T a, T b) {
  if (x < a) {
    return a;
  } else if (x > b) {
    return b;
  } else {
    return x;
  }
}

template <typename T, unsigned D>
T sphere(Vector<T,D> p, T r)
{
  return p.norm() - r;
}

template <typename T>
T box(Vector<T,3> p, Vector<T,3> b)
{
  auto q = p.abs() - b;
  return maxv(q, Vector<T,3>(0.0)).norm() + min(max(q[0], max(q[1],q[2])), 0.0);
}

template <typename T>
T torus(Vector<T,3> p, Vector<T,2> t)
{
  Vector<T,2> b{p[0], p[2]};
  Vector<T,2> q{b.norm()-t[0], p[1]};
  return q.norm() - t[1];
}

template <typename T>
T solid_angle(Vector<T,3> p, Vector<T,2> c, T r)
{
  Vector<T,2> q{Vector<T,2>{p[0], p[2]}.norm(), p[1]};
  T l = q.norm() - r;
  T m = (q - c*clamp(q * c, 0.0, r)).norm();
  return max(l, m*sign(c[1]*q[0]-c[0]*q[1]));
}

template <typename T, unsigned D>
Vector<T,D> translate(Vector<T,D> p, Vector<T,D> origin)
{
  return p - origin;
}

template <typename T>
Vector<T,3> flip(Vector<T,3> p)
{
  return {p[1], p[0], p[2]};
}

template <typename T>
T substract(T a, T b) {
  return max(-a, b);
}

template <typename T>
T unify(T a, T b)
{
  return min(a, b);
}

template <typename T>
T intersection(T a, T b)
{
  return max(a, b);
}

template <typename T>
T smooth_union(T a, T b, T k)
{
  T h = clamp(0.5 + 0.5*(b-a)/k, 0.0, 1.0);
  return mix(a, b, h) - k*h*(1.0-h);
}

template <typename T>
T smooth_subtraction(T a, T b, T k)
{
  T h = clamp(0.5 - 0.5*(b+a)/k, 0.0, 1.0);
  return mix(b, -a, h) + k*h*(1.0-h);
}

template <typename T>
T smooth_intersection(T a, T b, T k) {
  T h = clamp(0.5 - 0.5*(b-a)/k, 0.0, 1.0);
  return mix(b, a, h) + k*h*(1.0-h);
}

template <typename T>
T rounding(T a, T r)
{
  return a - r;
}

}

}

#endif
