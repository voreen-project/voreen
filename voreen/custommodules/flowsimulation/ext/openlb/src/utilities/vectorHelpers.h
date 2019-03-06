/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Lukas Baron, Mathias J. Krause
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

#ifndef VECTOR_HELPERS_H
#define VECTOR_HELPERS_H

#include<assert.h>
#include<cmath>
#include<vector>
#include<string>
#include<sstream>
#include<limits>

#include "io/ostreamManager.h"
#include "core/vector.h"

namespace olb {

template<typename T, unsigned Size > class Vector;

namespace util {


/// return true if a is close to zero
template <class T>
inline bool nearZero(const T& a)
{
  T EPSILON = std::numeric_limits<T>::epsilon();
  if (a > -EPSILON && a < EPSILON) {
    return true;
  } else {
    return false;
  }
}

template<typename T>
inline bool approxEqual(const T& a, const T& b)
{
    return nearZero<T>(a - b);
}

template <class T>
inline void copyN(T c[], const T a[], const unsigned dim)
{
  for (unsigned i=0; i<dim; i++) {
    c[i] = a[i];
  }
}

template <class T>
inline void copy3(T c[], const T a[])
{
  for (unsigned i=0; i<3; i++) {
    c[i] = a[i];
  }
}


template <typename T>
std::vector<T> fromVector3(const Vector<T,3>& vec)
{
  std::vector<T> v;
  v.push_back(vec[0]);
  v.push_back(vec[1]);
  v.push_back(vec[2]);
  return v;
}
template <typename T>
std::vector<T> fromVector2(const Vector<T,2>& vec)
{
  std::vector<T> v;
  v.push_back(vec[0]);
  v.push_back(vec[1]);
  return v;
}


/// l2 norm of a vector of arbitrary length
template <typename T>
T norm(const std::vector<T>& a)
{
  T v(0);
  for (unsigned iD=0; iD<a.size(); iD++) {
    v += a[iD]*a[iD];
  }
  v = sqrt(v);
  return v;
}

/// l2 norm to the power of 2 of a vector of arbitrary length
template <typename T>
T norm2(const std::vector<T>& a)
{
  T v = T();
  for (unsigned iD=0; iD<a.size(); iD++) {
    v += a[iD]*a[iD];
  }
  return v;
}

/// dot product, only valid in 3d
template <typename T>
T dotProduct3D(const Vector<T,3>& a, const Vector<T,3>& b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/// dot product, only valid in 2d
template <typename T>
T dotProduct2D(const Vector<T,2>& a, const Vector<T,2>& b)
{
  return a[0]*b[0] + a[1]*b[1];
}

/// returns a normalized vector, works for arbitrary lengths
template <typename T>
std::vector<T> normalize(const std::vector<T>& a)
{
  std::vector<T> out(a);
  T scale = norm(a);
  assert(scale>0);
  for (unsigned int iDim=0; iDim<a.size(); iDim++) {
    out[iDim] /= scale;
  }
  return out;
}

/*
/// algorithm by Möller–Trumbore (TODO add ref), implemented by Lucas Cruz and Mathias J. Krause
/// returns true if there is an intersection of a triangle given by (point0, point1, point1) and a ray given by its origin and direction and computes the distance
template <typename T>
bool triangleIntersectionWithNormalDirection(const std::vector<T>& point0,
    const std::vector<T>& point1, const std::vector<T>& point2,
    const std::vector<T>& origin, const std::vector<T>& normalDirection,
    T& distance)
{
  T EPSILON = std::numeric_limits<T>::epsilon();
  std::vector<T> e1, e2;
  std::vector<T> P, Q, TT;
  T det, inv_det;
  T t, u, v;
  e1 = point1 - point0;
  e2 = point2 - point0;
  P = crossProduct3D(normalDirection, e2);
  det = dotProduct3D(P, e1);
  if (det > -EPSILON && det < EPSILON) {
    return false;
  }
  inv_det = T(1) / det;
  TT = origin - point0;
  u = dotProduct3D(TT, P)*inv_det;
  if (u < T() || u > T(1)) {
    return false;
  }
  Q = crossProduct3D(TT, e1);
  v = dotProduct3D(normalDirection, Q) * inv_det;
  if (v < T() || u + v  > T(1)) {
    return false;
  }
  t = dotProduct3D(e2, Q)*inv_det;
  if (t > EPSILON) {
    distance = t;
    return true;
  }
  return false;
}

template <typename T>
bool triangleIntersection(const std::vector<T>& point0, const std::vector<T>& point1, const std::vector<T>& point2, const std::vector<T>& origin, const std::vector<T>& direction, T& distance)
{
  std::vector<T> normalDirection(normalize(direction) );
  return triangleIntersectionWithNormalDirection(point0, point1, point2, origin, normalDirection, distance );
}
*/
template <typename T>
std::vector<T> assign(T a, T b)
{
  std::vector<T> v1;
  v1.push_back(a);
  v1.push_back(b);
  return v1;
}

template <typename T>
std::vector<T> assign(T a, T b, T c)
{
  std::vector<T> v1;
  v1.push_back(a);
  v1.push_back(b);
  v1.push_back(c);
  return v1;
}

/// prints a vector of arbitrary length
template <typename T>
void print(const T& a, std::string name="", OstreamManager clout = OstreamManager(std::cout,"print"))
{
  if (name != "") {
    clout << name << "=";
  }
  clout << a << std::endl;
}

/// prints a vector of arbitrary length
template <typename T>
void print(const std::vector<T>& a, std::string name="", OstreamManager clout = OstreamManager(std::cout,"print"))
{
  if (name != "") {
    clout << name << "=";
  }
  clout << "(";
  for (unsigned iD=0; iD<a.size()-1; iD++) {
    clout << a[iD] << ",";
  }
  clout << a[a.size()-1] << ")" << std::endl;
}

/// prints a vector of arbitrary length
template <typename T>
void print(const T a[2], std::string name="", OstreamManager clout = OstreamManager(std::cout,"print"))
{
  if (name != "") {
    clout << name << "=";
  }
  unsigned size = 2;
  clout << "(";
  for (unsigned iD=0; iD<size-1; iD++) {
    clout << a[iD] << ",";
  }
  clout << a[size-1] << ")" << std::endl;
}

/// prints a vector of arbitrary length
template <typename T>
void print(const T a[3], const unsigned& size, std::string name="", OstreamManager clout = OstreamManager(std::cout,"print"))
{
  if (name != "") {
    clout << name << "=";
  }
  clout << "(";
  for (unsigned iD=0; iD<size-1; iD++) {
    clout << a[iD] << ",";
  }
  clout << a[size-1] << ")" << std::endl;
}

} // namespace util
} // namespace olb

#endif
