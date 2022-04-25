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

#ifndef SDF_H
#define SDF_H

#include "indicatorF2D.h"
#include "indicatorF3D.h"

#include "core/vector.h"
#include <functional>

namespace olb {

template <typename T>
class IndicatorSDF2D : public IndicatorF2D<T> {
private:
  std::function<T(Vector<T,2>)> _f;

public:
  IndicatorSDF2D(std::function<T(Vector<T,2>)> f);

  bool operator() (bool output[], const T input[]) override;
};

template <typename T>
class IndicatorSDF3D : public IndicatorF3D<T> {
private:
  std::function<T(Vector<T,3>)> _f;

public:
  IndicatorSDF3D(std::function<T(Vector<T,3>)> f);

  bool operator() (bool output[], const T input[]) override;
};

namespace sdf {

template <typename T, unsigned D>
T sphere(Vector<T,D> p, T r);

template <typename T>
T box(Vector<T,2> p, Vector<T,2> b);

template <typename T>
T box(Vector<T,3> p, Vector<T,3> extend);

template <typename T>
T triangle(Vector<T,2> p, Vector<T,2> a, Vector<T,2> b, Vector<T,2> c);

template <typename T>
T cylinder(Vector<T,3> p, Vector<T,3> a, Vector<T,3> ba, T baba, T r);

template <typename T>
T cylinder(Vector<T,3> p, Vector<T,3> a, Vector<T,3> b, T r);

template <typename T>
T cone(Vector<T,3> p, Vector<T,3> a, Vector<T,3> ba, T baba, T ra, T rb);

template <typename T>
T cone(Vector<T,3> p, Vector<T,3> a, Vector<T,3> b, T ra, T rb);

template <typename T>
T ellipsoid(Vector<T,3> p, Vector<T,3> r);

template <typename T>
T torus(Vector<T,3> p, Vector<T,2> t);

template <typename T>
T solidAngle(Vector<T,3> p, Vector<T,2> c, T r);

template <typename T, unsigned D>
Vector<T,D> translate(Vector<T,D> p, Vector<T,D> origin);

template <typename T>
Vector<T,3> flip(Vector<T,3> p);

template <typename T>
T subtraction(T a, T b);

template <typename T>
T unify(T a, T b);

template <typename T>
T intersection(T a, T b);

template <typename T>
T smooth_union(T d1, T d2, T k);

template <typename T>
T smooth_subtraction(T d1, T d2, T k);

template <typename T>
T smooth_intersection(T d1, T d2, T k);

template <typename T>
T rounding(T a, T r);

template <typename T>
T elongation(std::function<T(Vector<T,3>)> sdf, Vector<T,3> p, Vector<T,3> h, Vector<T,3> center=( T(0) ));

template <typename T>
T signedDistanceToPorosity(T signedDist, T eps);

template <typename T>
bool evalPorosity( T output[], T signedDist, T eps );

}

}

#endif
