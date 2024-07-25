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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI2
#define M_PI2 1.57079632679489661923
#endif

namespace olb {

template <typename T>
IndicatorSDF2D<T>::IndicatorSDF2D(std::function<T(Vector<T,2>)> f):
  _f(f) { }

template <typename T>
bool IndicatorSDF2D<T>::operator()(bool output[], const T input[])
{
  output[0] = _f(input) <= 0.0;
  return true;
}

template <typename T>
IndicatorSDF3D<T>::IndicatorSDF3D(std::function<T(Vector<T,3>)> f):
  _f(f) { }

template <typename T>
bool IndicatorSDF3D<T>::operator()(bool output[], const T input[])
{
  output[0] = _f(input) <= 0.0;
  return true;
}

namespace sdf {

template <typename T>
T mix(T a, T b, T h)
{
  return b*(1.0-h) + a*h;
}

template <typename T0, typename T1, typename T2>
decltype(T0{}*T1{}*T2{}) clamp(T0 x, T1 a, T2 b)
{
  if (x < a) {
    return a;
  }
  else if (x > b) {
    return b;
  }
  else {
    return x;
  }
}

/// A rough test for symmetry
template <typename T>
Vector<bool,3> isSymmetric(std::function<T(Vector<T,3>)> sdf, Vector<T,3> center)
{
  Vector<bool,3> isSymmetric;

  isSymmetric[0] = sdf(Vector<T,3> {1., 1., 1.} + center) == sdf(Vector<T,3> {-1., 1., 1.} + center) && sdf(Vector<T,3> {1., -1., 1.} + center) == sdf(Vector<T,3> {-1., -1., 1.} + center) &&
                   sdf(Vector<T,3> {1., 1., -1.} + center) == sdf(Vector<T,3> {-1., 1., -1.} + center) && sdf(Vector<T,3> {1., -1., -1.} + center) == sdf(Vector<T,3> {-1., -1., -1.} + center);
  isSymmetric[1] = sdf(Vector<T,3> {1., 1., 1.} + center) == sdf(Vector<T,3> {1., -1., 1.} + center) && sdf(Vector<T,3> {-1., 1., 1.} + center) == sdf(Vector<T,3> {-1., -1., 1.} + center) &&
                   sdf(Vector<T,3> {1., 1., -1.} + center) == sdf(Vector<T,3> {1., -1., -1.} + center) && sdf(Vector<T,3> {-1., 1., -1.} + center) == sdf(Vector<T,3> {-1., -1., -1.} + center);
  isSymmetric[2] = sdf(Vector<T,3> {1., 1., 1.} + center) == sdf(Vector<T,3> {1., 1., -1.} + center) && sdf(Vector<T,3> {-1., 1., 1.} + center) == sdf(Vector<T,3> {-1., 1., -1.} + center) &&
                   sdf(Vector<T,3> {1., -1., 1.} + center) == sdf(Vector<T,3> {1., -1., -1.} + center) && sdf(Vector<T,3> {-1., -1., 1.} + center) == sdf(Vector<T,3> {-1., -1., -1.} + center);

  return isSymmetric;
}


/** Exact signed distance to the surface of circle (in 2D) or sphere (in 3D).
 *
 * \param p         point for which the distance to the surface is calculated
 * \param r         radius describing the geometry
*/
template <typename T, unsigned D>
T sphere(Vector<T,D> p, T r)
{
  return norm(p) - r;
}

/** Exact signed distance to the surface of two-dimensional cuboid.
 *
 * \param p         point for which the distance to the surface is calculated
 * \param b         vector containing half side lengths to describe the cuboid
*/
template <typename T>
T box(Vector<T,2> p, Vector<T,2> b)
{
  Vector<T,2> q = abs(p) - b;
  return norm(maxv(q, Vector<T,2>(0.0))) + util::min(util::max({q[0],q[1]}), T());
}

/** Exact signed distance to the surface of three-dimensional cuboid.
 *
 * \param p         point for which the distance to the surface is calculated
 * \param b         vector containing half side lengths to describe the cuboid
*/
template <typename T>
T box(Vector<T,3> p, Vector<T,3> b)
{
  Vector<T,3> q = abs(p) - b;
  return norm(maxv(q, Vector<T,3>(0.0))) + util::min(util::max({q[0],q[1],q[2]}), T());
}

/** Exact signed distance to the surface of two-dimensional triangle.
 *
 * \param p         point for which the distance to the surface is calculated
 * \param a         vector containing first vertex
 * \param b         vector containing second vertex
 * \param c         vector containing third vertex
*/
template <typename T>
T triangle(Vector<T,2> p, Vector<T,2> a, Vector<T,2> b, Vector<T,2> c)
{
  Vector<T,2> e0 = b - a;
  Vector<T,2> v0 = p - a;
  Vector<T,2> e1 = c - b;
  Vector<T,2> v1 = p - b;
  Vector<T,2> e2 = a - c;
  Vector<T,2> v2 = p - c;
  Vector<T,2> pq0 = v0 - e0 * clamp( (v0 * e0) / (e0 * e0), 0., 1.);
  Vector<T,2> pq1 = v1 - e1 * clamp( (v1 * e1) / (e1 * e1), 0., 1.);
  Vector<T,2> pq2 = v2 - e2 * clamp( (v2 * e2) / (e2 * e2), 0., 1.);
  T s = util::sign(e0[0]*e2[1] - e0[1]*e2[0]);

  T dx = util::min( util::min( (pq0*pq0),
                               (pq1*pq1)),
                    (pq2*pq2));
  T dy = util::min( util::min( s*(v0[0]*e0[1]-v0[1]*e0[0]),
                               s*(v1[0]*e1[1]-v1[1]*e1[0])),
                    s*(v2[0]*e2[1]-v2[1]*e2[0]));
  return -util::sqrt(dx)*util::sign(dy);


}

/** Exact signed distance to the surface of three-dimensional cylinder.
 *
 * \param p         point for which the distance to the surface is calculated
 * \param a         vector containing center of the first circular end of the cylinder
 * \param ba        vector representing the axis of the cylinder
 * \param baba      length of the cylinder height squared
*/
template <typename T>
T cylinder(Vector <T,3> p, Vector <T,3> a, Vector<T,3> ba, T baba, T r)
{
  Vector<T,3> pa = p-a;
  T paba = pa * ba;
  T t = paba / baba;
  T x = norm(p - a - t*ba) - r;
  T y = (util::abs(t-0.5) - 0.5) * util::sqrt(baba);
  return ( util::min(util::max(x, y), T()) + util::sqrt ( util::max(x, T()) * util::max(x, T()) + util::max(y, T()) * util::max(y, T()) ) );
}

/** Calculate signed distance to the surface of three-dimensional cylinder defined by the centers of the two extremities and the radius.
 *
 * \param p         point for which the distance to the surface is calculated
 * \param a         vector containing center of the first circular end of the cylinder
 * \param b         vector containing center of the second circular end of the cylinder
 * \param r         radius of the cylinder
*/
template <typename T>
T cylinder(Vector<T,3> p, Vector<T,3> a, Vector<T,3> b, T r)
{
  Vector<T,3> ba = b - a;
  T baba = ba * ba;
  return sdf::cylinder(p, a, ba, baba, r);
}

/** Exact signed distance to the surface of three-dimensional (capped) cone.
 *
 * \param p         point for which the distance to the surface is calculated
 * \param a         vector containing center of the first circular end of the (capped) cone
 * \param ba        vector representing the axis of the (capped) cone
 * \param baba      length of the cone height squared
*/
template <typename T>
T cone(Vector<T,3> p, Vector<T,3> a, Vector<T,3> ba, T baba, T ra, T rb)
{
  T rba = rb - ra;
  T papa = (p-a) * (p-a);
  T paba = (p-a) * (ba) / baba;
  // x distance to the axis of the cone
  // note: abs() prevents from negative values, e.g. caused by rounding errors
  T x = util::sqrt(util::abs(papa - paba*paba*baba));
  // cax and cay for the distances to the caps
  T cax = util::max(0.0,x-((paba<0.5)?ra:rb));
  T cay = util::abs(paba-0.5)-0.5;
  // cbx and cby for the distance to the side wall
  T k = rba*rba + baba;
  T f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );
  T cbx = x-ra - f*rba;
  T cby = paba - f;
  T s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
  return s*util::sqrt( util::min(cax*cax + cay*cay*baba, cbx*cbx + cby*cby*baba) );
}

/** Calculate signed distance to the surface of three-dimensional capped cone defined by the centers of the two extremities and the corresponding radius.
 * For the calculation of a cone, the second radius is set to zero
 * \param p         point for which the distance to the surface is calculated
 * \param a         vector containing center of the first circular end of the capped cone
 * \param b         vector containing center of the second circular end of the capped cone or the apex of the cone
 * \param ra        radius of the first circular end of the capped cone
 * \param rb        radius of the second circular end of the capped cone
*/
template <typename T>
T cone(Vector<T,3> p, Vector<T,3> a, Vector<T,3> b, T ra, T rb)
{
  Vector<T,3> ba = b-a;
  T baba = ba*ba;
  return sdf::cone(p, a, ba, baba, ra, rb);
}

/** Calculate the NOT EXACT (!) signed distance to the surface of three-dimensional ellipsoid.
 * \param p         point for which the distance to the surface is calculated
 * \param r         vector containing half the length of the principal axes
*/
template <typename T>
T ellipsoid(Vector<T,3> p, Vector<T,3> r)
{
  const Vector<T,3> a(p[0]/r[0], p[1]/r[1], p[2]/r[2]);
  T k0 = norm(a);
  const Vector<T,3> r2(r[0]*r[0], r[1]*r[1], r[2]*r[2]);
  const Vector<T,3> b(p[0]/r2[0], p[1]/r2[1], p[2]/r2[2]);
  T k1 = norm(b);
  return (k0<1.0) ? (k0-1.0)*util::min(util::min(r[0],r[1]),r[2]) : k0*(k0-1.0)/k1;
}

/** Exact signed distance to the surface of a torus placed in the XZ-plane of the coordinate system.
 *
 * \param p         point for which the distance to the surface is calculated
 * \param t[0]      distance center of the tube to center of the torus
 * \param t[1]      radius of the tube
*/
template <typename T>
T torus(Vector<T,3> p, Vector<T,2> t)
{
  Vector<T,2> b{p[0], p[2]};
  Vector<T,2> q{norm(b)-t[0], p[1]};
  return norm(q) - t[1];
}

/** Exact signed distance to the surface of a solid angle which represents a combination of a cone (axis y=0 to y=r) and a sphere
 *
 * \param p         point for which the distance to the surface is calculated
 * \param c[0]      sin of the angle
 * \param c[1]      cos of the angle
 * \param r         radius of the sphere
*/
template <typename T>
T solidAngle(Vector<T,3> p, Vector<T,2> c, T r)
{
  Vector<T,2> q{norm(Vector<T,2>{p[0], p[2]}), p[1]};
  T l = norm(q) - r;
  T m = norm(q - c*clamp(q * c, 0.0, r));
  return util::max(l, m*util::sign(c[1]*q[0]-c[0]*q[1]));
}

/** Translation: The output of this function is used for the calculation of the signed distance to translated/shifted objects
 * \param p         point for which the distance to the surface is calculated
 * \param origin    point to which the object is shifted to
*/
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

/** Volume of a is subtracted from b.
 * The returned distance is not exact if the signed Distance of both objects
 * points to the part of the surface which is removed for creating the combined object
 * \param a         signed Distance of the first object
 * \param b         signed Distance of the second object
*/
template <typename T>
T subtraction(T a, T b)
{
  return util::max(-a, b);
}

/** Volume of a and volume of b are combined as a new object.
 * The returned distance is not exact if the signed Distance of both objects
 * points to the part of the surface which is removed for creating the combined object
 * \param a         signed Distance of the first object
 * \param b         signed Distance of the second object
*/
template <typename T>
T unify(T a, T b)
{
  return util::min(a, b);
}

/** Volume which is shared by a and b creates a new object.
 * The returned distance is not exact if the signed Distance of both objects
 * points to the part of the surface which is removed for creating the combined object
 * \param a         signed Distance of the first object
 * \param b         signed Distance of the second object
*/
template <typename T>
T intersection(T a, T b)
{
  return util::max(a, b);
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
T smooth_intersection(T a, T b, T k)
{
  T h = clamp(0.5 - 0.5*(b-a)/k, 0.0, 1.0);
  return mix(b, a, h) + k*h*(1.0-h);
}

/** Computes a layer of a constant thickness around the surface
 * \param a         signed Distance of the object
 * \param r         layer thickness
*/
template <typename T>
T rounding(T a, T r)
{
  return a - r;
}

/** Elongation splits the object in 2 (4 or 8) parts, moves them apart and connects them again
 * The object has to be placed in the origin of the coodinate system.
 * Symmetry is required so that the splitted parts are symmetric
 * \param sdf       signed Distance Function of the object which is altered
 * \param p         point for which the distance to the surface is calculated
 * \param h     vector containing the information how far the splitted parts are moved apart
 * \param center  vector containing the center, to shift the object to the origin of the coordinate system for the calculation of the elongation
*/
template <typename T>
T elongation(std::function<T(Vector<T,3>)> sdf, Vector<T,3> p, Vector<T,3> h, Vector<T,3> center)
{
  Vector<T,3> q = abs(p - center) - h;
  Vector<T,3> p_elong = p - center;
  // elongation requires symmetry e.g. for elongation in x --> plane of symmetry = YZ-plane
  Vector<bool,3> sdfIsSymmetric = isSymmetric(sdf, center);

  if (h[0] != 0) {
    p_elong[0] = util::max(q[0], 0.);
    if ( !sdfIsSymmetric[0] ) {
      std::cout << "Warning: symmetry in x is not met" << std::endl;
    }
  }

  if (h[1] != 0) {
    p_elong[1] = util::max(q[1], 0.);
    if ( !sdfIsSymmetric[1] ) {
      std::cout << "Warning: symmetry in y is not met" << std::endl;
    }
  }

  if (h[2] != 0) {
    p_elong[2] = util::max(q[2], 0.);
    if ( !sdfIsSymmetric[2] ) {
      std::cout << "Warning: symmetry in z is not met"<< std::endl;
    }
  }

  // Second term needed in case of 3D-Elongation if all parts of Vector q are negative
  return sdf(p_elong + center)
         + util::min(util::max({q[0], q[1], q[2]}), 0.);
}

/// Converts signed distance to values for the smooth epsilon boundary
// TODO: Rename function
template <typename T>
T signedDistanceToPorosity(T signedDist, T eps)
{
  const T d = signedDist+.5*eps;
  return util::pow(util::cos(M_PI2*d/eps), 2);
}

template <typename T>
bool evalPorosity( T output[], T signedDist, T eps )
{
  T const halfEps = .5*eps;
  if (signedDist <= -halfEps) {
    output[0] = 1.;
    return true;
  }
  else if (signedDist < halfEps) {
    output[0] = signedDistanceToPorosity(signedDist, eps);
    return true;
  }
  output[0] = 0.;
  return false;
}

}

}

#endif
