/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender, Jan E. Marquardt
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

#include "sdf.h"
#include "sdf.hh"

namespace olb {

// indicators
template class IndicatorSDF2D<double>;
template class IndicatorSDF3D<double>;

namespace sdf{

// help functions
template double mix(double a, double b, double h);
template double clamp(double x, double a, double b);

// primitive geometries
template double sphere(Vector<double,2> p, double r);
template double sphere(Vector<double,3> p, double r);
template double box(Vector<double,2>p, Vector<double,2> b);
template double box(Vector<double,3>p, Vector<double,3> b);
template double cylinder(Vector <double,3> p, Vector <double,3> a, Vector<double,3> ba, double baba, double r);
template double cylinder(Vector<double,3> p, Vector<double,3> a, Vector<double,3> b, double r);
template double cone(Vector<double,3> p, Vector<double,3> a, Vector<double,3> ba, double baba, double ra, double rb);
template double cone(Vector<double,3> p, Vector<double,3> a, Vector<double,3> b, double ra, double rb);
template double torus(Vector<double,3> p, Vector<double,2> t);
template double solidAngle(Vector<double,3> p, Vector<double,2> c, double r);

// combinations
template double subtraction(double a, double b);
template double unify(double a, double b);
template double intersection(double a, double b);
template double smooth_union(double d1, double d2, double k);
template double smooth_subtraction(double d1, double d2, double k);
template double smooth_intersection(double d1, double d2, double k);

// positioning
template Vector<double,2> translate(Vector<double,2> p, Vector<double,2> origin);
template Vector<double,3> translate(Vector<double,3> p, Vector<double,3> origin);
template Vector<double,3> flip(Vector<double,3> p);

// primitive alterations
template double rounding(double a, double r);
template double elongation(std::function<double(Vector<double,3>)> sdf, Vector<double,3> p, Vector<double,3> h, Vector<double,3> center);
template bool evalPorosity( double output[], double signedDist, double eps );

} // namespace sdf
} // namespace olb
