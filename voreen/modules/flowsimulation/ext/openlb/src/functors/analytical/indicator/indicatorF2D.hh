/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Cyril Masquelier, Jan Marquardt, Mathias J. Krause
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

#ifndef INDICATOR_F_2D_HH
#define INDICATOR_F_2D_HH

#include "utilities/omath.h"

#include "indicatorF2D.h"
#include "utilities/vectorHelpers.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define M_PI2 1.57079632679489661923

namespace olb {

// Warning : the cuboid is only defined parallel to the plans x=0 and y=0 !!!
// For other cuboids, please use the parallelepiped version
template <typename S>
IndicatorF2DfromIndicatorF3D<S>::IndicatorF2DfromIndicatorF3D(IndicatorF3D<S>& indicator3D)
  : _indicator3D(indicator3D)
{
  this->_myMin[0] = _indicator3D.getMin()[0];
  this->_myMin[1] = _indicator3D.getMin()[1];
  this->_myMax[0] = _indicator3D.getMax()[0];
  this->_myMax[1] = _indicator3D.getMax()[1];
}

template <typename S>
bool IndicatorF2DfromIndicatorF3D<S>::operator()(bool output[], const S input[])
{
  S input3D[3];
  input3D[0] = input[0];
  input3D[1] = input[1];
  input3D[2] = (_indicator3D.getMax()[2] - _indicator3D.getMin()[2]) * 0.5 + _indicator3D.getMin()[2];
  _indicator3D(output, input3D);
  return true;
}



// Warning : the cuboid is only defined parallel to the plans x=0 and y=0 !!!
// For other cuboids, please use the parallelepiped version
template <typename S>
IndicatorCuboid2D<S>::IndicatorCuboid2D(Vector<S,2> extend, Vector<S,2> origin, S theta)
  : _center(origin + S(.5)*extend), _xLength(extend[0]), _yLength(extend[1]), _theta(theta)

{
  this->_myMin = origin;
  this->_myMax = origin + extend;
}

template <typename S>
IndicatorCuboid2D<S>::IndicatorCuboid2D(S xLength, S yLength, Vector<S,2> center, S theta )
  : _center(center), _xLength(xLength), _yLength(yLength), _theta(-theta)
{
  this->_myMin = {_center[0] - _xLength/2., _center[1] - _yLength/2.};
  this->_myMax = {_center[0] + _xLength/2., _center[1] + _yLength/2.};
}


// returns true if x is inside the cuboid
template <typename S>
bool IndicatorCuboid2D<S>::operator()(bool output[], const S input[])
{
  S x, y;
  if ( !util::nearZero(_theta) ) {
    x = _center[0] + (input[0] - _center[0])*util::cos(_theta) - (input[1] - _center[1])*util::sin(_theta);
    y = _center[1] + (input[0] - _center[0])*util::sin(_theta) + (input[1] - _center[1])*util::cos(_theta);
  }
  else {
    x = input[0];
    y = input[1];
  }

  output[0] = (  (util::fabs(_center[0] - x) < _xLength/2. || util::approxEqual(util::fabs(_center[0] - x),_xLength/2.) )
                 && (util::fabs(_center[1] - y) < _yLength/2. || util::approxEqual(util::fabs(_center[1] - y), _yLength/2.) ) );
  return true;
}

template <typename S>
Vector<S,2> const& IndicatorCuboid2D<S>::getCenter() const
{
  return _center;
}

template <typename S>
S const IndicatorCuboid2D<S>::getxLength() const
{
  return _xLength;
}

template <typename S>
S const IndicatorCuboid2D<S>::getyLength() const
{
  return _yLength;
}

template <typename S>
S IndicatorCuboid2D<S>::signedDistance( const Vector<S,2>& input )
{
  // TODO: Implementation should be analogous to other indicators.

  Vector<S,2> ptransl = {input[0]-_center[0], input[1]-_center[1]};
  Vector<S,2> prot = {util::cos(-_theta)*ptransl[0]+util::sin(-_theta)*ptransl[1], util::cos(-_theta)*ptransl[1]-util::sin(-_theta)*ptransl[0]};

  return sdf::box(prot, Vector<S,2>(.5*_xLength, .5*_yLength));
}


// creator function
template <typename S>
IndicatorCuboid2D<S>* createIndicatorCuboid2D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCuboid2D");
  params.setWarningsOn(verbose);

  Vector<S,2> center;
  S xLength;
  S yLength;

  std::stringstream xmlCenter( params.getAttribute("center") );
  xmlCenter >> center[0] >> center[1];
  std::stringstream xmlRadius( params.getAttribute("length") );
  xmlRadius >> xLength >> yLength;

  return new IndicatorCuboid2D<S>(xLength, yLength, center);
}


template <typename S>
IndicatorCircle2D<S>::IndicatorCircle2D(Vector<S,2> center, S radius)
  :  _center(center),
     _radius(radius),
     _radius2(radius*radius)
{
  this->_myMin = _center - radius;
  this->_myMax = _center + radius;
}

template <typename S>
Vector<S,2> const& IndicatorCircle2D<S>::getCenter() const
{
  return _center;
}

template <typename S>
S const IndicatorCircle2D<S>::getRadius() const
{
  return _radius;
}

// returns true if x is inside the circle
template <typename S>
bool IndicatorCircle2D<S>::operator()(bool output[], const S input[])
{
  output[0] = ( util::pow(_center[0] - input[0],2) + util::pow(_center[1] - input[1], 2) <= _radius2 );
  return output[0];
}


template <typename S>
bool IndicatorCircle2D<S>::distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC)
{
  S a = direction[0]*direction[0] + direction[1]*direction[1];

  // returns 0 if point is at the boundary of the sphere
  if ( util::approxEqual(a,_radius2) ) {
    distance = S(0);
    return true;
  }
  // norm of direction
  a = util::sqrt(a);

  S b = 2.*((origin[0] - _center[0])*direction[0] +
            (origin[1] - _center[1])*direction[1])/a;
  S c = -_radius2 + (origin[0] - _center[0])*(origin[0] - _center[0])
        + (origin[1] - _center[1])*(origin[1] - _center[1]);

  // discriminant
  S d = b*b - 4.*c;
  if (d < 0) {
    return false;
  }

  S x1 = (- b + util::sqrt(d)) *0.5;
  S x2 = (- b - util::sqrt(d)) *0.5;

  // case if origin is inside the sphere
  if ((x1<0.) || (x2<0.)) {
    if (x1>0.) {
      distance = x1;
      return true;
    }
    if (x2>0.) {
      distance = x2;
      return true;
    }
  }
  // case if origin is ouside the sphere
  else {
    distance = util::min(x1,x2);
    return true;
  }

  return false;
}


template <typename S>
bool IndicatorCircle2D<S>::normal(Vector<S,2>& normal, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC)
{
  S dist;
  if (!(distance(dist, origin, direction, iC)) ) {
    return false;
  }

  Vector<S,2> intresection(origin + dist*direction);

  normal = intresection - _center;

  return true;
}

template <typename S>
S IndicatorCircle2D<S>::signedDistance(const Vector<S,2>& input)
{
  Vector<S,2> p = input - _center;
  return sdf::sphere(p, _radius);
}

template <typename S>
IndicatorCircle2D<S>* createIndicatorCircle2D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCircle2D");
  params.setWarningsOn(verbose);

  Vector<S,2> center;
  S radius = 1;

  std::stringstream xmlCenter( params.getAttribute("center") );
  xmlCenter >> center[0] >> center[1];
  std::stringstream xmlRadius( params.getAttribute("radius") );
  xmlRadius >> radius;

  return new IndicatorCircle2D<S>(center, radius);
}

template <typename S>
IndicatorTriangle2D<S>::IndicatorTriangle2D(Vector<S,2> a, Vector<S,2> b, Vector<S,2> c)
  :  _a(a),
     _b(b),
     _c(c)
{

  this->_myMin = {util::min(_a[0], util::min(_b[0], _c[0])), util::min(_a[1], util::min(_b[1], _c[1]))};
  this->_myMax = {util::max(_a[0], util::max(_b[0], _c[0])), util::max(_a[1], util::max(_b[1], _c[1]))};
}


template <typename S>
Vector<S,2> const& IndicatorTriangle2D<S>::getVertexA() const
{
  return _a;
}

template <typename S>
Vector<S,2> const& IndicatorTriangle2D<S>::getVertexB() const
{
  return _b;
}

template <typename S>
Vector<S,2> const& IndicatorTriangle2D<S>::getVertexC() const
{
  return _c;
}

// returns true if x is inside the triangle
template <typename S>
bool IndicatorTriangle2D<S>::operator()(bool output[], const S input[])
{
  output[0] = ( sdf::triangle({input[0], input[1]}, _a, _b, _c) <= 0 );
  return output[0];
}

template <typename S>
S IndicatorTriangle2D<S>::signedDistance(const Vector<S,2>& input)
{

  return sdf::triangle(input, _a, _b, _c);
}

template <typename S>
IndicatorTriangle2D<S>* createIndicatorTriangle2D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorTriangle2D");
  params.setWarningsOn(verbose);

  Vector<S,2> a;
  Vector<S,2> b;
  Vector<S,2> c;

  std::stringstream xmla( params.getAttribute("a") );
  xmla>> a[0] >> a[1];
  std::stringstream xmlb( params.getAttribute("b") );
  xmlb >> b[0] >> b[1];
  std::stringstream xmlc( params.getAttribute("c") );
  xmlc >> c[0] >> c[1];

  return new IndicatorTriangle2D<S>(a, b, c);
}

template <typename S>
IndicatorEquiTriangle2D<S>::IndicatorEquiTriangle2D(Vector<S,2> center, S radius)
  :  _center(center),
     _radius(radius),
     _a({center[0], center[1]+radius}),
_b({center[0]-util::sqrt(3)/2.*radius, center[1]-0.5*radius}),
_c({center[0]+util::sqrt(3)/2.*radius, center[1]-0.5*radius})
{

  this->_myMin = {util::min(_a[0], util::min(_b[0], _c[0])), util::min(_a[1], util::min(_b[1], _c[1]))};
  this->_myMax = {util::max(_a[0], util::max(_b[0], _c[0])), util::max(_a[1], util::max(_b[1], _c[1]))};
}


template <typename S>
Vector<S,2> const& IndicatorEquiTriangle2D<S>::getVertexA() const
{
  return _a;
}

template <typename S>
Vector<S,2> const& IndicatorEquiTriangle2D<S>::getVertexB() const
{
  return _b;
}

template <typename S>
Vector<S,2> const& IndicatorEquiTriangle2D<S>::getVertexC() const
{
  return _c;
}

template <typename S>
Vector<S,2> const& IndicatorEquiTriangle2D<S>::getCenter() const
{
  return _center;
}

template <typename S>
S const IndicatorEquiTriangle2D<S>::getRadius() const
{
  return _radius;
}

// returns true if x is inside the triangle
template <typename S>
bool IndicatorEquiTriangle2D<S>::operator()(bool output[], const S input[])
{
  output[0] = ( sdf::triangle({input[0], input[1]}, _a, _b, _c) <= 0 );
  return output[0];
}

template <typename S>
S IndicatorEquiTriangle2D<S>::signedDistance(const Vector<S,2>& input)
{

  return sdf::triangle( input, _a, _b, _c);
}

template <typename S>
IndicatorEquiTriangle2D<S>* createIndicatorEquiTriangle2D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorEquiTriangle2D");
  params.setWarningsOn(verbose);

  Vector<S,2> center;
  S radius = 1;

  std::stringstream xmlCenter( params.getAttribute("center") );
  xmlCenter>> center[0] >> center[1];
  std::stringstream xmlRadius( params.getAttribute("radius") );
  xmlRadius >> radius;

  return new IndicatorEquiTriangle2D<S>(center, radius);
}

template <typename S>
IndicatorLayer2D<S>::IndicatorLayer2D(IndicatorF2D<S>& indicatorF, S layerSize)
  :  _indicatorF(indicatorF), _layerSize(layerSize)
{
  this->_myMin = indicatorF.getMin() - layerSize;
  this->_myMax = indicatorF.getMax() + layerSize;
  OLB_ASSERT( (this->_myMax[0]-this->_myMin[0]) > std::numeric_limits<S>::epsilon(),"Indicator reduced to zero-set in x direction");
  OLB_ASSERT( (this->_myMax[1]-this->_myMin[1]) > std::numeric_limits<S>::epsilon(),"Indicator reduced to zero-set in y direction");
  _isPositive = std::signbit(layerSize);
}

// returns true if x is inside the layer
template <typename S>
bool IndicatorLayer2D<S>::operator()(bool output[], const S input[])
{
  output[0] = !_isPositive;
  S r[2];
  for (int iX =- 1; iX < 2; ++iX) {
    for (int iY =- 1; iY < 2; ++iY) {
      r[0] = input[0] + iX*_layerSize;
      r[1] = input[1] + iY*_layerSize;
      _indicatorF(output,r);
      if (output[0] == !_isPositive) {
        return true;
      }
    }
  }
  return true;
}

} // namespace olb

#endif
