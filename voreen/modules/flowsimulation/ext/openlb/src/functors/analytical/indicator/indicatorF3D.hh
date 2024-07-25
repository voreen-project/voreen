/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Albert Mink
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

#ifndef INDICATOR_F_3D_HH
#define INDICATOR_F_3D_HH

#include <vector>
#include "utilities/omath.h"
#include <cassert>
#include <sstream>

#include "indicatorF3D.h"
#include "indicComb3D.h"
#include "utilities/vectorHelpers.h"


namespace olb {

template <typename S>
IndicatorTranslate3D<S>::IndicatorTranslate3D(std::array<S,3> translate, IndicatorF3D<S>& indicator)
  : _translate(translate), _indicator(indicator)
{
}

template< typename S>
bool IndicatorTranslate3D<S>::operator() (bool output[], const S input[] )
{
  const S inputTranslated[3] = {
    input[0] - _translate[0],
    input[1] - _translate[1],
    input[2] - _translate[2]
  };

  _indicator.operator()(output, inputTranslated);
  return output[0];
}

template <typename S>
S IndicatorTranslate3D<S>::signedDistance( const Vector<S,3>& input )
{
  Vector<S,3> translate = sdf::translate(input, {_translate[0], _translate[1], _translate[2]});
  return _indicator.signedDistance(translate);
}



template <typename S>
IndicatorCircle3D<S>::IndicatorCircle3D(Vector<S,3> center, Vector<S,3> normal, S radius)
  :  _center(center), _normal(normal), _radius2(radius*radius)
{
  _normal = normalize(_normal);
  this->_myMin = _center - radius;
  this->_myMax = _center + radius;
}

template <typename S>
IndicatorCircle3D<S>::IndicatorCircle3D(S center0, S center1, S center2,
                                        S normal0, S normal1, S normal2, S radius)
  :  _radius2(radius*radius)
{
  _center = {center0, center1, center2};
  _normal = {normal0, normal1, normal2};
  _normal = normalize(_normal);
  this->_myMin = _center - radius;
  this->_myMax = _center + radius;
}

// returns true if x is inside the circle
template <typename S>
bool IndicatorCircle3D<S>::operator()(bool output[], const S input[])
{
  S eps = std::numeric_limits<S>::epsilon();
  IndicatorCylinder3D<S> cylinder(_center, _normal, getRadius(), 5*eps);
  return cylinder(output,input);
}

template <typename S>
Vector<S,3> const& IndicatorCircle3D<S>::getCenter() const
{
  return _center;
}

template <typename S>
Vector<S,3> const& IndicatorCircle3D<S>::getNormal() const
{
  return _normal;
}

template <typename S>
S IndicatorCircle3D<S>::getRadius() const
{
  return util::sqrt(_radius2);
}



template <typename S>
IndicatorSphere3D<S>::IndicatorSphere3D(Vector<S,3> center, S radius)
  :  _center(center), _radius(radius), _radius2(_radius*_radius)
{
  this->_myMin = _center - radius;
  this->_myMax = _center + radius;
}

template <typename S>
IndicatorSphere3D<S>::IndicatorSphere3D(const IndicatorSphere3D& sphere)
{
  this->_myMin = sphere._myMin;
  this->_myMax = sphere._myMax;
  _center = sphere._center;
  _radius = sphere._radius;
  _radius2 = sphere._radius2;
}

template <typename S>
Vector<S,3> const& IndicatorSphere3D<S>::getCenter() const
{
  return _center;
}

template <typename S>
S const IndicatorSphere3D<S>::getRadius() const
{
  return _radius;
}

template <typename S>
S IndicatorSphere3D<S>::signedDistance( const Vector<S,3>& input )
{
  Vector<S,3> p = input - _center;
  return sdf::sphere(p, _radius);
}

template <typename S>
bool IndicatorSphere3D<S>::distance(S& distance, const Vector<S,3>& origin,
                                    const Vector<S,3>& direction, int iC)
{
  // computes pos. distance by solving quadratic equation by a-b-c-formula
  S a = direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2];

  // returns 0 if point is at the boundary of the sphere
  if ( util::approxEqual(a,_radius2) ) {
    distance = S();
    return true;
  }
  // norm of direction
  a = util::sqrt(a);

  S b = 2.*((origin[0] - _center[0])*direction[0] +
            (origin[1] - _center[1])*direction[1] +
            (origin[2] - _center[2])*direction[2])/a;
  S c = -_radius2 + (origin[0] - _center[0])*(origin[0] - _center[0])
        + (origin[1] - _center[1])*(origin[1] - _center[1])
        + (origin[2] - _center[2])*(origin[2] - _center[2]);

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
  // case if origin is outside the sphere
  else {
    distance = util::min(x1,x2);
    return true;
  }

  return false;
}


template <typename S>
IndicatorLayer3D<S>::IndicatorLayer3D(IndicatorF3D<S>& indicatorF, S layerSize)
  :  _indicatorF(indicatorF), _layerSize(layerSize)
{
  this->_myMin = indicatorF.getMin() - layerSize;
  this->_myMax = indicatorF.getMax() + layerSize;
}

// returns true if x is inside the layer
template <typename S>
bool IndicatorLayer3D<S>::operator()(bool output[], const S input[])
{
  output[0] = false;
  S r[3];
  for (int iX =- 1; iX < 2; ++iX) {
    for (int iY =- 1; iY < 2; ++iY) {
      for (int iZ =- 1; iZ < 2; ++iZ) {
        r[0] = input[0] + iX*_layerSize;
        r[1] = input[1] + iY*_layerSize;
        r[2] = input[2] + iZ*_layerSize;
        _indicatorF(output,r);
        if (output[0]) {
          return true;
        }
      }
    }
  }
  return true;
}

template <typename S>
S IndicatorLayer3D<S>::signedDistance( const Vector<S,3>& input )

{

// Rounding: Creates a Layer with a constant thickness --> edges are rounded
  return sdf::rounding(_indicatorF.signedDistance(input), _layerSize);


}

template <typename S>
IndicatorInternal3D<S>::IndicatorInternal3D(IndicatorF3D<S>& indicatorF, S layerSize)
  :  _indicatorF(indicatorF), _layerSize(layerSize)
{
  this->_myMin = indicatorF.getMin() + layerSize;
  this->_myMax = indicatorF.getMax() - layerSize;
}

// returns true if x is inside the layer
template <typename S>
bool IndicatorInternal3D<S>::operator()(bool output[], const S input[])
{
  output[0] = false;
  _indicatorF(output,input);
  if (!output[0]) {
    return true;
  }

  S r[3];
  for (int iX =- 1; iX < 2; ++iX) {
    for (int iY =- 1; iY < 2; ++iY) {
      for (int iZ =- 1; iZ < 2; ++iZ) {
        r[0] = input[0] + iX*_layerSize;
        r[1] = input[1] + iY*_layerSize;
        r[2] = input[2] + iZ*_layerSize;
        _indicatorF(output,r);
        if (!output[0]) {
          return true;
        }
      }
    }
  }
  return true;
}


/// Indicator function for a cylinder
template <typename S>
IndicatorCylinder3D<S>::IndicatorCylinder3D(Vector<S,3> center1,
    Vector<S,3> center2, S radius)
  :  _center1(center1), _center2(center2),
     _ba(_center2 - _center1),
     _baba(_ba[0]*_ba[0] + _ba[1]*_ba[1] + _ba[2]*_ba[2]),
     _radius2(radius*radius), _length(util::sqrt(_baba))
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
  init();
}

/// Indicator function for a cylinder
template <typename S>
IndicatorCylinder3D<S>::IndicatorCylinder3D(Vector<S,3> center1,
    Vector<S,3> normal, S radius, S eps)
  :  _center1(center1-.5*eps*normalize(normal)),
     _center2(_center1 + eps*normalize(normal)),
     _ba(_center2 - _center1),
     _baba(_ba[0]*_ba[0] + _ba[1]*_ba[1] + _ba[2]*_ba[2]),
     _radius2(radius*radius), _length(util::sqrt(_baba))
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
  init();
}

/// Indicator function for a cylinder
template <typename S>
IndicatorCylinder3D<S>::IndicatorCylinder3D(IndicatorCircle3D<S> const& circleF, S eps)
  :  _center1(circleF.getCenter() - .5*eps*circleF.getNormal()),
     _center2(_center1 + eps*circleF.getNormal()),
     _ba(_center2 - _center1),
     _baba(_ba[0]*_ba[0] + _ba[1]*_ba[1] + _ba[2]*_ba[2]),
     _radius2(circleF.getRadius()*circleF.getRadius()), _length(util::sqrt(_baba))
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
  init();
}

// returns true if x is inside the cylinder
template <typename S>
bool IndicatorCylinder3D<S>::operator()(bool output[], const S input[])
{
  Vector<S,3> pa(Vector<S,3>(input) - _center1);

  S X = _I[0]*pa[0] + _I[1]*pa[1] + _I[2]*pa[2];
  S Y = _J[0]*pa[0] + _J[1]*pa[1] + _J[2]*pa[2];
  S Z = _K[0]*pa[0] + _K[1]*pa[1] + _K[2]*pa[2];

  // X^2 + Y^2 <= _radius2
  output[0] = ( Z <= _length && Z >= 0 && X*X + Y*Y <= _radius2 );
  return output[0];
}

template <typename S>
void IndicatorCylinder3D<S>::init()
{
  _K = _ba / _length;

  // _I and _J form an orthonormal base with _K
  if ( util::approxEqual(_center2[1],_center1[1]) && util::approxEqual(_center2[0],_center1[0]) ) {
    if ( util::approxEqual(_center2[2],_center1[2]) ) {
      OstreamManager clout = OstreamManager(std::cout,"IndicatorCylinder3D");
      clout << "Warning: in the cylinder, the two centers have the same coordinates" << std::endl;
      clout << _center1 << std::endl;
      clout << _center2 << std::endl;
    }
    _I = {1,0,0};
    _J = {0,1,0};
  }
  else {
    S normi = util::sqrt (_K[1]*_K[1]+_K[0]*_K[0]);
    _I = {-_K[1]/normi, _K[0]/normi,0};
    _J = {_K[1]*_I[2] - _K[2]*_I[1], _K[2]*_I[0] - _K[0]*_I[2], _K[0]*_I[1] - _K[1]*_I[0]};
  }

  double r = util::sqrt(_radius2);

  S maxx= _center1[0] + util::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + util::max(_K[0]*_length, 0.);
  S minx= _center1[0] - util::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + util::min(_K[0]*_length, 0.);

  S maxy= _center1[1] + util::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + util::max(_K[1]*_length, 0.);
  S miny= _center1[1] - util::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + util::min(_K[1]*_length, 0.);

  S maxz= _center1[2] + util::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + util::max(_K[2]*_length, 0.);
  S minz= _center1[2] - util::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + util::min(_K[2]*_length, 0.);

  this->_myMin = {minx, miny, minz};
  this->_myMax = {maxx, maxy, maxz};
}



template <typename S>
Vector<S,3> const& IndicatorCylinder3D<S>::getCenter1() const
{
  return _center1;
}

template <typename S>
Vector<S,3> const& IndicatorCylinder3D<S>::getCenter2() const
{
  return _center2;
}

template <typename S>
S IndicatorCylinder3D<S>::getRadius() const
{
  return util::sqrt(_radius2);
}

template <typename S>
S IndicatorCylinder3D<S>::signedDistance( const Vector<S,3>& input )
{
  return sdf::cylinder(input, _center1, _ba, _baba, util::sqrt(_radius2));
}



//TODO: Rename to IndicatorCappedCone3D and create a real IndicatorCone3D ??

// cone defined by the centers of the two extremities and the radiuses of the two extremities
// the 2nd radius is optional: if it is not defined, the 2nd center is the vertex of the cone
template <typename S>
IndicatorCone3D<S>::IndicatorCone3D(Vector<S,3> center1, Vector<S,3> center2,
                                    S radius1, S radius2)
  :  _center1(center1), _center2(center2),
     _ba(center2 - center1), _baba(_ba*_ba),
     _radius1(radius1), _radius2(radius2),
     _length(util::sqrt((_center2 - _center1)*(_center2 - _center1)))
{
  // _I,_J,_K is the new base where _K is the axe of the cone
  _K = (_center2 - _center1)/_length;

  // _I and _J form an orthonormal base with _K
  if ( util::approxEqual(_center2[1],_center1[1]) && util::approxEqual(_center2[0],_center1[0]) ) {
    if ( util::approxEqual(_center2[2],_center1[2]) ) {
      OstreamManager clout = OstreamManager(std::cout,"IndicatorCone3D");
      clout << "Warning: in the cone, the two centers have the same coordinates" << std::endl;
      clout << _center1 << std::endl;
      clout << _center2 << std::endl;
    }
    _I = {1,0,0};
    _J = {0,1,0};
  }
  else {
    S normi = util::sqrt(_K[1]*_K[1] + _K[0]*_K[0]);
    _I = {-_K[1]/normi, _K[0]/normi,0};
    _J = {_K[1]*_I[2] - _K[2]*_I[1], _K[2]*_I[0] - _K[0]*_I[2], _K[0]*_I[1] - _K[1]*_I[0]};
  }

  S maxx= _center1[0] + util::max( util::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                                   util::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);
  S minx= _center1[0] + util::min(-util::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                                  -util::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);

  S maxy= _center1[1] + util::max( util::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                                   util::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);
  S miny= _center1[1] + util::min(-util::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                                  -util::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);

  S maxz= _center1[2] + util::max( util::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                                   util::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);
  S minz= _center1[2] + util::min(-util::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                                  -util::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);

  this->_myMin = {minx, miny, minz};
  this->_myMax = {maxx, maxy, maxz};
}

// returns true if x is inside the cone(Vector<S,3> center1, Vector<S,3> center2, S radius1
template <typename S>
bool IndicatorCone3D<S>::operator()(bool output[], const S input[])
{
  // radius: the radius of the cone at the point x
  Vector<S,3> pa(Vector<S,3>(input) - _center1);
  S X = _I[0]*pa[0] + _I[1]*pa[1] + _I[2]*pa[2];
  S Y = _J[0]*pa[0] + _J[1]*pa[1] + _J[2]*pa[2];
  S Z = _K[0]*pa[0] + _K[1]*pa[1] + _K[2]*pa[2];
  S radius = _radius1 + (_radius2 - _radius1)*Z / _length;

  output[0] = ( Z <= _length && Z >= 0 && X*X + Y*Y <= radius*radius );
  return true;
}

template <typename S>
Vector<S,3> const& IndicatorCone3D<S>::getCenter1() const
{
  return _center1;
}

template <typename S>
Vector<S,3> const& IndicatorCone3D<S>::getCenter2() const
{
  return _center2;
}

template <typename S>
S IndicatorCone3D<S>::getRadius1() const
{
  return _radius1;
}

template <typename S>
S IndicatorCone3D<S>::getRadius2() const
{
  return _radius2;
}

template <typename S>
S IndicatorCone3D<S>::signedDistance( const Vector<S,3>& input )
{
  return sdf::cone(input, _center1, _ba, _baba, _radius1, _radius2);
}

template <typename S>
IndicatorEllipsoid3D<S>::IndicatorEllipsoid3D(Vector<S,3> center, Vector<S,3> radius)
  :  _center(center), _radius(radius)
{
#ifdef OLB_DEBUG
  OstreamManager clout(std::cout, "IndicatorEllipsoid3D");
  clout << "WARNING: The indicator doesn't provide an exact signed distance function." << std::endl;
#endif

  assert(_radius[0]>0 && _radius[1]>0 && _radius[2]>0);

  this->_myMin = center - radius;
  this->_myMax = center + radius;
}

template <typename S>
Vector<S,3> const& IndicatorEllipsoid3D<S>::getRadius() const
{
  return _radius;
}

template <typename S>
Vector<S,3> const& IndicatorEllipsoid3D<S>::getCenter() const
{
  return _center;
}

template <typename S>
S IndicatorEllipsoid3D<S>::signedDistance( const Vector<S,3>& input )
{
  Vector<S,3> p = input - _center;
  return sdf::ellipsoid(p, _radius);
}

template <typename S>
IndicatorSuperEllipsoid3D<S>::IndicatorSuperEllipsoid3D(Vector<S,3> center, S xHalfAxis, S yHalfAxis, S zHalfAxis, S exponent1, S exponent2)
  : _center(center), _xHalfAxis(xHalfAxis), _yHalfAxis(yHalfAxis), _zHalfAxis(zHalfAxis), _exp1(exponent1), _exp2(exponent2)
{
  assert(_xHalfAxis>0 && _yHalfAxis>0 && _zHalfAxis>0);

  S max_axis = util::max( xHalfAxis, util::max(yHalfAxis, zHalfAxis) );
  this->_myMin = {
    center[0] - util::sqrt(2.)*max_axis,
    center[1] - util::sqrt(2.)*max_axis,
    center[2] - util::sqrt(2.)*max_axis
  };
  this->_myMax = {
    center[0] + util::sqrt(2.)*max_axis,
    center[1] + util::sqrt(2.)*max_axis,
    center[2] + util::sqrt(2.)*max_axis
  };
}

template <typename S>
Vector<S,3> const& IndicatorSuperEllipsoid3D<S>::getCenter() const
{
  return _center;
}

template <typename S>
S IndicatorSuperEllipsoid3D<S>::getXHalfAxis() const
{
  return _xHalfAxis;
}

template <typename S>
S IndicatorSuperEllipsoid3D<S>::getYHalfAxis() const
{
  return _yHalfAxis;
}

template <typename S>
S IndicatorSuperEllipsoid3D<S>::getZHalfAxis() const
{
  return _zHalfAxis;
}

template <typename S>
S IndicatorSuperEllipsoid3D<S>::getExponent1() const
{
  return _exp1;
}

template <typename S>
S IndicatorSuperEllipsoid3D<S>::getExponent2() const
{
  return _exp2;
}

template <typename S>
bool IndicatorSuperEllipsoid3D<S>::operator()(bool output[], const S input[])
{
  S a = util::pow ( util::abs( (input[0] - _center[0]) / _xHalfAxis ), _exp1 );
  S b = util::pow ( util::abs( (input[1] - _center[1]) / _yHalfAxis ), _exp1 );
  S c = util::pow ( util::abs( (input[2] - _center[2]) / _zHalfAxis ), _exp2 );
  S ab = util::pow( a+b, _exp2/_exp1 );

  if ( (ab+c) <= 1. ) {
    output[0] = 1.;
    return true;
  }

  output[0] = 0.;
  return false;
}


// Warning : the cuboid is only defined parallel to the plans x=0, y=0 and z=0 !!!
// For other cuboids, please use the parallelepiped version
template <typename S>
IndicatorCuboid3D<S>::IndicatorCuboid3D(Vector<S,3> extend, Vector<S,3> origin)
  : _center(origin+.5*extend),_xLength(extend[0]), _yLength(extend[1]), _zLength(extend[2])
{
  assert(_xLength>0 && _yLength>0 && _zLength>0);

  this->_myMin = origin;
  this->_myMax = origin + extend;
}

template <typename S>
IndicatorCuboid3D<S>::IndicatorCuboid3D(S xLength, S yLength, S zLength, Vector<S,3> center)
  : _center(center), _xLength(xLength), _yLength(yLength), _zLength(zLength)
{
  assert(_xLength>0 && _yLength>0 && _zLength>0);

  this->_myMin = {_center[0] - S{0.5}*_xLength, _center[1] - S{0.5}*_yLength, _center[2] - S{0.5}*_zLength};
  this->_myMax = {_center[0] + S{0.5}*_xLength, _center[1] + S{0.5}*_yLength, _center[2] + S{0.5}*_zLength};
}

template <typename S>
Vector<S,3> const& IndicatorCuboid3D<S>::getCenter() const
{
  return _center;
}

template <typename S>
S const IndicatorCuboid3D<S>::getxLength() const
{
  return _xLength;
}

template <typename S>
S const IndicatorCuboid3D<S>::getyLength() const
{
  return _yLength;
}

template <typename S>
S const IndicatorCuboid3D<S>::getzLength() const
{
  return _zLength;
}

template <typename S>
bool IndicatorCuboid3D<S>::operator()(bool output[], const S input[])
{
  // returns true if x is inside the cuboid
  Vector<S,3> q = distanceXYZ(input);
  output[0] = ( (q[0] < std::numeric_limits<S>::epsilon())
                && (q[1] < std::numeric_limits<S>::epsilon())
                && (q[2] < std::numeric_limits<S>::epsilon()) );
  return output[0];
}

template <typename S>
Vector<S,3> IndicatorCuboid3D<S>::distanceXYZ( Vector<S,3> input )
{
  return  { S(util::abs(S(input[0] - _center[0])) - 0.5 * _xLength),
            S(util::abs(S(input[1] - _center[1])) - 0.5 * _yLength),
            S(util::abs(S(input[2] - _center[2])) - 0.5 * _zLength)
          };
}

template <typename S>
S IndicatorCuboid3D<S>::signedDistance( const Vector<S,3>& input )
{
  return sdf::box(input - _center, Vector<S,3>(.5*_xLength, .5*_yLength, .5*_zLength));
}



template <typename S>
IndicatorCuboidRotate3D<S>::IndicatorCuboidRotate3D(Vector<S,3> extend, Vector<S,3> origin, S theta, int plane, Vector<S,3> centerRotation)
  : IndicatorCuboid3D<S>(extend, origin), _theta(theta), _plane(plane), _centerRotation(centerRotation)
{
  assert(plane==0 || plane==1 || plane==2);
}

template <typename S>
IndicatorCuboidRotate3D<S>::IndicatorCuboidRotate3D(S xLength, S yLength, S zLength, Vector<S,3> origin, S theta, int plane, Vector<S,3> centerRotation)
  : IndicatorCuboid3D<S>(xLength, yLength, zLength, origin), _theta(theta), _plane(plane), _centerRotation(centerRotation)
{
  assert(plane==0 || plane==1 || plane==2);
}

template <typename S>
void IndicatorCuboidRotate3D<S>::transformInput(const S input[3], S newInput[3])
{
  //initialize for _plane == 2
  int i=0;
  int j=1;
  if (_plane == 1) { // rotation around y axis
    i=0;
    j=2;
  }
  else if (_plane == 0) {    // rotation around x axis
    i=1;
    j=2;
  }
  // translation to _centerRotation
  S x = input[i] - _centerRotation[i];
  S y = input[j] - _centerRotation[j];
  // rotation of _theta in rad
  S xx = x*util::cos(_theta) - y*util::sin(_theta);
  S yy = x*util::sin(_theta) + y*util::cos(_theta);
  // change back to standard coordinate system
  x = xx + _centerRotation[i];
  y = yy + _centerRotation[j];
  newInput[_plane] = input[_plane];
  newInput[i] = x;
  newInput[j] = y;
}

// do transformation to axis aligned cuboid, then call operator() of basic cuboid.
template <typename S>
bool IndicatorCuboidRotate3D<S>::operator()(bool output[], const S input[])
{
  S newInput[3];
  transformInput(input, newInput);
  IndicatorCuboid3D<S>::operator()(output, newInput);
  return output[0];
}

template <typename S>
S IndicatorCuboidRotate3D<S>::signedDistance( const Vector<S,3>& input )
{
  S newInput[3], tmp[3] = {input[0], input[1], input[2]};
  transformInput(tmp, newInput);
  return IndicatorCuboid3D<S>::signedDistance(Vector<S,3>(newInput));
}


////// creator functions /////
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCircle3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCircle3D");

  Vector<S,3> center;
  Vector<S,3> normal;
  S radius = 1;

  //  params.setWarningsOn(false);
  params.setWarningsOn(true);

  std::stringstream xmlCenter1( params.getAttribute("center") );
  xmlCenter1 >> center[0] >> center[1] >> center[2];
  std::stringstream xmlCenter2( params.getAttribute("normal") );
  xmlCenter2 >> normal[0] >> normal[1] >> normal[2];
  std::stringstream xmlRadius( params.getAttribute("radius") );
  xmlRadius >> radius;

  return std::make_shared<IndicatorCircle3D<S>>(center, normal, radius);
}

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorSphere3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorSphere3D");

  Vector<S,3> center;
  S radius = 1;

  //  params.setWarningsOn(false);
  params.setWarningsOn(true);

  std::stringstream xmlCenter1( params.getAttribute("center") );
  xmlCenter1 >> center[0] >> center[1] >> center[2];
  std::stringstream xmlRadius( params.getAttribute("radius") );
  xmlRadius >> radius;

  return std::make_shared<IndicatorSphere3D<S>>(center, radius);
}

// creator function for a cylinder3d
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCylinder3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCylinder3D");

  Vector<S,3> center1;
  Vector<S,3> center2(S(1),S(1),S(1));
  S radius = 1;

  //  params.setWarningsOn(false);
  //  params.setWarningsOn(true);

  std::stringstream xmlCenter1( (params).getAttribute("center1") );
  xmlCenter1 >> center1[0] >> center1[1] >> center1[2];
  std::stringstream xmlCenter2( (params).getAttribute("center2") );
  xmlCenter2 >> center2[0] >> center2[1] >> center2[2];
  std::stringstream xmlRadius( (params).getAttribute("radius") );
  xmlRadius >> radius;

  /// for debugging purpose
//  print(center1, "center1: ");
//  print(center2, "center2: ");
//  print(radius, "radius: ");

  return std::make_shared<IndicatorCylinder3D<S>>(center1, center2, radius);
}

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCone3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCone3D");

  Vector<S,3> center1;
  Vector<S,3> center2(S(1), S(1), S(1));
  S radius1 = S(0);
  S radius2 = S(1);

  //  params.setWarningsOn(false);
  params.setWarningsOn(true);

  std::stringstream xmlCenter1( params.getAttribute("center1") );
  xmlCenter1 >> center1[0] >> center1[1] >> center1[2];
  std::stringstream xmlCenter2( params.getAttribute("center2") );
  xmlCenter2 >> center2[0] >> center2[1] >> center2[2];
  std::stringstream xmlRadius1( params.getAttribute("radius1") );
  xmlRadius1 >> radius1;
  std::stringstream xmlRadius2( params.getAttribute("radius2") );
  xmlRadius2 >> radius2;

  return std::make_shared<IndicatorCone3D<S>>(center1, center2, radius1, radius2);
}

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCuboid3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorCuboid3D");

  Vector<S,3> origin;
  Vector<S,3> extend(S(1),S(1),S(1));

  //  params.setWarningsOn(false);
  params.setWarningsOn(true);

  std::stringstream xmlOrigin( params.getAttribute("origin") );
  xmlOrigin >> origin[0] >> origin[1] >> origin[2];
  std::stringstream xmlExtend( params.getAttribute("extend") );
  xmlExtend >> extend[0] >> extend[1] >> extend[2];

  return std::make_shared<IndicatorCuboid3D<S>>(extend, origin);
}

// Create Union with XML - file
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorUnion3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorUnion3D");

  //  clout << "xml position: " << params.getName() << std::endl;
  //  params.print(2);

  std::shared_ptr<IndicatorF3D<S>> output = createIndicatorF3D<S>(**params.begin());
  for (auto it = params.begin()+1; it != params.end(); ++it) {
    output = output + createIndicatorF3D<S>(**it);
  }
  return output;
}

// Create Without with XML - file
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorWithout3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorWithout3D");

  //  clout << "xml position: " << params.getName() << std::endl;
  //  params.print(2);

  std::shared_ptr<IndicatorF3D<S>> output = createIndicatorF3D<S>(**params.begin());
  for (auto it = params.begin()+1; it != params.end(); ++it) {
    output = output - createIndicatorF3D<S>(**it);
  }
  return output;
}

// Create Intersection with XML - file
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorIntersection3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorIntersection3D");

  //  clout << "xml position: " << params.getName() << std::endl;
  //  params.print(2);

  std::shared_ptr<IndicatorF3D<S>> output = createIndicatorF3D<S>(**params.begin());
  for (auto it = params.begin()+1; it != params.end(); ++it) {
    output = output * createIndicatorF3D<S>(**it);
  }
  return output;
}

// Create Geometry
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorF3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorF3D");

  //  clout << "XML element: "<< params.getName() << std::endl;
  //  params.print(2);

  std::string actualName = params.getName();
  if ( actualName == "IndicatorCircle3D" ) {
    return createIndicatorCircle3D<S>(params);
  }
  else if ( actualName == "IndicatorSphere3D" ) {
    return createIndicatorSphere3D<S>(params);
  }
  else if ( actualName == "IndicatorCylinder3D" ) {
    return createIndicatorCylinder3D<S>(params);
  }
  else if ( actualName == "IndicatorCone3D" ) {
    return createIndicatorCone3D<S>(params);
  }
  else if ( actualName == "IndicatorCuboid3D" ) {
    return createIndicatorCuboid3D<S>(params);
  }
  else if ( actualName == "IndicatorUnion3D" ) {
    return createIndicatorUnion3D<S>(params);
  }
  else if ( actualName == "IndicatorWithout3D" ) {
    return createIndicatorWithout3D<S>(params);
  }
  else if ( actualName == "IndicatorIntersection3D" ) {
    return createIndicatorIntersection3D<S>(params);
  }
  else {
    auto firstChild = params.begin(); // get iterator of childTree
    return createIndicatorF3D<S>( **firstChild );
  }
}


} // namespace olb

#endif