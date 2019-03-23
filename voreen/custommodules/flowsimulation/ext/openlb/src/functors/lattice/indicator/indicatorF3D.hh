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
#include <cmath>
#include <cassert>
#include <sstream>

#include "indicatorF3D.h"
#include "indicCalcF3D.h"
#include "utilities/vectorHelpers.h"
#include "functors/analytical/interpolationF3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace olb {



template <typename S>
IndicatorCircle3D<S>::IndicatorCircle3D(Vector<S,3> center, Vector<S,3> normal, S radius)
  :  _center(center), _normal(normal), _radius2(radius*radius)
{
  _normal.normalize();
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
  _normal.normalize();
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
  return sqrt(_radius2);
}


template <typename S>
IndicatorCircle3D<S>* createIndicatorCircle3D(XMLreader const& params, bool verbose)
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

  return new IndicatorCircle3D<S>(center, normal, radius);
}

template <typename S>
IndicatorSphere3D<S>::IndicatorSphere3D(Vector<S,3> center, S radius)
  :  _center(center), _radius2(radius*radius)
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
  _radius2 = sphere._radius2;
}

// returns true if x is inside the sphere
template <typename S>
bool IndicatorSphere3D<S>::operator()(bool output[], const S input[])
{
  output[0] = (  (_center[0] - input[0]) * (_center[0]-input[0])
                 +(_center[1] - input[1]) * (_center[1]-input[1])
                 +(_center[2] - input[2]) * (_center[2]-input[2]) <= _radius2 );
  return true;
}

template <typename S>
bool IndicatorSphere3D<S>::distance(S& distance, const Vector<S,3>& origin,
                                    const Vector<S,3>& direction, int iC)
{
  // computes pos. distance by solving quadratic equation by a-b-c-formula
  S a = direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2];

  // returns 0 if point is at the boundary of the sphere
  if ( util::nearZero(a-_radius2) ) {
    distance = S();
    return true;
  }
  // norm of direction
  a = std::sqrt(a);

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

  S x1 = (- b + std::sqrt(d)) *0.5;
  S x2 = (- b - std::sqrt(d)) *0.5;

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
    distance = std::min(x1,x2);
    return true;
  }

  return false;
}

template <typename S>
IndicatorSphere3D<S>* createIndicatorSphere3D(XMLreader const& params, bool verbose)
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

  return new IndicatorSphere3D<S>(center, radius);
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


/// Indicator function for a cylinder
template <typename S>
IndicatorCylinder3D<S>::IndicatorCylinder3D(Vector<S,3> center1,
    Vector<S,3> center2, S radius)
  :  _center1(center1), _center2(center2), _radius2(radius*radius)
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
  init();
}

/// Indicator function for a cylinder
template <typename S>
IndicatorCylinder3D<S>::IndicatorCylinder3D(Vector<S,3> center1,
    Vector<S,3> normal, S radius, S eps)
  :  _center1(center1), _center2(center1), _radius2(radius*radius)
{
  normal.normalize();
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
  _center1 -= .5*eps*normal;
  _center2 = _center1 + eps*normal;

  init();
}

/// Indicator function for a cylinder
template <typename S>
IndicatorCylinder3D<S>::IndicatorCylinder3D(IndicatorCircle3D<S> const& circleF, S eps)
  :  _center1(circleF.getCenter()), _center2(circleF.getCenter()),
     _radius2(circleF.getRadius()*circleF.getRadius())
{
  // cylinder defined by the centers of the two extremities and the radius
  // _I,_J,_K is the new base where _K is the axe of the cylinder0
  _center1 -= .5*eps*circleF.getNormal();
  _center2 = _center1 + eps*circleF.getNormal();

  init();
}

// returns true if x is inside the cylinder
template <typename S>
bool IndicatorCylinder3D<S>::operator()(bool output[], const S input[])
{
  S X = _I[0]*(input[0]-_center1[0]) + _I[1]*(input[1]-_center1[1]) + _I[2]*(input[2]-_center1[2]);
  S Y = _J[0]*(input[0]-_center1[0]) + _J[1]*(input[1]-_center1[1]) + _J[2]*(input[2]-_center1[2]);
  S Z = _K[0]*(input[0]-_center1[0]) + _K[1]*(input[1]-_center1[1]) + _K[2]*(input[2]-_center1[2]);

  // X^2 + Y^2 <= _radius2
  output[0] = ( Z <= _length && Z >= 0 && X*X + Y*Y <= _radius2 );
  return output[0];
}

template <typename S>
void IndicatorCylinder3D<S>::init()
{
  _length = std::sqrt( (_center2[0]-_center1[0]) * (_center2[0]-_center1[0])
                  +(_center2[1]-_center1[1]) * (_center2[1]-_center1[1])
                  +(_center2[2]-_center1[2]) * (_center2[2]-_center1[2]) );

  // _K = centre2 - centre1 (normalized)
  _K = {(_center2[0]-_center1[0]) / _length, (_center2[1]-_center1[1]) / _length,
        (_center2[2]-_center1[2]) / _length
       };

  // _I and _J form an orthonormal base with _K
  if ( util::nearZero(_center2[1]-_center1[1]) && util::nearZero(_center2[0]-_center1[0]) ) {
    if ( util::nearZero(_center2[2]-_center1[2]) ) {
      std::cout << "Warning: in the cylinder, the two centers have the same coordinates" << std::endl;
      util::print(_center1);
      util::print(_center2);
    }
    _I = {1,0,0};
    _J = {0,1,0};
  } else {
    S normi = std::sqrt (_K[1]*_K[1]+_K[0]*_K[0]);
    _I = {-_K[1]/normi, _K[0]/normi,0};
    _J = {_K[1]*_I[2] - _K[2]*_I[1], _K[2]*_I[0] - _K[0]*_I[2], _K[0]*_I[1] - _K[1]*_I[0]};
  }

  double r = std::sqrt(_radius2);
  S maxx, maxy, maxz;
  S minx, miny, minz;

  maxx= _center1[0] + std::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + std::max(_K[0]*_length, 0.);
  minx= _center1[0] - std::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*r + std::min(_K[0]*_length, 0.);

  maxy= _center1[1] + std::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + std::max(_K[1]*_length, 0.);
  miny= _center1[1] - std::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*r + std::min(_K[1]*_length, 0.);

  maxz= _center1[2] + std::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + std::max(_K[2]*_length, 0.);
  minz= _center1[2] - std::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*r + std::min(_K[2]*_length, 0.);

  this->_myMin = {minx, miny, minz};
  this->_myMax = {maxx, maxy, maxz};
}

// creator function for a cylinder3d
template <typename S>
IndicatorCylinder3D<S>* createIndicatorCylinder3D(XMLreader const& params, bool verbose)
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

  return new IndicatorCylinder3D<S>(center1, center2, radius);
}

// cone defined by the centers of the two extremities and the radiuses of the two extremities
// the 2nd radius is optional: if it is not defined, the 2nd center is the vertex of the cone
template <typename S>
IndicatorCone3D<S>::IndicatorCone3D(Vector<S,3> center1, Vector<S,3> center2,
                                    S radius1, S radius2)
  :  _center1(center1), _center2(center2),
     _radius1(radius1), _radius2(radius2)
{
  // _I,_J,_K is the new base where _K is the axe of the cone
  // _K = centre2 - centre1 (normalized)
  _length = std::sqrt( (_center2[0]-_center1[0]) * (_center2[0]-_center1[0])
                  +(_center2[1]-_center1[1]) * (_center2[1]-_center1[1])
                  +(_center2[2]-_center1[2]) * (_center2[2]-_center1[2]) );
  // _K = centre2 - centre1 (normalized)
  _K = {(_center2[0]-_center1[0]) / _length, (_center2[1]-_center1[1]) / _length,
        (_center2[2]-_center1[2]) / _length
       };

  // _I and _J form an orthonormal base with _K
  if ( util::nearZero(_center2[1]-_center1[1]) && util::nearZero(_center2[0]-_center1[0]) ) {
    if ( util::nearZero(_center2[2]-_center1[2]) ) {
      std::cout << "Warning: in the cone, the two center have the same coordinates" << std::endl;
      util::print(_center1);
      util::print(_center2);
    }
    _I = {1,0,0};
    _J = {0,1,0};
  } else {
    S normi = std::sqrt(_K[1]*_K[1] + _K[0]*_K[0]);
    _I = {-_K[1]/normi, _K[0]/normi,0};
    _J = {_K[1]*_I[2] - _K[2]*_I[1], _K[2]*_I[0] - _K[0]*_I[2], _K[0]*_I[1] - _K[1]*_I[0]};
  }

  S maxx, maxy, maxz;
  S minx, miny, minz;

  maxx= _center1[0] + std::max( std::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                                std::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);
  minx= _center1[0] + std::min(-std::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius1,
                               -std::sqrt(_I[0]*_I[0]+_J[0]*_J[0])*_radius2 + _K[0]*_length);

  maxy= _center1[1] + std::max( std::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                                std::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);
  miny= _center1[1] + std::min(-std::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius1,
                               -std::sqrt(_I[1]*_I[1]+_J[1]*_J[1])*_radius2 + _K[1]*_length);

  maxz= _center1[2] + std::max( std::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                                std::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);
  minz= _center1[2] + std::min(-std::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius1,
                               -std::sqrt(_I[2]*_I[2]+_J[2]*_J[2])*_radius2 + _K[2]*_length);


  this->_myMin = {minx, miny, minz};
  this->_myMax = {maxx, maxy, maxz};
}

// returns true if x is inside the cone(Vector<S,3> center1, Vector<S,3> center2, S radius1
template <typename S>
bool IndicatorCone3D<S>::operator()(bool output[], const S input[])
{
  // radius: the radius of the cone at the point x
  S X = _I[0]*(input[0]-_center1[0]) + _I[1]*(input[1]-_center1[1]) + _I[2]*(input[2]-_center1[2]);
  S Y = _J[0]*(input[0]-_center1[0]) + _J[1]*(input[1]-_center1[1]) + _J[2]*(input[2]-_center1[2]);
  S Z = _K[0]*(input[0]-_center1[0]) + _K[1]*(input[1]-_center1[1]) + _K[2]*(input[2]-_center1[2]);
  S radius = _radius1 + (_radius2 - _radius1)*Z / _length;

  output[0] = ( Z <= _length && Z >= 0 && X*X + Y*Y <= radius*radius );
  return true;
}

template <typename S>
IndicatorCone3D<S>* createIndicatorCone3D(XMLreader const& params, bool verbose)
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

  return new IndicatorCone3D<S>(center1, center2, radius1, radius2);
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

  this->_myMin = {_center[0] - _xLength/2., _center[1] - _yLength/2., _center[2] - _zLength/2.};
  this->_myMax = {_center[0] + _xLength/2., _center[1] + _yLength/2., _center[2] + _zLength/2.};
}


template <typename S>
bool IndicatorCuboid3D<S>::operator()(bool output[], const S input[])
{
  // returns true if x is inside the cuboid
  output[0] = ( (fabs(_center[0] - input[0]) < _xLength/2. || nearZero(fabs(_center[0] - input[0]) - _xLength/2.))
                && (fabs(_center[1] - input[1]) < _yLength/2. || nearZero(fabs(_center[1] - input[1]) - _yLength/2.))
                && (fabs(_center[2] - input[2]) < _zLength/2. || nearZero(fabs(_center[2] - input[2]) - _zLength/2.)) );
  return output[0];
}


template <typename S>
IndicatorCuboid3D<S>* createIndicatorCuboid3D(XMLreader const& params, bool verbose)
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

  return new IndicatorCuboid3D<S>(extend, origin);
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

// do transformation to axis aligned cuboid, then call operator() of basic cuboid.
template <typename S>
bool IndicatorCuboidRotate3D<S>::operator()(bool output[], const S input[])
{
  //initialize for _plane == 2
  int i=0;
  int j=1;
   if (_plane == 1) { // rotation around y axis
    i=0;
    j=2;
  } else if (_plane == 0) {  // rotation around x axis
    i=1;
    j=2;
  }
  // translation to _centerRotation
  S x = input[i] - _centerRotation[i];
  S y = input[j] - _centerRotation[j];
  // rotation of _theta in rad
  S xx = x*cos(_theta) - y*sin(_theta);
  S yy = x*sin(_theta) + y*cos(_theta);
  // change back to standard coordinate system
  x = xx + _centerRotation[i];
  y = yy + _centerRotation[j];
  S newInput[3] = {x, y, input[_plane]};
  IndicatorCuboid3D<S>::operator()(output, newInput);
  return output[0];
}

// Create Union with XML - file
template <typename S>
IndicatorF3D<S>* createIndicatorUnion3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorUnion3D");

  //  clout << "xml position: " << params.getName() << std::endl;
  //  params.print(2);

  std::vector< IndicatorF3D<S>* > vecPlus;
  for ( std::vector<XMLreader*>::const_iterator it = params.begin(); it != params.end(); ++it ) {
    vecPlus.push_back( createIndicatorF3D<S>(**it) );
  }

  IndicatorF3D<S>* output = vecPlus.front();
  typename std::vector< IndicatorF3D<S>* >::iterator it_vector;
  for (it_vector = vecPlus.begin()+1; it_vector != vecPlus.end(); ++it_vector) {
    output = new IndicPlus3D<S>( *output, **it_vector );
  }
  return output;
}

// Create Without with XML - file
template <typename S>
IndicatorF3D<S>* createIndicatorWithout3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorWithout3D");

  //  clout << "xml position: " << params.getName() << std::endl;
  //  params.print(2);

  std::vector< IndicatorF3D<S>* > vecMinus;
  for ( std::vector<XMLreader*>::const_iterator it = params.begin(); it != params.end(); ++it ) {
    vecMinus.push_back( createIndicatorF3D<S>(**it) );
  }

  IndicatorF3D<S>* output = vecMinus.front();
  typename std::vector< IndicatorF3D<S>* >::iterator it_vector;
  for (it_vector = vecMinus.begin()+1; it_vector != vecMinus.end(); ++it_vector) {
    output = new IndicMinus3D<S>( *output, **it_vector );
  }
  return output;
}

// Create Intersection with XML - file
template <typename S>
IndicatorF3D<S>* createIndicatorIntersection3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorIntersection3D");

  //  clout << "xml position: " << params.getName() << std::endl;
  //  params.print(2);

  std::vector< IndicatorF3D<S>* > vecIntersection;
  for ( std::vector<XMLreader*>::const_iterator it = params.begin(); it != params.end(); ++it ) {
    vecIntersection.push_back( createIndicatorF3D<S>(**it) );
  }

  IndicatorF3D<S>* output = vecIntersection.front();
  typename std::vector< IndicatorF3D<S>* >::iterator it_indi;
  for (it_indi = vecIntersection.begin()+1; it_indi != vecIntersection.end(); ++it_indi) {
    output = new IndicMultiplication3D<S>( *output, **it_indi );
  }
  return output;
}

// Create Geometry
template <typename S>
IndicatorF3D<S>* createIndicatorF3D(XMLreader const& params, bool verbose)
{
  OstreamManager clout(std::cout,"createIndicatorF3D");
  IndicatorF3D<S>* output = nullptr;

  //  clout << "XML element: "<< params.getName() << std::endl;
  //  params.print(2);

  std::string actualName = params.getName();
  if ( actualName == "IndicatorCircle3D" ) {
    output = createIndicatorCircle3D<S>(params);
    return output;
  } else if ( actualName == "IndicatorSphere3D" ) {
    output = createIndicatorSphere3D<S>(params);
    return output;
  } else if ( actualName == "IndicatorCylinder3D" ) {
    output = createIndicatorCylinder3D<S>(params);
    return output;
  } else if ( actualName == "IndicatorCone3D" ) {
    output = createIndicatorCone3D<S>(params);
    return output;
  } else if ( actualName == "IndicatorCuboid3D" ) {
    output = createIndicatorCuboid3D<S>(params);
    return output;
  } else if ( actualName == "IndicatorUnion3D" ) {
    output = createIndicatorUnion3D<S>(params);
    return output;
  } else if ( actualName == "IndicatorWithout3D" ) {
    output = createIndicatorWithout3D<S>(params);
    return output;
  } else if ( actualName == "IndicatorIntersection3D" ) {
    output = createIndicatorIntersection3D<S>(params);
    return output;
  }

  // descend a level in the xml tree
  std::vector<XMLreader*>::const_iterator it = params.begin();
  // \todo This is how we should iterate through all children and make a union
//  for(XMLreader* child : params) {
//    clout << "iterator to xml-child: " << child->getName() << std::endl;
//  }
  std::string childName = (**it).getName();
  //  clout << "iterator to xml-child: " << (**it).getName() << std::endl;
  if ( childName == "IndicatorUnion3D" ) {
    output = createIndicatorUnion3D<S>(**it);
    return output;
  } else if ( childName == "IndicatorWithout3D" ) {
    output = createIndicatorWithout3D<S>(**it);
    return output;
  } else if ( childName == "IndicatorIntersection3D" ) {
    output = createIndicatorIntersection3D<S>(**it);
    return output;
    // simple building geometry without using union,intersection
    // i dont like this
  } else if ( childName == "IndicatorCircle3D" ) {
    output = createIndicatorCircle3D<S>(**it);
    return output;
  } else if ( childName == "IndicatorSphere3D" ) {
    output = createIndicatorSphere3D<S>(**it);
    return output;
  } else if ( childName == "IndicatorCylinder3D" ) {
    output = createIndicatorCylinder3D<S>(**it);
    return output;
  } else if ( childName == "IndicatorCone3D" ) {
    output = createIndicatorCone3D<S>(**it);
    return output;
  } else if ( childName == "IndicatorCuboid3D" ) {
    output = createIndicatorCuboid3D<S>(**it);
    return output;
  } else if ( childName == "IndicatorUnion3D" ) {
    output = createIndicatorUnion3D<S>(**it);
    return output;
  } else if ( childName == "IndicatorWithout3D" ) {
    output = createIndicatorWithout3D<S>(**it);
    return output;
  } else if ( childName == "IndicatorIntersection3D" ) {
    output = createIndicatorIntersection3D<S>(**it);
    return output;
  } else {
    clout << "createIndicatorF3D throughs error" << std::endl;
    exit(-1);
  }

  //  return output;
}


} // namespace olb

#endif
