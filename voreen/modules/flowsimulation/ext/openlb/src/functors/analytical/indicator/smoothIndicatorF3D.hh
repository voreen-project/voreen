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

#ifndef SMOOTH_INDICATOR_F_3D_HH
#define SMOOTH_INDICATOR_F_3D_HH

#include <vector>
#include "utilities/omath.h"
#include <sstream>

#include "smoothIndicatorF3D.h"
#include "smoothIndicatorBaseF3D.h"
#include "smoothIndicatorCalcF3D.h"
#include "utilities/vectorHelpers.h"
#include "dynamics/descriptorAlias.h"
#include "functors/analytical/interpolationF3D.h"
#include "functors/lattice/reductionF3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI2
#define M_PI2 1.57079632679489661923
#endif

namespace olb {

//Constructor: SmoothIndicatorCuboid3D
template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCuboid3D<T,S,PARTICLE>::SmoothIndicatorCuboid3D(IndicatorCuboid3D<S>& ind,
    S epsilon, Vector<S,3> theta)
  :SmoothIndicatorCuboid3D(ind.getxLength(), ind.getyLength(), ind.getzLength(), ind.getCenter(), epsilon, theta)
{ }

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCuboid3D<T,S,PARTICLE>::SmoothIndicatorCuboid3D(S xLength, S yLength, S zLength,
    Vector<S,3> center, S epsilon, Vector<S,3> theta)
  :_ind(xLength, yLength, zLength, center)
{
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos = _ind.getCenter();
    this->_theta = {
      theta[0] * (M_PI/180.),
      theta[1] * (M_PI/180.),
      theta[2] * (M_PI/180.)
    };
  }

  this->_circumRadius = .5*(util::sqrt(xLength*xLength+yLength*yLength+zLength*zLength))+0.5*epsilon;
  if constexpr (!PARTICLE) {
    this->_myMin = {
      this->_pos[0] - this->getCircumRadius(),
      this->_pos[1] - this->getCircumRadius(),
      this->_pos[2] - this->getCircumRadius()
    };

    this->_myMax = {
      this->_pos[0] + this->getCircumRadius(),
      this->_pos[1] + this->getCircumRadius(),
      this->_pos[2] + this->getCircumRadius()
    };
    this->init(this->_theta);
  }
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorCuboid3D<T, S, PARTICLE>::getVolume( )
{
  return _ind.getxLength()*_ind.getyLength()*_ind.getzLength();
}

template <typename T, typename S, bool PARTICLE>
Vector<S,4> SmoothIndicatorCuboid3D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  T const xLength = _ind.getxLength();
  T const yLength = _ind.getyLength();
  T const zLength = _ind.getzLength();
  T const mass = getVolume()*density;
  T const xLength2 = xLength*xLength;
  T const yLength2 = yLength*yLength;
  T const zLength2 = zLength*zLength;
  Vector<S,3> mofi;
  mofi[0] = 1./12.*mass*(yLength2+zLength2);
  mofi[1] = 1./12.*mass*(xLength2+zLength2);
  mofi[2] = 1./12.*mass*(yLength2+xLength2);
  return Vector<S,4>(mofi[0], mofi[1], mofi[2], mass);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,3> SmoothIndicatorCuboid3D<T, S, PARTICLE>::surfaceNormal( const Vector<S,3>& pos, const S meshSize )
{
  return _ind.surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorCuboid3D<T, S, PARTICLE>::signedDistance( const PhysR<S,3> input )
{
  Vector<S,3> p;
  if constexpr(!PARTICLE) {
    // counter-clockwise rotation by _theta=-theta around the current position of the center of mass & translation
    p = body::motion::rotation<3,S,true>::execute(input, this->_rotMat, this->getPos());
  }
  else {
    p = input;
  }
  return _ind.signedDistance(p + _ind.getCenter());
}


//Constructor: SmoothIndicatorEllipsoid3D
template <typename T, typename S, bool PARTICLE>
SmoothIndicatorEllipsoid3D<T,S,PARTICLE>::SmoothIndicatorEllipsoid3D(IndicatorEllipsoid3D<S>& ind,
    S epsilon, Vector<S,3> theta)
  :SmoothIndicatorEllipsoid3D(ind.getCenter(), ind.getRadius(), epsilon, theta)
{ }

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorEllipsoid3D<T,S,PARTICLE>::SmoothIndicatorEllipsoid3D(Vector<S,3> center, Vector<S,3> radius,
    S epsilon, Vector<S,3> theta)
  :_ind(center, radius)
{
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos = center;
    this->_theta = {
      theta[0] * (M_PI/180.),
      theta[1] * (M_PI/180.),
      theta[2] * (M_PI/180.)
    };
  }

  T const max_axis = util::max( radius[0], util::max(radius[1], radius[2]) );
  this->_circumRadius = max_axis+0.5*epsilon;
  if constexpr (!PARTICLE) {
    this->_myMin = {
      this->_pos[0] - this->getCircumRadius(),
      this->_pos[1] - this->getCircumRadius(),
      this->_pos[2] - this->getCircumRadius()
    };
    this->_myMax = {
      this->_pos[0] + this->getCircumRadius(),
      this->_pos[1] + this->getCircumRadius(),
      this->_pos[2] + this->getCircumRadius()
    };
    this->init(this->_theta);
  }
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorEllipsoid3D<T, S, PARTICLE>::getVolume( )
{
  Vector<S,3> const radius = _ind.getRadius();
  return 4./3.*M_PI*radius[0]*radius[1]*radius[2];
}

template <typename T, typename S, bool PARTICLE>
Vector<S,4> SmoothIndicatorEllipsoid3D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  Vector<S,3> const radius = _ind.getRadius();
  T const mass = getVolume()*density;
  T const xHalfAxis2 = radius[0]*radius[0];
  T const yHalfAxis2 = radius[1]*radius[1];
  T const zHalfAxis2 = radius[2]*radius[2];
  Vector<S,3> mofi;
  mofi[0] = 0.2*mass*(yHalfAxis2+zHalfAxis2);
  mofi[1] = 0.2*mass*(xHalfAxis2+zHalfAxis2);
  mofi[2] = 0.2*mass*(yHalfAxis2+xHalfAxis2);
  return Vector<S,4>(mofi[0], mofi[1], mofi[2], mass);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,3> SmoothIndicatorEllipsoid3D<T, S, PARTICLE>::surfaceNormal( const Vector<S,3>& pos, const S meshSize )
{
  return _ind.surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorEllipsoid3D<T, S, PARTICLE>::signedDistance( const PhysR<S,3> input )
{
  // counter-clockwise rotation by _theta=-theta around the center of mass
  Vector<S,3> p;
  if constexpr(!PARTICLE) {
    // counter-clockwise rotation by _theta=-theta around the current position of the center of mass & translation
    p = body::motion::rotation<3,S,true>::execute(input, this->_rotMat, this->getPos());
  }
  else {
    p = input;
  }
  return _ind.signedDistance(p + _ind.getCenter());
}

template <typename T, typename S, bool PARTICLE>
bool SmoothIndicatorEllipsoid3D<T,S,PARTICLE>::operator()(T output[], const S input[])
{
  Vector<T,3> pos(0.);
  Vector<T,9> rotMatrix = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  if constexpr (!PARTICLE) {
    pos = this->getPos();
    rotMatrix = this->getRotationMatrix();
  }

  //1.Calculate distance between point and center of unrotated indicator
  T xDist = input[0] - pos[0];
  T yDist = input[1] - pos[1];
  T zDist = input[2] - pos[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x = pos[0] + rotMatrix[0]*xDist + rotMatrix[3]*yDist + rotMatrix[6]*zDist;
  T y = pos[1] + rotMatrix[1]*xDist + rotMatrix[4]*yDist + rotMatrix[7]*zDist;
  T z = pos[2] + rotMatrix[2]*xDist + rotMatrix[5]*yDist + rotMatrix[8]*zDist;

  T a = (x - pos[0]) / (_ind.getRadius()[0] - 0.5*this->getEpsilon() );
  T b = (y - pos[1]) / (_ind.getRadius()[1] - 0.5*this->getEpsilon() );
  T c = (z - pos[2]) / (_ind.getRadius()[2] - 0.5*this->getEpsilon() );
  T aEps = (x - pos[0]) / (_ind.getRadius()[0] + 0.5*this->getEpsilon() );
  T bEps = (y - pos[1]) / (_ind.getRadius()[1] + 0.5*this->getEpsilon() );
  T cEps = (z - pos[2]) / (_ind.getRadius()[2] + 0.5*this->getEpsilon() );

  if ( (a*a+b*b+c*c) <= 1. ) {
    output[0] = 1.;
    return true;
  }
  if ( (aEps*aEps+bEps*bEps+cEps*cEps) <= 1. ) {
    // TODO: Here the correct distance to the ellipsoid has to be calculated for smooth transition
    //       For now the epsilon region is taken to be 0.5
    output[0] = .5;
    return true;
  }
  output[0] = 0.;
  return false;
}


//Constructor: SmoothIndicatorSuperEllipsoid3D
template <typename T, typename S, bool PARTICLE>
SmoothIndicatorSuperEllipsoid3D<T,S,PARTICLE>::SmoothIndicatorSuperEllipsoid3D(IndicatorSuperEllipsoid3D<S>& ind,
    S epsilon, Vector<S,3> theta)
  :SmoothIndicatorSuperEllipsoid3D(ind.getCenter(), ind.getXHalfAxis(), ind.getYHalfAxis(), ind.getZHalfAxis(), ind.getExponent1(), ind.getExponent2(), epsilon, theta)
{ }

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorSuperEllipsoid3D<T,S,PARTICLE>::SmoothIndicatorSuperEllipsoid3D(Vector<S,3> center,
    S xHalfAxis, S yHalfAxis, S zHalfAxis, S exponent1, S exponent2,
    S epsilon, Vector<S,3> theta)
  :_ind(center, xHalfAxis, yHalfAxis, zHalfAxis, exponent1, exponent2)
{
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos = _ind.getCenter();
    this->_theta = {
      theta[0] * (M_PI/180.),
      theta[1] * (M_PI/180.),
      theta[2] * (M_PI/180.)
    };
  }

  T const max_axis = util::max( xHalfAxis, util::max(yHalfAxis, zHalfAxis) );
  this->_circumRadius = util::sqrt(2.)*max_axis+0.5*epsilon;
  if constexpr (!PARTICLE) {
    this->_myMin = {
      this->_pos[0] - this->getCircumRadius(),
      this->_pos[1] - this->getCircumRadius(),
      this->_pos[2] - this->getCircumRadius()
    };
    this->_myMax = {
      this->_pos[0] + this->getCircumRadius(),
      this->_pos[1] + this->getCircumRadius(),
      this->_pos[2] + this->getCircumRadius()
    };
    this->init(this->_theta);
  }
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorSuperEllipsoid3D<T, S, PARTICLE>::getVolume( )
{
  return moments(0., 0., 0.);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,4> SmoothIndicatorSuperEllipsoid3D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  T const mass = getVolume() * density;
  Vector<S,3> mofi;
  mofi[0] = ( moments(0., 2., 0.) + moments(0., 0., 2.) ) * density;
  mofi[1] = ( moments(2., 0., 0.) + moments(0., 0., 2.) ) * density;
  mofi[2] = ( moments(0., 2., 0.) + moments(2., 0., 0.) ) * density;
  return Vector<S,4>(mofi[0], mofi[1], mofi[2], mass);
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorSuperEllipsoid3D<T,S,PARTICLE>::beta(S arg1, S arg2)
{
  return (std::tgamma(arg1)*std::tgamma(arg2)) / std::tgamma(arg1+arg2);
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorSuperEllipsoid3D<T,S,PARTICLE>::moments(S p, S q, S r)
{
  S ex1 = 2./_ind.getExponent1();
  S ex2 = 2./_ind.getExponent2();
  S tmp1 = 2./(p+q+2.);
  S tmp2 = util::pow(_ind.getXHalfAxis(), p+1.) * util::pow(_ind.getYHalfAxis(), q+1.) * util::pow(_ind.getZHalfAxis(), r+1.) * ex1 * ex2;
  S tmp3 = beta( (r+1.)*(ex1/2.), (p+q+2.)*(ex2/2.)+1. ) * beta( (q+1.)*(ex2/2.), (p+1.)*(ex2/2.) );
  return tmp1 * tmp2 * tmp3;
}

template <typename T, typename S, bool PARTICLE>
Vector<S,3> SmoothIndicatorSuperEllipsoid3D<T, S, PARTICLE>::surfaceNormal( const Vector<S,3>& pos, const S meshSize )
{
  return _ind.surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorSuperEllipsoid3D<T, S, PARTICLE>::signedDistance( const PhysR<S,3> input )
{
  // counter-clockwise rotation by _theta=-theta around the center of mass
  Vector<S,3> p;
  if constexpr(!PARTICLE) {
    // counter-clockwise rotation by _theta=-theta around the current position of the center of mass & translation
    p = body::motion::rotation<3,S,true>::execute(input, this->_rotMat, this->getPos());
  }
  else {
    p = input;
  }
  return _ind.signedDistance(p + _ind.getCenter());
}

template <typename T, typename S, bool PARTICLE>
bool SmoothIndicatorSuperEllipsoid3D<T,S,PARTICLE>::operator()(T output[], const S input[])
{
  Vector<T,3> pos(0.);
  Vector<T,9> rotMatrix = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  if constexpr (!PARTICLE) {
    pos = this->getPos();
    rotMatrix = this->getRotationMatrix();
  }

  //1.Calculate distance between point and center of unrotated indicator
  T xDist = input[0] - pos[0];
  T yDist = input[1] - pos[1];
  T zDist = input[2] - pos[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x = pos[0] + rotMatrix[0]*xDist + rotMatrix[3]*yDist + rotMatrix[6]*zDist;
  T y = pos[1] + rotMatrix[1]*xDist + rotMatrix[4]*yDist + rotMatrix[7]*zDist;
  T z = pos[2] + rotMatrix[2]*xDist + rotMatrix[5]*yDist + rotMatrix[8]*zDist;

  T a = util::pow ( util::abs( (x-pos[0]) / _ind.getXHalfAxis() ), _ind.getExponent1() );
  T b = util::pow ( util::abs( (y-pos[1]) / _ind.getYHalfAxis() ), _ind.getExponent1() );
  T c = util::pow ( util::abs( (z-pos[2]) / _ind.getZHalfAxis() ), _ind.getExponent2() );
  T ab = util::pow( a+b, _ind.getExponent2()/_ind.getExponent1() );

  if ( (ab+c) <= 1. ) {
    output[0] = 1.;
    return true;
  }

  output[0] = 0.;
  return false;
}


//Constructor: SmoothIndicatorSphere3D
template <typename T, typename S, bool PARTICLE>
SmoothIndicatorSphere3D<T,S,PARTICLE>::SmoothIndicatorSphere3D(IndicatorSphere3D<S>& ind, S epsilon)
  : SmoothIndicatorSphere3D(ind.getCenter(), ind.getRadius(), epsilon)
{ }

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorSphere3D<T,S,PARTICLE>::SmoothIndicatorSphere3D(Vector<S,3> center, S radius, S epsilon)
  : _ind(center, radius)
{
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos = center;
  }

  this->_circumRadius = radius + 0.5*epsilon;
  if constexpr (!PARTICLE) {
    this->_myMin = {
      this->_pos[0] - this->getCircumRadius(),
      this->_pos[1] - this->getCircumRadius(),
      this->_pos[2] - this->getCircumRadius()
    };
    this->_myMax = {
      this->_pos[0] + this->getCircumRadius(),
      this->_pos[1] + this->getCircumRadius(),
      this->_pos[2] + this->getCircumRadius()
    };
    this->init({0.,0.,0.});
  }
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorSphere3D<T, S, PARTICLE>::getVolume( )
{
  return 4./3.*M_PI*util::pow(_ind.getRadius(), 3);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,4> SmoothIndicatorSphere3D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  T const radius = _ind.getRadius();
  T const radius2 = radius*radius;
  T const mass = getVolume()*density;
  Vector<S,3> mofi;
  mofi[0] = 2./5.*mass*radius2;
  mofi[1] = 2./5.*mass*radius2;
  mofi[2] = 2./5.*mass*radius2;
  return Vector<S,4>(mofi[0], mofi[1], mofi[2], mass);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,3> SmoothIndicatorSphere3D<T, S, PARTICLE>::surfaceNormal( const Vector<S,3>& pos, const S meshSize )
{
  return _ind.surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorSphere3D<T, S, PARTICLE>::signedDistance( const PhysR<S,3> input )
{
  // only translation necessary
  Vector<S,3> dist = input + _ind.getCenter();
  if constexpr (!PARTICLE) {
    dist -= this->_pos;
  }
  return _ind.signedDistance(dist);
}


template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCylinder3D<T,S,PARTICLE>::initIndicatorCylinder3D(Vector<S,3> theta, S length)
{
  this->_circumRadius = util::sqrt(_ind.getRadius()*_ind.getRadius()+(0.5*length)*(0.5*length))+0.5*this->getEpsilon();
  if constexpr (!PARTICLE) {
    this->_myMin = {
      this->_pos[0] - this->getCircumRadius(),
      this->_pos[1] - this->getCircumRadius(),
      this->_pos[2] - this->getCircumRadius()
    };
    this->_myMax = {
      this->_pos[0] + this->getCircumRadius(),
      this->_pos[1] + this->getCircumRadius(),
      this->_pos[2] + this->getCircumRadius()
    };
    this->init(this->_theta);
  }
}

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCylinder3D<T,S,PARTICLE>::SmoothIndicatorCylinder3D(IndicatorCylinder3D<S>& ind, S epsilon, Vector<S,3> theta)
  : SmoothIndicatorCylinder3D(ind.getCenter1(), ind.getCenter2(), ind.getRadius(), epsilon, theta)
{ }

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCylinder3D<T,S,PARTICLE>::SmoothIndicatorCylinder3D(Vector<S,3> center1, Vector<S,3> center2, S radius, S epsilon, Vector<S,3> theta)
  : _ind(center1, center2, radius)
{
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos = 0.5 * (_ind.getCenter1() + _ind.getCenter2());
    this->_theta = {
      theta[0] * (M_PI/180.),
      theta[1] * (M_PI/180.),
      theta[2] * (M_PI/180.)
    };
  }
  Vector<T,3> dist = _ind.getCenter2() - _ind.getCenter1();
  S const length = norm(dist);
  if constexpr (!PARTICLE) {
    initIndicatorCylinder3D(this->_theta, length);
  }
  else {
    initIndicatorCylinder3D(Vector<T,3>(0.), length);
  }

}

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCylinder3D<T,S,PARTICLE>::SmoothIndicatorCylinder3D(Vector<S,3> center, Vector<S,3> normal, S radius, S height, S epsilon, Vector<S,3> theta)
  : SmoothIndicatorCylinder3D(center-.5*height*normalize(normal), center+.5*height*normalize(normal), radius, epsilon, theta)
{ }

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorCylinder3D<T, S, PARTICLE>::getVolume( )
{
  const Vector<T,3> dist = _ind.getCenter2() - _ind.getCenter1();
  S const length = norm(dist);
  return M_PI*_ind.getRadius()*_ind.getRadius()*length;
}

template <typename T, typename S, bool PARTICLE>
Vector<S,4> SmoothIndicatorCylinder3D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  Vector<T,3> dist = _ind.getCenter2() - _ind.getCenter1();
  S const length = norm(dist);
  T const radius2 = _ind.getRadius() * _ind.getRadius();
  T const mass = getVolume()*density;
  Vector<S,3> mofi;
  mofi[0] = 0.5*mass*radius2;
  mofi[1] = 1/12.*mass*(length*length+3.*radius2);
  mofi[2] = 1/12.*mass*(length*length+3.*radius2);
  return Vector<S,4>(mofi[0], mofi[1], mofi[2], mass);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,3> SmoothIndicatorCylinder3D<T, S, PARTICLE>::surfaceNormal( const Vector<S,3>& pos, const S meshSize )
{
  return _ind.surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorCylinder3D<T, S, PARTICLE>::signedDistance( const PhysR<S,3> input )
{
  Vector<S,3> p;
  if constexpr(!PARTICLE) {
    // counter-clockwise rotation by _theta=-theta around the current position of the center of mass & translation
    p = body::motion::rotation<3,S,true>::execute(input, this->_rotMat, this->getPos());
  }
  else {
    p = input;
  }

  // moving to indicator coordinate system by adding start position coordinate
  return _ind.signedDistance(p + 0.5 * (_ind.getCenter2() + _ind.getCenter1()));
}


template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCone3D<T,S,PARTICLE>::SmoothIndicatorCone3D(IndicatorCone3D<S>& ind, S epsilon, Vector<S,3> theta)
  :SmoothIndicatorCone3D(ind.getCenter1(), ind.getCenter2(), ind.getRadius1(), ind.getRadius2(), epsilon, theta)
{ }

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCone3D<T,S,PARTICLE>::SmoothIndicatorCone3D(Vector<S,3> center1, Vector<S,3> center2, S radius1, S radius2, S epsilon, Vector<S,3> theta)
  :_ind(center1, center2, radius1, radius2)
{
  this->_epsilon = epsilon;
  _startPos = this->calcCenterOfMass();

  if constexpr (!PARTICLE) {
    this->_pos = _startPos;
    this->_theta = {
      theta[0] * (M_PI/180.),
      theta[1] * (M_PI/180.),
      theta[2] * (M_PI/180.)
    };
  }

  Vector<T,3> dist = _ind.getCenter2() - _ind.getCenter1();
  S const length = norm(dist);
  initIndicatorCone3D(theta, length);
}

// TODO: Add Moment of inertia for truncated cone
template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCone3D<T,S,PARTICLE>::initIndicatorCone3D(Vector<S,3> theta, S length)
{
  const T _radiusA = _ind.getRadius1();
  const T _radiusB = _ind.getRadius2();

  if (_radiusA >= _radiusB) {
    this->_circumRadius = util::sqrt(_radiusA*_radiusA+(0.5*length)*(0.5*length))+0.5*this->getEpsilon();
  }
  else {
    this->_circumRadius = util::sqrt(_radiusB*_radiusB+(0.5*length)*(0.5*length))+0.5*this->getEpsilon();
  }

  if constexpr (!PARTICLE) {
    this->_myMin = {
      this->_pos[0] - this->getCircumRadius(),
      this->_pos[1] - this->getCircumRadius(),
      this->_pos[2] - this->getCircumRadius()
    };
    this->_myMax = {
      this->_pos[0] + this->getCircumRadius(),
      this->_pos[1] + this->getCircumRadius(),
      this->_pos[2] + this->getCircumRadius()
    };

    this->init(this->_theta);
  }
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorCone3D<T, S, PARTICLE>::getVolume( )
{
  const Vector<T,3> dist = _ind.getCenter2() - _ind.getCenter1();
  S const length = norm(dist);
  return (1/3.)*M_PI*length* (_ind.getRadius1()*_ind.getRadius1()
                              +_ind.getRadius1()*_ind.getRadius2()
                              +_ind.getRadius2()*_ind.getRadius2());
}

template <typename T, typename S, bool PARTICLE>
Vector<S,4> SmoothIndicatorCone3D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  const T _radiusA = _ind.getRadius1();
  const T _radiusB = _ind.getRadius2();

  Vector<T,3> dist = _ind.getCenter2() - _ind.getCenter1();
  S const length = norm(dist);

  const T mass = getVolume()*density;
  Vector<S,3> mofi;

  //TODO: check for correctness
  //TODO: Add check that orientation fits or give a warning
  // This is only valid if the same orientation is chosen. If A and B are not chosen according to this definiton, the calculation will be wrong.
  T topRadius = _radiusA;
  T baseRadius = _radiusB;
  if (_radiusA > _radiusB) {
    topRadius = _radiusB;
    baseRadius = _radiusA;
  }
  mofi[0] = 3/10.*mass*(util::pow(baseRadius,5) - util::pow(topRadius,5))/(util::pow(baseRadius,3) - util::pow(topRadius,3));
  mofi[1] = 3/20.*mass*( (util::pow(baseRadius-topRadius, 2) + 4*length*length) * (util::pow(baseRadius,5) - util::pow(topRadius,5)) )
            / (util::pow(baseRadius-topRadius,3) * (baseRadius*baseRadius+baseRadius*topRadius+topRadius*topRadius));
  mofi[2] = mofi[1];
  return Vector<S,4>(mofi[0], mofi[1], mofi[2], mass);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,3> SmoothIndicatorCone3D<T, S, PARTICLE>::calcCenterOfMass( )
{
  T topRadius, baseRadius;
  Vector<T,3> centerTop, centerBase;

  if (_ind.getRadius1() < _ind.getRadius2()) {
    topRadius = _ind.getRadius1();
    baseRadius =  _ind.getRadius2();
    centerTop = _ind.getCenter1();
    centerBase = _ind.getCenter2();
  }
  else if (_ind.getRadius1() > _ind.getRadius2()) {
    topRadius = _ind.getRadius2();
    baseRadius =  _ind.getRadius1();
    centerTop = _ind.getCenter2();
    centerBase = _ind.getCenter1();
  }
  else {
    std::cerr << "Error calculating a cone's center of mass." << std::endl;
    assert(false);
  }

  const T height = norm(centerTop - centerBase);
  const T centerOfMassHeight = 0.25*height * (baseRadius*baseRadius + 2*baseRadius*topRadius + 3*topRadius*topRadius)
                               / (baseRadius*baseRadius + baseRadius*topRadius + topRadius*topRadius);
  return centerBase + normalize(centerTop-centerBase) * centerOfMassHeight;
}

template <typename T, typename S, bool PARTICLE>
Vector<S,3> SmoothIndicatorCone3D<T, S, PARTICLE>::surfaceNormal( const Vector<S,3>& pos, const S meshSize )
{
  return _ind.surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorCone3D<T, S, PARTICLE>::signedDistance( const PhysR<S,3> input )
{
  // counter-clockwise rotation by _theta=-theta around the current position of the center of mass
  Vector<S,3> p;
  if constexpr(!PARTICLE) {
    // counter-clockwise rotation by _theta=-theta around the current position of the center of mass & translation
    p = body::motion::rotation<3,S,true>::execute(input, this->_rotMat, this->getPos());
  }
  else {
    p = input;
  }

  // moving to indicator coordinate system by adding start position coordinate
  return _ind.signedDistance(p + _startPos);
}


//TODO: TO Be Repaired
//TODO: Check for consitency
template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCustom3D<T,S,PARTICLE>::SmoothIndicatorCustom3D(T latticeSpacing,
    std::shared_ptr<IndicatorF3D<T>> indPtr,
    Vector<T,3> pos,
    T epsilon,
    Vector<T,3> theta)
  : _indPtr(indPtr),
    _latticeSpacing(latticeSpacing),
    _center(3)
{
  OstreamManager clout(std::cout,"createIndicatorCustom3D");
  this->_name = "custom3D";
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos = pos;         // global position of the local center
    this->_theta = (M_PI/180.) * theta;
  }

  initRotationMatrix();
  initBlockData(*indPtr);

  // calculate mass and centerpoint for rotation
  calcCenter();
  // calculate min and max from circumRadius
  calcCircumRadius();
}

template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCustom3D<T,S,PARTICLE>::initRotationMatrix()
{
  if constexpr(!PARTICLE) {
    // initialize rotation matrix
    this->_rotMat = particles::dynamics::rotation_matrix<3,T>::calculate( this->_theta );
  }
}

template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCustom3D<T,S,PARTICLE>::initBlockData(IndicatorF3D<T>& ind)
{
  OstreamManager clout(std::cout,"createIndicatorCustom3D");

  // initialize temporary values
  int blockDataSize[3];
  int blockDataPadding[3];
  for (unsigned iD=0; iD<3; ++iD) {
    blockDataSize[iD] = util::ceil( (ind.getMax()[iD] - ind.getMin()[iD]) / _latticeSpacing );
    // Add a padding so that the distance can be cached in the vicinity of the geometry
    blockDataPadding[iD] = 2*util::ceil(0.05*blockDataSize[iD]);
    blockDataSize[iD] += blockDataPadding[iD];
  }

  // create blockData containing signed distance information
  this->_blockData.reset(new BlockData<3,T,BaseType<T>>({blockDataSize, 0}));
  int iX[3];
  for (iX[0]=0; iX[0] < blockDataSize[0]; iX[0]++) {
    for (iX[1]=0; iX[1] < blockDataSize[1]; iX[1]++) {
      for (iX[2]=0; iX[2] < blockDataSize[2]; iX[2]++) {
        Vector<T,3> input;
        for (unsigned iD=0; iD<3; ++iD) {
          input[iD] = (iX[iD]-blockDataPadding[iD]/2)*_latticeSpacing+ind.getMin()[iD];
        }
        this->_blockData->get(iX) = ind.signedDistance(input);
      }
    }
  }
}

template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCustom3D<T,S,PARTICLE>::calcCenter()
{
  // TODO check again for correctness of center due to smooth boundary and coordinate system
  unsigned nCells = 0;
  int input[3];
  this->_center = PhysR<T,3>(0.);
  for (input[0] = 0; input[0] < this->_blockData->getNx(); ++input[0]) {
    for (input[1] = 0; input[1] < this->_blockData->getNy(); ++input[1]) {
      for (input[2] = 0; input[2] < this->_blockData->getNz(); ++input[2]) {
        if (regardCell(input)) {
          // Always use real boundary as in other geometries too
          this->_center[0] += this->_latticeSpacing*input[0];
          this->_center[1] += this->_latticeSpacing*input[1];
          this->_center[2] += this->_latticeSpacing*input[2];
          ++nCells;
        }
      }
    }
  }
  this->_center *= 1./nCells;
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorCustom3D<T,S,PARTICLE>::getVolume( )
{
  return _volume;
}

template <typename T, typename S, bool PARTICLE>
Vector<T,4> SmoothIndicatorCustom3D<T,S,PARTICLE>::calcMofiAndMass(T rhoP)
{
  // TODO - calculation
  T cuboidMofi = util::pow(this->_latticeSpacing, 2)/ 6.0; // Single cuboid mofi at center of gravity
  unsigned nCells = 0;
  Vector<T,3> mofi = {T(0), T(0), T(0)};
  T dx, dy, dz;
  int input[3];
  for (input[0] = 0; input[0] < this->_blockData->getNx(); ++input[0]) {
    dx = util::abs(this->_latticeSpacing*input[0] - this->_center[0]);
    for (input[1] = 0; input[1] < this->_blockData->getNy(); ++input[1]) {
      dy = util::abs(this->_latticeSpacing*input[1] - this->_center[1]);
      for (input[2] = 0; input[2] < this->_blockData->getNz(); ++input[2]) {
        if (regardCell(input)) {
          dz = util::abs(this->_latticeSpacing*input[2] - this->_center[2]);
          mofi[0] += (dy*dy+dz*dz+cuboidMofi);
          mofi[1] += (dx*dx+dz*dz+cuboidMofi);
          mofi[2] += (dx*dx+dy*dy+cuboidMofi);
          ++nCells;
        }
      }
    }
  }
  _volume = nCells * util::pow(_latticeSpacing,3);
  const T mass = rhoP * _volume;
  const T cuboidMass = mass/nCells;
  mofi *= cuboidMass;

  return Vector<T,4>(mofi[0], mofi[1], mofi[2], mass);
}


template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCustom3D<T,S,PARTICLE>::calcCircumRadius()
{
  Vector<T,3> min(std::numeric_limits<olb::BaseType<T>>::max());
  Vector<T,3> max(-std::numeric_limits<olb::BaseType<T>>::max());
  Vector<T,3> distance;
  T maxDistance = -std::numeric_limits<olb::BaseType<T>>::max();

  int input[3];
  for (input[0] = 0; input[0] < this->_blockData->getNx(); ++input[0]) {
    distance[0] = this->_center[0] - this->_latticeSpacing*input[0];
    for (input[1] = 0; input[1] < this->_blockData->getNy(); ++input[1]) {
      distance[1] = this->_center[1] - this->_latticeSpacing*input[1];
      for (input[2] = 0; input[2] < this->_blockData->getNz(); ++input[2]) {
        distance[2] = this->_center[2] - this->_latticeSpacing*input[2];
        if (regardCell(input)) {
          if constexpr (!PARTICLE) {
            for (unsigned iD=0; iD<3; ++iD) {
              min[iD] = util::min(distance[iD], min[iD]);
              max[iD] = util::max(distance[iD], max[iD]);
            }
          }
          maxDistance = util::max(norm(distance), maxDistance);
        }
      }
    }
  }

  if constexpr (!PARTICLE) {
    min -= Vector<T,3>(-this->_epsilon);
    max += Vector<T,3>(-this->_epsilon);
    this->_myMin = this->_pos+min;
    this->_myMax = this->_pos+max;
  }

  this->_circumRadius = maxDistance + .5*this->_epsilon;
}

template <typename T, typename S, bool PARTICLE>
Vector<T,3> SmoothIndicatorCustom3D<T,S,PARTICLE>::getLocalCenter()
{
  return this->_center;
}

template <typename T, typename S, bool PARTICLE>
Vector<S,3> SmoothIndicatorCustom3D<T,S,PARTICLE>::surfaceNormal( const Vector<S,3>& pos, const S meshSize )
{
  return _indPtr->surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,3> SmoothIndicatorCustom3D<T,S,PARTICLE>::surfaceNormal( const Vector<S,3>& pos, const S meshSize,
    std::function<Vector<S,3>(const Vector<S,3>&)> transformPos )
{
  return _indPtr->surfaceNormal( pos, meshSize, transformPos );
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorCustom3D<T,S,PARTICLE>::signedDistance( const PhysR<S,3> input )
{
  // Translation
  T xDist = input[0];
  T yDist = input[1];
  T zDist = input[2];
  if constexpr (!PARTICLE) {
    xDist -= this->getPos()[0];
    yDist -= this->getPos()[1];
    zDist -= this->getPos()[2];
  }

  // counter-clockwise rotation by _theta=-theta around (0/0) and movement from rotation center to local center
  int x;
  int y;
  int z;
  if constexpr (!PARTICLE) {
    x = util::round((this->_center[0] + this->_rotMat[0]*xDist + this->_rotMat[3]*yDist + this->_rotMat[6]*zDist)/this->_latticeSpacing);
    y = util::round((this->_center[1] + this->_rotMat[1]*xDist + this->_rotMat[4]*yDist + this->_rotMat[7]*zDist)/this->_latticeSpacing);
    z = util::round((this->_center[2] + this->_rotMat[2]*xDist + this->_rotMat[5]*yDist + this->_rotMat[8]*zDist)/this->_latticeSpacing);
  }
  else {
    x = util::round((this->_center[0] + xDist)/this->_latticeSpacing);
    y = util::round((this->_center[1] + yDist)/this->_latticeSpacing);
    z = util::round((this->_center[2] + zDist)/this->_latticeSpacing);
  }

  if (x >= 0 && x < _blockData->getNx() && y >= 0 && y < _blockData->getNy() && z >= 0 && z < _blockData->getNz()) {
    return this->_blockData->get({x, y, z});
  }
  return 1;
}

template <typename T, typename S, bool PARTICLE>
bool SmoothIndicatorCustom3D<T,S,PARTICLE>::regardCell(int input[3])
{
  return this->_blockData->get(input) < std::numeric_limits<T>::epsilon();
}


} // namespace olb

#endif
