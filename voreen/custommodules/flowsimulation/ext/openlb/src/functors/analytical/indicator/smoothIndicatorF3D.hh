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
#include <cmath>
#include <sstream>

#include "smoothIndicatorF3D.h"
#include "smoothIndicatorBaseF3D.h"
#include "smoothIndicatorCalcF3D.h"
#include "utilities/vectorHelpers.h"
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
template <typename T, typename S, bool HLBM>
SmoothIndicatorCuboid3D<T,S,HLBM>::SmoothIndicatorCuboid3D(Vector<S,3> center, S xLength, S yLength, S zLength, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
  : _xLength(xLength),_yLength(yLength),_zLength(zLength)
{
  this->_pos = center;
  this->_circumRadius = .5*(sqrt(xLength*xLength+yLength*yLength+zLength*zLength))+0.5*epsilon;
  this->_myMin = {
    center[0] - this->getCircumRadius(),
    center[1] - this->getCircumRadius(),
    center[2] - this->getCircumRadius()
  };
  this->_myMax = {
    center[0] + this->getCircumRadius(),
    center[1] + this->getCircumRadius(),
    center[2] + this->getCircumRadius()
  };
  this->_epsilon = epsilon;
  this->_theta = {
    theta[0] * (M_PI/180.),
    theta[1] * (M_PI/180.),
    theta[2] * (M_PI/180.)
  };
  T mass = xLength*yLength*zLength*density;
  T xLength2 = xLength*xLength;
  T yLength2 = yLength*yLength;
  T zLength2 = zLength*zLength;
  Vector<S,3> mofi;
  mofi[0] = 1./12.*mass*(yLength2+zLength2);
  mofi[1] = 1./12.*mass*(xLength2+zLength2);
  mofi[2] = 1./12.*mass*(yLength2+xLength2);
  this->init(this->_theta, vel, mass, mofi);
}

template <typename T, typename S, bool HLBM>
bool SmoothIndicatorCuboid3D<T,S,HLBM>::operator()(T output[], const S input[])
{
  //1.Calculate distance between point and center of unrotated indicator
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T zDist = input[2] - this->getPos()[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x = this->getPos()[0] + this->getRotationMatrix()[0]*xDist + this->getRotationMatrix()[3]*yDist + this->getRotationMatrix()[6]*zDist;
  T y = this->getPos()[1] + this->getRotationMatrix()[1]*xDist + this->getRotationMatrix()[4]*yDist + this->getRotationMatrix()[7]*zDist;
  T z = this->getPos()[2] + this->getRotationMatrix()[2]*xDist + this->getRotationMatrix()[5]*yDist + this->getRotationMatrix()[8]*zDist;

  //3.Calculate distance between projected point and rotated indicator bounds
  xDist = fabs(x-this->getPos()[0]) - 0.5*(_xLength-this->getEpsilon());
  yDist = fabs(y-this->getPos()[1]) - 0.5*(_yLength-this->getEpsilon());
  zDist = fabs(z-this->getPos()[2]) - 0.5*(_zLength-this->getEpsilon());

  //4.Evaluate distance
  if (xDist <= 0 && yDist <= 0 && zDist <= 0) {
    output[0] = 1.;
    return true;
  }
  if (xDist >= this->getEpsilon() || yDist >= this->getEpsilon() || zDist >= this->getEpsilon()) {
    output[0] = 0.;
    return false;
  }
  //Evaluate epsilon on edges and borders
  T dist2 = 0.;
  T epsilon2 = this->getEpsilon()*this->getEpsilon();
  if (xDist < this->getEpsilon() && xDist > 0) {
    dist2 += xDist*xDist;
  }
  if (yDist < this->getEpsilon() && yDist > 0) {
    dist2 += yDist*yDist;
  }
  if (zDist < this->getEpsilon() && zDist > 0) {
    dist2 += zDist*zDist;
  }
  if (dist2 > 0 && dist2 <= epsilon2) {
    output[0] = T(cos(M_PI2*sqrt(dist2)/this->getEpsilon())*cos(M_PI2*sqrt(dist2)/this->getEpsilon()));
    return true;
  }
  output[0] = 0.;
  return false;
}


//Constructor: SmoothIndicatorEllipsoid3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorEllipsoid3D<T,S,HLBM>::SmoothIndicatorEllipsoid3D(Vector<S,3> center, S xHalfAxis, S yHalfAxis, S zHalfAxis, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
  : _xHalfAxis(xHalfAxis),_yHalfAxis(yHalfAxis),_zHalfAxis(zHalfAxis)
{
  this->_pos = center;
  T max_axis = std::max( xHalfAxis, std::max(yHalfAxis, zHalfAxis) );
  this->_circumRadius = max_axis+0.5*epsilon;
  this->_myMin = {
    center[0] - this->getCircumRadius(),
    center[1] - this->getCircumRadius(),
    center[2] - this->getCircumRadius()
  };
  this->_myMax = {
    center[0] + this->getCircumRadius(),
    center[1] + this->getCircumRadius(),
    center[2] + this->getCircumRadius()
  };
  this->_epsilon = epsilon;
  this->_theta = {
    theta[0] * (M_PI/180.),
    theta[1] * (M_PI/180.),
    theta[2] * (M_PI/180.)
  };
  T mass = 4./3.*M_PI*xHalfAxis*yHalfAxis*zHalfAxis*density;
  T xHalfAxis2 = xHalfAxis*xHalfAxis;
  T yHalfAxis2 = yHalfAxis*yHalfAxis;
  T zHalfAxis2 = zHalfAxis*zHalfAxis;
  Vector<S,3> mofi;
  mofi[0] = 0.2*mass*(yHalfAxis2+zHalfAxis2);
  mofi[1] = 0.2*mass*(xHalfAxis2+zHalfAxis2);
  mofi[2] = 0.2*mass*(yHalfAxis2+xHalfAxis2);
  this->init(this->_theta, vel, mass, mofi);
}

template <typename T, typename S, bool HLBM>
bool SmoothIndicatorEllipsoid3D<T,S,HLBM>::operator()(T output[], const S input[])
{
  //1.Calculate distance between point and center of unrotated indicator
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T zDist = input[2] - this->getPos()[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x = this->getPos()[0] + this->getRotationMatrix()[0]*xDist + this->getRotationMatrix()[3]*yDist + this->getRotationMatrix()[6]*zDist;
  T y = this->getPos()[1] + this->getRotationMatrix()[1]*xDist + this->getRotationMatrix()[4]*yDist + this->getRotationMatrix()[7]*zDist;
  T z = this->getPos()[2] + this->getRotationMatrix()[2]*xDist + this->getRotationMatrix()[5]*yDist + this->getRotationMatrix()[8]*zDist;

  T a = (x - this->getPos()[0]) / (_xHalfAxis - 0.5*this->getEpsilon() );
  T b = (y - this->getPos()[1]) / (_yHalfAxis - 0.5*this->getEpsilon() );
  T c = (z - this->getPos()[2]) / (_zHalfAxis - 0.5*this->getEpsilon() );
  T aEps = (x - this->getPos()[0]) / (_xHalfAxis + 0.5*this->getEpsilon() );
  T bEps = (y - this->getPos()[1]) / (_yHalfAxis + 0.5*this->getEpsilon() );
  T cEps = (z - this->getPos()[2]) / (_zHalfAxis + 0.5*this->getEpsilon() );

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
template <typename T, typename S, bool HLBM>
SmoothIndicatorSuperEllipsoid3D<T,S,HLBM>::SmoothIndicatorSuperEllipsoid3D(Vector<S,3> center, S xHalfAxis, S yHalfAxis, S zHalfAxis, S exponent1, S exponent2, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
  : _xHalfAxis(xHalfAxis), _yHalfAxis(yHalfAxis), _zHalfAxis(zHalfAxis), _exp1(exponent1), _exp2(exponent2)
{
  this->_pos = center;
  T max_axis = std::max( xHalfAxis, std::max(yHalfAxis, zHalfAxis) );
  this->_circumRadius = std::sqrt(2.)*max_axis+0.5*epsilon;
  this->_myMin = {
    center[0] - this->getCircumRadius(),
    center[1] - this->getCircumRadius(),
    center[2] - this->getCircumRadius()
  };
  this->_myMax = {
    center[0] + this->getCircumRadius(),
    center[1] + this->getCircumRadius(),
    center[2] + this->getCircumRadius()
  };
  this->_epsilon = epsilon;
  this->_theta = {
    theta[0] * (M_PI/180.),
    theta[1] * (M_PI/180.),
    theta[2] * (M_PI/180.)
  };
  //T prod_axes = xHalfAxis * yHalfAxis * zHalfAxis;
  //T prod_exponents = exponent1 * exponent2;
  T mass = moments(0., 0., 0.) * density;
  Vector<S,3> mofi;
  mofi[0] = ( moments(0., 2., 0.) + moments(0., 0., 2.) ) * density;
  mofi[1] = ( moments(2., 0., 0.) + moments(0., 0., 2.) ) * density;
  mofi[2] = ( moments(0., 2., 0.) + moments(2., 0., 0.) ) * density;
  this->init(this->_theta, vel, mass, mofi);
}

template <typename T, typename S, bool HLBM>
S SmoothIndicatorSuperEllipsoid3D<T,S,HLBM>::beta(S arg1, S arg2)
{
  return (std::tgamma(arg1)*std::tgamma(arg2)) / std::tgamma(arg1+arg2);
}

template <typename T, typename S, bool HLBM>
S SmoothIndicatorSuperEllipsoid3D<T,S,HLBM>::moments(S p, S q, S r)
{
  S ex1 = 2./_exp1;
  S ex2 = 2./_exp2;
  S tmp1 = 2./(p+q+2.);
  S tmp2 = std::pow(_xHalfAxis, p+1.) * std::pow(_yHalfAxis, q+1.) * std::pow(_zHalfAxis, r+1.) * ex1 * ex2;
  S tmp3 = beta( (r+1.)*(ex1/2.), (p+q+2.)*(ex2/2.)+1. ) * beta( (q+1.)*(ex2/2.), (p+1.)*(ex2/2.) );
  return tmp1 * tmp2 * tmp3;
}


template <typename T, typename S, bool HLBM>
bool SmoothIndicatorSuperEllipsoid3D<T,S,HLBM>::operator()(T output[], const S input[])
{

  //1.Calculate distance between point and center of unrotated indicator
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T zDist = input[2] - this->getPos()[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x = this->getPos()[0] + this->getRotationMatrix()[0]*xDist + this->getRotationMatrix()[3]*yDist + this->getRotationMatrix()[6]*zDist;
  T y = this->getPos()[1] + this->getRotationMatrix()[1]*xDist + this->getRotationMatrix()[4]*yDist + this->getRotationMatrix()[7]*zDist;
  T z = this->getPos()[2] + this->getRotationMatrix()[2]*xDist + this->getRotationMatrix()[5]*yDist + this->getRotationMatrix()[8]*zDist;

  T a = std::pow ( std::abs( (x-this->getPos()[0]) / _xHalfAxis ), _exp1 );
  T b = std::pow ( std::abs( (y-this->getPos()[1]) / _yHalfAxis ), _exp1 );
  T c = std::pow ( std::abs( (z-this->getPos()[2]) / _zHalfAxis ), _exp2 );
  T ab = std::pow( a+b, _exp2/_exp1 );

  if ( (ab+c) <= 1. ) {
    output[0] = 1.;
    return true;
  }

  output[0] = 0.;
  return false;
}


//Constructor: SmoothIndicatorSphere3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorSphere3D<T,S,HLBM>::SmoothIndicatorSphere3D(Vector<S,3> center,
    S radius, S epsilon, S density, Vector<S,3> vel)
  : _radius(radius)
{
  this->_pos = center;
  this->_circumRadius = radius + 0.5*epsilon;
  this->_myMin = {
    center[0] - this->getCircumRadius(),
    center[1] - this->getCircumRadius(),
    center[2] - this->getCircumRadius()
  };
  this->_myMax = {
    center[0] + this->getCircumRadius(),
    center[1] + this->getCircumRadius(),
    center[2] + this->getCircumRadius()
  };
  this->_epsilon = epsilon;

  T mass = 4./3.*M_PI*std::pow(radius, 3.)*density;
  T radius2 = radius * radius;
  Vector<S,3> mofi;
  mofi[0] = 2./5.*mass*radius2;
  mofi[1] = 2./5.*mass*radius2;
  mofi[2] = 2./5.*mass*radius2;
  Vector<S,3> theta(0.,0.,0.);
  this->init(theta, vel, mass, mofi);
}

template <typename T, typename S, bool HLBM>
bool SmoothIndicatorSphere3D<T, S, HLBM>::operator()(T output[], const S input[])
{
  //1.Calculate distance between point and center of indicator
  T distToCenter2 = (this->getPos()[0]-input[0])*(this->getPos()[0]-input[0]) +
                    (this->getPos()[1]-input[1])*(this->getPos()[1]-input[1]) +
                    (this->getPos()[2]-input[2])*(this->getPos()[2]-input[2]);

  //3.Calculate distance between point and indicator bounds
  T rDist = sqrt(distToCenter2) - (_radius-this->getEpsilon()*0.5);

  //4. Evaluate distance
  if (rDist <= 0) {
    output[0] = 1.;
    return true;
  }
  else if (rDist > 0 && rDist < this->getEpsilon()) {
    output[0] = T(cos(M_PI2*rDist/this->getEpsilon())*cos(M_PI2*rDist/this->getEpsilon()));
    return true;
  }
  output[0] = 0.;
  return false;
}


template <typename T, typename S, bool HLBM>
void SmoothIndicatorCylinder3D<T,S,HLBM>::initIndicatorCylinder3D(Vector<S,3> normal, Vector<S,3> theta, S density, Vector<S,3> vel)
{
  T normLength = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  T dx = normal[0]*(_length/normLength);
  T dy = normal[1]*(_length/normLength);
  T dz = normal[2]*(_length/normLength);

  //Rotate according to normal orientation
  theta[0] -= asin(dy/_length);
  if (dz >= 0) {
    theta[1] += asin(dx/_length);
  }
  else {
    theta[1] += M_PI;
    theta[1] -= asin(dx/_length);
  }

  this->_circumRadius = std::sqrt(_radius*_radius+(0.5*_length)*(0.5*_length))+0.5*this->getEpsilon();
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

  T radius2 = _radius * _radius;
  T mass = M_PI*radius2*_length*density;
  Vector<S,3> mofi;
  mofi[0] = 0.5*mass*radius2;
  mofi[1] = 1/12.*mass*(_length*_length+3.*radius2);
  mofi[2] = 1/12.*mass*(_length*_length+3.*radius2);

  this->init(theta, vel, mass, mofi);
}

//Constructor 1: SmoothIndicatorCylinder3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorCylinder3D<T,S,HLBM>::SmoothIndicatorCylinder3D(Vector<S,3> pointA, Vector<S,3> pointB, S radius, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
  : _radius(radius)
{
  this->_epsilon = epsilon;
  this->_pos[0] = 0.5*(pointA[0]+pointB[0]);
  this->_pos[1] = 0.5*(pointA[1]+pointB[1]);
  this->_pos[2] = 0.5*(pointA[2]+pointB[2]);
  T dx = pointB[0]-pointA[0];
  T dy = pointB[1]-pointA[1];
  T dz = pointB[2]-pointA[2];
  theta[0] *= M_PI/180.;
  theta[1] *= M_PI/180.;
  theta[2] *= M_PI/180.;
  this->_theta = theta;
  _length = sqrt( dx*dx + dy*dy + dz*dz );
  Vector<S,3> normal (dx/_length, dy/_length, dz/_length);
  initIndicatorCylinder3D(normal, theta, density, vel);
}

//Constructor 2: SmoothIndicatorCylinder3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorCylinder3D<T,S,HLBM>::SmoothIndicatorCylinder3D(Vector<S,3> center, Vector<S,3> normal, S radius, S length, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
  : _radius(radius), _length(length)
{
  this->_pos = center;
  this->_epsilon = epsilon;
  theta[0] *= M_PI/180.;
  theta[1] *= M_PI/180.;
  theta[2] *= M_PI/180.;
  this->_theta = theta;
  initIndicatorCylinder3D(normal, theta, density, vel);
}


template <typename T, typename S, bool HLBM>
bool SmoothIndicatorCylinder3D<T,S,HLBM>::operator()(T output[],const S input[])
{
  //1.Calculate distance between point and center of unrotated (z aligned) indicator
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T zDist = input[2] - this->getPos()[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x= this->getPos()[0] + this->getRotationMatrix()[0]*xDist + this->getRotationMatrix()[3]*yDist + this->getRotationMatrix()[6]*zDist;
  T y= this->getPos()[1] + this->getRotationMatrix()[1]*xDist + this->getRotationMatrix()[4]*yDist + this->getRotationMatrix()[7]*zDist;
  T z= this->getPos()[2] + this->getRotationMatrix()[2]*xDist + this->getRotationMatrix()[5]*yDist + this->getRotationMatrix()[8]*zDist;

  //3.Calculate distance between projected point and indicator bounds
  T xyDistToCenter2 = (this->getPos()[0]-x)*(this->getPos()[0]-x)
                      + (this->getPos()[1]-y)*(this->getPos()[1]-y);
  T rDist = sqrt(xyDistToCenter2) - (_radius-this->getEpsilon()*0.5);
  zDist = fabs(z -this-> getPos()[2]) - 0.5*(_length-this->getEpsilon());

  //4.Evaluate distance
  if ( zDist <= 0 && rDist <= 0) {
    output[0] = 1.;
    return true;
  }
  if (zDist >= this->getEpsilon() || rDist >= this->getEpsilon()) {
    output[0] = 0.;
    return false;
  }
  //Evaluate epsilon on edges and borders
  T dist2 = 0.;
  T epsilon2 = this->getEpsilon()*this->getEpsilon();
  if (zDist < this->getEpsilon() && zDist > 0) {
    dist2 += zDist*zDist;
  }
  if (rDist < this->getEpsilon() && rDist > 0) {
    dist2 += rDist*rDist;
  }
  if (dist2 > 0 && dist2 < epsilon2) {
    output[0] = T(cos(M_PI2*sqrt(dist2)/this->getEpsilon())*cos(M_PI2*sqrt(dist2)/this->getEpsilon()));
    return true;
  }
  output[0] = 0.;
  return false;
}


// TODO: Add Moment of inertia for truncated cone
template <typename T, typename S, bool HLBM>
void SmoothIndicatorCone3D<T,S,HLBM>::initIndicatorCone3D(Vector<S,3> normal, Vector<S,3> theta, S density, Vector<S,3> vel)
{
  T normLength = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
  T dx = normal[0]*(_length/normLength);
  T dy = normal[1]*(_length/normLength);
  T dz = normal[2]*(_length/normLength);

  //Rotate according to normal orientation
  theta[0] -= asin(dy/_length);
  if (dz >= 0) {
    theta[1] += asin(dx/_length);
  }
  else {
    theta[1] += M_PI;
    theta[1] -= asin(dx/_length);
  }

  T radiusA2 = _radiusA *_radiusA;
  T radiusB2 = _radiusB *_radiusB;
  T mass = (1/3.)*M_PI*(radiusA2+_radiusA*_radiusB+radiusB2)*_length*density;
  Vector<S,3> mofi;

  //TODO: only valid for simple cone. Not for truncated one
  // Use the superposition principle
  // This is only valid if the same orientation is chosen. If A and B are not chosen according to this definiton, the calculation will be wrong.
  mofi[0] = 3/10.*mass*radiusA2;
  mofi[1] = 3/20.*mass*(4*_length*_length+radiusA2);
  mofi[2] = 3/20.*mass*(4*_length*_length+radiusA2);

  if (_radiusA >= _radiusB) {
    this->_circumRadius = std::sqrt(_radiusA*_radiusA+(0.5*_length)*(0.5*_length))+0.5*this->getEpsilon();
  }
  else {
    this->_circumRadius = std::sqrt(_radiusB*_radiusB+(0.5*_length)*(0.5*_length))+0.5*this->getEpsilon();
  }

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

  this->init(theta, vel, mass, mofi);
}

//Constructor 1: SmoothIndicatorCone3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorCone3D<T,S,HLBM>::SmoothIndicatorCone3D(Vector<S,3> pointA, Vector<S,3> pointB, S radiusA, S radiusB, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
  : _radiusA(radiusA), _radiusB(radiusB)
{
  this->_epsilon = epsilon;
  // TODO: Calculation of center of mass (here it is just the center of the bounding box)
  this->_pos[0] = 0.5*(pointA[0]+pointB[0]);
  this->_pos[1] = 0.5*(pointA[1]+pointB[1]);
  this->_pos[2] = 0.5*(pointA[2]+pointB[2]);
  T dx = pointB[0]-pointA[0];
  T dy = pointB[1]-pointA[1];
  T dz = pointB[2]-pointA[2];
  theta[0] *= M_PI/180.;
  theta[1] *= M_PI/180.;
  theta[2] *= M_PI/180.;
  this->_theta = theta;
  _length = sqrt( dx*dx + dy*dy + dz*dz );
  Vector<S,3> normal (dx/_length, dy/_length, dz/_length);
  initIndicatorCone3D(normal, theta, density, vel);
}


//Constructor 2: SmoothIndicatorCone3D
template <typename T, typename S, bool HLBM>
SmoothIndicatorCone3D<T,S,HLBM>::SmoothIndicatorCone3D(Vector<S,3> center, Vector<S,3> normal, S radiusA, S radiusB, S length, S epsilon, Vector<S,3> theta, S density, Vector<S,3> vel)
  : _radiusA(radiusA), _radiusB(radiusB), _length(length)
{
  this->_pos = center;
  this->_epsilon = epsilon;
  theta[0] *= M_PI/180.;
  theta[1] *= M_PI/180.;
  theta[2] *= M_PI/180.;
  this->_theta = theta;
  initIndicatorCone3D(normal, theta, density, vel);
}

template <typename T, typename S, bool HLBM>
bool SmoothIndicatorCone3D<T,S,HLBM>::operator()(T output[],const S input[])
{
  //1.Calculate distance between point and center of unrotated (z aligned) indicator
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T zDist = input[2] - this->getPos()[2];

  //2.Calculate point projected to rotated indicator
  // counter-clockwise rotation by _theta=-theta around center
  T x= this->getPos()[0] + this->getRotationMatrix()[0]*xDist + this->getRotationMatrix()[3]*yDist + this->getRotationMatrix()[6]*zDist;
  T y= this->getPos()[1] + this->getRotationMatrix()[1]*xDist + this->getRotationMatrix()[4]*yDist + this->getRotationMatrix()[7]*zDist;
  T z= this->getPos()[2] + this->getRotationMatrix()[2]*xDist + this->getRotationMatrix()[5]*yDist + this->getRotationMatrix()[8]*zDist;

  //3.Calculate distance between projected point and indicator bounds
  zDist = fabs(z -this-> getPos()[2]) - 0.5*(_length-this->getEpsilon());
  T axDist = z -this-> getPos()[2];
  T radiusAtZ = _radiusA - (axDist/_length+0.5)*(_radiusA-_radiusB);
  T xyDistToCenter2 = (this->getPos()[0]-x)*(this->getPos()[0]-x)
                      + (this->getPos()[1]-y)*(this->getPos()[1]-y);
  T rDist = sqrt(xyDistToCenter2) - (radiusAtZ-this->getEpsilon()*0.5);


  //4.Evaluate distance
  if ( zDist <= 0 && rDist <= 0) {
    output[0] = 1.;
    return true;
  }
  if (zDist >= this->getEpsilon() || rDist >= this->getEpsilon()) {
    output[0] = 0.;
    return false;
  }
  //Evaluate epsilon on edges and borders
  T dist2 = 0.;
  T epsilon2 = this->getEpsilon()*this->getEpsilon();
  if (zDist < this->getEpsilon() && zDist > 0) {
    dist2 += zDist*zDist;
  }
  if (rDist < this->getEpsilon() && rDist > 0) {
    dist2 += rDist*rDist;
  }
  if (dist2 > 0 && dist2 < epsilon2) {
    output[0] = T(cos(M_PI2*sqrt(dist2)/this->getEpsilon())*cos(M_PI2*sqrt(dist2)/this->    getEpsilon()));
    return true;
  }
  output[0] = 0.;
  return false;
}


//TODO: TO Be Repaired
//TODO: Check for consitency
template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::SmoothIndicatorCustom3D(T latticeSpacing,
    IndicatorF3D<T>& ind,
    Vector<T,3> pos,
    T rhoP,
    T epsilon,
    Vector<T,3> theta,
    Vector<T,3> vel,
    T sigma,
    bool verbose,
    bool useRealBoundary )
  : _verbose(verbose),
    _useRealBoudnary(useRealBoundary),
    _latticeSpacing(latticeSpacing),
    _sigma(sigma),
    _center(3)
{
  OstreamManager clout(std::cout,"createIndicatorCustom3D");
  this->_name = "custom3D";
  this->_pos = pos;         // global position of the local center
  this->_vel = vel;
  this->_epsilon = epsilon;
  this->_theta = {
    theta[0] * (M_PI/180.),
    theta[1] * (M_PI/180.),
    theta[2] * (M_PI/180.)
  };

  initRotationMatrix();
  initBlockData(ind);
  //const T tmpNcells = initBlockData(ind);
  //const T invNcells = 1./tmpNcells;

  // calculate mass and centerpoint for rotation
  calcCenter();
  // calculate moment of inertia and particle mass
  calcMofi(rhoP);
  // calculate min and max from circumRadius
  calcCircumRadius();

  if (this->_verbose) {
    clout << "----->>>>> Grid cell count: " << this->_blockData.getNx() << " // " << this->_blockData.getNy() << " // " << this->_blockData.getNz() << '\n';
    clout << "----->>>>> Cell dimension of object: " << this->_myMax[0]-this->_myMin[0] << " // "
          << this->_myMax[1]-this->_myMin[1] << " // "
          << this->_myMax[2]-this->_myMin[2] << '\n';
    clout << "----->>>>> Local Center: " << this->_center[0] << " // " << this->_center[1] << " // " << this->_center[2] << '\n';
    clout << "----->>>>> Local Lattice Center: " << round(this->_center[0]/this->_latticeSpacing) << " // "
          << round(this->_center[1]/this->_latticeSpacing) << " // "
          << round(this->_center[2]/this->_latticeSpacing) << '\n';
    clout << "----->>>>> Mofi: " << this->_mofi[0] << " // " << this->_mofi[1] << " // " << this->_mofi[2] << '\n';
    clout << "----->>>>> Mass: " << this->_mass << '\n';
    clout << "----->>>>> Lenght: " << this->_latticeSpacing*this->_blockData.getNx() << " // " << this->_latticeSpacing*this->_blockData.getNy() << " // " << this->_latticeSpacing*this->_blockData.getNz() << '\n';
  }
}

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::SmoothIndicatorCustom3D(UnitConverter<T,DESCRIPTOR> const& converter,
    IndicatorF3D<T>& ind,
    Vector<T,3> pos,
    T rhoP,
    T epsilon,
    Vector<T,3> theta,
    Vector<T,3> vel,
    T sigma,
    bool verbose,
    bool useRealBoundary )
  : SmoothIndicatorCustom3D(converter.getConversionFactorLength(),ind,pos,rhoP,epsilon,theta,vel,sigma,verbose,useRealBoundary)
{ }

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
void SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::initRotationMatrix()
{
  // initialize rotation matrix
  Vector<int,3> ct;
  Vector<int,3> st;
  ct[0] = std::cos(this->_theta[0]);
  ct[1] = std::cos(this->_theta[1]);
  ct[2] = std::cos(this->_theta[2]);
  st[0] = std::sin(this->_theta[0]);
  st[1] = std::sin(this->_theta[1]);
  st[2] = std::sin(this->_theta[2]);

  this->_rotMat[0] = ct[1]*ct[2];
  this->_rotMat[1] = st[0]*st[1]*ct[2] - ct[0]*st[2];
  this->_rotMat[2] = ct[0]*st[1]*ct[2] + st[0]*st[2];
  this->_rotMat[3] = ct[1]*st[2];
  this->_rotMat[4] = st[0]*st[1]*st[2] + ct[0]*ct[2];
  this->_rotMat[5] = ct[0]*st[1]*st[2] - st[0]*ct[2];
  this->_rotMat[6] = -st[1];
  this->_rotMat[7] = st[0]*ct[1];
  this->_rotMat[8] = ct[0]*ct[1];
}

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
void SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::initBlockData(IndicatorF3D<T>& ind)
{
  OstreamManager clout(std::cout,"createIndicatorCustom3D");

  // initialize temporary values
  SmoothBlockIndicator3D<T,olb::descriptors::D3Q19<>> smoothBlock(ind, this->_latticeSpacing, this->_epsilon/this->_latticeSpacing,this->_sigma);
  const int _nX = smoothBlock.getBlockData().getNx();
  const int _nY = smoothBlock.getBlockData().getNy();
  const int _nZ = smoothBlock.getBlockData().getNz();
  T tmpNcells = 0.0;

  this->_myMax[0] = 0;
  this->_myMin[0] = _nX;
  this->_myMax[1] = 0;
  this->_myMin[1] = _nY;
  this->_myMax[2] = 0;
  this->_myMin[2] = _nZ;
  // create smoothed blockData
  BlockData3D<T,BaseType<T>> block_tmp(_nX, _nY, _nZ);
  for (int iX=0; iX < _nX; iX++) {
    for (int iY=0; iY < _nY; iY++) {
      for (int iZ=0; iZ < _nZ; iZ++) {
        block_tmp.get(iX, iY, iZ) = smoothBlock.getBlockData().get(iX, iY, iZ);
        if (iX>this->_myMax[0]) {
          this->_myMax[0] = iX;
        }
        if (iX<this->_myMin[0]) {
          this->_myMin[0] = iX;
        }
        if (iY>this->_myMax[1]) {
          this->_myMax[1] = iY;
        }
        if (iY<this->_myMin[1]) {
          this->_myMin[1] = iY;
        }
        if (iZ>this->_myMax[2]) {
          this->_myMax[2] = iZ;
        }
        if (iZ<this->_myMin[2]) {
          this->_myMin[2] = iZ;
        }
        // real boundary is at 0.5 due to smoothing
        if (regardCell(block_tmp,iX,iY,iZ)) {
          tmpNcells += block_tmp.get(iX, iY, iZ);
        }
      }
    }
  }

  if (this->_verbose) {
    clout << "------MIN/MAX----------- " << this->_myMin[0] << " " << this->_myMax[0] << " " << this->_myMin[1] << " " << this->_myMax[1] << " " << this->_myMin[2] << " " << this->_myMax[2] << '\n';
  }

  this->_blockData = block_tmp;
  //return tmpNcells;
}

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
void SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::calcCenter()
{
  // TODO check again for correctness of center due to smooth boundary and coordinate system
  T invNcells = 0;
  std::fill(this->_center.begin(), this->_center.end(), 0.0);
  for (int iX = 0; iX < this->_blockData.getNx(); ++iX) {
    for (int iY = 0; iY < this->_blockData.getNy(); ++iY) {
      for (int iZ = 0; iZ < this->_blockData.getNz(); ++iZ) {
        if (regardCell(this->_blockData,iX,iY,iZ)) {
          this->_center[0] += this->_blockData.get(iX,iY,iZ)*this->_latticeSpacing*iX;
          this->_center[1] += this->_blockData.get(iX,iY,iZ)*this->_latticeSpacing*iY;
          this->_center[2] += this->_blockData.get(iX,iY,iZ)*this->_latticeSpacing*iZ;
          invNcells += this->_blockData.get(iX,iY,iZ);
        }
      }
    }
  }
  invNcells = 1./invNcells;
  std::transform(this->_center.begin(), this->_center.end(), this->_center.begin(), [&invNcells](auto& c) {
    return c*invNcells;
  });
}

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
void SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::calcMofi(T rhoP)
{
  // TODO - calculation
  T cuboidMofi = pow(this->_latticeSpacing, 2)/ 6.0; // Single cuboid mofi at center of gravity
  T invNcells = 0;
  this->_mofi[0] = 0; // x
  this->_mofi[1] = 0; // y
  this->_mofi[2] = 0; // z
  T dx, dz, dy;
  for (int iX = 0; iX < this->_blockData.getNx(); ++iX) {
    dx = std::abs(this->_latticeSpacing*iX - this->_center[0]);
    for (int iY = 0; iY < this->_blockData.getNy(); ++iY) {
      dy = std::abs(this->_latticeSpacing*iY - this->_center[1]);
      for (int iZ = 0; iZ < this->_blockData.getNz(); ++iZ) {
        if (regardCell(this->_blockData,iX,iY,iZ)) {
          dz = std::abs(this->_latticeSpacing*iZ - this->_center[2]);
          this->_mofi[0] += (dy*dy+dz*dz+cuboidMofi)*this->_blockData.get(iX,iY,iZ);
          this->_mofi[1] += (dx*dx+dz*dz+cuboidMofi)*this->_blockData.get(iX,iY,iZ);
          this->_mofi[2] += (dx*dx+dy*dy+cuboidMofi)*this->_blockData.get(iX,iY,iZ);
          invNcells += this->_blockData.get(iX,iY,iZ);
        }
      }
    }
  }
  this->_mass = rhoP * invNcells * pow(_latticeSpacing,3);
  invNcells = 1./invNcells;
  const T cuboidMass = this->_mass*invNcells;
  this->_mofi[0] *= cuboidMass;
  this->_mofi[1] *= cuboidMass;
  this->_mofi[2] *= cuboidMass;
}


template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
void SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::calcCircumRadius()
{
  T distance = 0.;
  for (int iX = 0; iX < this->_blockData.getNx(); ++iX) {
    T x = this->_latticeSpacing*iX;
    for (int iY = 0; iY < this->_blockData.getNy(); ++iY) {
      T y = this->_latticeSpacing*iY;
      for (int iZ = 0; iZ < this->_blockData.getNz(); ++iZ) {
        T z = this->_latticeSpacing*iZ;
        if (this->_blockData.get(iX,iY,iZ) > std::numeric_limits<T>::epsilon()) {
          T tmpDist = sqrt(pow(this->_center[0]-x,2)+pow(this->_center[1]-y,2)+pow(this->_center[2]-z,2));
          if (tmpDist > distance) {
            distance = tmpDist;
          }
        }
      }
    }
  }
  this->_myMin[0] = -this->_epsilon - distance;
  this->_myMin[1] = -this->_epsilon - distance;
  this->_myMin[2] = -this->_epsilon - distance;
  this->_myMax[0] = this->_epsilon + distance;
  this->_myMax[1] = this->_epsilon + distance;
  this->_myMax[2] = this->_epsilon + distance;

  this->_circumRadius = .5*(sqrt(pow(this->_myMax[0]-this->_myMin[0],2)
                                 +pow(this->_myMax[1]-this->_myMin[1],2)
                                 +pow(this->_myMax[2]-this->_myMin[2],2)))+0.5*this->_epsilon;
}

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
Vector<T,3> SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::getLocalCenter()
{
  return this->_center;
}

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
Vector<T,3> SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::getMofi()
{
  return this->_mofi;
}

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
bool SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::regardCell(BlockData3D<T,T>& blockData, int x, int y, int z)
{
  return (!this->_useRealBoudnary || this->_blockData.get(x,y,z) > 0.5-std::numeric_limits<T>::epsilon());
}

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
bool SmoothIndicatorCustom3D<T,S,DESCRIPTOR,HLBM>::operator() (T output[], const S input[])
{
  // Translation
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];
  T zDist = input[2] - this->getPos()[2];

  // counter-clockwise rotation by _theta=-theta around (0/0) and movement from rotation center to local center
  int x = round((this->_center[0] + this->_rotMat[0]*xDist + this->_rotMat[3]*yDist + this->_rotMat[6]*zDist)/this->_latticeSpacing);
  int y = round((this->_center[1] + this->_rotMat[1]*xDist + this->_rotMat[4]*yDist + this->_rotMat[7]*zDist)/this->_latticeSpacing);
  int z = round((this->_center[2] + this->_rotMat[2]*xDist + this->_rotMat[5]*yDist + this->_rotMat[8]*zDist)/this->_latticeSpacing);

  if (x >= 0 && x < _blockData.getNx() && y >= 0 && y < _blockData.getNy() && z >= 0 && z < _blockData.getNz()) {
    if (this->_blockData.get(x, y, z) > std::numeric_limits<T>::epsilon()) {
      output[0] = T(this->_blockData.get(x, y, z));
      return true;
    }
  }
  output[0] = T(0);
  return false;
}


} // namespace olb

#endif
