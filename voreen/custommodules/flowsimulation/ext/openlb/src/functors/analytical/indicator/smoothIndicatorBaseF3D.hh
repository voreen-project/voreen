/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Albert Mink, Mathias J. Krause, Benjamin Förster
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

#ifndef SMOOTH_INDICATOR_BASE_F_3D_HH
#define SMOOTH_INDICATOR_BASE_F_3D_HH

#include "utilities/omath.h"

#include "smoothIndicatorBaseF3D.h"
#include "utilities/vectorHelpers.h"

namespace olb {


template <typename T, typename S>
SmoothIndicatorF3D<T, S, false>::SmoothIndicatorF3D()
  : AnalyticalF3D<T,S>(1),
    _myMin(S()), _myMax(S()), _pos(S()), _rotMat(S()), _circumRadius(S()),
    _theta(S()), _epsilon(S())
{ }

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,false>::init(Vector<S,3> theta)
{
  this->_rotMat[0] = util::cos(theta[1])*util::cos(theta[2]);
  this->_rotMat[1] = util::sin(theta[0])*util::sin(theta[1])*util::cos(theta[2]) - util::cos(theta[0])*util::sin(theta[2]);
  this->_rotMat[2] = util::cos(theta[0])*util::sin(theta[1])*util::cos(theta[2]) + util::sin(theta[0])*util::sin(theta[2]);
  this->_rotMat[3] = util::cos(theta[1])*util::sin(theta[2]);
  this->_rotMat[4] = util::sin(theta[0])*util::sin(theta[1])*util::sin(theta[2]) + util::cos(theta[0])*util::cos(theta[2]);
  this->_rotMat[5] = util::cos(theta[0])*util::sin(theta[1])*util::sin(theta[2]) - util::sin(theta[0])*util::cos(theta[2]);
  this->_rotMat[6] = -util::sin(theta[1]);
  this->_rotMat[7] = util::sin(theta[0])*util::cos(theta[1]);
  this->_rotMat[8] = util::cos(theta[0])*util::cos(theta[1]);
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, false>::getMin() const
{
  return _myMin;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T, S, false>::getMax() const
{
  return _myMax;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,false>::getPos() const
{
  return _pos;
}

template <typename T, typename S>
const S& SmoothIndicatorF3D<T,S,false>::getCircumRadius() const
{
  return _circumRadius;
}

template <typename T, typename S>
const Vector<S,9>& SmoothIndicatorF3D<T,S,false>::getRotationMatrix() const
{
  return _rotMat;
}

template <typename T, typename S>
const Vector<S,3>& SmoothIndicatorF3D<T,S,false>::getTheta() const
{
  return _theta;
}

template <typename T, typename S>
const S& SmoothIndicatorF3D<T,S,false>::getEpsilon() const
{
  return _epsilon;
}

template <typename T, typename S>
std::string SmoothIndicatorF3D<T,S,false>::name()
{
  return _name;
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,false>::setPos(Vector<S,3> pos)
{
  _pos = pos;
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,false>::setTheta(Vector<S,3> theta)
{
  _theta[0] = theta[0];
  _theta[1] = theta[1];
  _theta[2] = theta[2];

  T const cos0 = util::cos(theta[0]);
  T const cos1 = util::cos(theta[1]);
  T const cos2 = util::cos(theta[2]);
  T const sin0 = util::sin(theta[0]);
  T const sin1 = util::sin(theta[1]);
  T const sin2 = util::sin(theta[2]);

  this->_rotMat[0] = cos1 * cos2;
  this->_rotMat[1] = sin0*sin1*cos2 - cos0*sin2;
  this->_rotMat[2] = cos0*sin1*cos2 + sin0*sin2;
  this->_rotMat[3] = cos1*sin2;
  this->_rotMat[4] = sin0*sin1*sin2 + cos0*cos2;
  this->_rotMat[5] = cos0*sin1*sin2 - sin0*cos2;
  this->_rotMat[6] = -sin1;
  this->_rotMat[7] = sin0*cos1;
  this->_rotMat[8] = cos0*cos1;
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,false>::setEpsilon(S epsilon)
{
  _epsilon = epsilon;
}

template <typename T, typename S>
S SmoothIndicatorF3D<T, S, false>::getVolume( )
{
  // TODO: Fallback
  assert(false);
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename T, typename S>
Vector<S,4> SmoothIndicatorF3D<T, S, false>::calcMofiAndMass( S density )
{
  // TODO: Fallback
  assert(false);
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename T, typename S>
Vector<S,3> SmoothIndicatorF3D<T, S, false>::calcCenterOfMass( )
{
  // TODO: Fallback
  assert(false);
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename T, typename S>
bool SmoothIndicatorF3D<T, S, false>::operator()( T output[], const S input[] )
{
  T const signedDist = this->signedDistance(input);
  return sdf::evalPorosity(output, signedDist, this->getEpsilon());
}

template <typename T, typename S>
bool SmoothIndicatorF3D<T, S, false>::isInsideCircumRadius( const PhysR<S,3>& input )
{
  return norm(input-this->getPos()) <= this->getCircumRadius();
}

template <typename T, typename S>
Vector<S,3> SmoothIndicatorF3D<T, S, false>::surfaceNormal( const Vector<S,3>& pos, const S meshSize )
{
  return surfaceNormal(pos, meshSize, [&](const Vector<S,3>& pos) {
    return pos;
  });
}

template <typename T, typename S>
Vector<S,3> SmoothIndicatorF3D<T, S, false>::surfaceNormal( const Vector<S,3>& pos, const S meshSize,
    std::function<Vector<S,3>(const Vector<S,3>&)> transformPos )
{
  return util::surfaceNormal(pos, meshSize, [&](const Vector<S,3>& pos) {
    return this->signedDistance( transformPos(pos) );
  });
}

template <typename T, typename S>
const S SmoothIndicatorF3D<T, S, false>::signedDistance( const PhysR<T,3> input )
{
  // TODO: Raymarching as fallback
  assert(false);
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename T, typename S>
bool SmoothIndicatorF3D<T, S, false>::distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, S precision, S pitch)
{
  S const halfEps = this->getEpsilon() * 0.5;
  bool originValue;
  bool currentValue;

  // start at origin and move into given direction
  PhysR<S,3> currentPoint(origin);

  originValue = this->signedDistance(origin) <= halfEps;
  currentValue = this->signedDistance(currentPoint) <= halfEps;

  while (currentValue == originValue && this->isInsideCircumRadius(currentPoint)) {
    currentPoint += direction;
    // update currentValue until the first point on the other side (inside/outside) is found
    currentValue = this->signedDistance(currentPoint) <= halfEps;
  }

  // return false if no point was found in given direction
  if (!this->isInsideCircumRadius(currentPoint) && !originValue) {
    return false;
  }


  while (pitch >= precision) {
    if (!this->isInsideCircumRadius(currentPoint) && originValue) {
      currentPoint -= pitch * direction;
      pitch /= 2.;
    }
    else {
      currentValue = this->signedDistance(currentPoint) <= halfEps;
      if (currentValue == originValue) {
        currentPoint += pitch * direction;
        pitch /= 2.;
      }
      else {
        currentPoint -= pitch * direction;
        pitch /= 2.;
      }
    }
  }

  distance = norm(currentPoint - origin);
  return true;
}


// identity to "store results"
template <typename T, typename S>
SmoothIndicatorIdentity3D<T,S>::SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T,S,false>& f)
  : _f(f)
{
  this->_myMin = _f.getMin();
  this->_myMax = _f.getMax();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool SmoothIndicatorIdentity3D<T,S>::operator() (T output[], const S input[])
{
  _f(output, input);
  return true;
}

///////////////////////// for template specialisation HLBM=true

template <typename T, typename S>
SmoothIndicatorF3D<T,S,true>::SmoothIndicatorF3D()
  : AnalyticalF3D<T,S>(1),
    _circumRadius(S()), _epsilon(S())
{ }

template <typename T, typename S>
bool SmoothIndicatorF3D<T, S, true>::operator()( T output[], const S input[] )
{
#ifdef OLB_DEBUG
  OstreamManager clout(std::cout, "SmoothIndicator3D");
  clout << "WARNING: SmoothIndicatorF3D::operator() a particle (= true) SmoothIndicator does not consider the current position of the particle. Please use the evalPorosity method for this." << std::endl;
#endif
  T const signedDist = this->signedDistance(input);
  return sdf::evalPorosity(output, signedDist, this->getEpsilon());
}

template <typename T, typename S>
bool SmoothIndicatorF3D<T, S, true>::isInsideCircumRadius( const PhysR<S,3>& input )
{
  return norm(input) <= this->getCircumRadius();
}

template <typename T, typename S>
Vector<S,3> SmoothIndicatorF3D<T, S, true>::surfaceNormal( const Vector<S,3>& pos, const S meshSize )
{
  return surfaceNormal(pos, meshSize, [&](const Vector<S,3>& pos) {
    return pos;
  });
}

template <typename T, typename S>
Vector<S,3> SmoothIndicatorF3D<T, S, true>::surfaceNormal( const Vector<S,3>& pos, const S meshSize,
    std::function<Vector<S,3>(const Vector<S,3>&)> transformPos )
{
  return util::surfaceNormal(pos, meshSize, [&](const Vector<S,3>& pos) {
    return this->signedDistance( transformPos(pos) );
  });
}

template <typename T, typename S>
const S SmoothIndicatorF3D<T, S, true>::signedDistance( const PhysR<T,3> input )
{
  // TODO: Raymarching as fallback
  assert(false);
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename T, typename S>
bool SmoothIndicatorF3D<T, S, true>::distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, S precision, S pitch)
{
  S const halfEps = this->getEpsilon() * 0.5;
  bool originValue;
  bool currentValue;

  // start at origin and move into given direction
  PhysR<S,3> currentPoint(origin);

  originValue = this->signedDistance(origin) <= halfEps;
  currentValue = this->signedDistance(currentPoint) <= halfEps;

  while (currentValue == originValue && this->isInsideCircumRadius(currentPoint)) {
    currentPoint += direction;
    // update currentValue until the first point on the other side (inside/outside) is found
    currentValue = this->signedDistance(currentPoint) <= halfEps;
  }

  // return false if no point was found in given direction
  if (!this->isInsideCircumRadius(currentPoint) && !originValue) {
    return false;
  }


  while (pitch >= precision) {
    if (!this->isInsideCircumRadius(currentPoint) && originValue) {
      currentPoint -= pitch * direction;
      pitch /= 2.;
    }
    else {
      currentValue = this->signedDistance(currentPoint) <= halfEps;
      if (currentValue == originValue) {
        currentPoint += pitch * direction;
        pitch /= 2.;
      }
      else {
        currentPoint -= pitch * direction;
        pitch /= 2.;
      }
    }
  }

  distance = norm(currentPoint - origin);
  return true;
}

template <typename T, typename S>
const S& SmoothIndicatorF3D<T,S,true>::getCircumRadius() const
{
  return _circumRadius;
}

template <typename T, typename S>
const S& SmoothIndicatorF3D<T, S, true>::getEpsilon() const
{
  return _epsilon;
}

template <typename T, typename S>
std::string SmoothIndicatorF3D<T,S,true>::name()
{
  return _name;
}

template <typename T, typename S>
void SmoothIndicatorF3D<T,S,true>::setEpsilon(S epsilon)
{
  _epsilon = epsilon;
}

template <typename T, typename S>
S SmoothIndicatorF3D<T, S, true>::getVolume( )
{
  // TODO: Fallback
  assert(false);
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename T, typename S>
Vector<S,4> SmoothIndicatorF3D<T, S, true>::calcMofiAndMass( S density )
{
  // TODO: Fallback
  assert(false);
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename T, typename S>
Vector<S,3> SmoothIndicatorF3D<T, S, true>::calcCenterOfMass( )
{
  // TODO: Fallback
  assert(false);
  return std::numeric_limits<double>::quiet_NaN();
}

} // namespace olb


#endif
