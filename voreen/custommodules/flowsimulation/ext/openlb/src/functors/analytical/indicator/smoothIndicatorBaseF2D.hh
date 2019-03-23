/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Benjamin Förster
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

#ifndef SMOOTH_INDICATOR_BASE_F_2D_HH
#define SMOOTH_INDICATOR_BASE_F_2D_HH

#include <cmath>

#include "smoothIndicatorBaseF2D.h"
#include "math.h"

namespace olb {


template <typename T, typename S>
SmoothIndicatorF2D<T,S>::SmoothIndicatorF2D()
  : AnalyticalF2D<T,S>(1),
    _myMin(S()), _myMax(S()), _center(S()),
    _vel(S()), _acc(S()), _theta(S()), _omega(S()), _alpha(S()), _mass(S()), _mofi(S())
{ }

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorF2D<T,S>::getMin()
{
  return _myMin;
}

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorF2D<T,S>::getMax()
{
  return _myMax;
}

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorF2D<T,S>::getCenter()
{
  return _center;
};

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorF2D<T,S>::getVel()
{
  return _vel;
}

template <typename T, typename S>
Vector<S,2>& SmoothIndicatorF2D<T,S>::getAcc()
{
  return _acc;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF2D<T,S>::getAcc2()
{
  return _acc2;
}

template <typename T, typename S>
S& SmoothIndicatorF2D<T,S>::getTheta()
{
  return _theta;
}

template <typename T, typename S>
S& SmoothIndicatorF2D<T,S>::getOmega()
{
  return _omega;
}

template <typename T, typename S>
S& SmoothIndicatorF2D<T,S>::getAlpha()
{
  return _alpha;
}

template <typename T, typename S>
S& SmoothIndicatorF2D<T,S>::getMass()
{
  return _mass;
}

template <typename T, typename S>
S& SmoothIndicatorF2D<T,S>::getMofi()
{
  return _mofi;
}

template <typename T, typename S>
S SmoothIndicatorF2D<T,S>::getDiam()
{
  return 2*_radius;
};

template <typename T, typename S>
S SmoothIndicatorF2D<T,S>::getRadius()
{
  return _radius+_epsilon;
};


// identity to "store results"
template <typename T, typename S>
SmoothIndicatorIdentity2D<T,S>::SmoothIndicatorIdentity2D(SmoothIndicatorF2D<T,S>& f)
  : _f(f)
{
  this->_myMin = _f.getMin();
  this->_myMax = _f.getMax();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool SmoothIndicatorIdentity2D<T,S>::operator() (T output[], const S input[])
{
  _f(output, input);
  return true;
}


template <typename T, typename S>
ParticleIndicatorF2D<T,S>::ParticleIndicatorF2D()
  : AnalyticalF2D<T,S>(1),
    _pos(S()),
    _vel(S()), _acc(S()), _acc2(S()), _theta(S()), _omega(S()), _alpha(S()), _alpha2(S()), _mass(S()), _mofi(S()), _circumradius(S())
{ }

// identity to "store results"
template <typename T, typename S>
ParticleIndicatorIdentity2D<T,S>::ParticleIndicatorIdentity2D(ParticleIndicatorF2D<T,S>& f)
  : _f(f)
{
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool ParticleIndicatorIdentity2D<T,S>::operator() (T output[], const S input[])
{
  _f(output, input);
  return true;
}


} // namespace olb

#endif
