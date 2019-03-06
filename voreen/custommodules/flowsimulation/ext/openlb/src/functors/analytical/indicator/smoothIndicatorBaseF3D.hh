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

#include <cmath>

#include "smoothIndicatorBaseF3D.h"
#include "utilities/vectorHelpers.h"
#include "math.h"

namespace olb {


template <typename T, typename S>
SmoothIndicatorF3D<T,S>::SmoothIndicatorF3D()
  : AnalyticalF3D<T,S>(1),
    _myMin(S()), _myMax(S()), _center(S()),
    _vel(S()), _acc(S()), _theta(S()), _omega(S()), _alpha(S()), _mass(S()), _mofi(S())
{ }

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getMin()
{
  return _myMin;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getMax()
{
  return _myMax;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getCenter()
{
  return _center;
};

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getVel()
{
  return _vel;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getAcc()
{
  return _acc;
}

template <typename T, typename S>
Vector<S,6>& SmoothIndicatorF3D<T,S>::getAcc2()
{
  return _acc2;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getTheta()
{
  return _theta;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getOmega()
{
  return _omega;
}

template <typename T, typename S>
Vector<S,3>& SmoothIndicatorF3D<T,S>::getAlpha()
{
  return _alpha;
}

template <typename T, typename S>
S& SmoothIndicatorF3D<T,S>::getMass()
{
  return _mass;
}

template <typename T, typename S>
S& SmoothIndicatorF3D<T,S>::getMofi()
{
  return _mofi;
}

template <typename T, typename S>
S SmoothIndicatorF3D<T,S>::getDiam()
{
  return 2*_radius;
};

template <typename T, typename S>
S SmoothIndicatorF3D<T,S>::getRadius()
{
  return _radius+_epsilon;
};

template <typename T, typename S>
void SmoothIndicatorF3D<T,S>::setCenter(S centerX, S centerY, S centerZ)
{
  _center[0] = centerX;
  _center[1] = centerY;
  _center[2] = centerZ;
};

template <typename T, typename S>
void SmoothIndicatorF3D<T,S>::setTheta(S thetaX, S thetaY, S thetaZ)
{
  _theta[0] = thetaX;
  _theta[1] = thetaY;
  _theta[2] = thetaZ;
};

// identity to "store results"
template <typename T, typename S>
SmoothIndicatorIdentity3D<T,S>::SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T,S>& f)
  : _f(f)
{
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = _f.getMin()[i];
    this->_myMax[i] = _f.getMax()[i];
  }
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool SmoothIndicatorIdentity3D<T,S>::operator() (T output[], const S input[])
{
  _f(output, input);
  return true;
}


///////////////////////////////////////
/////     ParticleIndicator3D     /////
///////////////////////////////////////

template <typename T, typename S>
ParticleIndicatorF3D<T,S>::ParticleIndicatorF3D()
  : AnalyticalF3D<T,S>(1),
    _myMin(S()), _myMax(S()), _pos(S()),
    _vel(S()), _acc(S()), _acc2(S()), _theta(S()), _omega(S()), _alpha(S()), _alpha2(S()), _mofi(S()), _rotMat(S()), _mass(S()), _radius(S())
{ }

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getMin()
{
  return _myMin;
};


template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getMax()
{
  return _myMax;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getPos()
{
  return _pos;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getVel()
{
  return _vel;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getAcc()
{
  return _acc;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getAcc2()
{
  return _acc2;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getTheta()
{
  return _theta;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getOmega()
{
  return _omega;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getAlpha()
{
  return _alpha;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getAlpha2()
{
  return _alpha2;
};

template <typename T, typename S>
Vector<S,3>& ParticleIndicatorF3D<T,S>::getMofi()
{
  return _mofi;
};

template <typename T, typename S>
Vector<S,9>& ParticleIndicatorF3D<T,S>::getRotationMat()
{
  return _rotMat;
};

template <typename T, typename S>
S& ParticleIndicatorF3D<T,S>::getMass()
{
  return _mass;
};


template <typename T, typename S>
S& ParticleIndicatorF3D<T,S>::getRadius()
{
  return _radius;
};


template <typename T, typename S>
S ParticleIndicatorF3D<T,S>::getDiam()
{
  return 2.*_radius;
};


template <typename T, typename S>
ParticleIndicatorIdentity3D<T,S>::ParticleIndicatorIdentity3D(ParticleIndicatorF3D<T,S>& f)
  : _f(f)
{
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = _f.getMin()[i];
    this->_myMax[i] = _f.getMax()[i];
  }
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool ParticleIndicatorIdentity3D<T,S>::operator() (T output[], const S input[])
{
  _f(output, input);
  return true;
}


} // namespace olb

#endif
