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

#ifndef SMOOTH_INDICATOR_BASE_F_2D_H
#define SMOOTH_INDICATOR_BASE_F_2D_H

#include <vector>

#include "core/vector.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalBaseF.h"

namespace olb {


/** SmoothIndicatorF2D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _center
  * \param _diam
  */
template <typename T, typename S>
class SmoothIndicatorF2D : public AnalyticalF2D<T,S> {
protected:
  SmoothIndicatorF2D();
  Vector<S,2> _myMin;
  Vector<S,2> _myMax;
  Vector<S,2> _center;
  Vector<S,2> _vel;
  Vector<S,2> _acc;
  Vector<S,3> _acc2;
  S _theta;
  S _omega;
  S _alpha;

  S _mass;
  S _mofi; //Moment of Inertia

  S _epsilon;
  S _radius;
public:
  virtual Vector<S,2>& getMin();
  virtual Vector<S,2>& getMax();
  virtual Vector<S,2>& getCenter();
  virtual Vector<S,2>& getVel();
  virtual Vector<S,2>& getAcc();
  virtual Vector<S,3>& getAcc2();
  virtual S& getTheta();
  virtual S& getOmega();
  virtual S& getAlpha();
  virtual S& getMass();
  virtual S& getMofi();
  virtual S getDiam();
  virtual S getRadius();

  SmoothIndicatorF2D<T,S>& operator+(SmoothIndicatorF2D<T,S>& rhs);
};


template <typename T, typename S>
class SmoothIndicatorIdentity2D : public SmoothIndicatorF2D<T,S> {
protected:
  SmoothIndicatorF2D<T,S>& _f;
public:
  SmoothIndicatorIdentity2D(SmoothIndicatorF2D<T,S>& f);
  bool operator() (T output[], const S input[]) override;
};

/** ParticleIndicatorF2D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _center
  * \param _diam
  */
template <typename T, typename S>
class ParticleIndicatorF2D : public AnalyticalF2D<T,S> {
protected:
  ParticleIndicatorF2D();
  Vector<S,2> _pos;
  Vector<S,2> _vel;
  Vector<S,2> _acc;
  Vector<S,2> _acc2;
  Vector<S,4> _rotMat;  //saved values of rotation matrix
  S _theta;
  S _omega;
  S _alpha;
  S _alpha2;
  S _mass;
  S _mofi; //Moment of Inertia
  S _epsilon;
  S _circumradius;

public:
  Vector<S,2>& getVel()
  {
    return _vel;
  };
  Vector<S,2>& getAcc()
  {
    return _acc;
  };
  Vector<S,2>& getAcc2()
  {
    return _acc2;
  };
  Vector<S,2>& getPos()
  {
    return _pos;
  };
  S& getTheta()
  {
    return _theta;
  };
  S& getOmega()
  {
    return _omega;
  };
  S& getAlpha()
  {
    return _alpha;
  };
  S& getAlpha2()
  {
    return _alpha2;
  };
  S& getMass()
  {
    return _mass;
  };
  S& getMofi()
  {
    return _mofi;
  };
  S& getCircumRadius()
  {
    return _circumradius;
  };
  Vector<S,4>& getRotationMat()
  {
    return _rotMat;
  };

  ParticleIndicatorF2D<T,S>& operator+(ParticleIndicatorF2D<T,S>& rhs);
};


template <typename T, typename S>
class ParticleIndicatorIdentity2D : public ParticleIndicatorF2D<T,S> {
protected:
  ParticleIndicatorF2D<T,S>& _f;
public:
  ParticleIndicatorIdentity2D(ParticleIndicatorF2D<T,S>& f);
  bool operator() (T output[], const S input[]) override;
};


}

#endif
