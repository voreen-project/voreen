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

#ifndef SMOOTH_INDICATOR_BASE_F_3D_H
#define SMOOTH_INDICATOR_BASE_F_3D_H

#include <vector>

#include "core/vector.h"
#include "core/blockStructure.h"
#include "functors/analytical/analyticalBaseF.h"
#include "functors/genericF.h"
#include "sdf.h"

namespace olb {


template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorF3D;

/** SmoothIndicatorF3D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _center
  * \param _diam_
  */
template <typename T, typename S>
class SmoothIndicatorF3D<T,S,false> : public AnalyticalF3D<T,S> {
protected:
  SmoothIndicatorF3D();
  //~SmoothIndicatorF3D() override {};
  Vector<S,3> _myMin;
  Vector<S,3> _myMax;
  Vector<S,3> _pos;
  Vector<S,9> _rotMat;  //saved values of rotation matrix
  S _circumRadius;
  Vector<S,3> _theta;
  //Vector<S,3> _scale;
  S _epsilon;
  std::string _name = "smoothIndicator3D";
public:
  void init(Vector<S,3> theta);
  const Vector<S,3>& getMin() const;
  const Vector<S,3>& getMax() const;
  const Vector<S,3>& getPos() const;
  const Vector<S,9>& getRotationMatrix() const;
  const Vector<S,3>& getTheta() const;
  const S& getCircumRadius() const;
  const S& getEpsilon() const;
  std::string name();
  void setPos(Vector<S,3> pos);
  void setTheta(Vector<S,3> theta);   //set angle in radian
  void setEpsilon(S epsilon);
  virtual S getVolume( );
  virtual Vector<S,4> calcMofiAndMass(S density);
  virtual Vector<S,3> calcCenterOfMass();
  virtual Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize);
  virtual Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize,
                                    std::function<Vector<S,3>(const Vector<S,3>&)> transformPos);
  virtual const S signedDistance(const PhysR<T,3> input);
  virtual bool distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, S precision, S pitch);
  virtual bool operator()(T output[], const S input[]);
  bool isInsideCircumRadius(const PhysR<S,3>& input);

  SmoothIndicatorF3D<T,S,false>& operator+(SmoothIndicatorF3D<T,S,false>& rhs);
};


template <typename T, typename S>
class SmoothIndicatorIdentity3D : public SmoothIndicatorF3D<T, S, false> {
protected:
  SmoothIndicatorF3D<T, S, false>& _f;
public:
  SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T, S, false>& f);
  bool operator() (T output[], const S input[]) override;
};

/** SmoothIndicatorF3D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * Base class for specific SmoothIndicator implementation providing common data..
  */
template <typename T, typename S>
class SmoothIndicatorF3D<T, S, true> : public AnalyticalF3D<T,S> {
protected:
  SmoothIndicatorF3D();
  S _circumRadius;
  S _epsilon;
  std::string _name="3D-Particle surface";

public:
  const S& getCircumRadius() const;
  const S& getEpsilon() const;
  std::string name();
  void setEpsilon(S epsilon);
  virtual S getVolume( );
  virtual Vector<S,4> calcMofiAndMass(S density);
  virtual Vector<S,3> calcCenterOfMass();
  virtual Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize);
  virtual Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize,
                                    std::function<Vector<S,3>(const Vector<S,3>&)> transformPos);
  virtual const S signedDistance(const PhysR<T,3> input);
  virtual bool distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, S precision, S pitch);
  virtual bool operator()(T output[], const S input[]);
  bool isInsideCircumRadius(const PhysR<S,3>& input);
};

}

#endif
