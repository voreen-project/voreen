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

#ifndef INDICATOR_BASE_F_3D_H
#define INDICATOR_BASE_F_3D_H

#include <functional>
#include <vector>

#include "core/vector.h"
#include "functors/genericF.h"
#include "indicatorBase.h"

namespace olb {


/** IndicatorF3D is an application from \f$ \Omega \subset R^3 \to \{0,1\} \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  */
template <typename S>
class IndicatorF3D : public GenericF<bool,S> {
protected:
  IndicatorF3D();
  Vector<S,3> _myMin;
  Vector<S,3> _myMax;
public:
  virtual Vector<S,3>& getMin();
  virtual Vector<S,3>& getMax();
  /** \returns false or true and pos. distance if there was one found for a given origin and direction.
   * Mind that the default computation is done by a numerical approximation which searches .. [TODO: CYRIL]
   */
  virtual bool distance(S& distance, const Vector<S,3>& origin, S precision, const Vector<S,3>& direction);
  virtual bool distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, S precision, S pitch);
  virtual bool distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC=-1);
  virtual bool distance(S& distance, const Vector<S,3>& origin);
  virtual bool distance(S& distance, const S input[]);
  /// returns true and the normal if there was one found for an given origin and direction
  /**
   * (mind that the default computation is done by a numerical approximation which searches .. [TODO])
   */
  virtual bool normal(Vector<S,3>& normal, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC=-1);
  ///Rotate vector around axis by angle theta
  virtual bool rotOnAxis(Vector<S,3>& vec_rot, const Vector<S,3>& vec, const Vector<S,3>& axis, S& theta);
  /// Returns true if input is inside the indicator
  virtual bool operator() (bool output[1], const S input[3]);
  /// Returns signed distance to the nearest point on the indicator surface
  virtual S signedDistance(const Vector<S,3>& input);
  /// Return surface normal
  virtual Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize);
  /// Return surface normal after possible translation and rotation
  Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize,
                            std::function<Vector<S,3>(const Vector<S,3>&)> transformPos);
  /// Returns true if `point` is inside a cube with corners `_myMin` and `_myMax`
  bool isInsideBox(Vector<S,3> point);
};


template <typename S>
class IndicatorIdentity3D : public IndicatorF3D<S> {
public:
  std::shared_ptr<IndicatorF3D<S>> _f;

  IndicatorIdentity3D(std::shared_ptr<IndicatorF3D<S>> f);
  bool operator() (bool output[1], const S input[3]);
};


}

#endif
