/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Albert Mink, Jan E. Marquardt, Anna Husfeldt
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

#ifndef INDIC_ALT_3D_HH
#define INDIC_ALT_3D_HH

#include "indicAlt3D.h"

namespace olb {

template <typename S>
IndicElongation3D<S>::IndicElongation3D(FunctorPtr<IndicatorF3D<S>> f, Vector<S, 3> h)
  : _f(std::move(f))
{
  setElongation(h);
}

template <typename S>
Vector<S,3> IndicElongation3D<S>::getEstimatedCenter()
{
  return 0.5 * (_f->getMin() + _f->getMax());
}

template <typename S>
S IndicElongation3D<S>::signedDistance(const Vector<S, 3>& input)
{
  std::function<S(Vector<S, 3>)> sdf = [this](Vector<S, 3> input) {
    return this->_f->signedDistance(input);
  };
  return sdf::elongation(sdf, Vector<S, 3>(input), _h, getEstimatedCenter());
}

template <typename S>
void IndicElongation3D<S>::setElongation(Vector<S, 3> h)
{
  _h = h;
  for (int i = 0; i < 3; i++) {
    this->_myMin[i] = _f->getMin()[i] - h[i];
    this->_myMax[i] = _f->getMax()[i] + h[i];
  }
}

template <typename S>
Vector<S, 3> IndicElongation3D<S>::getElongation()
{
  return _h;
}

}
#endif