/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Jan E. Marquardt
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

#ifndef INDIC_MOD_HH
#define INDIC_MOD_HH

#include "indicMod.h"

namespace olb {

template <typename S, unsigned D>
IndicInverse<S,D>::IndicInverse( FunctorPtr<IndicatorF<S,D>> f, PhysR<S,D> min, PhysR<S,D> max )
  : _f(std::move(f))
{
  this->_myMin = min;
  this->_myMax = max;
}

template <typename S, unsigned D>
S IndicInverse<S,D>::signedDistance(const Vector<S, D>& input)
{
  return -this->_f->signedDistance(input);
}

/*
template <typename S, unsigned D>
bool IndicInverse<S,D>::operator()(bool output[], const S input[])
{
  this->_f->operator()(output, input);
  output[0] = !output[0];
  return output[0];
}
*/

}
#endif