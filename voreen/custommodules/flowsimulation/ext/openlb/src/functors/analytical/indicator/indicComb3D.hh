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

#ifndef INDIC_COMB_3D_HH
#define INDIC_COMB_3D_HH

#include "indicComb3D.h"

namespace olb {


template <typename S, template<typename U> class F>
IndicComb3D<S,F>::IndicComb3D(std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g)
  : _f(f), _g(g)
{ }

template <typename S, template<typename U> class F>
bool IndicComb3D<S,F>::operator()( bool output[], const S input[3])
{
  // componentwise operation on equidimensional functors
  bool* outputF = output;
  _f->operator()(outputF, input);

  bool outputG[this->getTargetDim()];
  _g->operator()(outputG, input);

  for (int i = 0; i < this->getTargetDim(); i++) {
    output[i] = F<S>()(outputF[i], outputG[i]);
  }
  return output;
}

template <typename S>
IndicPlus3D<S>::IndicPlus3D(std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g)
  : IndicComb3D<S,util::plus>(f, g)
{
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = util::min(this->_f->getMin()[i], this->_g->getMin()[i]);
    this->_myMax[i] = util::max(this->_f->getMax()[i], this->_g->getMax()[i]);
  }
}

template <typename S>
S IndicPlus3D<S>::signedDistance( const Vector<S,3>& input )
{
  return sdf::unify(this->_f->signedDistance(input), this->_g->signedDistance(input));
}


template <typename S>
IndicMinus3D<S>::IndicMinus3D(std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g)
  : IndicComb3D<S,util::minus>(f, g)
{
  // TODO: Improve
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = this->_f->getMin()[i];
    this->_myMax[i] = this->_f->getMax()[i];
  }
}

template <typename S>
S IndicMinus3D<S>::signedDistance( const Vector<S,3>& input )
{
  return sdf::subtraction(this->_f->signedDistance(input), this->_g->signedDistance(input));
}


template <typename S>
IndicMultiplication3D<S>::IndicMultiplication3D(std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g)
  : IndicComb3D<S,util::multiplies>(f, g)
{
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = util::max(this->_f->getMin()[i], this->_g->getMin()[i]);
    this->_myMax[i] = util::min(this->_f->getMax()[i], this->_g->getMax()[i]);
  }
}

template <typename S>
S IndicMultiplication3D<S>::signedDistance( const Vector<S,3>& input )
{
  return sdf::intersection(this->_f->signedDistance(input), this->_g->signedDistance(input));
}


//// no association to a operator+ from a class is needed, thus we have these free functions
template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator+(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs)
{
  return std::make_shared<IndicPlus3D<S>>(lhs, rhs);
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator-(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs)
{
  return std::make_shared<IndicMinus3D<S>>(rhs, lhs);
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator*(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs)
{
  return std::make_shared<IndicMultiplication3D<S>>(lhs, rhs);
}

// template specialization for indicatorIdentity
template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator+(F1<S> & lhs, std::shared_ptr<F2<S>> rhs)
{
  return lhs._f + rhs;
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator-(F1<S> & lhs, std::shared_ptr<F2<S>> rhs)
{
  return lhs._f - rhs;
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator*(F1<S> & lhs, std::shared_ptr<F2<S>> rhs)
{
  return lhs._f * rhs;
}


} // end namespace olb

#endif