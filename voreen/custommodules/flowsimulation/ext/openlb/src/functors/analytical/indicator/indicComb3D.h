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

#ifndef INDIC_COMB_3D_H
#define INDIC_COMB_3D_H

#include "indicatorBaseF3D.h"
#include "utilities/arithmetic.h"
#include "sdf.h"

namespace olb {


/*
 *  arithmetic helper classes for IndicatorF3D, smoothIndicator3D
 *  UNION         +
 *  WITHOUT       -
 *  INTERSECTION  *
*/

//////////////////////////////// IndicComb3D ////////////////////////////////
/// arithmetic helper class for Indicator 3d functors
template <typename S, template<typename U> class F>
class IndicComb3D : public IndicatorF3D<S> {
protected:
  IndicComb3D( std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g );
  std::shared_ptr<IndicatorF3D<S>> _f;
  std::shared_ptr<IndicatorF3D<S>> _g;
public:
  virtual S signedDistance( const Vector<S,3>& input )=0;
  bool operator() (bool output[], const S input[3]);
};

/// Union
template <typename S>
class IndicPlus3D : public IndicComb3D<S,util::plus> {
public:
  IndicPlus3D( std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g );
  S signedDistance( const Vector<S,3>& input ) override;
};

/// Subtraction
template <typename S>
class IndicMinus3D : public IndicComb3D<S,util::minus> {
public:
  IndicMinus3D( std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g );
  S signedDistance( const Vector<S,3>& input ) override;
};

/// Intersection
template <typename S>
class IndicMultiplication3D : public IndicComb3D<S,util::multiplies> {
public:
  IndicMultiplication3D( std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g );
  S signedDistance( const Vector<S,3>& input ) override;
};

/** Free function implements lhs+rhs, only for IndicaotrsF3D types through enable_if and is_base_of
 *
 * \tparam S usual type for source dimension of the functor
 * \tparam F1 lhs has to be derived from IndicatorF3D, otherwise function is disabled
 * \tparam F2 rhs
 */
template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorF3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator+(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorF3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator-(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorF3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator*(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorIdentity3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator+(F1<S> & lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorIdentity3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator-(F1<S> & lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorIdentity3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator*(F1<S> & lhs, std::shared_ptr<F2<S>> rhs);
} // end namespace olb

#endif