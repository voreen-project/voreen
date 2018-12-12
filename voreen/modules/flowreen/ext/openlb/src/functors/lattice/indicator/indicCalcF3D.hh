/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Mathias J. Krause, Cyril Masquelier,
 *  Benjamin Förster, Albert Mink
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

#ifndef INDIC_CALC_F_3D_HH
#define INDIC_CALC_F_3D_HH

#include "indicCalcF3D.h"

namespace olb {


template <typename S>
IndicCalc3D<S>::IndicCalc3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g)
  : _f(f), _g(g)
{
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = std::min(f.getMin()[i], g.getMin()[i]);
    this->_myMax[i] = std::max(f.getMax()[i], g.getMax()[i]);
  }
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}

template <typename S>
IndicPlus3D<S>::IndicPlus3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g)
  : IndicCalc3D<S>(f, g)
{}

// returns 1 if( f==1 || g==1 ) UNION
template <typename S>
bool IndicPlus3D<S>::operator() (bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] |= tmp;
  return true;
}


template <typename S>
IndicMinus3D<S>::IndicMinus3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g)
  : IndicCalc3D<S>(f, g)
{}

// returns 1 if( f==1 && g==0 ) WITHOUT
template <typename S>
bool IndicMinus3D<S>::operator()(bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] &= !tmp;
  return true;
}


// returns 1 if( f==1 && g==1 ) INTERSECTION
template <typename S>
IndicMultiplication3D<S>::IndicMultiplication3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g)
  : IndicCalc3D<S>(f, g)
{}

template <typename S>
bool IndicMultiplication3D<S>::operator() (bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] &= tmp;
  return true;
}


template <typename S>
IndicatorF3D<S>& IndicatorF3D<S>::operator+(IndicatorF3D<S>& rhs)
{
  auto tmp = std::make_shared< IndicPlus3D<S> >(*this, rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename S>
IndicatorF3D<S>& IndicatorF3D<S>::operator-(IndicatorF3D<S>& rhs)
{
  auto tmp = std::make_shared< IndicMinus3D<S> >(*this, rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename S>
IndicatorF3D<S>& IndicatorF3D<S>::operator*(IndicatorF3D<S>& rhs)
{
  auto tmp = std::make_shared< IndicMultiplication3D<S> >(*this, rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}


} // end namespace olb

#endif
