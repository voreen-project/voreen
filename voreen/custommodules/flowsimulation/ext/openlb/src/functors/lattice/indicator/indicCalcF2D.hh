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

#ifndef INDIC_CALC_F_2D_HH
#define INDIC_CALC_F_2D_HH

#include "indicCalcF2D.h"

namespace olb {


//////////////////////////////// IndicCalc1D ////////////////////////////////
template <typename S>
IndicCalc1D<S>::IndicCalc1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g)
  : _f(f), _g(g)
{
  this->_myMin[0] = std::min(f.getMin()[0], g.getMin()[0]);
  this->_myMax[0] = std::max(f.getMax()[0], g.getMax()[0]);
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}


template <typename S>
IndicPlus1D<S>::IndicPlus1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g)
  : IndicCalc1D<S>(f, g)
{}

// returns 1 if( f==1 || g==1 ) UNION
template <typename S>
bool IndicPlus1D<S>::operator() (bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] |= tmp;
  return true;
}


template <typename S>
IndicMinus1D<S>::IndicMinus1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g)
  : IndicCalc1D<S>(f, g)
{}

// returns 1 if( f==1 && g==0 ) WITHOUT
template <typename S>
bool IndicMinus1D<S>::operator()(bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] &= !tmp;
  return true;
}



template <typename S>
IndicMultiplication1D<S>::IndicMultiplication1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g)
  : IndicCalc1D<S>(f, g)
{}

// returns 1 if( f==1 && g==1 ) INTERSECTION
template <typename S>
bool IndicMultiplication1D<S>::operator() (bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] &= tmp;
  return true;
}



template <typename S>
IndicatorF1D<S>& IndicatorF1D<S>::operator+(IndicatorF1D<S>& rhs)
{
  auto tmp = std::make_shared< IndicPlus1D<S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename S>
IndicatorF1D<S>& IndicatorF1D<S>::operator-(IndicatorF1D<S>& rhs)
{
  auto tmp = std::make_shared< IndicMinus1D<S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename S>
IndicatorF1D<S>& IndicatorF1D<S>::operator*(IndicatorF1D<S>& rhs)
{
  auto tmp = std::make_shared< IndicMultiplication1D<S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}




//////////////////////////////// IndicCalc2D ////////////////////////////////
template <typename S>
IndicCalc2D<S>::IndicCalc2D(IndicatorF2D<S>& f, IndicatorF2D<S>& g)
  : _f(f), _g(g)
{
  for ( int i=0; i<2; i++) {
    this->_myMin[i] = std::min(f.getMin()[i], g.getMin()[i]);
    this->_myMax[i] = std::max(f.getMax()[i], g.getMax()[i]);
  }
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}


template <typename S>
IndicPlus2D<S>::IndicPlus2D(IndicatorF2D<S>& f, IndicatorF2D<S>& g)
  : IndicCalc2D<S>(f, g)
{}

// returns 1 if( f==1 || g==1 ) UNION
template <typename S>
bool IndicPlus2D<S>::operator() (bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] |= tmp;
  return true;
}


template <typename S>
IndicMinus2D<S>::IndicMinus2D(IndicatorF2D<S>& f, IndicatorF2D<S>& g)
  : IndicCalc2D<S>(f, g) {}

// returns 1 if( f==1 && g==0 ) WITHOUT
template <typename S>
bool IndicMinus2D<S>::operator()(bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] &= !tmp;
  return true;
}



template <typename S>
IndicMultiplication2D<S>::IndicMultiplication2D(IndicatorF2D<S>& f,
    IndicatorF2D<S>& g) : IndicCalc2D<S>(f, g) {}

// returns 1 if( f==1 && g==1 ) INTERSECTION
template <typename S>
bool IndicMultiplication2D<S>::operator() (bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] &= tmp;
  return true;
}



template <typename S>
IndicatorF2D<S>& IndicatorF2D<S>::operator+(IndicatorF2D<S>& rhs)
{
  auto tmp = std::make_shared< IndicPlus2D<S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename S>
IndicatorF2D<S>& IndicatorF2D<S>::operator-(IndicatorF2D<S>& rhs)
{
  auto tmp = std::make_shared< IndicMinus2D<S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename S>
IndicatorF2D<S>& IndicatorF2D<S>::operator*(IndicatorF2D<S>& rhs)
{
  auto tmp = std::make_shared< IndicMultiplication2D<S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}


} // end namespace olb

#endif
