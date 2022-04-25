/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_BASE_F_HH
#define ANALYTICAL_BASE_F_HH

#include "analyticalBaseF.h"

namespace olb {

// identity to "store results"
template<unsigned D, typename T, typename S>
AnalyticalIdentity<D,T,S>::AnalyticalIdentity(AnalyticalF<D,T,S>& f)
  : AnalyticalF<D,T,S>(f.getTargetDim()), _f(f)
{
  this->getName() = _f.getName();
  // pass through the shared_ptr from _f, e.g. an arithemticClass, to the identity
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template<unsigned D, typename T, typename S>
bool AnalyticalIdentity<D,T,S>::operator()(T output[], const S input[])
{
  _f(output,input);
  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
template <typename T, typename S>
AnalyticalFfromIndicatorF3D<T, S>::AnalyticalFfromIndicatorF3D(IndicatorF3D<T>& indicatorF)
  : AnalyticalF3D<T,S>(1), _indicatorF(indicatorF)
{
  this->getName() = "IndicatorFfrom" + _indicatorF.getName();
}

template <typename T, typename S>
bool AnalyticalFfromIndicatorF3D<T, S>::operator() (T output[], const S input[])
{
  bool tmp = false;
  _indicatorF(&tmp, input);
  if ( tmp ) {
    output[0] = T(1);
  }
  else {
    output[0] = T(0);
  }
  return tmp;
}



} // end namespace olb

#endif
