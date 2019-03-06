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

#ifndef INDIC_CALC_F_2D_H
#define INDIC_CALC_F_2D_H

#include "indicatorBaseF2D.h"

namespace olb {


/*
 *  arithmetic helper classes for IndicatorF1D, IndicatorF2D, smoothIndicator2D
 *  UNION         +
 *  WITHOUT       -
 *  INTERSECTION  *
*/

//////////////////////////////// IndicCalc1D ////////////////////////////////
/// arithmetic helper class for Indicator 1d functors
template <typename S>
class IndicCalc1D : public IndicatorF1D<S> {
protected:
  IndicatorF1D<S>& _f;
  IndicatorF1D<S>& _g;
public:
  // set image/target dimensions of IndicCalc1D as well
  IndicCalc1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g);
};

/// addition functor acts as union
template <typename S>
class IndicPlus1D : public IndicCalc1D<S> {
public:
  IndicPlus1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};

/// subtraction functor acts as without
template <typename S>
class IndicMinus1D : public IndicCalc1D<S> {
public:
  IndicMinus1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};

/// multiplication functor acts as intersection
template <typename S>
class IndicMultiplication1D : public IndicCalc1D<S> {
public:
  IndicMultiplication1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};



//////////////////////////////// IndicCalc2D ////////////////////////////////
/// arithmetic helper class for Indicator 2d functors
template <typename S>
class IndicCalc2D : public IndicatorF2D<S> {
protected:
  IndicatorF2D<S>& _f;
  IndicatorF2D<S>& _g;
public:
  IndicCalc2D(IndicatorF2D<S>& f, IndicatorF2D<S>& g);
};

/// addition functor acts as union
template <typename S>
class IndicPlus2D : public IndicCalc2D<S> {
public:
  IndicPlus2D(IndicatorF2D<S>& f, IndicatorF2D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};

/// subtraction functor acts as without
template <typename S>
class IndicMinus2D : public IndicCalc2D<S> {
public:
  IndicMinus2D(IndicatorF2D<S>& f, IndicatorF2D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};

/// multiplication functor acts as intersection
template <typename S>
class IndicMultiplication2D : public IndicCalc2D<S> {
public:
  IndicMultiplication2D(IndicatorF2D<S>& f, IndicatorF2D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};


} // end namespace olb

#endif
