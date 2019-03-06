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

#ifndef INDIC_CALC_F_3D_H
#define INDIC_CALC_F_3D_H

#include "indicatorBaseF3D.h"

namespace olb {


/*
 *  arithmetic helper classes for IndicatorF3D, smoothIndicator3D
 *  UNION         +
 *  WITHOUT       -
 *  INTERSECTION  *
*/

//////////////////////////////// IndicCalc3D ////////////////////////////////
/// arithmetic helper class for Indicator 3d functors
template <typename S>
class IndicCalc3D : public IndicatorF3D<S> {
protected:
  IndicatorF3D<S>& _f;
  IndicatorF3D<S>& _g;
public:
  IndicCalc3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g);
};

/// addition functor acts as union
template <typename S>
class IndicPlus3D : public IndicCalc3D<S> {
public:
  IndicPlus3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};

/// subtraction functor acts as without
template <typename S>
class IndicMinus3D : public IndicCalc3D<S> {
public:
  IndicMinus3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};

/// multiplication functor acts as intersection
template <typename S>
class IndicMultiplication3D : public IndicCalc3D<S> {
public:
  IndicMultiplication3D(IndicatorF3D<S>& f, IndicatorF3D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};


} // end namespace olb

#endif
