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

#ifndef ANALYTICAL_CALC_F_H
#define ANALYTICAL_CALC_F_H


#include "analyticalBaseF.h"


namespace olb {

/*
    arithmetic helper classes for AnalyticalF1D, AnalyticalF3D, AnalyticalF3D

    pointwise: difference, plus, multiplication, division

*/

//////////////////////////////// AnalyticCalc1D ////////////////////////////////

/// arithmetic helper class for analytical 1D functors
template <typename T, typename S>
class AnalyticCalc1D : public AnalyticalF1D<T,S> {
protected:
  AnalyticalF1D<T,S>& _f;
  AnalyticalF1D<T,S>& _g;
public:
  AnalyticCalc1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g);
};

/// addition functor
template <typename T, typename S>
class AnalyticPlus1D : public AnalyticCalc1D<T,S> {
public:
  AnalyticPlus1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};

/// subtraction functor
template <typename T, typename S>
class AnalyticMinus1D : public AnalyticCalc1D<T,S> {
public:
  AnalyticMinus1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};

/// multiplication functor
template <typename T, typename S>
class AnalyticMultiplication1D : public AnalyticCalc1D<T,S> {
public:
  AnalyticMultiplication1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};

/// division functor
template <typename T, typename S>
class AnalyticDivision1D : public AnalyticCalc1D<T,S> {
public:
  AnalyticDivision1D(AnalyticalF1D<T,S>& f, AnalyticalF1D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};



//////////////////////////////// AnalyticCalc2D ////////////////////////////////

/// arithmetic helper class for analytical 2D functors
template <typename T, typename S>
class AnalyticCalc2D : public AnalyticalF2D<T,S> {
protected:
  AnalyticalF2D<T,S>& _f;
  AnalyticalF2D<T,S>& _g;
public:
  AnalyticCalc2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g);
};

/// addition functor
template <typename T, typename S>
class AnalyticPlus2D : public AnalyticCalc2D<T,S> {
public:
  AnalyticPlus2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};

/// subtraction functor
template <typename T, typename S>
class AnalyticMinus2D : public AnalyticCalc2D<T,S> {
public:
  AnalyticMinus2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};

/// multiplication functor
template <typename T, typename S>
class AnalyticMultiplication2D : public AnalyticCalc2D<T,S> {
public:
  AnalyticMultiplication2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};

/// division functor
template <typename T, typename S>
class AnalyticDivision2D : public AnalyticCalc2D<T,S> {
public:
  AnalyticDivision2D(AnalyticalF2D<T,S>& f, AnalyticalF2D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};



//////////////////////////////// AnalyticCalc3D ////////////////////////////////

/// arithmetic helper class for analytical 3D functors
template <typename T, typename S>
class AnalyticCalc3D : public AnalyticalF3D<T,S> {
protected:
  AnalyticalF3D<T,S>& _f;
  AnalyticalF3D<T,S>& _g;
public:
  AnalyticCalc3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g);
};

/// addition functor
template <typename T, typename S>
class AnalyticPlus3D : public AnalyticCalc3D<T,S> {
public:
  AnalyticPlus3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};

/// subtraction functor
template <typename T, typename S>
class AnalyticMinus3D : public AnalyticCalc3D<T,S> {
public:
  AnalyticMinus3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};

/// multiplication functor
template <typename T, typename S>
class AnalyticMultiplication3D : public AnalyticCalc3D<T,S> {
public:
  AnalyticMultiplication3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};

/// division functor
template <typename T, typename S>
class AnalyticDivision3D : public AnalyticCalc3D<T,S> {
public:
  AnalyticDivision3D(AnalyticalF3D<T,S>& f, AnalyticalF3D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};


} // end namespace olb

#endif
