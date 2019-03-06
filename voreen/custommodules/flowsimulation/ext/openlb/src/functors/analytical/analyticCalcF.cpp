/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011-2013 Lukas Baron, Tim Dornieden, Mathias J. Krause
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

#include "analyticCalcF.h"
#include "analyticCalcF.hh"

namespace olb {


// arithmetic helper class for analytical 1d functors

template class AnalyticCalc1D<double,int>;
template class AnalyticCalc1D<double,double>;

template class AnalyticPlus1D<double,int>;
template class AnalyticPlus1D<double,double>;

template class AnalyticMinus1D<double,int>;
template class AnalyticMinus1D<double,double>;

template class AnalyticMultiplication1D<double,int>;
template class AnalyticMultiplication1D<double,double>;

template class AnalyticDivision1D<double,int>;
template class AnalyticDivision1D<double,double>;

// arithmetic helper class for analytical 2d functors

template class AnalyticCalc2D<double,int>;
template class AnalyticCalc2D<double,double>;

template class AnalyticPlus2D<double,int>;
template class AnalyticPlus2D<double,double>;

template class AnalyticMinus2D<double,int>;
template class AnalyticMinus2D<double,double>;

template class AnalyticMultiplication2D<double,int>;
template class AnalyticMultiplication2D<double,double>;

template class AnalyticDivision2D<double,int>;
template class AnalyticDivision2D<double,double>;


// arithmetic helper class for analytical 3d functors

template class AnalyticCalc3D<double,int>;
template class AnalyticCalc3D<double,double>;

template class AnalyticPlus3D<double,int>;
template class AnalyticPlus3D<double,double>;

template class AnalyticMinus3D<double,int>;
template class AnalyticMinus3D<double,double>;

template class AnalyticMultiplication3D<double,int>;
template class AnalyticMultiplication3D<double,double>;

template class AnalyticDivision3D<double,int>;
template class AnalyticDivision3D<double,double>;

}

