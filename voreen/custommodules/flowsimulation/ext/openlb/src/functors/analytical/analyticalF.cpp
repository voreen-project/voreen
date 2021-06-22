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

#include "analyticalF.h"
#include "analyticalF.hh"

namespace olb {

template class AnalyticalComposed<2,double,int>;
template class AnalyticalComposed<3,double,int>;
template class AnalyticalComposed<2,double,double>;
template class AnalyticalComposed<3,double,double>;

template class AnalyticalConst<1,double,int>;
template class AnalyticalConst<2,double,int>;
template class AnalyticalConst<3,double,int>;
template class AnalyticalConst<1,double,double>;
template class AnalyticalConst<2,double,double>;
template class AnalyticalConst<3,double,double>;

template class AnalyticalNormal<1,double,int>;
template class AnalyticalNormal<2,double,int>;
template class AnalyticalNormal<3,double,int>;
template class AnalyticalNormal<1,double,double>;
template class AnalyticalNormal<2,double,double>;
template class AnalyticalNormal<3,double,double>;

template class AnalyticalRandomBase<1,double,int>;
template class AnalyticalRandomBase<2,double,int>;
template class AnalyticalRandomBase<3,double,int>;
template class AnalyticalRandomBase<1,double,double>;
template class AnalyticalRandomBase<2,double,double>;
template class AnalyticalRandomBase<3,double,double>;

template class AnalyticalRandomUniform<1,double,int>;
template class AnalyticalRandomUniform<2,double,int>;
template class AnalyticalRandomUniform<3,double,int>;
template class AnalyticalRandomUniform<1,double,double>;
template class AnalyticalRandomUniform<2,double,double>;
template class AnalyticalRandomUniform<3,double,double>;

template class AnalyticalRandomNormal<1,double,int>;
template class AnalyticalRandomNormal<2,double,int>;
template class AnalyticalRandomNormal<3,double,int>;
template class AnalyticalRandomNormal<1,double,double>;
template class AnalyticalRandomNormal<2,double,double>;
template class AnalyticalRandomNormal<3,double,double>;

template class AnalyticalRandomOld<1,double,int>;
template class AnalyticalRandomOld<2,double,int>;
template class AnalyticalRandomOld<3,double,int>;
template class AnalyticalRandomOld<1,double,double>;
template class AnalyticalRandomOld<2,double,double>;
template class AnalyticalRandomOld<3,double,double>;

////////////// OLD IMPLEMENTATION //////////////////////////

template class AnalyticalLinear1D<double,int>;
template class AnalyticalLinear1D<double,double>;

template class AnalyticalSquare1D<double,int>;
template class AnalyticalSquare1D<double,double>;

template class PolynomialStartScale<double,int>;
template class PolynomialStartScale<double,double>;

template class SinusStartScale<double,int>;
template class SinusStartScale<double,double>;

template class AnalyticalDiffFD1D<double>;

template class Cosinus<double,int>;
template class Cosinus<double,double>;

template class CosinusComposite<double, int>;
template class CosinusComposite<double, double>;




template class AnalyticalLinear2D<double,int>;
template class AnalyticalLinear2D<double,double>;


template class AnalyticalLinear3D<double,int>;
template class AnalyticalLinear3D<double,double>;

template class AnalyticalScaled3D<double,int>;
template class AnalyticalScaled3D<double,double>;

}

