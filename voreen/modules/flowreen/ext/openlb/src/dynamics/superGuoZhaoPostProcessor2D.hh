/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Davide Dapelo, Mathias J. Krause
 *  OpenLB e-mail contact: info@openlb.net
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

/** \file
 * Class to define the external fields involved
 * in the porous modelling (i.e., porosity, porous conductivity
 * with the addition of a generic body force)
 * -- generic file
 */

#ifndef SUPER_GUO_ZAO_POST_PROCESSOR_2D_HH
#define SUPER_GUO_ZAO_POST_PROCESSOR_2D_HH

#include "dynamics/guoZhaoLbHelpers.h"

namespace olb {

template<typename T, template<typename U> class Lattice, class dynamicsManager>
SuperGuoZhaoInstantiator2D<T, Lattice, dynamicsManager>::SuperGuoZhaoInstantiator2D (
  SuperLattice2D<T, Lattice>& sLattice_) :
  sLattice(sLattice_)
{}

template<typename T, template<typename U> class Lattice, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, Lattice, dynamicsManager>::definePorousFields (
  AnalyticalF2D<T,T>& epsilon_, AnalyticalF2D<T,T>& K_)
{

}

template<typename T, template<typename U> class Lattice, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, Lattice, dynamicsManager>::defineEpsilon (
  SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& epsilon)
{

  sLattice.defineExternalField(sGeometry, material,
                               Lattice<T>::ExternalField::epsilonAt, 1, epsilon );
}

template<typename T, template<typename U> class Lattice, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, Lattice, dynamicsManager>::defineK (
  UnitConverter<T,Lattice> const& converter, SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& K)
{

  AnalyticalConst2D<T,T> normFactor(converter.getConversionFactorLength()*converter.getConversionFactorLength());
  AnalyticalIdentity2D<T,T> KLb(K / normFactor);
  sLattice.defineExternalField(sGeometry, material, Lattice<T>::ExternalField::KAt, 1, KLb);
}

template<typename T, template<typename U> class Lattice, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, Lattice, dynamicsManager>::defineNu (
  UnitConverter<T,Lattice> const& converter, SuperGeometry2D<T>& sGeometry, int material)
{

  AnalyticalConst2D<T,T> nu(converter.getLatticeViscosity());
  sLattice.defineExternalField(sGeometry, material, Lattice<T>::ExternalField::nuAt, 1, nu);
}

template<typename T, template<typename U> class Lattice, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, Lattice, dynamicsManager>::defineBodyForce (
  UnitConverter<T,Lattice> const& converter, SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& BodyForce)
{

  std::vector<T> normFactorValue ( 2,
                                   converter.getConversionFactorLength() / (converter.getConversionFactorTime()*converter.getConversionFactorTime()) );
  AnalyticalConst2D<T,T> normFactor(normFactorValue);
  AnalyticalIdentity2D<T,T> BodyForceLb(BodyForce / normFactor);
  sLattice.defineExternalField(sGeometry, material,
                               Lattice<T>::ExternalField::bodyForceBeginsAt, 2, BodyForceLb);
}

}

#endif
