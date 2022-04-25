/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Max Gaedtke, Albert Mink, Davide Dapelo
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

/** \file
 * Function to create a unitConverter with fractional timestep (and lattice size).
 */

#ifndef FRACTIONAL_UNITCONVERTER_H
#define FRACTIONAL_UNITCONVERTER_H

#include "unitConverter.h"

namespace olb {

template <typename T, typename DESCRIPTOR, typename DESCRIPTOR_AD>
UnitConverter<T,DESCRIPTOR> createADfractionalUnitConverter (
  const UnitConverter<T,DESCRIPTOR>& converter,
  int fraction, T targetLatticeRelaxationTime )
{
  if (fraction <= 0) {
    throw std::out_of_range("fracton must be positive.");
  }

  T conversionViscosityOrDiffusivity = converter.getPhysDeltaX() * converter.getPhysDeltaX() / (converter.getPhysDeltaT() * fraction);
  T physViscosiyOrDiffusivity = (targetLatticeRelaxationTime - 0.5) * conversionViscosityOrDiffusivity / descriptors::invCs2<T,DESCRIPTOR_AD>();

  return UnitConverter<T,DESCRIPTOR> (
           converter.getPhysDeltaX() / fraction,
           converter.getPhysDeltaT() / fraction,
           converter.getCharPhysLength(),
           converter.getCharPhysVelocity(),
           physViscosiyOrDiffusivity,
           converter.getPhysDensity(),
           converter.getCharPhysPressure()
         );
}

template <typename T, typename DESCRIPTOR>
UnitConverter<T,DESCRIPTOR> createFractionalUnitConverter (
  const UnitConverter<T,DESCRIPTOR>& converter,
  int fraction, T targetLatticeRelaxationTime )
{
  return createADfractionalUnitConverter<T,DESCRIPTOR,DESCRIPTOR> (
           converter, fraction, targetLatticeRelaxationTime );
}

template <typename T, typename DESCRIPTOR>
T residualPhysDiffusivity ( const UnitConverter<T,DESCRIPTOR>& converterFractional, T physDiffusivity )
{
  return physDiffusivity - converterFractional.getPhysViscosity();
}

} // olb

#endif
