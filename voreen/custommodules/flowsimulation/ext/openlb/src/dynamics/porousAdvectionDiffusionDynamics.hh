/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef POROUS_ADVECTION_DIFFUSION_DYNAMICS_HH
#define POROUS_ADVECTION_DIFFUSION_DYNAMICS_HH

#include <algorithm>
#include <limits>

#include "porousAdvectionDiffusionDynamics.h"
#include "core/util.h"
#include "lbm.h"

namespace olb {


//==================================================================//
//========== BGK Model for porous Advection diffusion=======//
//==================================================================//

template<typename T, typename DESCRIPTOR, typename MOMENTA>
PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR, MOMENTA>::PorousAdvectionDiffusionBGKdynamics (
  T omega, T tSolid )
  : legacy::BasicDynamics<T, DESCRIPTOR, MOMENTA>( ),
    _omega(omega), _tSolid(tSolid)
{
  this->getName() = "PorousAdvectionDiffusionBGKdynamics";
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR, MOMENTA>::PorousAdvectionDiffusionBGKdynamics (
  const UnitConverter<T,DESCRIPTOR>& converter, T tSolid )
  : legacy::BasicDynamics<T, DESCRIPTOR, MOMENTA>( ),
    _omega(converter.getLatticeRelaxationFrequency()), _tSolid(tSolid)
{
  this->getName() = "PorousAdvectionDiffusionBGKdynamics";
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR, MOMENTA>::computeEquilibrium( int iPop, T rho,
    const T u[DESCRIPTOR::d] ) const
{
  // does temperature need to be considered here?
  return equilibrium<DESCRIPTOR>::firstOrder( iPop, rho, u );
}


template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR, MOMENTA>::collide( Cell<T,DESCRIPTOR>& cell )
{
  T temperature = MomentaF().computeRho( cell );

  // apply temperature scaling
  auto porosity = cell.template getField<descriptors::POROSITY>();
  temperature = scaleTemp(temperature, porosity);

  auto u = cell.template getField<descriptors::VELOCITY>();

  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, temperature, u, _omega);

  return {temperature, uSqr};
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR, MOMENTA>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR, MOMENTA>::setOmega( T omega )
{
  _omega = omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR, MOMENTA>::scaleTemp(const T temp, const T porosity) const
{
  T rho = porosity * temp + ( T(1) - porosity ) * _tSolid;
  return rho;
}


} // namespace olb

#endif
