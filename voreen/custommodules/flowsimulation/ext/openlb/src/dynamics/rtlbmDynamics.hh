/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink, Christopher McHardy
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
 * A collection of radiative transport dynamics classes -- generic implementation.
 */

#ifndef RTLBM_DYNAMICS_HH
#define RTLBM_DYNAMICS_HH

#include "rtlbmDynamics.h"
#include "lbHelpers.h"

using namespace olb::descriptors;

namespace olb {



//==================================================================//
//============= BGK Model for Advection diffusion plus sink term ===//
//==================================================================//

template<typename T, template<typename U> class Lattice>
RTLBMdynamicsMink<T,Lattice>::RTLBMdynamicsMink
( T omega, Momenta<T,Lattice>& momenta, T latticeAbsorption, T latticeScattering )
  : BasicDynamics<T, Lattice>( momenta ), _omega(omega), _sink( 3.0*latticeAbsorption*(latticeAbsorption+latticeScattering) / 8.0)
{
  static_assert( std::is_base_of<D3Q7DescriptorBaseRTLBM<T>, Lattice<T> >::value, "Descriptor not derived from D3Q7DescriptorBase.");
}

template<typename T, template<typename U> class Lattice>
T RTLBMdynamicsMink<T, Lattice>::computeEquilibrium
( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const
{
  return lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, template<typename U> class Lattice>
void RTLBMdynamicsMink<T, Lattice>::collide
( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics )
{
  T intensity = this->_momenta.computeRho( cell );
  T uSqr = lbHelpers<T,Lattice>::sinkCollision( cell, intensity, _omega, _sink );
  statistics.incrementStats( intensity, uSqr );
}

template<typename T, template<typename U> class Lattice>
void RTLBMdynamicsMink<T, Lattice>::staticCollide( Cell<T, Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  assert( false );
}

template<typename T, template<typename U> class Lattice>
T RTLBMdynamicsMink<T, Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void RTLBMdynamicsMink<T, Lattice>::setOmega( T omega )
{
  _omega = omega;
}

template<typename T, template<typename U> class Lattice>
T RTLBMdynamicsMink<T, Lattice>::getSink() const
{
  return _sink;
}


template<typename T, template<typename U> class Lattice>
RTLBMconstDynamicsMink<T,Lattice>::RTLBMconstDynamicsMink
( Momenta<T,Lattice>& momenta, T latticeAbsorption, T latticeScattering )
  : BasicDynamics<T, Lattice>( momenta ), _sink( 3.0*latticeAbsorption*(latticeAbsorption+latticeScattering) / 8.0)
{
  constexpr bool is_d3q7_descriptor = std::is_base_of<D3Q7DescriptorBaseRTLBM<T>, Lattice<T> >::value
                               || std::is_base_of<D3Q7DescriptorBase<T>, Lattice<T> >::value;

  static_assert( is_d3q7_descriptor, "Descriptor not derived from D3Q7DescriptorBase.");
}

template<typename T, template<typename U> class Lattice>
T RTLBMconstDynamicsMink<T, Lattice>::computeEquilibrium
( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const
{
  return lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, template<typename U> class Lattice>
void RTLBMconstDynamicsMink<T, Lattice>::collide
( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics )
{
  T intensity = this->_momenta.computeRho( cell );
  T uSqr = lbHelpers<T,Lattice>::sinkCollision( cell, intensity, 1, _sink );
  statistics.incrementStats( intensity, uSqr );
}

template<typename T, template<typename U> class Lattice>
void RTLBMconstDynamicsMink<T, Lattice>::staticCollide( Cell<T, Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  assert( false );
}

template<typename T, template<typename U> class Lattice>
T RTLBMconstDynamicsMink<T, Lattice>::getOmega() const
{
  return 1;
}

template<typename T, template<typename U> class Lattice>
void RTLBMconstDynamicsMink<T, Lattice>::setOmega( T omega )
{}

//==================================================================//
//============= BGK Model for Advection diffusion anisotropic ===//
//==================================================================//

template<typename T, template<typename U> class Lattice>
RTLBMdynamicsMcHardy<T, Lattice>::RTLBMdynamicsMcHardy
(Momenta<T, Lattice>& momenta, T latticeAbsorption, T latticeScattering)
  : BasicDynamics<T, Lattice>(momenta), _absorption(latticeAbsorption), _scattering(latticeScattering)
{ }

template<typename T, template<typename U> class Lattice>
T RTLBMdynamicsMcHardy<T, Lattice>::computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const
{
  return lbHelpers<T,Lattice>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, template<typename U> class Lattice>
void RTLBMdynamicsMcHardy<T, Lattice>::collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics )
{
  T temperature = this->_momenta.computeRho(cell );
//  T uSqr = advectionDiffusionLbHelpers<T,Lattice>::sinkCollision( cell, temperature, omega, 0. );
  T uSqr = lbHelpers<T, Lattice>::
           anisoCollision( cell, temperature, _absorption, _scattering );
  statistics.incrementStats( temperature, uSqr );
}

template<typename T, template<typename U> class Lattice>
void RTLBMdynamicsMcHardy<T, Lattice>::staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  assert( false );
}

template<typename T, template<typename U> class Lattice>
T RTLBMdynamicsMcHardy<T, Lattice>::getOmega() const
{
  return -1;
}

template<typename T, template<typename U> class Lattice>
void RTLBMdynamicsMcHardy<T, Lattice>::setOmega( T omega )
{
}

//==================================================================================//
template<typename T, template<typename U> class Lattice>
RTLBMdynamicsMcHardyWH<T, Lattice>::RTLBMdynamicsMcHardyWH
(Momenta<T, Lattice>& momenta, T latticeAbsorption, T latticeScattering)
  : RTLBMdynamicsMcHardy<T, Lattice>(momenta, latticeAbsorption, latticeScattering)
{ }

template<typename T, template<typename U> class Lattice>
T RTLBMdynamicsMcHardyWH<T, Lattice>::computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const
{
  return lbHelpers<T,Lattice>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, template<typename U> class Lattice>
void RTLBMdynamicsMcHardyWH<T, Lattice>::collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics )
{
  T temperature = this->_momenta.computeRho(cell );
//  T uSqr = advectionDiffusionLbHelpers<T,Lattice>::sinkCollision( cell, temperature, omega, 0. );
  T uSqr = lbHelpers<T, Lattice>::
           anisoCollisionWH( cell, temperature, this->_absorption, this->_scattering );
  statistics.incrementStats( temperature, uSqr );
}

template<typename T, template<typename U> class Lattice>
void RTLBMdynamicsMcHardyWH<T, Lattice>::staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  assert( false );
}

template<typename T, template<typename U> class Lattice>
T RTLBMdynamicsMcHardyWH<T, Lattice>::getOmega() const
{
  return -1;
}

template<typename T, template<typename U> class Lattice>
void RTLBMdynamicsMcHardyWH<T, Lattice>::setOmega( T omega )
{
}


} // namespace olb


#endif
