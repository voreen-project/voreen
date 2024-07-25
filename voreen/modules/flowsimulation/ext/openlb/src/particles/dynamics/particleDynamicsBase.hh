/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef PARTICLE_DYNAMICS_BASE_HH
#define PARTICLE_DYNAMICS_BASE_HH


namespace olb {

namespace particles {

namespace dynamics {

template <typename T, typename PARTICLETYPE>
std::string& ParticleDynamics<T,PARTICLETYPE>::getName()
{
  return _name;
}

template <typename T, typename PARTICLETYPE>
std::string const& ParticleDynamics<T,PARTICLETYPE>::getName() const
{
  return _name;
}


template<typename T, typename PARTICLETYPE>
NoParticleDynamics<T,PARTICLETYPE>::NoParticleDynamics( T rhoDummy )
{
  this->getName() = "NoParticleDynamics";
}

template<typename T, typename PARTICLETYPE>
void NoParticleDynamics<T,PARTICLETYPE>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize)
{ }



template<typename T, typename PARTICLETYPE>
VerletParticleDynamics<T,PARTICLETYPE>::VerletParticleDynamics ( )
{
  this->getName() = "VerletParticleDynamics";
}

template<typename T, typename PARTICLETYPE>
void VerletParticleDynamics<T,PARTICLETYPE>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Calculate acceleration
  auto acceleration = getAcceleration<T,PARTICLETYPE>( particle );
  //Check for angular components
  if constexpr ( providesAngle<PARTICLETYPE>() ) {
    //Calculate angular acceleration
    auto angularAcceleration = getAngAcceleration<T,PARTICLETYPE>( particle );
    //Verlet algorithm
    particles::dynamics::velocityVerletIntegration<T, PARTICLETYPE>(
      particle, timeStepSize, acceleration, angularAcceleration );
    //Check if rotation matrix provided
    if constexpr ( providesRotationMatrix<PARTICLETYPE>() ) {
      //Update rotation matrix
      updateRotationMatrix<T,PARTICLETYPE>( particle );
    }
  }
  else {
    //Verlet algorithm without rotation
    particles::dynamics::velocityVerletIntegration<T, PARTICLETYPE>(
      particle, timeStepSize, acceleration );
  }
}



template<typename T, typename PARTICLETYPE>
VerletParticleDynamicsCubicBoundsCheck<T,PARTICLETYPE>::VerletParticleDynamicsCubicBoundsCheck (
  PhysR<T,PARTICLETYPE::d>& domainMin,
  PhysR<T,PARTICLETYPE::d>& domainMax
)
  : _domainMin(domainMin), _domainMax(domainMax)
{
  this->getName() = "VerletParticleDynamicsCubicBoundsCheck";
}

template<typename T, typename PARTICLETYPE>
void VerletParticleDynamicsCubicBoundsCheck<T,PARTICLETYPE>::process (
  Particle<T,PARTICLETYPE>& particle, T timeStepSize )
{
  //Calculate acceleration
  auto acceleration = getAcceleration<T,PARTICLETYPE>( particle );
  //Note position before calculating movement
  auto positionPre = getPosition<T,PARTICLETYPE>( particle );
  //Check for angular components
  if constexpr ( providesAngle<PARTICLETYPE>() ) {
    //Calculate angular acceleration
    auto angularAcceleration = getAngAcceleration<T,PARTICLETYPE>( particle );
    //Verlet algorithm
    particles::dynamics::velocityVerletIntegration<T, PARTICLETYPE>(
      particle, timeStepSize, acceleration, angularAcceleration );
    //Update rotation matrix
    updateRotationMatrix<T,PARTICLETYPE>( particle );
  }
  else {
    //Verlet algorithm without angle
    particles::dynamics::velocityVerletIntegration<T, PARTICLETYPE>(
      particle, timeStepSize, acceleration );
  }

  //Check domain contact and potentially reset to positionPre and velocity zero
  doAtCubicBoundPenetration<T,PARTICLETYPE>( particle, _domainMin, _domainMax,
  [&](unsigned iDim, Vector<T,PARTICLETYPE::d>& normal, T distToBound ) {
    //Reset (only!) penetration direction
    resetDirection( particle, positionPre, iDim );
  });
}

} //namespace dynamics
  
} //namespace particles

} //namespace olb

#endif
