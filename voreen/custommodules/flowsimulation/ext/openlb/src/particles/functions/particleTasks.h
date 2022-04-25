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


/* This file containes particle tasks that can be passed to the particle manager.
 * Those need to provide an execute method, and two booleans specifying the coupling
 * and the necessity for beeing looped over all particles.
*/

#ifndef PARTICLE_TASKS_H
#define PARTICLE_TASKS_H


#include <cassert>

namespace olb {

namespace particles {

/// Couple lattice to particles
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE,
         typename FORCEFUNCTOR=SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>>
struct couple_lattice_to_particles{
  auto static execute(
    ParticleSystem<T,PARTICLETYPE>& particleSystem,
    SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
    SuperLattice<T,DESCRIPTOR>& sLattice,
    UnitConverter<T,DESCRIPTOR> const& converter,
    Vector<bool,DESCRIPTOR::d> periodicity = Vector<bool,DESCRIPTOR::d> (false),
    std::size_t iP0=0 )
  {
    //Momentum Exchange/Loss
    FORCEFUNCTOR forceFunctor( sLattice, sGeometry,
                               particleSystem, converter, periodicity, iP0 );
    //Store boundary force in particles
    dynamics::applyParticleForce<T,PARTICLETYPE>( forceFunctor, particleSystem, iP0 );
  }
  static constexpr bool latticeCoupling = true;
  static constexpr bool particleLoop = false;
};

/// Apply gravity
template<typename T, typename PARTICLETYPE>
struct apply_gravity{
  auto static execute(
    ParticleSystem<T,PARTICLETYPE>& particleSystem,
    Vector<T,PARTICLETYPE::d> externalAcceleration,
    int iP, T timeStepSize)
  {
    using namespace descriptors;
    auto particle = particleSystem.get(iP);
    //Apply acceleration
    Vector<T,PARTICLETYPE::d> force = particle.template getField<FORCING,FORCE>();
    T mass = particle.template getField<PHYSPROPERTIES,MASS>();
    force += externalAcceleration * mass;
    particle.template setField<FORCING,FORCE>( force );
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = true;
};

/// Process particle dynamics
template<typename T, typename PARTICLETYPE>
struct process_dynamics{
  auto static execute(
    ParticleSystem<T,PARTICLETYPE>& particleSystem,
    Vector<T,PARTICLETYPE::d>& externalAcceleration,
    int iP, T timeStepSize )
  {
    auto particle = particleSystem.get(iP);
    //Call process on particle dynamics
    particle.process(timeStepSize);
  }
  static constexpr bool latticeCoupling = false;
  static constexpr bool particleLoop = true;
};

/// Couple particles to lattice
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
struct couple_particles_to_lattice{
  auto static execute(
    ParticleSystem<T,PARTICLETYPE>& particleSystem,
    SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
    SuperLattice<T,DESCRIPTOR>& sLattice,
    UnitConverter<T,DESCRIPTOR> const& converter,
    int iP )
  {
    using namespace descriptors;
    //Retrieve particle
    auto particle = particleSystem.get(iP);
    //Create eccentric velocity field
    EccentricLatticeVelocityField<T,T,DESCRIPTOR> eccentricVelocityField(
      particle.template getField<GENERAL,POSITION>(),
      particle.template getField<MOBILITY,VELOCITY>(),
      particle.template getField<MOBILITY,ANG_VELOCITY>(),
      converter
    );
    //Write particle field
    setSuperParticleField( sGeometry, eccentricVelocityField, sLattice, particle );
  }
  static constexpr bool latticeCoupling = true;
  static constexpr bool particleLoop = true;
};


} //namespace particles

} //namespace olb


#endif
