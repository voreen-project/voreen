/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen
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


#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H


namespace olb {

namespace particles {


template<typename T, typename PARTICLETYPE>
class ParticleSystem {
private:
  /// Main data storage
  Container<T,PARTICLETYPE,DynamicFieldGroupsD<T, typename PARTICLETYPE::fieldList>> _fieldGroupsContainer;
  /// Danamics Map
  Container<T,PARTICLETYPE,OldBlockDynamicsMap<T,PARTICLETYPE,dynamics::ParticleDynamics<T,PARTICLETYPE>>> _dynamicsMapContainer;
  //Additional storage for multiple complex data types linked to particle (extendable by user...)
  using associatedTypes = meta::list<
    std::vector<std::shared_ptr<SmoothIndicatorF2D<T,T,true>>>,
    std::vector<std::shared_ptr<SmoothIndicatorF3D<T,T,true>>>
  >;
  typename associatedTypes::template decompose_into<std::tuple> _associatedData; //Storage for complex types

public:
  static_assert(PARTICLETYPE::d == 2 || PARTICLETYPE::d == 3, "Only D=2 and D=3 are supported");

  /// Default constructor
  ParticleSystem();

  /// Constructor with initial particle size
  ParticleSystem( std::size_t count );

  /// Constructor with initial particle size and associated Data
  ParticleSystem( std::size_t count,
                  typename associatedTypes::template decompose_into<std::tuple> associatedData );

  /// Size of ParticleSystem
  constexpr std::size_t size();

  /// Get the dynamics for specific particle
  dynamics::ParticleDynamics<T,PARTICLETYPE>* getDynamics(int iP);

  /// Define the dynamics for specific particle
  void defineDynamics(int iP, dynamics::ParticleDynamics<T,PARTICLETYPE>* dynamics);

  /// Extend particle system by one particle
  void extend();
 
  /// Swap particles by index 
  void swapParticles(std::size_t iP, std::size_t jP);

  /// Upate/process particle quantities on a subset (mostly local equivalent to collide)
  void process( std::size_t p0, std::size_t p1, T timeStepSize );

  /// Upate/process particle quantities on all particles
  void process(T timeStepSize);

  /// Update particles (equivalent to stream())
  void update();

  /// Expose container
  Container<T,PARTICLETYPE,DynamicFieldGroupsD<T,typename PARTICLETYPE::fieldList>>& get();

  /// Create and return particle (analogous to cell)
  Particle<T,PARTICLETYPE> get(std::size_t iParticle);

  /// Get whole Field by GROUP (or base of GROUP) and FIELD
  template <typename GROUP, typename FIELD>
  auto& getFieldD();

  /// Get FieldPointer by GROUP (or base of GROUP) and FIELD for specific iParticle
  template <typename GROUP, typename FIELD>
  auto getFieldPointer( std::size_t iParticle );

  /// Get associated data by specifying data type
  template<typename TYPE>
  auto& getAssociatedData();

};


} //namespace particles

} //namespace olb


#endif
