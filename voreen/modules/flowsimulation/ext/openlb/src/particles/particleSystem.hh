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


#ifndef PARTICLE_SYSTEM_HH
#define PARTICLE_SYSTEM_HH



namespace olb {

namespace particles{

template<typename T, typename PARTICLETYPE>
ParticleSystem<T, PARTICLETYPE>::ParticleSystem()
  : _fieldGroupsContainer( Container<T,PARTICLETYPE,DynamicFieldGroupsD<T, typename PARTICLETYPE::fieldList>>(0) ),
    _dynamicsMapContainer( Container<T,PARTICLETYPE,OldBlockDynamicsMap<T,PARTICLETYPE,dynamics::ParticleDynamics<T,PARTICLETYPE>>>(0) )
{}

template<typename T, typename PARTICLETYPE>
ParticleSystem<T, PARTICLETYPE>::ParticleSystem( std::size_t count )
  : _fieldGroupsContainer( Container<T,PARTICLETYPE,DynamicFieldGroupsD<T, typename PARTICLETYPE::fieldList>>(count) ),
    _dynamicsMapContainer( Container<T,PARTICLETYPE,OldBlockDynamicsMap<T,PARTICLETYPE,dynamics::ParticleDynamics<T,PARTICLETYPE>>>(count) )
{}

template<typename T, typename PARTICLETYPE>
ParticleSystem<T, PARTICLETYPE>::ParticleSystem( std::size_t count,
    typename associatedTypes::template decompose_into<std::tuple> associatedData )
  : _associatedData( associatedData ),
    _fieldGroupsContainer( Container<T,PARTICLETYPE,DynamicFieldGroupsD<T, typename PARTICLETYPE::fieldList>>(count) ),
    _dynamicsMapContainer( Container<T,PARTICLETYPE,OldBlockDynamicsMap<T,PARTICLETYPE,dynamics::ParticleDynamics<T,PARTICLETYPE>>>(count) )
{}

template<typename T, typename PARTICLETYPE>
constexpr std::size_t ParticleSystem<T, PARTICLETYPE>::size()
{
  return _fieldGroupsContainer.size();
}

template<typename T, typename PARTICLETYPE>
dynamics::ParticleDynamics<T,PARTICLETYPE>* ParticleSystem<T,PARTICLETYPE>::getDynamics (
  int iP)
{
  return get(iP).getDynamics();
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::defineDynamics (
  int iP, dynamics::ParticleDynamics<T,PARTICLETYPE>* dynamics )
{
  get(iP).defineDynamics(dynamics);
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::extend()
{
  _fieldGroupsContainer.push_back();
  _dynamicsMapContainer.push_back();
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::swapParticles(std::size_t iP, std::size_t jP)
{
  _fieldGroupsContainer.swapElements(iP,jP);
  _dynamicsMapContainer.swapElements(iP,jP);
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::process (
  std::size_t iP0, std::size_t iPmax, T timeStepSize )
{
  auto particle = get(iP0);
  for (std::size_t iP=iP0; iP<iPmax; ++iP) {
    particle.process(timeStepSize);
    particle.advanceId();
  }
}

template<typename T, typename PARTICLETYPE>
void ParticleSystem<T,PARTICLETYPE>::process(T timeStepSize)
{
  process(0, size(), timeStepSize);
}


template<typename T, typename PARTICLETYPE>
Container<T,PARTICLETYPE,DynamicFieldGroupsD<T,
          typename PARTICLETYPE::fieldList>>& ParticleSystem<T,PARTICLETYPE>::get()
{
  return _fieldGroupsContainer;
}


template<typename T, typename PARTICLETYPE>
Particle<T,PARTICLETYPE> ParticleSystem<T, PARTICLETYPE>::get(std::size_t iParticle)
{
  return Particle<T,PARTICLETYPE>( _fieldGroupsContainer.data(),
                                 _dynamicsMapContainer.data(),
                                 iParticle );
}

template<typename T, typename PARTICLETYPE>
template <typename GROUP, typename FIELD>
auto& ParticleSystem<T, PARTICLETYPE>::getFieldD()
{
  using GROUP_EVAL = typename PARTICLETYPE::template derivedField<GROUP>;
  using FIELD_EVAL = typename GROUP_EVAL::template derivedField<FIELD>;
  return _fieldGroupsContainer.data().template get<GROUP_EVAL>()
         .template get<FIELD_EVAL>();
}

template<typename T, typename PARTICLETYPE>
template <typename GROUP, typename FIELD>
auto ParticleSystem<T, PARTICLETYPE>::getFieldPointer( std::size_t iParticle )
{
  using GROUP_EVAL = typename PARTICLETYPE::template derivedField<GROUP>;
  using FIELD_EVAL = typename GROUP_EVAL::template derivedField<FIELD>;
  return _fieldGroupsContainer.data().template get<GROUP_EVAL>()
         .template get<FIELD_EVAL>()
         .getRowPointer( iParticle );
}

template<typename T, typename PARTICLETYPE>
template<typename TYPE>
auto& ParticleSystem<T, PARTICLETYPE>::getAssociatedData()
{
  const std::size_t idx = associatedTypes::template index<TYPE>();
  return std::get<idx>(_associatedData);
}



} //namespace particles

} //namespace olb


#endif
