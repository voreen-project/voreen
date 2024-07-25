/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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

#ifndef PARTICLE_MANAGER_HH
#define PARTICLE_MANAGER_HH


namespace olb {

namespace particles {

namespace dynamics {

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::ParticleManager( 
   ParticleSystem<T,PARTICLETYPE>& particleSystem,
   SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
   SuperLattice<T,DESCRIPTOR>& sLattice,
   UnitConverter<T,DESCRIPTOR> const& converter,
   Vector<T,PARTICLETYPE::d> externalAcceleration )
   : _particleSystem(particleSystem), _sGeometry(sGeometry),
     _sLattice(sLattice), _converter(converter), _externalAcceleration(externalAcceleration)
 {}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
template<typename taskList>
void ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::unpackTasksLooped(int iP, T timeStepSize){
  auto executeTask = [&](auto task){
    //Derive type of task
    using TASK = typename decltype(task)::type; 
    //Evaluate coupling method 
    if constexpr(TASK::latticeCoupling){
      TASK::execute( _particleSystem, _sGeometry, _sLattice, _converter, iP );
    } else {
      TASK::execute( _particleSystem, _externalAcceleration, iP, timeStepSize );
    }
  };
  //Execute single tasks
  using taskListReversed = typename taskList::template decompose_into<meta::reverse_t>;
  taskListReversed::for_each(executeTask); //reverse, as for_each runs reversed
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
template<typename taskList>
void ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::unpackTasks(){
  auto executeTask = [&](auto task){
    //Derive type of task
    using TASK = typename decltype(task)::type; 
    TASK::execute( _particleSystem, _sGeometry, _sLattice, _converter );
  };
  //Execute single tasks
  using taskListReversed = typename taskList::template decompose_into<meta::reverse_t>;
  taskListReversed::for_each(executeTask); //reverse, as for_each runs reversed
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
template<typename ...TASKLIST>
void ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::execute(T timeStepSize){
  //Unpack task not requiring loop
  using tasksWithoutLoop = meta::filter_t<requires_no_loop, TASKLIST...>;
  unpackTasks<tasksWithoutLoop>();
  //Unpack tasks requiring loop
  using tasksRequiringLoop = meta::filter_t<requires_loop, TASKLIST...>;
  for (std::size_t iP=0; iP<_particleSystem.size(); ++iP) {
    unpackTasksLooped<tasksRequiringLoop>( iP, timeStepSize );
  }
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
template<typename ...TASKLIST>
void ParticleManager<T,DESCRIPTOR,PARTICLETYPE>::execute(){ execute<TASKLIST...>(_converter.getPhysDeltaT()); }


} //namespace dynamics

} //namespace particles

} //namespace olb


#endif
