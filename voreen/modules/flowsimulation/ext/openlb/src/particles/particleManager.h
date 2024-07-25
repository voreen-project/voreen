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


/* The particle manager is intended to provide a semi-generic wrapper for dynamic functions in 
 * struct form. Its main purpose is the reduction of loops over all particles. Fitting structs
 * can both be administered by the particle manager or by calling them directly it, if desired.
*/

#ifndef PARTICLE_MANAGER_H
#define PARTICLE_MANAGER_H

namespace olb {

namespace particles {

namespace dynamics {

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
class ParticleManager{
private:
  ParticleSystem<T,PARTICLETYPE>& _particleSystem;
  SuperGeometry<T,DESCRIPTOR::d>& _sGeometry;
  SuperLattice<T,DESCRIPTOR>& _sLattice;
  UnitConverter<T,DESCRIPTOR> const& _converter;
  Vector<T,PARTICLETYPE::d> _externalAcceleration;

  //Condition for TASKs requiring a particle loop
  template <typename TASK>
  using requires_loop = std::integral_constant<bool, TASK::particleLoop>;
  template <typename TASK>
  using requires_no_loop = std::integral_constant<bool, !TASK::particleLoop>;

  //Unpack tasks requiring loop
  template<typename taskList>
  void unpackTasksLooped(int iP, T timeStepSize);

  //Unpack tasks not requiring loop
  template<typename taskList>
  void unpackTasks();

public:
  //Constructor
  ParticleManager( 
    ParticleSystem<T,PARTICLETYPE>& particleSystem,
    SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
    SuperLattice<T,DESCRIPTOR>& sLattice,
    UnitConverter<T,DESCRIPTOR> const& converter,
    Vector<T,PARTICLETYPE::d> externalAcceleration = Vector<T,PARTICLETYPE::d>(0.) );

  //Execute task and pass time step
  template<typename ...TASKLIST>
  void execute(T timeStepSize);

  //Execute task and use timestep provided by the converter 
  template<typename ...TASKLIST>
  void execute();
};


} //namespace dynamics

} //namespace particles

} //namespace olb


#endif
