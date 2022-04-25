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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "core/fieldArrayD.h"
#include "core/blockDynamicsMap.h"

namespace olb {

//Legacy danamics support for particle
template<typename T, typename PARTICLETYPE, typename DYNAMICS=Dynamics<T,PARTICLETYPE>>
class OldBlockDynamicsMap {
private:
  std::size_t _count;
  std::vector<DYNAMICS*> _dynamics_of_cell;

public:
  OldBlockDynamicsMap(std::size_t count, DYNAMICS* default_dynamics = nullptr):
    _count(count),
    _dynamics_of_cell(count, default_dynamics)
  { }

  OldBlockDynamicsMap(DYNAMICS* default_dynamics = nullptr):
    OldBlockDynamicsMap(1, default_dynamics)
  { }

  const DYNAMICS& get(std::size_t iCell) const
  {
    return *_dynamics_of_cell[iCell];
  }

  DYNAMICS& get(std::size_t iCell)
  {
    return *_dynamics_of_cell[iCell];
  }

  void set(std::size_t iCell, DYNAMICS* dynamics)
  {
    _dynamics_of_cell[iCell] = dynamics;
  }

  std::size_t count()
  {
    return _count;
  }

  void resize(std::size_t newCount) {
    _count=newCount;
    _dynamics_of_cell.resize(newCount);
  }

  void swap(std::size_t i, std::size_t j)
  {
    std::swap(_dynamics_of_cell[i], _dynamics_of_cell[j]);
  }

};

namespace particles{

//Forward declaration
namespace dynamics {
template<typename T, typename PARTICLETYPE> struct ParticleDynamics;
}


//Particle interface as pendent to cell interface
template <typename T, typename PARTICLETYPE>
class Particle {
private:
  DynamicFieldGroupsD<T, typename PARTICLETYPE::fieldList>& _dynamicFieldGroupsD;
  OldBlockDynamicsMap<T,PARTICLETYPE,dynamics::ParticleDynamics<T,PARTICLETYPE>>& _dynamicsMap;
  std::size_t _iParticle;
public:
  Particle( DynamicFieldGroupsD<T, typename PARTICLETYPE::fieldList>& dynamicFieldGroupsD,
            OldBlockDynamicsMap<T,PARTICLETYPE,dynamics::ParticleDynamics<T,PARTICLETYPE>>& _dynamicsMap,
            std::size_t iParticle );
  
  /// Initialize
  void init();

  ///Print 
  void print(std::size_t iParticle);
  void print();

  /// Return memory ID of the currently represented particle
  std::size_t getId() const;

  /// Jump to arbitrary particle memory ID
  void setId( std::size_t iParticle );

  /// Jump to next particle in linearization sequence
  void advanceId();

  /// Define dynamics for specific particle
  void defineDynamics(dynamics::ParticleDynamics<T,PARTICLETYPE>* dynamics);

  /// Get a pointer to the dynamics
  dynamics::ParticleDynamics<T,PARTICLETYPE>* getDynamics();

  /// Apply processing to the particle according to local dynamics
  void process(T timeStepSize);

  ////////// Get and Set functions

  /// Return copy of descriptor-declared FIELD as a vector
  template <typename GROUP, typename FIELD>
  std::enable_if_t<(
    PARTICLETYPE::template size<
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
      >() > 1),
      FieldD<T,PARTICLETYPE,
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>>>
          getField();

  /// Return copy of descriptor-declared FIELD as a scalar (more specifically as value type of FIELD)
  template <typename GROUP, typename FIELD>
  std::enable_if_t<(
    PARTICLETYPE::template size<
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
      >() == 1),
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>::template value_type<T>>
  getField();

  /// Set descriptor-declared FIELD
  template <typename GROUP, typename FIELD>
  std::enable_if_t<(
    PARTICLETYPE::template size<
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
      >() > 1),
      void>
      setField(const FieldD<T,PARTICLETYPE,
               typename PARTICLETYPE::template derivedField<GROUP>
               ::template derivedField<FIELD>>& v);

  /// Set descriptor-declared FIELD
  template <typename GROUP, typename FIELD>
  std::enable_if_t<(
    PARTICLETYPE::template size<
      typename PARTICLETYPE::template derivedField<GROUP>::template derivedField<FIELD>
      >() == 1),
      void>
      setField( typename PARTICLETYPE::template derivedField<GROUP>
                ::template derivedField<FIELD>
                ::template value_type<T> value );

};

} //namespace particles

} //namespace olb

#endif
