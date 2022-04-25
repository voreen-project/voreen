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

#ifndef PARTICLE_DYNAMICS_BASE_H
#define PARTICLE_DYNAMICS_BASE_H

namespace olb {

namespace particles {

namespace dynamics {


/// Basic particle dynamics
template<typename T, typename PARTICLETYPE>
struct ParticleDynamics {
  /// Destructor: virtual to enable inheritance
  virtual ~ParticleDynamics() { }
  /// Implementation of the processing step
  virtual void process(Particle<T,PARTICLETYPE>& particle, T timeStepSize) =0;
  /// read and write access to name
  std::string& getName();
  /// read only access to name
  std::string const& getName() const;
private:
  std::string _name;
};

/// No particle dynamics equivalent to no lattice dynamics
template<typename T, typename PARTICLETYPE>
class NoParticleDynamics : public ParticleDynamics<T,PARTICLETYPE> {
public:
  NoParticleDynamics( T rhoDummy );
  /// Processing step
  void process(Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
};

/// Standard dynamics for particles
template<typename T, typename PARTICLETYPE>
class VerletParticleDynamics : public ParticleDynamics<T,PARTICLETYPE> {
public:
  /// Constructor
  VerletParticleDynamics( );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
};

/// Velocity verlet particle dynamics with limitation of position and velocity by checking domain bounds
/// in cartesion direcion
template<typename T, typename PARTICLETYPE>
class VerletParticleDynamicsCubicBoundsCheck : public ParticleDynamics<T,PARTICLETYPE> {
public:
  /// Constructor
  VerletParticleDynamicsCubicBoundsCheck( PhysR<T,PARTICLETYPE::d>& domainMin,
                                          PhysR<T,PARTICLETYPE::d>& domainMax );
  /// Procesisng step
  void process (Particle<T,PARTICLETYPE>& particle, T timeStepSize) override;
private:
  PhysR<T,PARTICLETYPE::d> _domainMin;
  PhysR<T,PARTICLETYPE::d> _domainMax;
};

} //namespace dynamics

} //namespace particles

} //namespace olb

#endif
