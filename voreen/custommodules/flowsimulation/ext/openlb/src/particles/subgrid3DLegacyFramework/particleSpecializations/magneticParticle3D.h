/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause, Davide Dapelo
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

#ifndef MAGNETIC_PARTICLE_3D_H
#define MAGNETIC_PARTICLE_3D_H

#include <set>
#include <vector>
#include <list>
#include <deque>
#include <string>
#include <iostream>
#include "particles/subgrid3DLegacyFramework/particle3D.h"

namespace olb {


/////////////////////////////////////////// MagneticParticle3D ///////////////////////////////////////////

/*
 * Particles for magnetic force
 */
template<typename T>
class MagneticParticle3D: public Particle3D<T> {
private:
  std::vector<T> _dMoment;
  std::vector<T> _aVel;
  std::vector<T> _torque;
  T _magnetisation;
  // _sActivity values range from 1 to 3
  // Gives information about the position of the particle in the simulation area
  // One collision model:
  // 1: Particle is not deposited
  // 2: Not in use
  // 3: Particle is deposited on mag. wire.
  // Multiple collision models:
  // 1: Particle is in sufficiently large distance to the inhom. mag. field of the wire
  //    Particle collisions are managed by a collision model function
  // 2: Particle is within the inhom. mag. field of the wire
  //    Particle collisions are managed by a mechanic contact force
  //    The sActivity change (1 --> 2) is done in explicitEuler()
  // 3: Particle is deposited on mag. wire.
  //    Particle collisions are managed by a mechanic contact force
  int _sActivity;
  // Assigns an agglomerate to the particle
  typename std::deque<std::list<MagneticParticle3D<T>*>>::iterator _aggloItr;
  // Assigns the particle position in an Agglomerate to a Particle
  // typename std::list<MagneticParticle3D<T>*>::iterator _aggloItrPos;

public:
  MagneticParticle3D();
  MagneticParticle3D(std::vector<T> pos, T mas = 1., T rad = 1., int id = 0);
  MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1., int id = 0, T masAdd = 1.);
  MagneticParticle3D(const MagneticParticle3D<T>& p);
  MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id,
                     std::vector<T> dMoment, std::vector<T> aVel, std::vector<T> torque, T magnetisation);
  MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id,
                     std::vector<T> dMoment, std::vector<T> aVel, std::vector<T> torque, T magnetisation, int sActivity);
  static const int serialPartSize = 28;
  void serialize(T serial[]);
  void unserialize(T*);
  // Set torque zero
  inline void resetTorque();
  // Set and get orientation of magnetic dipolemoment
  inline void setMoment(std::vector<T> moment);
  inline std::vector<T>& getMoment();
  inline const std::vector<T>& getMoment() const;
  // Set and get angular velocity
  inline void setAVel(std::vector<T> aVel);
  inline std::vector<T>& getAVel();
  inline const std::vector<T>& getAVel() const;
  // Set and get torque
  inline void setTorque(std::vector<T> torque);
  inline std::vector<T>& getTorque();
  inline const std::vector<T>& getTorque() const;
  // Set and get magnetisation
  inline void setMagnetisation(T magnetisation);
  inline T& getMagnetisation();
  inline const T& getMagnetisation() const;
  // Set and get sActivity
  inline void setSActivity(int sActivity);
  inline int& getSActivity();
  // Set and get aggloItr
  inline void setAggloItr(typename std::deque<std::list<MagneticParticle3D<T>*>>::iterator aggloItr);
  inline typename std::deque<std::list<MagneticParticle3D<T>*>>::iterator& getAggloItr();
};

/////////////////////////////////////// SimulateParticles<T, MagneticParticle3D> ///////////////////////////////////////
template<typename T>
class SimulateParticles<T, MagneticParticle3D> {
public:
  SimulateParticles(ParticleSystem3D<T, MagneticParticle3D>* ps);
  inline void simulate(T dT, bool scale = false);
  inline void simulateWithTwoWayCoupling_Mathias ( T dT,
      ForwardCouplingModel<T,MagneticParticle3D>& forwardCoupling,
      BackCouplingModel<T,MagneticParticle3D>& backCoupling,
      int material, int subSteps = 1, bool scale = false );
  inline void simulateWithTwoWayCoupling_Davide ( T dT,
      ForwardCouplingModel<T,MagneticParticle3D>& forwardCoupling,
      BackCouplingModel<T,MagneticParticle3D>& backCoupling,
      int material, int subSteps = 1, bool scale = false );
//   multiple collision models
  inline void simulate(T dT, std::set<int> sActivityOfParticle, bool scale = false);

private:
  ParticleSystem3D<T, MagneticParticle3D>* _pSys;

};


}
#endif

