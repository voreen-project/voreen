/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Davide Dapelo
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

#ifndef PARTICLE_3D_HH
#define PARTICLE_3D_HH

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <deque>

#include "particle3D.h"

namespace olb {


  /////////////////////////////////////// Particle3D ///////////////////////////////////////

  template<typename T>
  Particle3D<T>::Particle3D()
    : _pos(3, 0.),
    _vel(3, 0.),
    _force(3, 0.),
    _effectiveMas(1.),
    _mas(1.),
    _masAdd(0.),
    _rad(0),
    _cuboid(0),
    _id(0),
    _active(false),
    _storeForce(3, 0.)
  { }

  /*
  template<typename T>
  Particle3D<T>::Particle3D(std::vector<T> pos, T mas, T rad, int id)
    : _pos(pos),
      _vel(3, 0.),
      _force(3, 0.),
      _rad(rad),
      _cuboid(0),
      _id(id),
      _active(true),
      _storeForce(3, 0.)
  {
    setMass(mas);
    setAddedMass(mas);
    // RK4
  //  _positions = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _velocities = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _forces = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  }
  */

  template<typename T>
  Particle3D<T>::Particle3D(std::vector<T> pos, T mas, T rad, int id, T masAdd)
    : _pos(pos),
    _vel(3, 0.),
    _force(3, 0.),
    _rad(rad),
    _cuboid(0),
    _id(id),
    _active(true),
    _storeForce(3, 0.)
  {
    setMass(mas);
    setAddedMass(masAdd);
    // RK4
  //  _positions = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _velocities = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _forces = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  }

  template<typename T>
  Particle3D<T>::Particle3D(const Particle3D<T>& p)
    : _pos(p._pos),
    _vel(p._vel),
    _force(p._force),
    _rad(p._rad),
    _cuboid(p._cuboid),
    _id(p._id),
    _active(p._active),
    _storeForce(p._storeForce)
  {
    setMass(p._mas);
    setAddedMass(p._masAdd);
    // RK4
  //  _positions = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _velocities = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _forces = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  }

  /*
  template<typename T>
  Particle3D<T>::Particle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id)
    : _pos(pos),
      _vel(vel),
      _force(12, 0.),
      _rad(rad),
      _cuboid(0),
      _id(id),
      _active(true),
      _storeForce(3, 0.)
  {
    setMass(mas);
    setAddedMass(mas);
    _vel.resize(12, 0.);
    // RK4
  //  _positions = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _velocities = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _forces = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  }
  */

  template<typename T>
  Particle3D<T>::Particle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id, T masAdd)
    : _pos(pos),
    _vel(vel),
    _force(12, 0.),
    _rad(rad),
    _cuboid(0),
    _id(id),
    _active(true),
    _storeForce(3, 0.)
  {
    setMass(mas);
    setAddedMass(masAdd);
    _vel.resize(12, 0.);
    // RK4
  //  _positions = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _velocities = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  //  _forces = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
  }

  template<typename T>
  inline void Particle3D<T>::setPos(std::vector<T> pos)
  {
    _pos = pos;
  }

  template<typename T>
  inline void Particle3D<T>::setStoredPos(std::vector<T> pos)
  {
    _storePos = pos;
  }

  template<typename T>
  inline std::vector<T>& Particle3D<T>::getStoredPos()
  {
    return _storePos;
  }

  template<typename T>
  inline std::vector<T>& Particle3D<T>::getPos()
  {
    return _pos;
  }

  template<typename T>
  inline const std::vector<T>& Particle3D<T>::getPos() const
  {
    return _pos;
  }

  template<typename T>
  inline void Particle3D<T>::setVel(std::vector<T> vel)
  {
    _vel = vel;
  }

  template<typename T>
  inline void Particle3D<T>::setStoredVel(std::vector<T> vel)
  {
    _storeVel = vel;
  }

  template<typename T>
  inline std::vector<T>& Particle3D<T>::getStoredVel()
  {
    return _storeVel;
  }

  template<typename T>
  inline int Particle3D<T>::getID()
  {
    return _id;
  }

  template<typename T>
  inline void Particle3D<T>::setID(int id)
  {
    _id = id;
  }

  template<typename T>
  inline std::vector<T>& Particle3D<T>::getVel()
  {
    return _vel;
  }

  template<typename T>
  inline const std::vector<T>& Particle3D<T>::getVel() const
  {
    return _vel;
  }

  template<typename T>
  inline void Particle3D<T>::addForce(std::vector<T>& force)
  {
    for (int i = 0; i < 3; i++) {
      _force[i] += force[i];
    }
  }
  // set and get force
  template<typename T>
  inline void Particle3D<T>::setForce(std::vector<T>& force)
  {
    _force = force;
  }
  template<typename T>
  inline void Particle3D<T>::resetForce()
  {
    for (int i = 0; i < 3; i++) {
      _force[i] = 0.;
    }
  }

  // set and get storedForce
  template<typename T>
  inline void Particle3D<T>::setStoreForce(std::vector<T>& storeForce)
  {
    for (int i = 0; i < 3; i++) {
      _storeForce[i] = storeForce[i];
    }
  }

  template<typename T>
  inline void Particle3D<T>::resetStoreForce()
  {
    for (int i = 0; i < 3; i++) {
      _storeForce[i] = T(0);
    }
  }

  template<typename T>
  inline std::vector<T>& Particle3D<T>::getForce()
  {
    return _force;
  }

  template<typename T>
  inline const std::vector<T>& Particle3D<T>::getForce() const
  {
    return _force;
  }

  template<typename T>
  inline std::vector<T>& Particle3D<T>::getStoreForce()
  {
    return _storeForce;
  }

  template<typename T>
  inline const std::vector<T>& Particle3D<T>::getStoreForce() const
  {
    return _storeForce;
  }

  // RK4
  //template<typename T>
  //inline void Particle3D<T>::setPos(std::vector<T> pos, int i)
  //{
  //  _positions[i] = pos;
  //}

  //template<typename T>
  //inline std::vector<T>& Particle3D<T>::getPos(int i)
  //{
  //  return _positions[i];
  //}

  //template<typename T>
  //inline const std::vector<T>& Particle3D<T>::getPos(int i) const
  //{
  //  return _positions[i];
  //}

  //template<typename T>
  //inline void Particle3D<T>::setVel(std::vector<T> vel, int i)
  //{
  //  _velocities[i] = vel;
  //}

  //template<typename T>
  //inline std::vector<T>& Particle3D<T>::getVel(int i)
  //{
  //  return _velocities[i];
  //}

  //template<typename T>
  //inline const std::vector<T>& Particle3D<T>::getVel(int i) const
  //{
  //  return _velocities[i];
  //}

  //template<typename T>
  //inline void Particle3D<T>::setForce(std::vector<T> force, int i)
  //{
  //  _forces[i] = force;
  //}

  //template<typename T>
  //inline std::vector<T>& Particle3D<T>::getForce(int i)
  //{
  //  return _forces[i];
  //}

  //template<typename T>
  //inline const std::vector<T>& Particle3D<T>::getForce(int i) const
  //{
  //  return _forces[i];
  //}

  template<typename T>
  void Particle3D<T>::serialize(T serial[])
  {
    for (int i = 0; i < 3; i++) {
      serial[i] = _pos[i];
      serial[i + 3] = _vel[i];
      serial[i + 6] = _force[i];
    }
    serial[9] = _mas;
    serial[10] = _masAdd;
    serial[11] = _effectiveMas;
    serial[12] = _rad;
    serial[13] = _cuboid;
    serial[14] = _active;
    serial[15] = _id;

    for (int i = 0; i < 3; i++) {
      serial[i + 16] = _storeForce[i];
    }

    //for (int i = 0; i < 18; i++) {
    //cout << "serialize " << i << ": " << serial[i]  << " tn: " << typeid(serial[i]).name() << std::endl;
    //}
  }

  template<typename T>
  void Particle3D<T>::unserialize(T* data)
  {
    for (int i = 0; i < 3; i++) {
      _pos[i] = data[i];
      _vel[i] = data[i + 3];
      _force[i] = data[i + 6];
    }
    _mas = data[9];
    _masAdd = data[10];
    _effectiveMas = data[11];
    _rad = data[12];
    _cuboid = int(data[13]);
    _active = data[14];
    _invMas = 1. / _mas;
    _invMasAdd = 1. / _masAdd;
    _invEffectiveMas = 1. / _effectiveMas;
    _id = data[15];

    for (int i = 0; i < 3; i++) {
      _storeForce[i] = data[i + 16];
    }
  }

  template<typename T>
  void Particle3D<T>::print()
  {
    std::cout << "Pos=(" << _pos[0] << ", " << _pos[1] << ", " << _pos[2] << ") "
      << "Vel=(" << _vel[0] << ", " << _vel[1] << ", " << _vel[2] << ") "
      << "Cub=" << _cuboid
      << std::endl;
  }

  template<typename T>
  void Particle3D<T>::printDeep(std::string message)
  {
    std::cout << message
      << " ID="          << this->getID()
//      << " rad="         << this->getRad()
//      << " mass="        << this->getMass()
//      << " invMass="     << this->getInvMass()
//      << " addedMass="   << this->getAddedMass()
//      << " force=(" << this->getForce()[0] << ", " << this->getForce()[1] << ", " << this->getForce()[2] << ")"
//      << " storeForce=(" << this->getStoreForce()[0] << " " << this->getStoreForce()[1] << ", " << this->getStoreForce()[2] << ")"
//      << " active="      << this->getActive()
//      << std::endl;
    << " ";
    this->print();
//    std::cout << std::endl;
  }

  template<typename T>
  inline const T& Particle3D<T>::getMass()
  {
    return _mas;
  }

  template<typename T>
  inline const T& Particle3D<T>::getAddedMass()
  {
    return _masAdd;
  }

  template<typename T>
  inline const T& Particle3D<T>::getInvMass()
  {
    return _invMas;
  }

  template<typename T>
  inline const T& Particle3D<T>::getInvEffectiveMass()
  {
    return _invEffectiveMas;
  }

  template<typename T>
  inline const T& Particle3D<T>::getEffectiveMass()
  {
    return _effectiveMas;
  }

  template<typename T>
  inline const T& Particle3D<T>::getMass() const
  {
    return _mas;
  }

  template<typename T>
  inline const T& Particle3D<T>::getAddedMass() const
  {
    return _masAdd;
  }

  template<typename T>
  inline const T& Particle3D<T>::getEffectiveMass() const
  {
    return _effectiveMas;
  }

  template<typename T>
  inline void Particle3D<T>::setMass(T m)
  {
    if (m <= T()) {
      throw std::invalid_argument("Exception called in particle3D::setMass(T). Input value must be > 0, but instead was " + std::to_string(m));
    }
    _mas = m;
    _invMas = 1. / _mas;
    _effectiveMas = _mas + _masAdd;
    _invEffectiveMas = 1. / _effectiveMas;
  }

  template<typename T>
  inline void Particle3D<T>::setAddedMass(T m)
  {
    if (m < T()) {
      throw std::invalid_argument("Exception called in particle3D::setAddedMass(T). Input value must be >= 0, but instead was " + std::to_string(m));
    }
    _masAdd = m;
    _effectiveMas = _mas + _masAdd;
    _invEffectiveMas = 1. / _effectiveMas;
  }

  template<typename T>
  inline const T& Particle3D<T>::getRad()
  {
    return _rad;
  }

  template<typename T>
  inline const T& Particle3D<T>::getRad() const
  {
    return _rad;
  }

  template<typename T>
  inline void Particle3D<T>::setRad(T r)
  {
    _rad = r;
  }

  template<typename T>
  inline const int& Particle3D<T>::getCuboid()
  {
    return _cuboid;
  }

  template<typename T>
  inline void Particle3D<T>::setCuboid(int c)
  {
    _cuboid = c;
  }

  template<typename T>
  inline const bool& Particle3D<T>::getActive()
  {
    return _active;
  }

  template<typename T>
  inline const bool& Particle3D<T>::getActive() const
  {
    return _active;
  }

  template<typename T>
  inline void Particle3D<T>::setActive(bool act)
  {
    _active = act;
    if (!act) {
      _vel[0] = 0;
      _vel[1] = 0;
      _vel[2] = 0;
    }
  }


  /////////////////////////////////////// SimulateParticles ///////////////////////////////////////

  template<typename T, template<typename U> class PARTICLETYPE>
  SimulateParticles<T,PARTICLETYPE>::SimulateParticles(ParticleSystem3D<T, PARTICLETYPE>* ps)
    : _pSys(ps)
  { }

  template<typename T, template<typename U> class PARTICLETYPE>
  inline void SimulateParticles<T,PARTICLETYPE>::simulate(T dT, bool scale)
  {
    _pSys->computeForce();
    _pSys->explicitEuler(dT, scale);
    //_pSys->rungeKutta4(dT);
  }

  template<typename T, template<typename U> class PARTICLETYPE>
  inline void SimulateParticles<T,PARTICLETYPE>::simulate(T dT, std::set<int> sActivityOfParticle, bool scale)
  {
    simulate(dT, scale);
  }

  template<typename T, template<typename U> class PARTICLETYPE>
  inline void SimulateParticles<T,PARTICLETYPE>::simulateWithTwoWayCoupling_Mathias ( T dT,
                                                                                     ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                                                                     BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                                                                     int material, int subSteps, bool scale )
  {
    for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
      if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
        std::cout << " on substep " << iSubStep << std::endl;
        singleton::exit(1);
      }
      _pSys->computeForce();
      _pSys->explicitEuler(dT/(T)(subSteps), scale);
      //_pSys->rungeKutta4(dT/(T)(subSteps));
    }
    _pSys->executeBackwardCoupling(backCoupling, material);
  }

  template<typename T, template<typename U> class PARTICLETYPE>
  inline void SimulateParticles<T,PARTICLETYPE>::simulateWithTwoWayCoupling_Davide ( T dT,
                                                                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                                                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                                                                    int material, int subSteps, bool scale )
  {
    for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
      if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
        std::cout << " on substep " << iSubStep << std::endl;
        singleton::exit(1);
      }
      _pSys->executeBackwardCoupling(backCoupling, material, subSteps);
      _pSys->computeForce();
      _pSys->explicitEuler(dT/(T)(subSteps), scale);
      //_pSys->rungeKutta4(dT/(T)(subSteps));
    }
  }


}

#endif /* PARTICLE_3D_HH */
