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

#ifndef PARTICLE_3D_H
#define PARTICLE_3D_H

#include <set>
#include <vector>
#include <list>
#include <deque>
#include <string>
#include <iostream>
#include "twoWayCouplings/twoWayCouplings3D.h"
//#include "../core/superLattice3D.h"

namespace olb {


/////////////////////////////////////////// Particle3D ///////////////////////////////////////////

template<typename T>
class Particle3D {
public:
  Particle3D();
  Particle3D(std::vector<T> pos, T mas = 1., T rad = 1., int id = 0);
  Particle3D(std::vector<T> pos, T mas = 1., T masAdd = 1., T rad = 1., int id = 0);
  Particle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1., int id = 0);
  Particle3D(const Particle3D<T>& p);
  inline void setPos(std::vector<T> pos);
  inline void setStoredPos(std::vector<T> pos);
  inline std::vector<T>& getStoredPos();
  std::vector<T>& getPos();
  const std::vector<T>& getPos() const;
  inline void setVel(std::vector<T> vel);
  inline void setStoredVel(std::vector<T> vel);
  inline std::vector<T>& getStoredVel();
  inline int getID();
  inline void setID(int id);
  inline std::vector<T>& getVel();
  inline const std::vector<T>& getVel() const;
  inline void addForce(std::vector<T>& frc);
  inline void setForce(std::vector<T>& frc);
  inline void resetForce();
  inline std::vector<T>& getForce();
  inline const std::vector<T>& getForce() const;
  inline void setStoreForce(std::vector<T>& storeForce);
  inline void resetStoreForce();
  inline std::vector<T>& getStoreForce();
  inline const std::vector<T>& getStoreForce() const;
  // RK4
//  inline void setPos(std::vector<T> pos, int i);
//  inline std::vector<T>& getPos(int i);
//  inline const std::vector<T>& getPos(int i) const;
//  inline void setVel(std::vector<T> vel, int i);
//  inline std::vector<T>& getVel(int i);
//  inline const std::vector<T>& getVel(int i) const;
//  inline void setForce(std::vector<T> force, int i);
//  inline std::vector<T>& getForce(int i);
//  inline const std::vector<T>& getForce(int i) const;
// write particle properties into serial
  void serialize(T serial[]);
  // write data into particle properties
  void unserialize(T*);
  void print();
  void printDeep(std::string message);
  inline const T& getMass();
  inline const T& getAddedMass();
  inline const T& getInvMass();
  inline const T& getInvAddedMass();
  inline const T& getMass() const;
  inline const T& getAddedMass() const;
  inline void setMass(T m);
  inline void setAddedMass(T m);
  inline const T& getRad();
  inline const T& getRad() const;
  inline void setRad(T r);
  inline const int& getCuboid();
  inline void setCuboid(int c);
  inline const bool& getActive();
  inline const bool& getActive() const;
  inline void setActive(bool act);

  static const int serialPartSize = 18; // pos, vel, force, mas, masAdd, rad, cuboid, id, active, storeForce
  std::vector<std::pair<size_t, T> > _verletList;

protected:
  std::vector<T> _pos;
  std::vector<T> _vel;
  std::vector<T> _force;
  T _invMas;
  T _invMasAdd;
  T _mas;
  T _masAdd;
  T _rad;
  ///globIC
  int _cuboid;
  int _id;
  bool _active;
  std::vector<T> _storePos;
  std::vector<T> _storeVel;
  std::vector<T> _storeForce;
  // RK4
//  std::vector<std::vector<T> > _positions;
//  std::vector<std::vector<T> > _velocities;
//  std::vector<std::vector<T> > _forces;

};

template<typename T, template<typename U> class PARTICLETYPE> class ParticleSystem3D;

template<typename T, template<typename U> class PARTICLETYPE>
class SimulateParticles {
public:
  SimulateParticles(ParticleSystem3D<T, PARTICLETYPE>* ps);
  inline void simulate(T dT, bool scale = false);
  inline void simulate(T dT, std::set<int> sActivityOfParticle, bool scale = false);
  inline void simulateWithTwoWayCoupling_Mathias ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps = 1, bool scale = false );
  inline void simulateWithTwoWayCoupling_Davide ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps = 1, bool scale = false );

private:
  ParticleSystem3D<T, PARTICLETYPE>* _pSys;
};

}
#endif /* PARTICLE_3D_H */
