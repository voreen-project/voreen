/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Particle Reynold number classes for Lagrangian two-way coupling methods -- header file.
 */

#ifndef LB_REYNOLDS_NUMBERS_3D_H
#define LB_REYNOLDS_NUMBERS_3D_H

#include "functors/lattice/reductionF3D.h"

namespace olb {

/** Abstract class for particle Reynolds number computation within drag model.
  * Its raison d'etre consists of not being templetized in Lattice.
  */
template<typename T, template<typename V> class Particle>
class ParticleReynoldsNumber {
public:
  /// Returns the particle Reynolds number. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T magU, int globicFull[])=0;
  /// Destructor
  virtual ~ParticleReynoldsNumber() {};
protected:
  T _RePmin = 0.01;
};

/** Abstract class for particle Reynolds number computation within drag model.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class ParticleReynoldsNumberBase : public ParticleReynoldsNumber<T, Particle> {
public:
  /// Destructor
  virtual ~ParticleReynoldsNumberBase() {};
protected:
  /// Constructor
  ParticleReynoldsNumberBase(UnitConverter<T, Lattice>& converter);
  UnitConverter<T, Lattice>& _converter; // reference to a UnitConverter
};

/** Class class for Newtonian particle Reynolds number computation within drag model.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class NewtonianParticleReynoldsNumber: public ParticleReynoldsNumberBase<T,Lattice,Particle> {
public:
  /// Constructor
  NewtonianParticleReynoldsNumber(UnitConverter<T, Lattice>& converter);
  /// Destructor
  ~NewtonianParticleReynoldsNumber() {};
  /// Returns the particle Reynolds number. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T magU, int globicFull[]) override;
};

/** Class class for power-law particle Reynolds number computation within drag model.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class PowerLawParticleReynoldsNumber: public ParticleReynoldsNumberBase<T,Lattice,Particle> {
public:
  /// Constructor
  PowerLawParticleReynoldsNumber (
    UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice );
  /// Destructor
  ~PowerLawParticleReynoldsNumber() {};
  /// Returns the particle Reynolds number. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T magU, int globicFull[]) override;
protected:
  SuperLattice<T, Lattice>& _sLattice; // reference to a lattice
};

}

#endif
