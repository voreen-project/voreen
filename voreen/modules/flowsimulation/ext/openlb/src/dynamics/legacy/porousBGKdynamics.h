/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Jonas Latt
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

#ifndef LEGACY_POROUS_BGK_DYNAMICS_H
#define LEGACY_POROUS_BGK_DYNAMICS_H

#include "dynamics.h"

namespace olb {

/// Implementation of the BGK collision step for subgridscale particles
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class SubgridParticleBGKdynamics : public legacy::BGKdynamics<T,DESCRIPTOR,MOMENTA> {
private:
  T omega;      ///< relaxation parameter
  T _fieldTmp[4];

public:
  template<typename M>
  using exchange_momenta = SubgridParticleBGKdynamics<T,DESCRIPTOR,M>;

  SubgridParticleBGKdynamics(T omega_);
  /// extended Collision step, computes local drag in a given direction
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell) override;

  /// get relaxation parameter
  T    getOmega() const;

};

/* Generic implementation for moving porous media (HLBM approach).
 * As this scheme requires additionla data stored in an external field,
 * it is meant to be used along with a PorousParticle descriptor.
 * \param omega Lattice relaxation frequency
 * \param momenta A standard object for the momenta computation
 */
template<typename T, typename DESCRIPTOR, bool isStatic=false>
class PorousParticleDynamics {
protected:
  template <bool isStatic_ = isStatic>
  static std::enable_if_t<isStatic_> calculate(ConstCell<T,DESCRIPTOR>& cell, T* pVelocity);
  template <bool isStatic_ = isStatic>
  static std::enable_if_t<!isStatic_> calculate(Cell<T,DESCRIPTOR>& cell, T* pVelocity);
};

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple, bool isStatic=false>
class PorousParticleBGK final : public legacy::BGKdynamics<T,DESCRIPTOR,MOMENTA>
                              , public PorousParticleDynamics<T,DESCRIPTOR,isStatic> {
public:
  PorousParticleBGK(T omega_);
protected:
  T porousParticleBgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], T omega);
};



/* Implementation of the BGK collision for Zeta-Field formulation (Geng2019).
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic=false>
class DBBParticleBGKdynamics : public legacy::BGKdynamics<T,DESCRIPTOR,MOMENTA>
                             , public PorousParticleDynamics<T,DESCRIPTOR,isStatic> {
public:
  template<typename M>
  using exchange_momenta = DBBParticleBGKdynamics<T,DESCRIPTOR,M,isStatic>;

  DBBParticleBGKdynamics(T omega_);
  /// extended Collision step, computes local drag in a given direction
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell) override;

protected:
  T dbbParticleBgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], T eta[DESCRIPTOR::d], T uPlus[DESCRIPTOR::d], T tmp_cell[(DESCRIPTOR::q+1)/2], T omega );

};


/// Implementation of the HBGK collision step for a porosity model enabling
/// drag computation for many particles
/// including the Krause turbulence modell
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class KrauseHBGKdynamics : public legacy::BGKdynamics<T,DESCRIPTOR,MOMENTA> {
private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  void computeOmega(T omega0_, Cell<T,DESCRIPTOR>& cell, T preFactor_, T rho_,
                    T u[DESCRIPTOR::d],
                    T newOmega[DESCRIPTOR::q] );

  T omega;      ///< relaxation parameter
  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;

  T _fieldTmp[4];

public:
  template<typename M>
  using exchange_momenta = KrauseHBGKdynamics<T,DESCRIPTOR,M>;

  KrauseHBGKdynamics(T omega_, T smagoConst_, T dx_ = 1, T dt_ = 1);
  /// extended Collision step, computes local drag in a given direction
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell) override;

  /// get relaxation parameter
  T    getOmega() const;

};

/// Implementation of the BGK collision step for a porosity model enabling
/// drag computation
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class ParticlePorousBGKdynamics : public legacy::BGKdynamics<T,DESCRIPTOR,MOMENTA> {
private:
  T omega;      ///< relaxation parameter

public:
  template<typename M>
  using exchange_momenta = ParticlePorousBGKdynamics<T,DESCRIPTOR,M>;

  ParticlePorousBGKdynamics(T omega_);
  /// extended Collision step, computes local drag in a given direction
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell) override;

  /// get relaxation parameter
  T    getOmega() const;

};

/// Implementation of the BGK collision step for a small particles enabling
/// two way coupling
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class SmallParticleBGKdynamics : public legacy::BGKdynamics<T,DESCRIPTOR,MOMENTA> {
private:
  T omega;      ///< relaxation parameter

public:
  template<typename M>
  using exchange_momenta = SmallParticleBGKdynamics<T,DESCRIPTOR,M>;

  SmallParticleBGKdynamics(T omega_);
  /// extended Collision step, computes local drag in a given direction
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell) override;

  /// get relaxation parameter
  T    getOmega() const;

};

} // olb

#endif

#include "porousBGKdynamics.hh"
