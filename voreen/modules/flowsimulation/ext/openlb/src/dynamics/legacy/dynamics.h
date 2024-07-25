/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef LB_LEGACY_DYNAMICS_H
#define LB_LEGACY_DYNAMICS_H

#include "dynamics/interface.h"

namespace olb {

namespace legacy {

template <typename T, typename DESCRIPTOR, typename MOMENTA>
struct BasicDynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override {
    return equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u);
  };

  std::type_index id() override {
    return typeid(BasicDynamics<T,DESCRIPTOR,MOMENTA>);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    throw std::bad_function_call();
  }

};

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class BGKdynamics : public BasicDynamics<T,DESCRIPTOR,MOMENTA> {
public:
  BGKdynamics(T omega);
  /// Collision step
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell) override;

  virtual T getOmega() const {
    return _omega;
  }

private:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

  T _omega;
};

template<typename T, typename DESCRIPTOR, typename MOMENTA>
BGKdynamics<T,DESCRIPTOR,MOMENTA>::BGKdynamics(T omega)
  : BasicDynamics<T,DESCRIPTOR,MOMENTA>(),
    _omega(omega)
{ }

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> BGKdynamics<T,DESCRIPTOR,MOMENTA>::collide(Cell<T,DESCRIPTOR>& cell)
{
  T rho, u[DESCRIPTOR::d];
  MomentaF().computeRhoU(cell, rho, u);
  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, _omega);
  return {rho, uSqr};
}

template<typename T, typename DESCRIPTOR>
class NoLatticeDynamics : public Dynamics<T,DESCRIPTOR> {
public:
  /// You may fix a fictitious density value on no dynamics node via this constructor.
  NoLatticeDynamics(T rho = T(1) );
  /// Yields 0;
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override any_platform;
  /// Collision step
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell) override;
  /// Yields 1;
  T computeRho(ConstCell<T,DESCRIPTOR>& cell) const override;
  /// Yields 0;
  void computeU (
    ConstCell<T,DESCRIPTOR>& cell,
    T u[DESCRIPTOR::d] ) const override;
  /// Yields 0;
  void computeJ (
    ConstCell<T,DESCRIPTOR>& cell,
    T j[DESCRIPTOR::d] ) const override;
  /// Yields NaN
  void computeStress (
    ConstCell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR >::n] ) const override;
  void computeRhoU (
    ConstCell<T,DESCRIPTOR>& cell,
    T& rho, T u[DESCRIPTOR::d]) const override;
  void computeAllMomenta (
    ConstCell<T,DESCRIPTOR>& cell,
    T& rho, T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR >::n] ) const override;
  /// Does nothing
  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override;
  /// Does nothing
  void defineU(Cell<T,DESCRIPTOR>& cell,
               const T u[DESCRIPTOR::d]) override;
  /// Does nothing
  void defineRhoU (
    Cell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d]) override;
  /// Does nothing
  void defineAllMomenta (
    Cell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d],
    const T pi[util::TensorVal<DESCRIPTOR >::n] ) override;

  std::type_index id() override {
    return typeid(NoLatticeDynamics<T,DESCRIPTOR>);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    throw std::bad_function_call();
  }

private:
  /// Default rho=1
  T _rho;
};

/// Dynamics for offLattice boundary conditions
/// OffDynamics are basically NoLatticeDynamics with the additional functionality
/// to store given velocities exactly at boundary links.
template<typename T, typename DESCRIPTOR>
class OffDynamics : public NoLatticeDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  OffDynamics(const T _location[DESCRIPTOR::d]);
  /// Constructor
  OffDynamics(const T _location[DESCRIPTOR::d], T _distances[DESCRIPTOR::q]);
  /// Returns local stored rho which is updated if the bc is used as velocity!=0 condition
  T computeRho(ConstCell<T,DESCRIPTOR>& cell) const override;
  /// Returns an average of the locally stored u
  void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const override;
  /// Set Intersection of the link and the boundary
  void setBoundaryIntersection(int iPop, T distance);
  /// Get Intersection of the link and the boundary
  bool getBoundaryIntersection(int iPop, T intersection[DESCRIPTOR::d]);
  /// Set particle density on the cell.
  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override;
  /// Set single velocity
  void defineRho(int iPop, T rho);
  /// Set fluid velocity on the cell.
  void defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) override;
  /// Set constant velocity
  void defineU(const T u[DESCRIPTOR::d]);
  /// Set single velocity
  void defineU(int iPop, const T u[DESCRIPTOR::d]);
  /// Get VelocitySummand for Bouzidi-Boundary Condition
  T getVelocityCoefficient(int iPop);

  std::type_index id() override {
    return typeid(OffDynamics<T,DESCRIPTOR>);
  }


private:
  T _rho;
  T _u[DESCRIPTOR::q][DESCRIPTOR::d];
  T location[DESCRIPTOR::d];
  T distances[DESCRIPTOR::q];
  T boundaryIntersection[DESCRIPTOR::q][DESCRIPTOR::d];
  T velocityCoefficient[DESCRIPTOR::q];
};

////////////////////// Class NoLatticeDynamics ///////////////////////////

template<typename T, typename DESCRIPTOR>
NoLatticeDynamics<T,DESCRIPTOR>::NoLatticeDynamics(T rho) :_rho(rho)
{
  this->getName() = "NoLatticeDynamics";
}

template<typename T, typename DESCRIPTOR>
T NoLatticeDynamics<T,DESCRIPTOR>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const
{
  return T();
}

template<typename T, typename DESCRIPTOR>
CellStatistic<T> NoLatticeDynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell) {
  return {-1, -1};
}

template<typename T, typename DESCRIPTOR>
T NoLatticeDynamics<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return _rho;
}

template<typename T, typename DESCRIPTOR>
void NoLatticeDynamics<T,DESCRIPTOR>::computeU (
  ConstCell<T,DESCRIPTOR>& cell,
  T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void NoLatticeDynamics<T,DESCRIPTOR>::computeJ (
  ConstCell<T,DESCRIPTOR>& cell,
  T j[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void NoLatticeDynamics<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<DESCRIPTOR >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, typename DESCRIPTOR>
void NoLatticeDynamics<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void NoLatticeDynamics<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void NoLatticeDynamics<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{ }

template<typename T, typename DESCRIPTOR>
void NoLatticeDynamics<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{ }

template<typename T, typename DESCRIPTOR>
void NoLatticeDynamics<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{ }

template<typename T, typename DESCRIPTOR>
void NoLatticeDynamics<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{ }

////////////////////// Class offDynamics ///////////////////////////

template<typename T, typename DESCRIPTOR>
OffDynamics<T,DESCRIPTOR>::OffDynamics(const T _location[DESCRIPTOR::d])
{
  this->getName() = "OffDynamics";  
  typedef DESCRIPTOR L;
  for (int iD = 0; iD < L::d; iD++) {
    location[iD] = _location[iD];
  }
  for (int iPop = 0; iPop < L::q; iPop++) {
    distances[iPop] = -1;
    velocityCoefficient[iPop] = 0;
    for (int iD = 0; iD < L::d; iD++) {
      boundaryIntersection[iPop][iD] = _location[iD];
      _u[iPop][iD] = T();
    }
  }
  _rho=T(1);
}

template<typename T, typename DESCRIPTOR>
OffDynamics<T,DESCRIPTOR>::OffDynamics(const T _location[DESCRIPTOR::d], T _distances[DESCRIPTOR::q])
{
  this->getName() = "OffDynamics";  
  typedef DESCRIPTOR L;
  for (int iD = 0; iD < L::d; iD++) {
    location[iD] = _location[iD];
  }
  for (int iPop = 0; iPop < L::q; iPop++) {
    distances[iPop] = _distances[iPop];
    velocityCoefficient[iPop] = 0;
    for (int iD = 0; iD < L::d; iD++) {
      boundaryIntersection[iPop][iD] = _location[iD] - _distances[iPop]*descriptors::c<L>(iPop,iD);
      _u[iPop][iD] = T();
    }
  }
  _rho=T(1);
}

template<typename T, typename DESCRIPTOR>
T OffDynamics<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  /*typedef DESCRIPTOR L;
  T rhoTmp = T();
  T counter = T();
  int counter2 = int();
  for (int iPop = 0; iPop < L::q; iPop++) {
    if (distances[iPop] != -1) {
      rhoTmp += (cell[iPop] + descriptors::t<T,L>(iPop))*descriptors::t<T,L>(iPop);
      counter += descriptors::t<T,L>(iPop);
      counter2++;
    }
  }
  //if (rhoTmp/counter + 1<0.1999) std::cout << rhoTmp/counter2 + 1 <<std::endl;
  //if (rhoTmp/counter + 1>1.001) std::cout << rhoTmp/counter2 + 1 <<std::endl;
  return rhoTmp/counter/counter;*/
  return _rho;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  typedef DESCRIPTOR L;
  for (int iD = 0; iD < L::d; iD++) {
    u[iD] = T();
  }
  int counter = 0;
  for (int iPop = 0; iPop < L::q; iPop++) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      for (int iD = 0; iD < L::d; iD++) {
        u[iD] += _u[iPop][iD];
      }
      counter++;
    }
  }
  if (counter!=0) {
    for (int iD = 0; iD < L::d; iD++) {
      u[iD] /= counter;
    }
  }
  return;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::setBoundaryIntersection(int iPop, T distance)
{
  /// direction points from the fluid node into the solid domain
  /// distance is the distance from the fluid node to the solid wall
  typedef DESCRIPTOR L;
  distances[iPop] = distance;
  for (int iD = 0; iD < L::d; iD++) {
    boundaryIntersection[iPop][iD] = location[iD] - distance*descriptors::c<L>(iPop,iD);
  }
}

template<typename T, typename DESCRIPTOR>
bool OffDynamics<T,DESCRIPTOR>::getBoundaryIntersection(int iPop, T intersection[DESCRIPTOR::d])
{
  typedef DESCRIPTOR L;
  if ( !util::nearZero(distances[iPop]+1) ) {
    for (int iD = 0; iD < L::d; iD++) {
      intersection[iD] = boundaryIntersection[iPop][iD];
    }
    return true;
  }
  return false;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  _rho=rho;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineRho(int iPop, T rho)
{
  _rho=rho;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  defineU(u);
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineU(const T u[DESCRIPTOR::d])
{
  typedef DESCRIPTOR L;
  for (int iPop = 0; iPop < L::q; iPop++) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      defineU(iPop, u);
    }
  }
}

/// Bouzidi velocity boundary condition formulas for the Coefficients:
/** 2*     invCs2*weight*(c,u)  for dist < 1/2
 *  1/dist*invCs2*weight*(c,u)  for dist >= 1/2
 */

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineU(
  int iPop, const T u[DESCRIPTOR::d])
{
  OLB_PRECONDITION(distances[iPop] != -1)
  typedef DESCRIPTOR L;
  velocityCoefficient[iPop] = 0;
  // scalar product of c(iPop) and u
  for (int sum = 0; sum < L::d; sum++) { // +/- problem because of first stream than postprocess
    velocityCoefficient[iPop] -= descriptors::c<L>(iPop,sum)*u[sum];
  }
  // compute summand for boundary condition
  velocityCoefficient[iPop] *= 2*descriptors::invCs2<T,L>() * descriptors::t<T,L>(iPop);

  for (int iD = 0; iD < L::d; iD++) {
    _u[iPop][iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
T OffDynamics<T,DESCRIPTOR>::getVelocityCoefficient(int iPop)
{
  return velocityCoefficient[iPop];
}

}

}

#endif
