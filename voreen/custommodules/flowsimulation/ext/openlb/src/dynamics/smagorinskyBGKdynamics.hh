/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2015 Mathias J. Krause, Jonas Latt, Patrick Nathen
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

/** \file
 * BGK Dynamics with adjusted omega -- generic implementation.
 */
#ifndef SMAGORINSKY_BGK_DYNAMICS_HH
#define SMAGORINSKY_BGK_DYNAMICS_HH

#include "smagorinskyBGKdynamics.h"
#include "core/cell.h"
#include "core/util.h"
#include "lbHelpers.h"
#include "math.h"

#include <complex> // For shear kalman Smagorinsky - Populations


namespace olb {

/// Smagorinsky Dynamics
template<typename T, template<typename U> class Lattice>
SmagorinskyDynamics<T,Lattice>::SmagorinskyDynamics(T smagoConst_)
  : smagoConst(smagoConst_), preFactor(computePreFactor())
{ }

template<typename T, template<typename U> class Lattice>
T SmagorinskyDynamics<T,Lattice>::computePreFactor()
{
  return (T)smagoConst*smagoConst*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyDynamics<T,Lattice>::getPreFactor()
{
  return preFactor;
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyDynamics<T,Lattice>::getSmagoConst()
{
  return smagoConst;
}

///////////////////////// ADM BGK /////////////////////////////

template<typename T, template<typename U> class Lattice>
ADMBGKdynamics<T,Lattice>::ADMBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_), omega(omega_)
{ }

template<typename T, template<typename U> class Lattice>
void ADMBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics )
{
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, *cell.getExternal(rhoIsAt), cell.getExternal(velocityBeginsAt), omega);
  statistics.incrementStats(*cell.getExternal(rhoIsAt), uSqr);
}

template<typename T, template<typename U> class Lattice>
void ADMBGKdynamics<T,Lattice>::staticCollide( Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, *cell.getExternal(rhoIsAt), u, omega);
  statistics.incrementStats(*cell.getExternal(rhoIsAt), uSqr);
}

///////////////////////// ForcedADM BGK /////////////////////////////

template<typename T, template<typename U> class Lattice>
ForcedADMBGKdynamics<T,Lattice>::ForcedADMBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, template<typename U> class Lattice>
void ForcedADMBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  OstreamManager clout(std::cout,"Forced ADM collide:");
  T rho, u[Lattice<T>::d], utst[Lattice<T>::d];

// this->momenta.computeAllMomenta(cell, rho, utst, pi);

  T* rho_fil = cell.getExternal(filRhoIsAt);
  T* u_filX = cell.getExternal(localFilVelXBeginsAt);
  T* u_filY = cell.getExternal(localFilVelYBeginsAt);
  T* u_filZ = cell.getExternal(localFilVelZBeginsAt);

  u[0] = *u_filX;/// *rho_fil;
  u[1] = *u_filY;/// *rho_fil;
  u[2] = *u_filZ;/// *rho_fil;

  T* force = cell.getExternal(forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }

  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, *rho_fil, u, omega);

  lbHelpers<T,Lattice>::addExternalForce(cell, u, omega);

  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void ForcedADMBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d];


  T* rho_fil = cell.getExternal(filRhoIsAt);
  T* u_filX = cell.getExternal(localFilVelXBeginsAt);
  T* u_filY = cell.getExternal(localFilVelYBeginsAt);
  T* u_filZ = cell.getExternal(localFilVelZBeginsAt);

  uTemp[0] = *u_filX;/// *rho_fil;
  uTemp[1] = *u_filY;/// *rho_fil;
  uTemp[2] = *u_filZ;/// *rho_fil;


  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, *rho_fil, uTemp, omega);
  statistics.incrementStats(rho, uSqr);
}

///////////////////// Class ConSmagorinskyBGKdynamics //////////////////////////

/////////////////// Consistent Smagorinsky BGK --> Malaspinas/Sagaut //////////////

template<typename T, template<typename U> class Lattice>
ConSmagorinskyBGKdynamics<T,Lattice>::ConSmagorinskyBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,Lattice>(omega_, momenta_, smagoConst_)
{ }

template<typename T, template<typename U> class Lattice>
T ConSmagorinskyBGKdynamics<T,Lattice>::computeEffectiveOmega (Cell<T,Lattice>& cell)
{

  T H[util::TensorVal<Lattice<T> >::n];
  T conSmagoR[Lattice<T>::q];
  T S[util::TensorVal<Lattice<T> >::n];
  T tau_mol = 1./this->getOmega();
  T cs2 = 1./Lattice<T>::invCs2;
  T smagoConst = this->getSmagoConst();

  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(2*PiNeqNormSqr);

  //Strain rate Tensor
  if (PiNeqNorm != 0) {
    T Phi = (-0.5*(-rho*tau_mol*cs2+sqrt(rho*rho*tau_mol*tau_mol*cs2*cs2+2.0*(smagoConst*smagoConst)*rho*PiNeqNorm))/(smagoConst*smagoConst*rho*PiNeqNorm));
    for (int n = 0; n < util::TensorVal<Lattice<T> >::n; ++n) {
      S[n] = Phi*pi[n];
    }
  } else {
    for (int n = 0; n < util::TensorVal<Lattice<T> >::n; ++n) {
      S[n] = 0;
    }
  }

  //Strain rate Tensor Norm
  T SNormSqr = S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    SNormSqr += S[2]*S[2] + S[3]*S[3] + 2.0*S[4]*S[4] + S[5]*S[5];
  }
  T SNorm    = sqrt(2*SNormSqr);

  //consistent Samagorinsky additional R term
  for (int q = 0; q < Lattice<T>::q; ++q) {
    T t = Lattice<T>::t[q]; //lattice weights

    //Hermite-Polynom H = c*c-cs^2*kronDelta
    H[0] = Lattice<T>::c[q][0]*Lattice<T>::c[q][0]-cs2;
    H[1] = Lattice<T>::c[q][0]*Lattice<T>::c[q][1];
    H[2] = Lattice<T>::c[q][1]*Lattice<T>::c[q][1]-cs2;//2D
    if (util::TensorVal<Lattice<T> >::n == 6) {
      H[2] = Lattice<T>::c[q][0]*Lattice<T>::c[q][2];//3D
      H[3] = Lattice<T>::c[q][1]*Lattice<T>::c[q][1]-cs2;
      H[4] = Lattice<T>::c[q][1]*Lattice<T>::c[q][2];
      H[5] = Lattice<T>::c[q][2]*Lattice<T>::c[q][2]-cs2;
    }

    //contraction or scalar product H*S
    T contractHS = H[0]*S[0] + 2.0*H[1]*S[1] + H[2]*S[2];
    if (util::TensorVal<Lattice<T> >::n == 6) {
      contractHS += H[2]*S[2] + H[3]*S[3] + 2.0*H[4]*S[4] + H[5]*S[5];
    }

    //additional term
    conSmagoR[q] = t*this->getPreFactor()*SNorm*contractHS;
  }
  return conSmagoR[0];
}

/////////////////// Consistent Strain Smagorinsky BGK --> Malaspinas/Sagaut //////////////

template<typename T, template<typename U> class Lattice>
ConStrainSmagorinskyBGKdynamics<T,Lattice>::ConStrainSmagorinskyBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,Lattice>(omega_, momenta_, smagoConst_)
{ }

template<typename T, template<typename U> class Lattice>
T ConStrainSmagorinskyBGKdynamics<T,Lattice>::computeEffectiveOmega(Cell<T,Lattice>& cell)
{
  T S[util::TensorVal<Lattice<T> >::n];
  T cs2 = 1./Lattice<T>::invCs2;
  T tau_mol = 1./this->getOmega();
  T smagoConst_ = this->getSmagoConst();

  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.0*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(2*PiNeqNormSqr);


  //Strain Tensor
  if ( !util::nearZero(PiNeqNorm) ) {
    T Phi = (-0.5*(-rho*tau_mol*cs2+sqrt(rho*rho*tau_mol*tau_mol*cs2*cs2+2.0*(smagoConst_*smagoConst_)*rho*PiNeqNorm))/(smagoConst_*smagoConst_*rho*PiNeqNorm));
    for (int n = 0; n < util::TensorVal<Lattice<T> >::n; ++n) {
      S[n] = Phi*pi[n];
    }
  } else {
    for (int n = 0; n < util::TensorVal<Lattice<T> >::n; ++n) {
      S[n] = 0;
    }
  }

  //Strain Tensor Norm
  T SNormSqr = S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    SNormSqr += S[2]*S[2] + S[3]*S[3] + 2.0*S[4]*S[4] + S[5]*S[5];
  }
  T SNorm    = sqrt(2*SNormSqr);

  /// Turbulent realaxation time
  T tau_turb = pow(smagoConst_,2)*SNorm/cs2;
  /// Effective realaxation time
  T tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

///////////////////////// DYNAMIC SMAGO BGK /////////////////////////////
template<typename T, template<typename U> class Lattice>
DynSmagorinskyBGKdynamics<T,Lattice>::DynSmagorinskyBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_)
  : SmagorinskyBGKdynamics<T,Lattice>(omega_, momenta_, T(0.))
{ }

template<typename T, template<typename U> class Lattice>
T DynSmagorinskyBGKdynamics<T,Lattice>::computeEffectiveOmega(Cell<T,Lattice>& cell)
{
  // computation of the relaxation time
  T v_t = 0;
  T* dynSmago = cell.getExternal(smagoConstIsAt);

  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  //v_t = *dynSmago*dx*dx*PiNeqNorm;
  v_t = *dynSmago*PiNeqNorm;
  T tau_t = 3*v_t;
  T tau_0 = 1/this->getOmega();
  T omega_new = 1/(tau_t+tau_0);
  return omega_new;
}

///////////////////SHEAR IMPROVED SMAGORINSKY//////////////////////////
template<typename T, template<typename U> class Lattice>
ShearSmagorinskyBGKdynamics<T,Lattice>::ShearSmagorinskyBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,Lattice>(omega_, momenta_, smagoConst_)
{ }

template<typename T, template<typename U> class Lattice>
void ShearSmagorinskyBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell,statistics.getTime());
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void ShearSmagorinskyBGKdynamics<T,Lattice>::staticCollide(Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell,statistics.getTime());
  T rho = this->_momenta.computeRho(cell);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T ShearSmagorinskyBGKdynamics<T,Lattice>::getEffectiveOmega(Cell<T,Lattice>& cell)
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.getExternal(avShearIsAt);
  T tau_0 = 1./this->getOmega();
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(this->getPreFactor()*PiNeqNorm_SISM/rho))-tau_0);

  T omega_new = 1./(tau_t+tau_0);

  return omega_new;
}

template<typename T, template<typename U> class Lattice>
T ShearSmagorinskyBGKdynamics<T,Lattice>::computeEffectiveOmega(Cell<T,Lattice>& cell, int iT)
{
  OstreamManager clout(std::cout,"shearImprovedCollide");

  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.getExternal(avShearIsAt);
  *avShear = (*avShear*iT + PiNeqNorm)/(iT+1);
  T tau_0 = 1./this->getOmega();
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(this->getPreFactor()*PiNeqNorm_SISM/rho))-tau_0);

  T omega_new = 1./(tau_t+tau_0);

  return omega_new;
}

///////////////////////// FORCED SHEAR SMAGO BGK /////////////////////////////
template<typename T, template<typename U> class Lattice>
ShearSmagorinskyForcedBGKdynamics<T,Lattice>::ShearSmagorinskyForcedBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_)
  : SmagorinskyForcedBGKdynamics<T,Lattice>(omega_, momenta_, smagoConst_)
{ }

template<typename T, template<typename U> class Lattice>
void ShearSmagorinskyForcedBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell, statistics.getTime());
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.getExternal(this->forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void ShearSmagorinskyForcedBGKdynamics<T,Lattice>::staticCollide( Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell, statistics.getTime());
  T rho = this->_momenta.computeRho(cell);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T ShearSmagorinskyForcedBGKdynamics<T,Lattice>::getEffectiveOmega(Cell<T,Lattice>& cell)
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.getExternal(avShearIsAt);
  T tau_0 = 1./this->getOmega();
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(this->getPreFactor()*PiNeqNorm_SISM/rho))-tau_0);

  T omega_new = 1./(tau_t+tau_0);

  return omega_new;
}

template<typename T, template<typename U> class Lattice>
T ShearSmagorinskyForcedBGKdynamics<T,Lattice>::computeEffectiveOmega(Cell<T,Lattice>& cell, int iT)
{
  OstreamManager clout(std::cout,"shearImprovedCollide");

  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.getExternal(avShearIsAt);
  *avShear = (*avShear*iT+PiNeqNorm)/(iT+1);

  T tau_0 = 1./this->getOmega();
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(this->getPreFactor()*PiNeqNorm_SISM/rho))-tau_0);

  T omega_new = 1./(tau_t+tau_0);

  return omega_new;
}

////////////////////// Class SmagorinskyBGKdynamics //////////////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
SmagorinskyBGKdynamics<T,Lattice>::SmagorinskyBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_)
  : SmagorinskyDynamics<T,Lattice>(smagoConst_),
    BGKdynamics<T,Lattice>(omega_,momenta_)
{ }

template<typename T, template<typename U> class Lattice>
void SmagorinskyBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell);
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyBGKdynamics<T,Lattice>::staticCollide(Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell);
  T rho = this->_momenta.computeRho(cell);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  BGKdynamics<T,Lattice>::setOmega(omega_);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyBGKdynamics<T,Lattice>::getEffectiveOmega(Cell<T,Lattice>& cell)
{
  T newOmega = computeEffectiveOmega(cell);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyBGKdynamics<T,Lattice>::computeEffectiveOmega(Cell<T,Lattice>& cell)
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /this->getOmega();
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + this->getPreFactor()/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  T tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;

}

///////////////////////// FORCED SMAGO BGK /////////////////////////////
template<typename T, template<typename U> class Lattice>
SmagorinskyForcedBGKdynamics<T,Lattice>::SmagorinskyForcedBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_)
  : SmagorinskyDynamics<T,Lattice>(smagoConst_),
    ForcedBGKdynamics<T,Lattice>(omega_,momenta_)
{ }

template<typename T, template<typename U> class Lattice>
void SmagorinskyForcedBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell);
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.getExternal(this->forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyForcedBGKdynamics<T,Lattice>::staticCollide( Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell);
  T rho = this->_momenta.computeRho(cell);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyForcedBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  ForcedBGKdynamics<T,Lattice>::setOmega(omega_);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyForcedBGKdynamics<T,Lattice>::getEffectiveOmega(Cell<T,Lattice>& cell)
{
  T newOmega = computeEffectiveOmega(cell);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyForcedBGKdynamics<T,Lattice>::computeEffectiveOmega(Cell<T,Lattice>& cell)
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /this->getOmega();
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + this->getPreFactor()/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  T tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

///////////////////////// FORCED Linear Velocity SMAGO BGK /////////////////////////////
template<typename T, template<typename U> class Lattice>
SmagorinskyLinearVelocityForcedBGKdynamics<T,Lattice>::SmagorinskyLinearVelocityForcedBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_)
  : SmagorinskyForcedBGKdynamics<T,Lattice>(omega_, momenta_, smagoConst_)
{ }

template<typename T, template<typename U> class Lattice>
void SmagorinskyLinearVelocityForcedBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeEffectiveOmega(cell);
  T* force = cell.getExternal(this->forceBeginsAt);
  int nDim = Lattice<T>::d;
  T forceSave[nDim];
  // adds a+Bu to force, where
  //   d=2: a1=v[0], a2=v[1], B11=v[2], B12=v[3], B21=v[4], B22=v[5]
  //   d=2: a1=v[0], a2=v[1], a3=v[2], B11=v[3], B12=v[4], B13=v[5], B21=v[6], B22=v[7], B23=v[8], B31=v[9], B32=v[10], B33=v[11]
  T* v = cell.getExternal(Lattice<T>::ExternalField::vBeginsAt);
  for (int iDim=0; iDim<nDim; ++iDim) {
    forceSave[iDim] = force[iDim];
    force[iDim] += v[iDim];
    for (int jDim=0; jDim<nDim; ++jDim) {
      force[iDim] += v[jDim + iDim*nDim + nDim]*u[jDim];
    }
  }
  for (int iVel=0; iVel<nDim; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }

  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
  // Writing back to froce fector
  for (int iVel=0; iVel<nDim; ++iVel) {
    force[iVel] = forceSave[iVel];
  }
}

////////////////////// Class KrauseBGKdynamics //////////////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
KrauseBGKdynamics<T,Lattice>::KrauseBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,Lattice>(omega_, momenta_, smagoConst_),
    preFactor(computePreFactor() )
{ }

template<typename T, template<typename U> class Lattice>
void KrauseBGKdynamics<T,Lattice>::collide(Cell<T,Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  T newOmega[Lattice<T>::q];
  this->_momenta.computeRhoU(cell, rho, u);
  computeEffectiveOmega(this->getOmega(), cell, preFactor, rho, u, newOmega);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void KrauseBGKdynamics<T,Lattice>::staticCollide(Cell<T,Lattice>& cell,
    const T u[Lattice<T>::d], LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  T newOmega[Lattice<T>::q];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  computeEffectiveOmega(this->getOmega(), cell, preFactor, rho, uTemp, newOmega);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T KrauseBGKdynamics<T,Lattice>::getEffectiveOmega(Cell<T,Lattice>& cell)
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  T newOmega[Lattice<T>::q];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  computeEffectiveOmega(this->getOmega(), cell, preFactor, rho, uTemp, newOmega);
  T newOmega_average = 0.;
  for (int iPop=0; iPop<Lattice<T>::q; iPop++) {
    newOmega_average += newOmega[iPop];
  }
  newOmega_average /= Lattice<T>::q;
  return newOmega_average;
}

template<typename T, template<typename U> class Lattice>
T KrauseBGKdynamics<T,Lattice>::computePreFactor()
{
  return (T)this->getSmagoConst()*this->getSmagoConst()*3*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}

template<typename T, template<typename U> class Lattice>
void KrauseBGKdynamics<T,Lattice>::computeEffectiveOmega(T omega0, Cell<T,Lattice>& cell, T preFactor_, T rho,
    T u[Lattice<T>::d], T newOmega[Lattice<T>::q])
{
  T uSqr = u[0]*u[0];
  for (int iDim=0; iDim<Lattice<T>::d; iDim++) {
    uSqr += u[iDim]*u[iDim];
  }
  /// Molecular realaxation time
  T tau_mol = 1./omega0;

  for (int iPop=0; iPop<Lattice<T>::q; iPop++) {
    T fNeq = std::fabs(cell[iPop] - lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr));
    /// Turbulent realaxation time
    T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor_/rho*fNeq) - tau_mol);
    /// Effective realaxation time
    T tau_eff = tau_mol + tau_turb;
    newOmega[iPop] = 1./tau_eff;
  }
}

////////////////////// Class WALEBGKdynamics //////////////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
WALEBGKdynamics<T,Lattice>::WALEBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,Lattice>(omega_, momenta_, smagoConst_)
{
  this->preFactor =  this->getSmagoConst()*this->getSmagoConst();
}

template<typename T, template<typename U> class Lattice>
T WALEBGKdynamics<T,Lattice>::computePreFactor()
{
  return (T)this->getSmagoConst()*this->getSmagoConst();
}

template<typename T, template<typename U> class Lattice>
T WALEBGKdynamics<T,Lattice>::computeEffectiveOmega(Cell<T,Lattice>& cell_)
{
  // velocity gradient tensor
  T g[3][3];
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      g[i][j] = *(cell_.getExternal(Lattice<T>::ExternalField::veloGradIsAt)+(i*3 + j));
    }
  }
  // strain rate tensor
  T s[3][3];
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      s[i][j] = (g[i][j] + g[j][i]) / 2.;
    }
  }
  // traceless symmetric part of the square of the velocity gradient tensor
  T G[3][3];
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      G[i][j] = 0.;
    }
  }

  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      for ( int k = 0; k < 3; k++) {
        G[i][j] += (g[i][k]*g[k][j] + g[j][k]*g[k][i]) / 2.;  // The change
      }
    }
  }

  T trace = 0.;
  for ( int i = 0; i < 3; i++) {
    trace += (1./3.) * g[i][i] * g[i][i];
  }

  for ( int i = 0; i < 3; i++) {
    G[i][i] -= trace;
  }


  // inner product of the traceless symmetric part of the square of the velocity gradient tensor
  T G_ip = 0;
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      G_ip = G[i][j] * G[i][j];
    }
  }

  // inner product of the strain rate
  T s_ip = 0;
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      s_ip = s[i][j] * s[i][j];
    }
  }

  // Turbulent relaxation time
  T tau_turb = 3. * this->getPreFactor() * (pow(G_ip,1.5) / (pow(s_ip,2.5) + pow(G_ip,1.25)));
  if ((pow(s_ip,2.5) + pow(G_ip,1.25)) == 0) {
    tau_turb = 0.;
  }

  // Physical turbulent viscosity must be equal or higher that zero
  if (tau_turb < 0.) {
    tau_turb = 0.;
  }

  /// Molecular relaxation time
  T tau_mol = 1. /this->getOmega();

  /// Effective relaxation time
  T tau_eff = tau_mol + tau_turb;
  T omega_new = 1. / tau_eff;

  return omega_new;

}

//////////////// Class ShearKalmanFDSmagorinskyBGKdynamics ///////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
FDKalmanShearSmagorinskyBGKdynamics<T,Lattice>::FDKalmanShearSmagorinskyBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_,  T u_char_lat, T f_char_lat)
  : SmagorinskyBGKdynamics<T,Lattice>(omega_, momenta_, smagoConst_),
    VarInVelKal(pow(2.0*(4.0 * std::atan(1.0))*(u_char_lat*f_char_lat)/sqrt(3),2)),
    UCharLat(u_char_lat)
{

  this->preFactor = this->getSmagoConst() * this->getSmagoConst() * Lattice<T>::invCs2;

}

template<typename T, template<typename U> class Lattice>
T FDKalmanShearSmagorinskyBGKdynamics<T,Lattice>::getEffectiveOmega(Cell<T,Lattice>& cell)
{
  T FNSR;
  computeNormStrainRate(cell, Lattice<T>::ExternalField::FilteredvelGradIsAt, FNSR);

  T INSR;
  computeNormStrainRate(cell, Lattice<T>::ExternalField::FilteredvelGradIsAt, INSR);

  // Turbulent relaxation time
  T tau_turb = 0.;
  if (INSR > FNSR) {
    tau_turb = this->getPreFactor() * (INSR - FNSR);
  }

  /// Molecular relaxation time
  T tau_mol = 1. /this->getOmega();

  /// Effective Omega
  T omega_new = 1. / tau_mol + tau_turb;

  return omega_new;
}

template<typename T, template<typename U> class Lattice>
T FDKalmanShearSmagorinskyBGKdynamics<T,Lattice>::computePreFactor()
{
  return this->getSmagoConst()*this->getSmagoConst()*Lattice<T>::invCs2;
}

template<typename T, template<typename U> class Lattice>
T FDKalmanShearSmagorinskyBGKdynamics<T,Lattice>::computeOmega(Cell<T,Lattice>& cell)
{
  OstreamManager clout(std::cout,"shearImprovedKalmanFDCollide");

  // Kalman procedure to update the filtered velocity
  KalmanStep(cell);

  // Norm of filtered Strain Rate
  T FNSR;
  computeNormStrainRate(cell, Lattice<T>::ExternalField::FilteredvelGradIsAt, FNSR);

  // Norm of Instantaneous Strain Rate
  T INSR;
  computeNormStrainRate(cell, Lattice<T>::ExternalField::FilteredvelGradIsAt, INSR);

  T tau_turb = 0.;
  if (INSR > FNSR) {
    tau_turb = this->getPreFactor() * (INSR - FNSR);
  }

  /// Molecular relaxation time
  T tau_mol = 1. /this->getOmega();

  /// Effective Omega
  T omega_new = 1. / tau_mol + tau_turb;

  return omega_new;

}

template<typename T, template<typename U> class Lattice>
void FDKalmanShearSmagorinskyBGKdynamics<T,Lattice>::computeNormStrainRate(Cell<T,Lattice>& cell, int PosVelGrad, T& NormStrainRate)
{
  int Dim = Lattice<T>::d;
  // Velocity gradient in 2D-3D
  T VG[Dim][Dim];
  for ( int i = 0; i < Dim; i++) {
    for ( int j = 0; j < Dim; j++) {
      VG[i][j] = *(cell.getExternal(PosVelGrad)+(i*Dim + j));
    }
  }
  // Strain rate tensor
  T S[Dim][Dim];
  for ( int i = 0; i < Dim; i++) {
    for ( int j = 0; j < Dim; j++) {
      S[i][j] = (VG[i][j] + VG[j][i]) / 2.;
    }
  }
  // inner product of the strain tensor
  T SIP = 0;
  for ( int i = 0; i < Dim; i++) {
    for ( int j = 0; j < Dim; j++) {
      SIP = S[i][j] * S[i][j];
    }
  }

  // Norm of the strain rate tensor
  NormStrainRate = sqrt(2. * SIP);

}

template<typename T, template<typename U> class Lattice>
void FDKalmanShearSmagorinskyBGKdynamics<T,Lattice>::KalmanStep(Cell<T,Lattice>& cell)
{
  // 1. Prediction Step
  T* ErrorCovariance = cell.getExternal(Lattice<T>::ExternalField::ErrorCovarianceIsAt);
  ErrorCovariance[0] += VarInVelKal;

  // 2. Update Step
  // 2.1. Smooothing Factor : K
  T* Variance = cell.getExternal(Lattice<T>::ExternalField::VarianceIsAt);
  T K = ErrorCovariance[0]/(ErrorCovariance[0] + Variance[0]);

  // 2.2. Kalman filtered Velocity -> Kalman filtered Populations
  //T* KalmanPopulation = cell.getExternal(FilteredPopulationIsAt);
  T j[Lattice<T>::d] = {0., 0., 0.};
  cell.getDynamics()->computeJ(cell, j);
  T rho = cell.computeRho();

  T* KalmanVel = cell.getExternal(Lattice<T>::ExternalField::velocityIsAt);
  for (int iVel=0; iVel<Lattice<T>::d; iVel++) {
    KalmanVel[iVel] = (KalmanVel[iVel] * (1-K)) + (K * (j[iVel]/rho));
  }

  // 2.3. Error covariance : P
  ErrorCovariance[0] *= (1-K);

  // 3. Adapt Step
  T epsilon = 0.1;
  T KalU_InstU[Lattice<T>::d] = {0., 0., 0.};
  for (int iVel=0; iVel < Lattice<T>::d; ++iVel) {
    KalU_InstU[iVel] = KalmanVel[iVel]-(j[iVel]/rho);
  }

  Variance[0] = std::max(UCharLat*util::normSqr<T,Lattice<T>::d>(KalU_InstU),epsilon*pow(UCharLat,2));
}



//////////////////////////////////////////////////////////////////////////////
//////////// Shear Improved - Kalman Filter - Smagorinsky BGK ////////////////
template<typename T, template<typename U> class Lattice>
ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::ShearKalmanSmagorinskyBGKdynamics(T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T u_char_lat, T f_char_lat)
  : SmagorinskyBGKdynamics<T,Lattice>(omega_, momenta_, smagoConst_),
    VarInVelKal(pow(2.0*(4.0 * std::atan(1.0))*(u_char_lat*f_char_lat)/sqrt(3),2)),
    UCharLat(u_char_lat)
{ }

template<typename T, template<typename U> class Lattice>
T ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::getEffectiveOmega(Cell<T,Lattice>& cell)
{
  OstreamManager clout(std::cout,"shearImprovedKalmanCollide");

  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [n+1][n+1]
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  cell.computeAllMomenta(rho, u, pi);
  T PiNeqNormSqrn1;
  computeNormSOM(pi, rho, PiNeqNormSqrn1);

  // Compute the norm of second moment of non-equilibrium filtered distribution <n+1><n+1>
  // Filtered Stress at time n+1
  T KalmanPiNeqNormSqrN1, KalmanUN1[Lattice<T>::d], KalmanPiNeqN1[util::TensorVal<Lattice<T> >::n];
  computeKalmanUStress(cell,KalmanUN1,KalmanPiNeqN1);
  computeNormSOM(KalmanPiNeqN1, KalmanPiNeqNormSqrN1);

  T tau_mol = 1./this->getOmega();
  T tau_sgs = T(0.);
  if (PiNeqNormSqrn1 > KalmanPiNeqNormSqrN1) {
    tau_sgs = 0.5*( ( pow(tau_mol,2.0) + ( this->getPreFactor()*( sqrt( PiNeqNormSqrn1)-sqrt(KalmanPiNeqNormSqrN1) ) ) ) - tau_mol );
  }

  T EffectiveOmega = 1.0/(tau_mol + tau_sgs);

  return EffectiveOmega;
}

template<typename T, template<typename U> class Lattice>
T ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeEffectiveOmega(Cell<T,Lattice>& cell)
{
  OstreamManager clout(std::cout,"shearImprovedKalmanCollide");

  // Update the filtered velocity wit a Kalman procedure
  KalmanStep(cell);

  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [n+1][n+1]
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  cell.computeAllMomenta(rho, u, pi);
  T PiNeqNormSqrn1;
  computeNormSOM(pi, rho, PiNeqNormSqrn1);

  // Compute the norm of second moment of non-equilibrium filtered distribution <n+1><n+1>
  // Filtered Stress at time n+1
  T KalmanPiNeqNormSqrN1, KalmanUN1[Lattice<T>::d], KalmanPiNeqN1[util::TensorVal<Lattice<T> >::n];
  computeKalmanUStress(cell,KalmanUN1,KalmanPiNeqN1);
  computeNormSOM(KalmanPiNeqN1, KalmanPiNeqNormSqrN1);

  T tau_mol = 1./this->getOmega();
  T tau_sgs = T(0.);
  if (PiNeqNormSqrn1 > KalmanPiNeqNormSqrN1) {
    tau_sgs = 0.5*( ( pow(tau_mol,2.0) + ( this->getPreFactor()*( sqrt( PiNeqNormSqrn1)-sqrt(KalmanPiNeqNormSqrN1) ) ) ) - tau_mol );
  }

  T EffectiveOmega = 1.0/(tau_mol + tau_sgs);

  return EffectiveOmega;
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::KalmanStep(Cell<T,Lattice>& cell)
{
  // The Kalman filter procedure //
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  cell.computeAllMomenta(rho, u, pi);

  T* KalmanPopulation = cell.getExternal(FilteredPopulationIsAt);
  if (KalmanPopulation[0] == (T)-1.0){
    for (int iPop=0; iPop<Lattice<T>::q; iPop++){
      KalmanPopulation[iPop] =  cell[iPop]/rho;
    }
  }

  // 1. Prediction Step
  T* ErrorCovariance = cell.getExternal(ErrorCovarianceIsAt);
  *ErrorCovariance += VarInVelKal;

  // 2. Update Step
  // 2.1. Smoothing Factor : K
  T* Variance = cell.getExternal(VarianceIsAt);
  T K = *ErrorCovariance/(*ErrorCovariance + *Variance);

  // 2.2. Kalman filtered Velocity -> Kalman filtered Populations
  for (int iPop=0; iPop<Lattice<T>::q; iPop++) {
    KalmanPopulation[iPop] = (KalmanPopulation[iPop] * (1-K)) + (K * cell[iPop]/rho);
  }

  // 2.3. Error covariance : P
  *ErrorCovariance *= (1-K);

  // 3. Adapt Step
  T epsilon = T(0.1);
  T KalU_InstU[Lattice<T>::d];
    // Filtered Stress at time n+1
    T KalmanUN1[Lattice<T>::d];
    computeKalmanU(cell, KalmanUN1);
  for (int iD=0; iD < Lattice<T>::d; ++iD) {
    KalU_InstU[iD] = KalmanUN1[iD]-u[iD];
  }

  *Variance = std::max(UCharLat*util::normSqr<T,Lattice<T>::d>(KalU_InstU),epsilon*pow(UCharLat,2));

  // Fin filtering procedure //
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeKalmanUStress(Cell<T,Lattice>& cell,
      T (&KalmanU)[Lattice<T>::d],T (&KalmanPi)[util::TensorVal<Lattice<T> >::n] )
{
  computeKalmanU(cell,KalmanU);
  computeKalmanStress(cell,KalmanU,KalmanPi);
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeKalmanU(Cell<T,Lattice>& cell, T (&KalmanU)[Lattice<T>::d])
{
  T* KalmanPopulation = cell.getExternal(FilteredPopulationIsAt);
  for (int iD=0; iD < Lattice<T>::d; ++iD) {
    KalmanU[iD] = T();
  }
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    for (int iD=0; iD < Lattice<T>::d; ++iD) {
      KalmanU[iD] += KalmanPopulation[iPop]*Lattice<T>::c[iPop][iD];
    }
  }
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeKalmanStress(Cell<T,Lattice>& cell,
      T (&KalmanU)[Lattice<T>::d],T (&KalmanPi)[util::TensorVal<Lattice<T> >::n] )
{
  T* KalmanPopulation = cell.getExternal(FilteredPopulationIsAt);

  T rhoRelative = T(0.);
  for (int iPop = 0; iPop < Lattice<T>::q; iPop++) {
    rhoRelative += KalmanPopulation[iPop];
  }

  int iPi = 0;
  for (int iAlpha=0; iAlpha < Lattice<T>::d; ++iAlpha) {
    for (int iBeta=iAlpha; iBeta < Lattice<T>::d; ++iBeta) {
      KalmanPi[iPi] = T();
      for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
        KalmanPi[iPi] += Lattice<T>::c[iPop][iAlpha]*Lattice<T>::c[iPop][iBeta] * KalmanPopulation[iPop];
      }
      // stripe off equilibrium contribution
      KalmanPi[iPi] -= KalmanU[iAlpha]*KalmanU[iBeta];
      if (iAlpha==iBeta) {
        KalmanPi[iPi] -= rhoRelative/Lattice<T>::invCs2;
      }
      ++iPi;
    }
  }
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeAndupdateTauSgs(Cell<T,Lattice>& cell,
      T rho, T pi[util::TensorVal<Lattice<T> >::n], T KalmanPiNeqN[util::TensorVal<Lattice<T> >::n],
      T KalmanPiNeqN1[util::TensorVal<Lattice<T> >::n], T K, T &tau_sgs)
{
  // Compute the norm of second moment of non-equilibrium filtered distribution <n><n>
  T KalmanPiNeqNormSqrN;
  computeNormSOM(KalmanPiNeqN, KalmanPiNeqNormSqrN);

  // Compute the norm of cross second moment of non-equilibrium filtered-Instantaneous distribution <n>[n+1]
  T KalmanInstPiNeqNormSqrN;
  computeNormSOM(KalmanPiNeqN, pi, rho, KalmanInstPiNeqNormSqrN);

  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [n+1][n+1]
  T PiNeqNormSqrn;
  computeNormSOM(pi, rho, PiNeqNormSqrn);

  computeTauSgs(cell, rho, KalmanPiNeqNormSqrN, KalmanInstPiNeqNormSqrN, PiNeqNormSqrn, K, tau_sgs);

  // Compute the norm of second moment of non-equilibrium filtered distribution <n+1><n+1>
  T KalmanPiNeqNormSqrN1;
  computeNormSOM(KalmanPiNeqN1, KalmanPiNeqNormSqrN1);

  updateTauSgsKalman(cell, KalmanPiNeqNormSqrN, KalmanInstPiNeqNormSqrN, PiNeqNormSqrn, KalmanPiNeqNormSqrN1, K, tau_sgs);
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeNormSOM(T pi[util::TensorVal<Lattice<T> >::n], T &piNorm)
{
  // Compute the norm of second moment of non-equilibrium filtered distribution <-><->
  piNorm = pow(pi[0],2) + 2.0*pow(pi[1],2) + pow(pi[2],2);
  if (util::TensorVal<Lattice<T> >::n == 6) {
    piNorm += pow(pi[2],2) + pow(pi[3],2) + 2*pow(pi[4],2) + pow(pi[5],2);
  }
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeNormSOM(T pi1[util::TensorVal<Lattice<T> >::n],
       T pi2[util::TensorVal<Lattice<T> >::n], T rho, T &piNorm)
{
  // Compute the norm of cross second moment of non-equilibrium filtered-Instantaneous distribution <->[-]
  piNorm = pi1[0]*pi2[0] + 2.0*pi1[1]*pi2[1] + pi1[2]*pi2[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    piNorm += pi1[2]*pi2[2] + pi1[3]*pi2[3] + 2*pi1[4]*pi2[4] + pi1[5]*pi2[5];
  }
  piNorm /= rho;
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeNormSOM(T pi[util::TensorVal<Lattice<T> >::n], T rho, T &piNorm)
{
  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [-][-]
  computeNormSOM(pi, piNorm);
  piNorm /= pow(rho,2.0);
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeTauSgs(Cell<T,Lattice>& cell,
      T rho, T KalmanPiNeqNormSqr, T KalmanInstPiNeqNormSqr, T PiNeqNormSqr, T K, T &tau_sgs)
{
  T tau_mol = this->getOmega();
  T tau_eff_n = tau_mol + *(cell.getExternal(TauSgsIsAt));
  T F = 1.0/(pow(this->getSmagoConst()/Lattice<T>::invCs2,2.0)/(2.0*tau_eff_n));

  // Coefficients of 4th polynomial : Ax⁴ + Bx³ + Cx² + Dx + E = 0
  T A = pow(F,2.0);
  T B = 2.0*A*tau_mol;

  T C = F*pow(tau_mol,2.0)
        - (2.0*F*sqrt(2.0*PiNeqNormSqr)*tau_eff_n)
        - (pow((1.0-K),2.0)*2.0*KalmanPiNeqNormSqr);

  T D = -((2.0*F*sqrt(2.0*PiNeqNormSqr)*tau_mol*tau_eff_n)
           + (2.0*pow((1.0-K),2.0)*tau_mol*2.0*KalmanPiNeqNormSqr)
           + (2.0*K*(1.0-K)*2.0*KalmanInstPiNeqNormSqr*tau_eff_n)
         );

  T E = ((1-pow(K,2.0))*2.0*PiNeqNormSqr*pow(tau_eff_n,2.0))
        + (pow((1.0-K)*tau_mol,2.0)*2.0*KalmanPiNeqNormSqr)
        - (2.0*K*(1.0-K)*2.0*KalmanInstPiNeqNormSqr*tau_mol*tau_eff_n);

  std::complex<T> Roots[4];
  computeRoots4thPoly(A, B, C, D, E, Roots);

  tau_sgs = 0.;
  for ( int i = 0; i < 4; i++){
    if (std::imag(Roots[i]) == T(0.)) {
      if (std::real(Roots[i]) > tau_sgs){
        tau_sgs = std::real(Roots[i]);
      }
    }
  }

  /// Update the value of instantaneous effective omega
  //T* EffectiveOmega = cell.getExternal(EffectiveOmegaIsAt);
  //*EffectiveOmega = 1.0/(tau_mol + tau_sgs);

}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::computeRoots4thPoly(T A, T B, T C, T D, T E, std::complex<T> (&Roots)[4])
{
  T p = T((8.*A*C - 3.*pow(B,2.0))/(8.0*pow(A,2.0)));
  T q = T((pow(B,3.0) - 4.0*A*B*C + 8.0*pow(A,2.0)*D)/(8.0*pow(A,3.0)));

  T Delta0 = T(pow(C,2.0) - 3.0*B*D + 12.0*A*E);
  T Delta1 = 2.0*pow(C,3.0) - 9.0*B*C*D + 27*pow(B,2.0)*E + 27.0*A*pow(D,2.0) - 72.0*A*C*E;

  // Discriminant
  std::complex<T> Dis = (pow(Delta1,2.0) - 4.0*pow(Delta0,3.0))/(-27.0);

  std::complex<T> Q = pow((Delta1+sqrt(Dis*(-27.0)))/2.0,1.0/3.0);
  std::complex<T> S = 0.5*sqrt((-2.0*p/3.0)+((1.0/(3.0*A))*(Q+(Delta0/Q))));

  std::complex<T> cas1, cas2;
  for ( int i = 0; i < 2; i++) {
    for ( int j = 0; j < 2; j++) {
      cas1 = T(2*i-1);
      cas2 = T(2*j-1);
      Roots[2*i+j] = (-B/4*A) + cas1*S + (cas2*0.5*sqrt((-4.0*S*S)+(2.0*p)+(q/S)));
    }
  }
}

template<typename T, template<typename U> class Lattice>
void ShearKalmanSmagorinskyBGKdynamics<T,Lattice>::updateTauSgsKalman(Cell<T,Lattice>& cell, T NN, T Nn1, T n1n1, T N1N1, T K, T tau_sgs_n1)
{
  T tau_mol = this->getOmega();
  T* tau_sgs_N = cell.getExternal(TauSgsIsAt);

  //T tau_eff_N = tau_mol + *tau_sgs_N;
  //T tau_eff_n1 = *(cell.getExternal(EffectiveOmegaIsAt));
  T tau_eff_n1 = *tau_sgs_N;
//  T A;
//  if (!util::nearZero(tau_sgs_n1)) {
//    A = sqrt( pow((1-K) * (*tau_sgs_N) / tau_eff_N,2.0) * sqrt(NN)
//              + 2.0 * K * (1.0-K) * (*tau_sgs_N * tau_sgs_n1) * sqrt(Nn1) / (tau_eff_N * tau_eff_n1)
//              + pow(K * tau_sgs_n1 / tau_eff_n1,2.0) * sqrt(n1n1)
//           );
//  } else {
//    A = sqrt( pow((1-K) * (*tau_sgs_N) / tau_eff_N,2.0) * sqrt(NN));
//  }
//
//  A /= sqrt(N1N1);
//
//  *tau_sgs_N = tau_mol/((1.0/A)-1.0);
  *tau_sgs_N = (sqrt(2.0*N1N1)
                /(sqrt(2.0*n1n1)
                  -(2.0*pow(1./(this->getSmagoConst()*Lattice<T>::invCs2),2.0)*tau_sgs_n1*tau_eff_n1)
                 )
               )*tau_eff_n1;

  if ((*tau_sgs_N - tau_mol) > T(0.)){
    *tau_sgs_N -= tau_mol;
  } else {
    *tau_sgs_N = T(0.);
  }

}


}

#endif
