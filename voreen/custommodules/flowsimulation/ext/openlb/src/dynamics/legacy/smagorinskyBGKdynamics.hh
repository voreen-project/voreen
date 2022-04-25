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
#ifndef LEGACY_SMAGORINSKY_BGK_DYNAMICS_HH
#define LEGACY_SMAGORINSKY_BGK_DYNAMICS_HH

#include "smagorinskyBGKdynamics.h"
#include "core/cell.h"
#include "core/util.h"
#include "dynamics/lbm.h"
#include "math.h"


namespace olb {

/// Smagorinsky Dynamics
template<typename T, typename DESCRIPTOR>
SmagorinskyDynamics<T,DESCRIPTOR>::SmagorinskyDynamics(T smagoConst_)
  : smagoConst(smagoConst_), preFactor(computePreFactor())
{ }

template<typename T, typename DESCRIPTOR>
T SmagorinskyDynamics<T,DESCRIPTOR>::computePreFactor()
{
  return (T)smagoConst*smagoConst*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*2*util::sqrt(2);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyDynamics<T,DESCRIPTOR>::getPreFactor()
{
  return preFactor;
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyDynamics<T,DESCRIPTOR>::getSmagoConst()
{
  return smagoConst;
}

///////////////////////// FORCED Linear Velocity SMAGO BGK /////////////////////////////
template<typename T, typename DESCRIPTOR, typename MOMENTA>
SmagorinskyLinearVelocityForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::SmagorinskyLinearVelocityForcedBGKdynamics(
  T omega_, T smagoConst_)
  : SmagorinskyForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>(omega_, smagoConst_)
{ }

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> SmagorinskyLinearVelocityForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::collide(Cell<T,DESCRIPTOR>& cell)
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  MOMENTA().computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeEffectiveOmega(cell);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  int nDim = DESCRIPTOR::d;
  T forceSave[nDim];
  // adds a+Bu to force, where
  //   d=2: a1=v[0], a2=v[1], B11=v[2], B12=v[3], B21=v[4], B22=v[5]
  //   d=2: a1=v[0], a2=v[1], a3=v[2], B11=v[3], B12=v[4], B13=v[5], B21=v[6], B22=v[7], B23=v[8], B31=v[9], B32=v[10], B33=v[11]
  auto v = cell.template getFieldPointer<descriptors::V12>();
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

  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  lbm<DESCRIPTOR>::addExternalForce(cell, u, newOmega, rho);
  //statistics.incrementStats(rho, uSqr);
  // Writing back to froce fector
  for (int iVel=0; iVel<nDim; ++iVel) {
    force[iVel] = forceSave[iVel];
  }
}

//////////////// Class ShearKalmanFDSmagorinskyBGKdynamics ///////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::FDKalmanShearSmagorinskyBGKdynamics(
  T omega_,  T smagoConst_,  T u_char_lat, T f_char_lat)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>(omega_, smagoConst_),
    VarInVelKal(util::pow(2.0*(4.0 * util::atan(1.0))*(u_char_lat*f_char_lat)/util::sqrt(3),2)),
    UCharLat(u_char_lat)
{

  this->preFactor = this->getSmagoConst() * this->getSmagoConst() * descriptors::invCs2<T,DESCRIPTOR>();

}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::getEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T FNSR;
  computeNormStrainRate(cell, FNSR);

  T INSR;
  computeNormStrainRate(cell, INSR);

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

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computePreFactor()
{
  return this->getSmagoConst()*this->getSmagoConst()*descriptors::invCs2<T,DESCRIPTOR>();
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeOmega(Cell<T,DESCRIPTOR>& cell)
{
  OstreamManager clout(std::cout,"shearImprovedKalmanFDCollide");

  // Kalman procedure to update the filtered velocity
  KalmanStep(cell);

  // Norm of filtered Strain Rate
  T FNSR;
  computeNormStrainRate(cell, FNSR);

  // Norm of Instantaneous Strain Rate
  T INSR;
  computeNormStrainRate(cell, INSR);

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

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeNormStrainRate(Cell<T,DESCRIPTOR>& cell, T& NormStrainRate)
{
  int Dim = DESCRIPTOR::d;
  // Velocity gradient in 2D-3D
  T VG[Dim][Dim];
  for ( int i = 0; i < Dim; i++) {
    for ( int j = 0; j < Dim; j++) {
      VG[i][j] = cell.template getFieldPointer<descriptors::FILTERED_VEL_GRAD>()[i*Dim + j];
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
  NormStrainRate = util::sqrt(2. * SIP);

}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::KalmanStep(Cell<T,DESCRIPTOR>& cell)
{
  // 1. Prediction Step
  auto ErrorCovariance = cell.template getFieldPointer<descriptors::ERROR_COVARIANCE>();
  ErrorCovariance[0] += VarInVelKal;

  // 2. Update Step
  // 2.1. Smooothing Factor : K
  auto Variance = cell.template getFieldPointer<descriptors::VARIANCE>();
  T K = ErrorCovariance[0]/(ErrorCovariance[0] + Variance[0]);

  // 2.2. Kalman filtered Velocity -> Kalman filtered Populations
  //T* KalmanPopulation = cell.template getFieldPointer<descriptors::FILTERED_POPULATION>();
  T u[DESCRIPTOR::d] = {0., 0., 0.};
  cell.computeU(u);

  auto KalmanVel = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iVel=0; iVel<DESCRIPTOR::d; iVel++) {
    KalmanVel[iVel] = (KalmanVel[iVel] * (1-K)) + (K * u[iVel]);
  }

  // 2.3. Error covariance : P
  ErrorCovariance[0] *= (1-K);

  // 3. Adapt Step
  T epsilon = 0.1;
  T KalU_InstU[DESCRIPTOR::d] = {0., 0., 0.};
  for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
    KalU_InstU[iVel] = KalmanVel[iVel]-u[iVel];
  }

  Variance[0] = util::max(UCharLat*util::normSqr<T,DESCRIPTOR::d>(KalU_InstU),epsilon*util::pow(UCharLat,2));
}



//////////////////////////////////////////////////////////////////////////////
//////////// Shear Improved - Kalman Filter - Smagorinsky BGK ////////////////
template<typename T, typename DESCRIPTOR, typename MOMENTA>
ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::ShearKalmanSmagorinskyBGKdynamics(
  T omega_, T smagoConst_, T u_char_lat, T f_char_lat)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>(omega_, smagoConst_),
    VarInVelKal(util::pow(2.0*(4.0 * util::atan(1.0))*(u_char_lat*f_char_lat)/util::sqrt(3),2)),
    UCharLat(u_char_lat)
{ }

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::getEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  OstreamManager clout(std::cout,"shearImprovedKalmanCollide");

  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [n+1][n+1]
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  cell.computeAllMomenta(rho, u, pi);
  T PiNeqNormSqrn1;
  computeNormSOM(pi, rho, PiNeqNormSqrn1);

  // Compute the norm of second moment of non-equilibrium filtered distribution <n+1><n+1>
  // Filtered Stress at time n+1
  T KalmanPiNeqNormSqrN1, KalmanUN1[DESCRIPTOR::d], KalmanPiNeqN1[util::TensorVal<DESCRIPTOR >::n];
  computeKalmanUStress(cell,KalmanUN1,KalmanPiNeqN1);
  computeNormSOM(KalmanPiNeqN1, KalmanPiNeqNormSqrN1);

  T tau_mol = 1./this->getOmega();
  T tau_sgs = T(0.);
  if (PiNeqNormSqrn1 > KalmanPiNeqNormSqrN1) {
    tau_sgs = 0.5*( ( util::pow(tau_mol,2.0) + ( this->getPreFactor()*( util::sqrt( PiNeqNormSqrn1)-util::sqrt(KalmanPiNeqNormSqrN1) ) ) ) - tau_mol );
  }

  T EffectiveOmega = 1.0/(tau_mol + tau_sgs);

  return EffectiveOmega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  OstreamManager clout(std::cout,"shearImprovedKalmanCollide");

  // Update the filtered velocity wit a Kalman procedure
  KalmanStep(cell);

  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [n+1][n+1]
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  cell.computeAllMomenta(rho, u, pi);
  T PiNeqNormSqrn1;
  computeNormSOM(pi, rho, PiNeqNormSqrn1);

  // Compute the norm of second moment of non-equilibrium filtered distribution <n+1><n+1>
  // Filtered Stress at time n+1
  T KalmanPiNeqNormSqrN1, KalmanUN1[DESCRIPTOR::d], KalmanPiNeqN1[util::TensorVal<DESCRIPTOR >::n];
  computeKalmanUStress(cell,KalmanUN1,KalmanPiNeqN1);
  computeNormSOM(KalmanPiNeqN1, KalmanPiNeqNormSqrN1);

  T tau_mol = 1./this->getOmega();
  T tau_sgs = T(0.);
  if (PiNeqNormSqrn1 > KalmanPiNeqNormSqrN1) {
    tau_sgs = 0.5*( ( util::pow(tau_mol,2.0) + ( this->getPreFactor()*( util::sqrt( PiNeqNormSqrn1)-util::sqrt(KalmanPiNeqNormSqrN1) ) ) ) - tau_mol );
  }

  T EffectiveOmega = 1.0/(tau_mol + tau_sgs);

  return EffectiveOmega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::KalmanStep(Cell<T,DESCRIPTOR>& cell)
{
  // The Kalman filter procedure //
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  cell.computeAllMomenta(rho, u, pi);

  auto KalmanPopulation = cell.template getFieldPointer<descriptors::FILTERED_POPULATION>();
  if (KalmanPopulation[0] == (T)-1.0) {
    for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
      KalmanPopulation[iPop] =  cell[iPop]/rho;
    }
  }

  // 1. Prediction Step
  auto ErrorCovariance = cell.template getFieldPointer<descriptors::ERROR_COVARIANCE>();
  ErrorCovariance[0] += VarInVelKal;

  // 2. Update Step
  // 2.1. Smoothing Factor : K
  auto Variance = cell.template getFieldPointer<descriptors::VARIANCE>();
  T K = ErrorCovariance[0]/(ErrorCovariance[0] + Variance[0]);

  // 2.2. Kalman filtered Velocity -> Kalman filtered Populations
  for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
    KalmanPopulation[iPop] = (KalmanPopulation[iPop] * (1-K)) + (K * cell[iPop]/rho);
  }

  // 2.3. Error covariance : P
  *ErrorCovariance *= (1-K);

  // 3. Adapt Step
  T epsilon = T(0.1);
  T KalU_InstU[DESCRIPTOR::d];
  // Filtered Stress at time n+1
  T KalmanUN1[DESCRIPTOR::d];
  computeKalmanU(cell, KalmanUN1);
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    KalU_InstU[iD] = KalmanUN1[iD]-u[iD];
  }

  *Variance = util::max(UCharLat*util::normSqr<T,DESCRIPTOR::d>(KalU_InstU),epsilon*util::pow(UCharLat,2));

  // Fin filtering procedure //
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeKalmanUStress(Cell<T,DESCRIPTOR>& cell,
    T (&KalmanU)[DESCRIPTOR::d],T (&KalmanPi)[util::TensorVal<DESCRIPTOR >::n] )
{
  computeKalmanU(cell,KalmanU);
  computeKalmanStress(cell,KalmanU,KalmanPi);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeKalmanU(
    Cell<T,DESCRIPTOR>& cell, T (&KalmanU)[DESCRIPTOR::d])
{
  auto KalmanPopulation = cell.template getFieldPointer<descriptors::FILTERED_POPULATION>();
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    KalmanU[iD] = T();
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      KalmanU[iD] += KalmanPopulation[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
    }
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeKalmanStress(Cell<T,DESCRIPTOR>& cell,
    T (&KalmanU)[DESCRIPTOR::d],T (&KalmanPi)[util::TensorVal<DESCRIPTOR >::n] )
{
  auto KalmanPopulation = cell.template getFieldPointer<descriptors::FILTERED_POPULATION>();

  T rhoRelative = T(0.);
  for (int iPop = 0; iPop < DESCRIPTOR::q; iPop++) {
    rhoRelative += KalmanPopulation[iPop];
  }

  int iPi = 0;
  for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
    for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
      KalmanPi[iPi] = T();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        KalmanPi[iPi] += descriptors::c<DESCRIPTOR>(iPop,iAlpha)*descriptors::c<DESCRIPTOR>(iPop,iBeta) * KalmanPopulation[iPop];
      }
      // stripe off equilibrium contribution
      KalmanPi[iPi] -= KalmanU[iAlpha]*KalmanU[iBeta];
      if (iAlpha==iBeta) {
        KalmanPi[iPi] -= rhoRelative/descriptors::invCs2<T,DESCRIPTOR>();
      }
      ++iPi;
    }
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeAndupdateTauSgs(Cell<T,DESCRIPTOR>& cell,
    T rho, T pi[util::TensorVal<DESCRIPTOR >::n], T KalmanPiNeqN[util::TensorVal<DESCRIPTOR >::n],
    T KalmanPiNeqN1[util::TensorVal<DESCRIPTOR >::n], T K, T &tau_sgs)
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

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeNormSOM(T pi[util::TensorVal<DESCRIPTOR >::n], T &piNorm)
{
  // Compute the norm of second moment of non-equilibrium filtered distribution <-><->
  piNorm = util::pow(pi[0],2) + 2.0*util::pow(pi[1],2) + util::pow(pi[2],2);
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    piNorm += util::pow(pi[2],2) + util::pow(pi[3],2) + 2*util::pow(pi[4],2) + util::pow(pi[5],2);
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeNormSOM(T pi1[util::TensorVal<DESCRIPTOR >::n],
    T pi2[util::TensorVal<DESCRIPTOR >::n], T rho, T &piNorm)
{
  // Compute the norm of cross second moment of non-equilibrium filtered-Instantaneous distribution <->[-]
  piNorm = pi1[0]*pi2[0] + 2.0*pi1[1]*pi2[1] + pi1[2]*pi2[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    piNorm += pi1[2]*pi2[2] + pi1[3]*pi2[3] + 2*pi1[4]*pi2[4] + pi1[5]*pi2[5];
  }
  piNorm /= rho;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeNormSOM(T pi[util::TensorVal<DESCRIPTOR >::n], T rho, T &piNorm)
{
  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [-][-]
  computeNormSOM(pi, piNorm);
  piNorm /= util::pow(rho,2.0);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeTauSgs(Cell<T,DESCRIPTOR>& cell,
    T rho, T KalmanPiNeqNormSqr, T KalmanInstPiNeqNormSqr, T PiNeqNormSqr, T K, T &tau_sgs)
{
  T tau_mol = this->getOmega();
  T tau_eff_n = tau_mol + cell.template getField<descriptors::TAU_SGS>();
  T F = 1.0/(util::pow(this->getSmagoConst()/descriptors::invCs2<T,DESCRIPTOR>(),2.0)/(2.0*tau_eff_n));

  // Coefficients of 4th polynomial : Ax^4 + Bx^3 + Cx^2 + Dx + E = 0
  T A = util::pow(F,2.0);
  T B = 2.0*A*tau_mol;

  T C = F*util::pow(tau_mol,2.0)
        - (2.0*F*util::sqrt(2.0*PiNeqNormSqr)*tau_eff_n)
        - (util::pow((1.0-K),2.0)*2.0*KalmanPiNeqNormSqr);

  T D = -((2.0*F*util::sqrt(2.0*PiNeqNormSqr)*tau_mol*tau_eff_n)
          + (2.0*util::pow((1.0-K),2.0)*tau_mol*2.0*KalmanPiNeqNormSqr)
          + (2.0*K*(1.0-K)*2.0*KalmanInstPiNeqNormSqr*tau_eff_n)
         );

  T E = ((1-util::pow(K,2.0))*2.0*PiNeqNormSqr*util::pow(tau_eff_n,2.0))
        + (util::pow((1.0-K)*tau_mol,2.0)*2.0*KalmanPiNeqNormSqr)
        - (2.0*K*(1.0-K)*2.0*KalmanInstPiNeqNormSqr*tau_mol*tau_eff_n);

  std::complex<T> Roots[4];
  computeRoots4thPoly(A, B, C, D, E, Roots);

  tau_sgs = 0.;
  for ( int i = 0; i < 4; i++) {
    if (std::imag(Roots[i]) == T(0.)) {
      if (std::real(Roots[i]) > tau_sgs) {
        tau_sgs = std::real(Roots[i]);
      }
    }
  }

  /// Update the value of instantaneous effective omega
  //T* EffectiveOmega = cell.template getFieldPointer<descriptors::EFFECTIVE_OMEGA>();
  //*EffectiveOmega = 1.0/(tau_mol + tau_sgs);

}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeRoots4thPoly(T A, T B, T C, T D, T E, std::complex<T> (&Roots)[4])
{
  T p = T((8.*A*C - 3.*util::pow(B,2.0))/(8.0*util::pow(A,2.0)));
  T q = T((util::pow(B,3.0) - 4.0*A*B*C + 8.0*util::pow(A,2.0)*D)/(8.0*util::pow(A,3.0)));

  T Delta0 = T(util::pow(C,2.0) - 3.0*B*D + 12.0*A*E);
  T Delta1 = 2.0*util::pow(C,3.0) - 9.0*B*C*D + 27*util::pow(B,2.0)*E + 27.0*A*util::pow(D,2.0) - 72.0*A*C*E;

  // Discriminant
  std::complex<T> Dis = (util::pow(Delta1,2.0) - 4.0*util::pow(Delta0,3.0))/(-27.0);

  std::complex<T> Q = util::pow((Delta1+util::sqrt(Dis*(-27.0)))/2.0,1.0/3.0);
  std::complex<T> S = 0.5*util::sqrt((-2.0*p/3.0)+((1.0/(3.0*A))*(Q+(Delta0/Q))));

  std::complex<T> cas1, cas2;
  for ( int i = 0; i < 2; i++) {
    for ( int j = 0; j < 2; j++) {
      cas1 = T(2*i-1);
      cas2 = T(2*j-1);
      Roots[2*i+j] = (-B/4*A) + cas1*S + (cas2*0.5*util::sqrt((-4.0*S*S)+(2.0*p)+(q/S)));
    }
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR,MOMENTA>::updateTauSgsKalman(Cell<T,DESCRIPTOR>& cell,
T NN, T Nn1, T n1n1, T N1N1, T K, T tau_sgs_n1)
{
  T tau_mol = this->getOmega();
  auto tau_sgs_N = cell.template getFieldPointer<descriptors::TAU_SGS>();

  //T tau_eff_N = tau_mol + *tau_sgs_N;
  //T tau_eff_n1 = *(cell[EffectiveOmegaIsAt]);
  T tau_eff_n1 = *tau_sgs_N;
//  T A;
//  if (!util::nearZero(tau_sgs_n1)) {
//    A = util::sqrt( util::pow((1-K) * (*tau_sgs_N) / tau_eff_N,2.0) * util::sqrt(NN)
//              + 2.0 * K * (1.0-K) * (*tau_sgs_N * tau_sgs_n1) * util::sqrt(Nn1) / (tau_eff_N * tau_eff_n1)
//              + util::pow(K * tau_sgs_n1 / tau_eff_n1,2.0) * util::sqrt(n1n1)
//           );
//  } else {
//    A = util::sqrt( util::pow((1-K) * (*tau_sgs_N) / tau_eff_N,2.0) * util::sqrt(NN));
//  }
//
//  A /= util::sqrt(N1N1);
//
//  *tau_sgs_N = tau_mol/((1.0/A)-1.0);
  *tau_sgs_N = (util::sqrt(2.0*N1N1)
                /(util::sqrt(2.0*n1n1)
                  -(2.0*util::pow(1./(this->getSmagoConst()*descriptors::invCs2<T,DESCRIPTOR>()),2.0)*tau_sgs_n1*tau_eff_n1)
                 )
               )*tau_eff_n1;

  if ((*tau_sgs_N - tau_mol) > T(0.)) {
    *tau_sgs_N -= tau_mol;
  }
  else {
    *tau_sgs_N = T(0.);
  }

}


}

#endif
