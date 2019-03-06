/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcekt
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
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 * Power Law similar to "Huidan Yu, Sharath S. Girimaji, Li-Shi Luo - DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method"
 */
#ifndef POWER_LAW_BGK_DYNAMICS_HH
#define POWER_LAW_BGK_DYNAMICS_HH

#include "dynOmegaLatticeDescriptors.h"
#include "powerLawBGKdynamics.h"
#include "core/cell.h"
#include "core/util.h"
#include "lbHelpers.h"
#include "latticeDescriptors.h"
#include "math.h"

namespace olb {

////////////////////// Class PowerLawBGKdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
PowerLawBGKdynamics<T,Lattice>::PowerLawBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_, T m_, T n_ , T dt_)
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    m(m_),
    n(n_),
    dt(dt_)
    //preFactor(computePreFactor(omega_,smagoConst_) )
{ }

template<typename T, template<typename U> class Lattice>
void PowerLawBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // load old omega from dyn. omega descriptor
  //T oldOmega = this->getOmega(); //compute with constant omega
  T oldOmega = cell.getExternal(Lattice<T>::ExternalField::omegaBeginsAt)[0]; //compute with dynamic omega
  T newOmega = computeOmega(
                 oldOmega, preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  // save new omega to dyn. omega descriptor
  cell.getExternal(Lattice<T>::ExternalField::omegaBeginsAt)[0] = newOmega; //compute with dynamic omega
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void PowerLawBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void PowerLawBGKdynamics<T,Lattice>::setOmega(T omega)
{
  this->BGKdynamics<T,Lattice>::setOmega(omega);
  preFactor = computePreFactor(omega, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T PowerLawBGKdynamics<T,Lattice>::getPowerLawOmega(Cell<T,Lattice>& cell )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T PowerLawBGKdynamics<T,Lattice>::computePreFactor(T omega, T smagoConst)
{
  return (T)0.5 * util::sqr(smagoConst*omega*Lattice<T>::invCs2);
}

template<typename T, template<typename U> class Lattice>
T PowerLawBGKdynamics<T,Lattice>::computeOmega(T omega0, T preFactor, T rho, T pi[util::TensorVal<Lattice<T> >::n] )
{

  // strain rate tensor without prefactor
  T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
  }

  T pre2 = pow(Lattice<T>::invCs2/2./dt* omega0/rho,2.); // prefactor to the strain rate tensor
  T D = pre2*PiNeqNormSqr; // Strain rate tensor
  T gamma = sqrt(2.*D); // shear rate

  T nuNew = m*pow(gamma,n-1.); //nu for non-Newtonian fluid
  T newOmega = 2./(nuNew*6.+1.);

  /*
     * problem if newOmega too small or too big is see e.g. "Daniel Conrad , Andreas Schneider, Martin Böhle:
     * A viscosity adaption method for Lattice Boltzmann simulations"
    */
  if (newOmega>1.965) {
    newOmega = 1.965;  //std::cout << newOmega << std::endl;
  }
  if (newOmega<0.1) {
    newOmega = 0.1;  //std::cout << newOmega << std::endl;
  }
  return newOmega;
  //return omega0;
}


////////////////////// Class ForcedPowerLawBGKdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
PowerLawForcedBGKdynamics<T,Lattice>::PowerLawForcedBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_, T m_, T n_, T dt_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    m(m_),
    n(n_),
    dt(dt_)
    //preFactor(computePreFactor(omega_,smagoConst_) )
{ }

template<typename T, template<typename U> class Lattice>
void PowerLawForcedBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // load old omega from dyn. omega descriptor
  //T oldOmega = this->getOmega(); //compute with constant omega
  T oldOmega = cell.getExternal(Lattice<T>::ExternalField::omegaBeginsAt)[0]; //compute with dynamic omega
  T newOmega = computeOmega(oldOmega, preFactor, rho, pi);
  T* force = cell.getExternal(Lattice<T>::ExternalField::forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  // save new omega to dyn. omega descriptor
  cell.getExternal(Lattice<T>::ExternalField::omegaBeginsAt)[0] = newOmega; //compute with dynamic omega
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void PowerLawForcedBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void PowerLawForcedBGKdynamics<T,Lattice>::setOmega(T omega)
{
  this->setOmega(omega);
  preFactor = computePreFactor(omega, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T PowerLawForcedBGKdynamics<T,Lattice>::getPowerLawOmega(Cell<T,Lattice>& cell )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T PowerLawForcedBGKdynamics<T,Lattice>::computePreFactor(T omega, T smagoConst)
{
  return (T)0.5 * util::sqr(smagoConst*omega*Lattice<T>::invCs2);
}

template<typename T, template<typename U> class Lattice>
T PowerLawForcedBGKdynamics<T,Lattice>::computeOmega(T omega0, T preFactor, T rho, T pi[util::TensorVal<Lattice<T> >::n] )
{

  /* T PiNeqNormSqr = (pi[0]+1./Lattice<T>::invCs2*(rho-(T)1))*(pi[0]+1./Lattice<T>::invCs2*(rho-(T)1)) + 2.0*pi[1]*pi[1] + (pi[2]+1./Lattice<T>::invCs2*(rho-(T)1))*(pi[2]+1./Lattice<T>::invCs2*(rho-(T)1));
   if (util::TensorVal<Lattice<T> >::n == 6)
     PiNeqNormSqr += + 2*pi[2]*pi[2] - (pi[2]+1./Lattice<T>::invCs2*(rho-(T)1))*(pi[2]+1./Lattice<T>::invCs2*(rho-(T)1)) + (pi[3]+1./Lattice<T>::invCs2*(rho-(T)1))*(pi[3]+1./Lattice<T>::invCs2*(rho-(T)1)) + 2*pi[4]*pi[4] + (pi[5]+1./Lattice<T>::invCs2*(rho-(T)1))*(pi[5]+1./Lattice<T>::invCs2*(rho-(T)1));
   T PiNeqNorm    = sqrt(PiNeqNormSqr)/rho; // TODO "*rho" or "/rho"?

   // Compute new omega0 by solving with newton's scheme
   T tau = 1./omega0;
   for (int i=0; i<10 ; i++) {
     T fTau = m*pow(3,n)*pow(tau,1-n)*pow(PiNeqNorm,n-1)+0.5-tau;
     T dfTau = m*pow(3,n)*(1.-n)*pow(tau,-n)*pow(PiNeqNorm,n-1)-1;
     tau=tau-fTau/dfTau;
     //std::cout << "Newton step=" << i << "; tau=" << tau << std::endl;
   }
   omega0=1./tau;
   return omega0;
  */

  // strain rate tensor without prefactor
  T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
  }

  T pre2 = pow(Lattice<T>::invCs2/2./dt* omega0/rho,2.); // prefactor to the strain rate tensor
  T D = pre2*PiNeqNormSqr; // Strain rate tensor
  T gamma = sqrt(2.*D); // shear rate

  T nuNew = m*pow(gamma,n-1.); //nu for non-Newtonian fluid
  T newOmega = 2./(nuNew*6.+1.);
  /*
   * problem if newOmega too small or too big is see e.g. "Daniel Conrad , Andreas Schneider, Martin Böhle:
   * A viscosity adaption method for Lattice Boltzmann simulations"
  */
  if (newOmega>1.965) {
    newOmega = 1.965;  //std::cout << newOmega << std::endl;
  }
  if (newOmega<0.5) {
    newOmega = 0.5;  //std::cout << newOmega << std::endl;
  }
  return newOmega;
  //return omega0;
}


}

#endif
