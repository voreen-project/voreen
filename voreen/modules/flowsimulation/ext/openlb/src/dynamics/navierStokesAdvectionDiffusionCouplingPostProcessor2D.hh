/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
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

#ifndef NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_HH
#define NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_HH

#include "latticeDescriptors.h"
#include "navierStokesAdvectionDiffusionCouplingPostProcessor2D.h"
#include "core/util.h"
#include "core/finiteDifference2D.h"


namespace olb {

//=====================================================================================
//==============  NavierStokesAdvectionDiffusionCouplingPostProcessor2D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<BlockStructureD<2>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), partners(partners_)
{
  this->getName() = "NavierStokesAdvectionDiffusionCouplingPostProcessor2D";
  // we normalize the direction of force vector
  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD]*dir[iD];
  }
  normDir = util::sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    forcePrefactor[iD] = gravity * dir[iD];
  }

  tPartner = static_cast<BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        auto cell = blockLattice.get(iX,iY);
        auto partnerCell = tPartner->get(iX,iY);
        // computation of the bousinessq force
        auto force = cell.template getFieldPointer<descriptors::FORCE>();
        T temperatureDifference = partnerCell.computeRho() - T0;
        for (unsigned iD = 0; iD < L::d; ++iD) {
          force[iD] = forcePrefactor[iD] * temperatureDifference;
        }
        // Velocity coupling
        auto u = partnerCell.template getField<descriptors::VELOCITY>();
        cell.computeU(u.data());
        partnerCell.template setField<descriptors::VELOCITY>(u);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionCouplingGenerator2D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionCouplingGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<2>* > partners) const
{
  return new NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionCouplingGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new NavierStokesAdvectionDiffusionCouplingGenerator2D<T,DESCRIPTOR>(*this);
}



//=====================================================================================
//==============  TotalEnthalpyPhaseChangeCouplingPostProcessor2D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
TotalEnthalpyPhaseChangeCouplingPostProcessor2D<T,DESCRIPTOR,DYNAMICS>::
TotalEnthalpyPhaseChangeCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<BlockStructureD<2>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), partners(partners_)
{
  this->getName() = "TotalEnthalpyPhaseChangeCouplingPostProcessor2D";
  // we normalize the direction of force vector
  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD]*dir[iD];
  }
  normDir = util::sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    forcePrefactor[iD] = gravity * dir[iD];
  }

  tPartner = static_cast<BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TEMPERATURE>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
void TotalEnthalpyPhaseChangeCouplingPostProcessor2D<T,DESCRIPTOR,DYNAMICS>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {
    auto* dynamics = static_cast<DYNAMICS*>(tPartner->template getDynamics<DYNAMICS>());
    typename DYNAMICS::ParametersD& parameters = static_cast<typename DYNAMICS::ParametersD&>(
      tPartner->template getData<DynamicsParameters<DYNAMICS>>());

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        auto cell = blockLattice.get(iX,iY);
        auto partnerCell = tPartner->get(iX,iY);

        T enthalpy = partnerCell.computeRho();

        cell.template setField<descriptors::POROSITY>(
          dynamics->template computeLiquidFraction<T>(parameters, enthalpy));
        auto temperature = partnerCell.template getFieldPointer<descriptors::TEMPERATURE>();
        temperature[0] = dynamics->template computeTemperature<T>(parameters, enthalpy);

        // computation of the bousinessq force
        auto force = cell.template getFieldPointer<descriptors::FORCE>();
        T temperatureDifference = temperature[0] - T0;
        for (unsigned iD = 0; iD < L::d; ++iD) {
          force[iD] = forcePrefactor[iD] * temperatureDifference;
        }
        // Velocity coupling
        T u[DESCRIPTOR::d] { };
        cell.computeU(u);
        partnerCell.template setField<descriptors::VELOCITY>(u);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
void TotalEnthalpyPhaseChangeCouplingPostProcessor2D<T,DESCRIPTOR,DYNAMICS>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
TotalEnthalpyPhaseChangeCouplingGenerator2D<T,DESCRIPTOR,DYNAMICS>::
TotalEnthalpyPhaseChangeCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
PostProcessor2D<T,DESCRIPTOR>* TotalEnthalpyPhaseChangeCouplingGenerator2D<T,DESCRIPTOR,DYNAMICS>::generate (
  std::vector<BlockStructureD<2>* > partners) const
{
  return new TotalEnthalpyPhaseChangeCouplingPostProcessor2D<T,DESCRIPTOR,DYNAMICS>(
           this->x0,this->x1,this->y0,this->y1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* TotalEnthalpyPhaseChangeCouplingGenerator2D<T,DESCRIPTOR,DYNAMICS>::clone() const
{
  return new TotalEnthalpyPhaseChangeCouplingGenerator2D<T,DESCRIPTOR,DYNAMICS>(*this);
}


//=====================================================================================
//==============  PhaseFieldCouplingPostProcessor2D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
PhaseFieldCouplingPostProcessor2D<T,DESCRIPTOR>::
PhaseFieldCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
                                  T rho_L, T rho_H, T mu_L, T mu_H, T surface_tension, T interface_thickness,
                                  std::vector<BlockStructureD<2>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
     _rho_L(rho_L), _rho_H(rho_H), _delta_rho(rho_H - rho_L), _mu_L(mu_L), _mu_H(mu_H), _surface_tension(surface_tension), _interface_thickness(interface_thickness),
     _beta(12.0 * surface_tension / interface_thickness), _kappa(1.5 * surface_tension * interface_thickness)
{
  this->getName() = "PhaseFieldCouplingPostProcessor2D";
  tPartner = static_cast<BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::INTERPHASE_NORMAL>> *>(partners_[0]);
}

template<typename T, typename DESCRIPTOR>
void PhaseFieldCouplingPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( x0, x1, y0, y1,
                         x0_, x1_, y0_, y1_,
                         newX0, newX1, newY0, newY1 ) ) {

    // generate phi cache
    auto& phi_cache = blockLattice.template getField<PHI_CACHE>()[0];
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        phi_cache[blockLattice.getCellId(iX,iY)] = util::max(util::min(tPartner->get(iX,iY).computeRho(), 1.0), 0.0);
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        auto cell = blockLattice.get(iX,iY);
        auto partnerCell = tPartner->get(iX,iY);

        T phi = phi_cache[blockLattice.getCellId(iX,iY)];

        // compute rho from phi
        T rho = _rho_L + phi * _delta_rho;

        // compute dynamic viscosity
        T viscosity = _mu_L + phi * (_mu_H - _mu_L);

        // get relaxation time
        T tau = cell.template getField<descriptors::TAU_EFF>();

        // compute grad phi and laplace phi
        Vector<T,L::d> grad_phi(0.0, 0.0);
        T laplace_phi = 0.0;
        for (int iPop = 1; iPop < L::q; ++iPop) {
          int nextX = iX + descriptors::c<L>(iPop,0);
          int nextY = iY + descriptors::c<L>(iPop,1);
          T neighbor_phi = phi_cache[blockLattice.getCellId(nextX,nextY)];

          laplace_phi += (neighbor_phi - phi) * descriptors::t<T,L>(iPop);

          neighbor_phi *= descriptors::t<T,L>(iPop);
          grad_phi += neighbor_phi * descriptors::c<L>(iPop);
        }
        grad_phi *= descriptors::invCs2<T,L>();
        laplace_phi *= 2.0 * descriptors::invCs2<T,L>();

        // compute grad rho
        Vector<T,L::d> grad_rho(_delta_rho, _delta_rho);
        grad_rho *= grad_phi;

        // compute interphase normal, save to external field
        T norm_grad_phi = norm(grad_phi);
        norm_grad_phi = util::max(norm_grad_phi, std::numeric_limits<T>::epsilon());
        partnerCell.template setField<descriptors::INTERPHASE_NORMAL>(grad_phi / norm_grad_phi);

        // compute forces (F_s, F_b, F_p, F_nu)
        // F_s (surface tension)
        T chemical_potential = (4.0 * _beta * (phi - 0.0) * (phi - 0.5) * (phi - 1.0)) - _kappa * laplace_phi;
        T surface_tension_force[] = {chemical_potential*grad_phi[0], chemical_potential*grad_phi[1]};

        // F_b (body force, e.g. bouyancy)
        T body_force[] = {0.0, 0.0};

        // F_p (pressure)
        T pressure = blockLattice.get(iX,iY).computeRho();
        T pressure_force[] = {-pressure / descriptors::invCs2<T,L>() * grad_rho[0], -pressure / descriptors::invCs2<T,L>() * grad_rho[1]};

        // F_nu (viscous)
        T viscous_force[] = {0.0, 0.0};
        T rho_tmp, u_tmp[2];
        cell.computeRhoU( rho_tmp, u_tmp );
        T p_tmp = rho_tmp / descriptors::invCs2<T,DESCRIPTOR>();
        T uSqr_tmp = util::normSqr<T,DESCRIPTOR::d>(u_tmp);
        for (int iPop = 0; iPop < L::q; ++iPop) {
          T fEq = cell.getDynamics()->computeEquilibrium( iPop, p_tmp, u_tmp);
          T fNeq = cell[iPop] - fEq;
          for (int iD = 0; iD < L::d; ++iD) {
            for (int jD = 0; jD < L::d; ++jD) {
              viscous_force[iD] += descriptors::c<L>(iPop,iD) * descriptors::c<L>(iPop,jD) * fNeq * grad_rho[jD];
            }
          }
        }
        for (int iD = 0; iD < L::d; ++iD) {
          viscous_force[iD] *= - viscosity / rho / tau * descriptors::invCs2<T,L>();
        }

        // save force/rho to external field
        auto force = cell.template getFieldPointer<descriptors::FORCE>();
        for (int iD = 0; iD < L::d; ++iD) {
          force[iD] = (surface_tension_force[iD] + body_force[iD] + pressure_force[iD] + viscous_force[iD]) / rho;
        }

        // compute u, save to external field
        Vector<T,DESCRIPTOR::template size<descriptors::VELOCITY>()> u;
        cell.computeU(u.data());
        partnerCell.template setField<descriptors::VELOCITY>(u);

        // compute relaxation time, save to external field

        tau = viscosity / rho * descriptors::invCs2<T,L>() + 0.5;
        cell.template setField<descriptors::TAU_EFF>(tau);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PhaseFieldCouplingPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
PhaseFieldCouplingGenerator2D<T,DESCRIPTOR>::
PhaseFieldCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
                              T rho_L, T rho_H, T mu_L, T mu_H, T surface_tension, T interface_thickness)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    _rho_L(rho_L), _rho_H(rho_H), _mu_L(mu_L), _mu_H(mu_H), _surface_tension(surface_tension), _interface_thickness(interface_thickness)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* PhaseFieldCouplingGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<2>* > partners) const
{
  return new PhaseFieldCouplingPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, _rho_L, _rho_H, _mu_L, _mu_H, _surface_tension, _interface_thickness, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* PhaseFieldCouplingGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new PhaseFieldCouplingGenerator2D<T,DESCRIPTOR>(*this);
}


//=====================================================================================
//==============  SmagorinskyBoussinesqCouplingPostProcessor2D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
SmagorinskyBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
SmagorinskyBoussinesqCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_, T smagoPrefactor_,
    std::vector<BlockStructureD<2>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), PrTurb(PrTurb_), smagoPrefactor(smagoPrefactor_), partners(partners_)
{
  this->getName() = "SmagorinskyBoussinesqCouplingPostProcessor2D";
  // we normalize the direction of force vector
  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD]*dir[iD];
  }
  normDir = util::sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    forcePrefactor[iD] = gravity * dir[iD];
  }

  tauTurbADPrefactor = descriptors::invCs2<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF>>() / descriptors::invCs2<T,DESCRIPTOR>() / PrTurb;
  tPartner = static_cast<BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void SmagorinskyBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {

        // computation of the bousinessq force
        auto force = blockLattice.get(iX,iY).template getFieldPointer<descriptors::FORCE>();
        T temperatureDifference = tPartner->get(iX,iY).computeRho() - T0;
        for (unsigned iD = 0; iD < L::d; ++iD) {
          force[iD] = forcePrefactor[iD] * temperatureDifference;
        }

        // Velocity coupling
        auto u = tPartner->get(iX,iY).template getField<descriptors::VELOCITY>();
        // tau coupling
        auto tauNS = blockLattice.get(iX,iY).template getFieldPointer<descriptors::TAU_EFF>();
        auto tauAD = tPartner->get(iX,iY).template getFieldPointer<descriptors::TAU_EFF>();

        T rho, pi[util::TensorVal<DESCRIPTOR >::n];
        blockLattice.get(iX,iY).computeAllMomenta(rho, u.data(), pi);
        tPartner->get(iX,iY).template setField<descriptors::VELOCITY>(u);
        T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
        if (util::TensorVal<DESCRIPTOR >::n == 6) {
          PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
        }
        T PiNeqNorm    = util::sqrt(PiNeqNormSqr);
        /// Molecular realaxation time
        T tau_mol_NS = 1. / blockLattice.get(iX,iY).getDynamics()->getParameters(blockLattice).template getOrFallback<descriptors::OMEGA>(0);
        T tau_mol_AD = 1. / tPartner->get(iX,iY).getDynamics()->getParameters(*tPartner).template getOrFallback<descriptors::OMEGA>(0);
        /// Turbulent realaxation time
        T tau_turb_NS = 0.5*(util::sqrt(tau_mol_NS*tau_mol_NS + smagoPrefactor/rho*PiNeqNorm) - tau_mol_NS);
        /// Effective realaxation time
        tauNS[0] = tau_mol_NS+tau_turb_NS;

        T tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
        tauAD[0] = tau_mol_AD+tau_turb_AD;
      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
void SmagorinskyBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
SmagorinskyBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::
SmagorinskyBoussinesqCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_, T smagoPrefactor_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_), PrTurb(PrTurb_), smagoPrefactor(smagoPrefactor_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* SmagorinskyBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<2>* > partners) const
{
  return new SmagorinskyBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, gravity, T0, deltaTemp, dir, PrTurb, smagoPrefactor, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* SmagorinskyBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new SmagorinskyBoussinesqCouplingGenerator2D<T,DESCRIPTOR>(*this);
}


//=====================================================================================
//==============  MixedScaleBoussinesqCouplingPostProcessor2D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
MixedScaleBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
MixedScaleBoussinesqCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
    std::vector<BlockStructureD<2>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), PrTurb(PrTurb_), partners(partners_)
{
  // we normalize the direction of force vector
  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD]*dir[iD];
  }
  normDir = util::sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    forcePrefactor[iD] = gravity * dir[iD];
  }

  tauTurbADPrefactor = descriptors::invCs2<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>>() / descriptors::invCs2<T,DESCRIPTOR>() / PrTurb;
  tPartner = static_cast<BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void MixedScaleBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  const T C_nu = 0.04;
  const T C_alpha = 0.5;
  const T deltaT = 1.0;

  const T invCs2_g = descriptors::invCs2<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>>();

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    auto& heatFluxCache = blockLattice.template getField<HEAT_FLUX_CACHE>()[0];
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        const T temperature = tPartner->get(iX,iY).computeRho();
        heatFluxCache[blockLattice.getCellId(iX, iY)] = temperature;

        // computation of the bousinessq force
        T temperatureDifference = temperature - T0;
        blockLattice.get(iX,iY).template setField<descriptors::FORCE>(temperatureDifference*forcePrefactor);

        FieldD<T,DESCRIPTOR,descriptors::VELOCITY> u;
        blockLattice.get(iX,iY).computeU(u.data());
        tPartner->get(iX,iY).template setField<descriptors::VELOCITY>(u);
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {

        auto u_pp = tPartner->get(iX+1, iY+1).template getField<descriptors::VELOCITY>();
        auto u_0p = tPartner->get(iX, iY+1).template getField<descriptors::VELOCITY>();
        auto u_np = tPartner->get(iX-1, iY+1).template getField<descriptors::VELOCITY>();
        auto u_p0 = tPartner->get(iX+1, iY  ).template getField<descriptors::VELOCITY>();
        auto u_pn = tPartner->get(iX+1, iY-1).template getField<descriptors::VELOCITY>();
        auto u_00 = tPartner->get(iX, iY  ).template getField<descriptors::VELOCITY>();
        auto u_0n = tPartner->get(iX, iY-1).template getField<descriptors::VELOCITY>();
        auto u_n0 = tPartner->get(iX-1, iY  ).template getField<descriptors::VELOCITY>();
        auto u_nn = tPartner->get(iX-1, iY-1).template getField<descriptors::VELOCITY>();

        const T *h_pp = & heatFluxCache[blockLattice.getCellId(iX+1, iY+1)];
        const T *h_0p = & heatFluxCache[blockLattice.getCellId(iX, iY+1)];
        const T *h_np = & heatFluxCache[blockLattice.getCellId(iX-1, iY+1)];
        const T *h_p0 = & heatFluxCache[blockLattice.getCellId(iX+1, iY  )];
        const T *h_pn = & heatFluxCache[blockLattice.getCellId(iX+1, iY-1)];
        const T *h_00 = & heatFluxCache[blockLattice.getCellId(iX, iY  )];
        const T *h_0n = & heatFluxCache[blockLattice.getCellId(iX, iY-1)];
        const T *h_n0 = & heatFluxCache[blockLattice.getCellId(iX-1, iY  )];
        const T *h_nn = & heatFluxCache[blockLattice.getCellId(iX-1, iY-1)];

        Vector<T, 2> filtered_u;
        T filtered_h;
        filtered_h =((h_pp[0] + 2.*h_0p[0] + h_np[0])
                     + 2.*(h_p0[0] + 2.*h_00[0] + h_n0[0])
                     + (h_pn[0] + 2.*h_0n[0] + h_nn[0]))*0.25*0.25;

        filtered_u =((u_pp + 2.*u_0p + u_np)
                     + 2.*(u_p0 + 2.*u_00 + u_n0)
                     + (u_pn + 2.*u_0n + u_nn))*0.25*0.25;

        Vector<T,2> filtered_u_reduced = u_00 - filtered_u;
        Vector<T,2> filtered_heatFlux_reduced = h_00[0]*u_00 - filtered_h*filtered_u;

        T cutoffKinEnergy = filtered_u_reduced[0]*filtered_u_reduced[0]
                            + filtered_u_reduced[1]*filtered_u_reduced[1];
        T cutoffHeatFlux = filtered_heatFlux_reduced[0]*filtered_heatFlux_reduced[0]
                           + filtered_heatFlux_reduced[1]*filtered_heatFlux_reduced[1];

        blockLattice.get(iX,iY).template setField<descriptors::CUTOFF_KIN_ENERGY>(util::pow(0.5*cutoffKinEnergy, 0.25));
        tPartner->get(iX,iY).template setField<descriptors::CUTOFF_HEAT_FLUX>(util::pow(0.5*cutoffHeatFlux, 0.25));
        // cout << cutoffKinEnergy << " " << u_00[0] << " " << u_00[1] << " " << cutoffKinEnergy_14[0] << std::endl;

        // tau coupling
        auto tauNS = blockLattice.get(iX,iY).template getField<descriptors::TAU_EFF>();
        auto tauAD = tPartner->get(iX,iY).template getField<descriptors::TAU_EFF>();

        /// Molecular realaxation time
        T tau_mol_NS = 1. / blockLattice.get(iX,iY).getDynamics()->getParameters(blockLattice).template getOrFallback<descriptors::OMEGA>(0);
        T tau_mol_AD = 1. / tPartner->get(iX,iY).getDynamics()->getParameters(*tPartner).template getOrFallback<descriptors::OMEGA>(0);

        const T temperature = tPartner->get(iX,iY).computeRho();

        // computation of the bousinessq force
        T temperatureDifference = temperature - T0;
        blockLattice.get(iX,iY).template
        setField<descriptors::FORCE>(temperatureDifference*forcePrefactor);

        auto u = tPartner->get(iX,iY).template getField<descriptors::VELOCITY>();
        T rho, pi[util::TensorVal<DESCRIPTOR>::n], j[DESCRIPTOR::d];
        // blockLattice.get(iX,iY).computeAllMomenta(rho, u, pi);
        rho = blockLattice.get(iX,iY).computeRho();
        blockLattice.get(iX,iY).computeStress(pi);

        auto force = blockLattice.get(iX,iY).template getField<descriptors::FORCE>();

        int iPi = 0;
        for (int Alpha=0; Alpha<DESCRIPTOR::d; ++Alpha) {
          for (int Beta=Alpha; Beta<DESCRIPTOR::d; ++Beta) {
            pi[iPi] += rho/2.*(force[Alpha]*u[Beta] + u[Alpha]*force[Beta]);
            ++iPi;
          }
        }
        const Vector<T,3> piSqr = {pi[0]*pi[0], pi[1]*pi[1], pi[2]*pi[2]};
        const T PiNeqNormSqr = piSqr[0] + 2.0*piSqr[1] + piSqr[2];
        const T PiNeqNorm    = util::sqrt(PiNeqNormSqr);

        tPartner->get(iX,iY).computeJ(j);
        const T tmp_preFactor = invCs2_g / rho / tauAD;
        const Vector<T,2> jNeq = {(j[0] - temperature * u[0]), (j[1] - temperature * u[1])};
        const Vector<T,2> jNeqSqr = {jNeq[0]*jNeq[0], jNeq[1]*jNeq[1]};
        const T jNeqSqr_prefacor = 2. * 0.25 * (jNeq[0] + jNeq[1]) * (jNeq[0] + jNeq[1]);

        const T TnormSqr = jNeqSqr_prefacor*PiNeqNormSqr;
        const T Tnorm    = util::sqrt(TnormSqr);

        /// Turbulent realaxation time
        // T tau_turb_NS = 0.5*(util::sqrt(tau_mol_NS*tau_mol_NS + dynamic_cast<SmagorinskyDynamics<T,DESCRIPTOR>*>(blockLattice.get(iX,iY).getDynamics())->getPreFactor()/rho*PiNeqNorm) - tau_mol_NS);

        // const T tmp_A = C_nu * util::sqrt(util::sqrt(2.)/2.) * descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>() * util::sqrt(PiNeqNorm / rho) * cutoffKinEnergy_14[0];
        // const T tmp_A_2 = tmp_A * tmp_A;
        // const T tmp_A_4 = tmp_A_2 * tmp_A_2;

        // const T tau_mol_NS_2 = tau_mol_NS * tau_mol_NS;
        // const T tau_mol_NS_3 = tau_mol_NS_2 * tau_mol_NS;

        // const T tmp_1_3 = 1./3.;
        // const T tmp_2_13 = util::pow(2., tmp_1_3);
        // const T tmp_3_3_12 = 3. * util::sqrt(3.);

        // const T tmp_sqrtA = util::sqrt(27.*tmp_A_4-4.*tmp_A_2*tau_mol_NS_3);

        //   // T tau_turb_NS = 1/3 ((27 A^2 + 3 util::sqrt(3) util::sqrt(27 A^4 - 4 A^2 b^3) - 2 b^3)^(1/3)/2^(1/3) + (2^(1/3) b^2)/(27 A^2 + 3 util::sqrt(3) util::sqrt(27 A^4 - 4 A^2 b^3) - 2 b^3)^(1/3) - b)
        // T tau_turb_NS = ( util::pow(27.*tmp_A_2 + tmp_3_3_12*util::sqrt(27.*tmp_A_4-4.*tmp_A_2*tau_mol_NS_3)-2.*tau_mol_NS_3, tmp_1_3) / tmp_2_13
        //                    + (tmp_2_13*tau_mol_NS_2) / util::pow(27.*tmp_A_2+tmp_3_3_12*util::sqrt(27.*tmp_A_4-4.*tmp_A_2*tau_mol_NS_3) - 2.*tau_mol_NS_3, tmp_1_3)
        //                    - tau_mol_NS
        //                    ) * tmp_1_3;

        // if ( tau_turb_NS != tau_turb_NS )
        //   tau_turb_NS = 0.;

        //cout << tau_turb_NS << " " << 27. * tmp_A_2 << " " << 4. * tau_mol_NS_3 << " " << PiNeqNorm << " " << " " << rho << std::endl;

        auto cutoffKinEnergy_14 = blockLattice.get(iX,iY).template getField<descriptors::CUTOFF_KIN_ENERGY>();
        auto cutoffHeatFlux_14 = tPartner->get(iX,iY).template getField<descriptors::CUTOFF_HEAT_FLUX>();

        const T tmp_A = C_nu * util::sqrt(util::sqrt(2.)/2.) * descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>() * util::sqrt(PiNeqNorm / rho / tauNS) * cutoffKinEnergy_14;
        const T tau_turb_NS = tmp_A;

        // T tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
        const T tmp_B = C_alpha * descriptors::invCs2<T,DESCRIPTOR>() / rho * util::sqrt(2.0 * Tnorm * invCs2_g / tauNS / tauAD) * cutoffHeatFlux_14;
        const T tau_turb_AD = tmp_B;
        // cout << jNeq[0] << " " << jNeq[1] << " " << util::sqrt(Tnorm * invCs2_g / tauNS / tauAD) << " " << TnormSqr << std::endl;

        /// Effective realaxation time
        blockLattice.get(iX,iY).template setField<descriptors::TAU_EFF>(tau_mol_NS+tau_turb_NS);
        tPartner->get(iX,iY).template setField<descriptors::TAU_EFF>(tau_mol_AD+tau_turb_AD);

      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
void MixedScaleBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
MixedScaleBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::
MixedScaleBoussinesqCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
                                        T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_), PrTurb(PrTurb_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* MixedScaleBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<2>* > partners) const
{
  return new MixedScaleBoussinesqCouplingPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, gravity, T0, deltaTemp, dir, PrTurb, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* MixedScaleBoussinesqCouplingGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new MixedScaleBoussinesqCouplingGenerator2D<T,DESCRIPTOR>(*this);
}

}  // namespace olb

#endif
