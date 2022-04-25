/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software you can redistribute it and/or
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

#ifndef NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_HH
#define NAVIER_STOKES_INTO_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_3D_HH

#include "latticeDescriptors.h"
#include "navierStokesAdvectionDiffusionCouplingPostProcessor3D.h"
#include "core/util.h"
#include "core/finiteDifference3D.h"
#include "advectionDiffusionForces.hh"
#include "advectionDiffusionForces.h"

#include "descriptorAlias.h"

namespace olb {

//=====================================================================================
//==============  NavierStokesAdvectionDiffusionCouplingPostProcessor3D ===========
//=====================================================================================

template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<BlockStructureD<3>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), partners(partners_)
{
  this->getName() = "NavierStokesAdvectionDiffusionCouplingPostProcessor3D";
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

  tPartner = static_cast<BlockLattice<T,descriptors::D3Q7<descriptors::VELOCITY>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          // computation of the bousinessq force
          auto force = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::FORCE>();
          T temperatureDifference = tPartner->get(iX,iY,iZ).computeRho() - T0;
          for (unsigned iD = 0; iD < L::d; ++iD) {
            force[iD] = forcePrefactor[iD] * temperatureDifference;
          }
          // Velocity coupling
          auto u = tPartner->get(iX,iY,iZ).template getField<descriptors::VELOCITY>();
          blockLattice.get(iX,iY,iZ).computeU(u.data());
          tPartner->get(iX,iY,iZ).template setField<descriptors::VELOCITY>(u);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void NavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}


//=====================================================================================
//==============  TotalEnthalpyPhaseChangeCouplingPostProcessor3D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
TotalEnthalpyPhaseChangeCouplingPostProcessor3D<T,DESCRIPTOR,DYNAMICS>::
TotalEnthalpyPhaseChangeCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<BlockStructureD<3>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), partners(partners_)
{
  this->getName() = "TotalEnthalpyPhaseChangeCouplingPostProcessor3D";
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

  tPartner = static_cast<BlockLattice<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::TEMPERATURE>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
void TotalEnthalpyPhaseChangeCouplingPostProcessor3D<T,DESCRIPTOR,DYNAMICS>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    auto* dynamics = static_cast<DYNAMICS*>(tPartner->template getDynamics<DYNAMICS>());
    typename DYNAMICS::ParametersD& parameters = static_cast<typename DYNAMICS::ParametersD&>(
      tPartner->template getData<DynamicsParameters<DYNAMICS>>());

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          auto cell = blockLattice.get(iX,iY,iZ);
          auto partnerCell = tPartner->get(iX,iY,iZ);

          T enthalpy = partnerCell.computeRho();

          cell.template setField<descriptors::POROSITY>(
            dynamics->template computeLiquidFraction<T>(parameters, enthalpy));
          auto temperature = partnerCell.template getFieldPointer<descriptors::TEMPERATURE>();
          temperature[0] = dynamics->template computeTemperature<T>(parameters, enthalpy);

          // computation of the bousinessq force
          auto force = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::FORCE>();
          T temperatureDifference = temperature[0] - T0;
          for (unsigned iD = 0; iD < L::d; ++iD) {
            force[iD] = forcePrefactor[iD] * temperatureDifference;
          }
          // Velocity coupling
          FieldD<T,DESCRIPTOR,descriptors::VELOCITY> u;
          blockLattice.get(iX,iY,iZ).computeU(u.data());
          tPartner->get(iX,iY,iZ).template setField<descriptors::VELOCITY>(u);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
void TotalEnthalpyPhaseChangeCouplingPostProcessor3D<T,DESCRIPTOR,DYNAMICS>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
TotalEnthalpyPhaseChangeCouplingGenerator3D<T,DESCRIPTOR,DYNAMICS>::
TotalEnthalpyPhaseChangeCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
PostProcessor3D<T,DESCRIPTOR>* TotalEnthalpyPhaseChangeCouplingGenerator3D<T,DESCRIPTOR,DYNAMICS>::generate (
  std::vector<BlockStructureD<3>* > partners) const
{
  return new TotalEnthalpyPhaseChangeCouplingPostProcessor3D<T,DESCRIPTOR,DYNAMICS>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* TotalEnthalpyPhaseChangeCouplingGenerator3D<T,DESCRIPTOR,DYNAMICS>::clone() const
{
  return new TotalEnthalpyPhaseChangeCouplingGenerator3D<T,DESCRIPTOR,DYNAMICS>(*this);
}


//=====================================================================================
//==============  PhaseFieldCouplingPostProcessor3D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
PhaseFieldCouplingPostProcessor3D<T,DESCRIPTOR>::
PhaseFieldCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                  T rho_L, T rho_H, T mu_L, T mu_H, T surface_tension, T interface_thickness,
                                  std::vector<BlockStructureD<3>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
     _rho_L(rho_L), _rho_H(rho_H), _delta_rho(rho_H - rho_L), _mu_L(mu_L), _mu_H(mu_H), _surface_tension(surface_tension), _interface_thickness(interface_thickness),
     _beta(12.0 * surface_tension / interface_thickness), _kappa(1.5 * surface_tension * interface_thickness)
{
  this->getName() = "PhaseFieldCouplingPostProcessor3D";
  tPartner = static_cast<BlockLattice<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::INTERPHASE_NORMAL>> *>(partners_[0]);
}

template<typename T, typename DESCRIPTOR>
void PhaseFieldCouplingPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    // generate phi cache
    auto& phi_cache = blockLattice.template getField<PHI_CACHE>()[0];
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        for (int iZ=newZ0-1; iZ<=newZ1+1; ++iZ) {
          phi_cache[blockLattice.getCellId(iX,iY,iZ)] = util::max(util::min(tPartner->get(iX,iY,iZ).computeRho(), 1.0), 0.0);
        }
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          T phi = phi_cache[blockLattice.getCellId(iX, iY, iZ)];

          // compute rho from phi
          T rho = _rho_L + phi * _delta_rho;

          // compute dynamic viscosity
          T viscosity = _mu_L + phi * (_mu_H - _mu_L);

          // get relaxation time
          auto tau = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::TAU_EFF>();

          // compute grad phi and laplace phi
          Vector<T,3> grad_phi;
          T laplace_phi = 0.0;
          for (int iPop = 1; iPop < L::q; ++iPop) {
            int nextX = iX + descriptors::c<L>(iPop,0);
            int nextY = iY + descriptors::c<L>(iPop,1);
            int nextZ = iZ + descriptors::c<L>(iPop,2);
            T neighbor_phi = phi_cache[blockLattice.getCellId(nextX, nextY, nextZ)];

            laplace_phi += (neighbor_phi - phi) * descriptors::t<T,L>(iPop);
            neighbor_phi *= descriptors::t<T,L>(iPop);
            grad_phi += neighbor_phi * descriptors::c<L>(iPop);
          }
          grad_phi *= descriptors::invCs2<T,L>();
          laplace_phi *= 2.0 * descriptors::invCs2<T,L>();

          // compute grad rho
          Vector<T,3> grad_rho = _delta_rho * grad_phi;

          // compute interphase normal, save to external field
          T norm_grad_phi = util::max(norm(grad_phi), std::numeric_limits<T>::epsilon());
          tPartner->get(iX,iY,iZ).template setField<descriptors::INTERPHASE_NORMAL>(grad_phi / norm_grad_phi);

          // compute _forces (F_s, F_b, F_p, F_nu)
          // F_s (surface tension)
          T chemical_potential = (4.0 * _beta * (phi - 0.0) * (phi - 0.5) * (phi - 1.0)) - _kappa * laplace_phi;
          Vector<T,3> surface_tension_force = chemical_potential * grad_phi;

          // F_b (body force, e.g. bouyancy)
          Vector<T,3> body_force;

          // F_p (pressure)
          T pressure = blockLattice.get(iX,iY,iZ).computeRho();
          Vector<T,3> pressure_force = -pressure / descriptors::invCs2<T,L>() * grad_rho;

          // F_nu (viscous)
          Vector<T,3> viscous_force;
          T rho_tmp, u_tmp[3];
          blockLattice.get(iX,iY,iZ).computeRhoU( rho_tmp, u_tmp );
          T p_tmp = rho_tmp / descriptors::invCs2<T,DESCRIPTOR>();
          T uSqr_tmp = util::normSqr<T,DESCRIPTOR::d>(u_tmp);
          for (int iPop = 0; iPop < L::q; ++iPop) {
            T fEq = blockLattice.get(iX,iY,iZ).getDynamics()->computeEquilibrium( iPop, p_tmp, u_tmp );
            T fNeq = blockLattice.get(iX,iY,iZ)[iPop] - fEq;
            for (int iD = 0; iD < L::d; ++iD) {
              for (int jD = 0; jD < L::d; ++jD) {
                viscous_force[iD] += descriptors::c<L>(iPop,iD) * descriptors::c<L>(iPop,jD) * fNeq * grad_rho[jD];
              }
            }
          }
          viscous_force *= -viscosity / rho / tau[0] * descriptors::invCs2<T,L>();

          // save force/rho to external field
          blockLattice.get(iX,iY,iZ).template setField<descriptors::FORCE>(
            (surface_tension_force + body_force + pressure_force + viscous_force) / rho
          );

          // compute u, save to external field
          auto u = tPartner->get(iX,iY,iZ).template getField<descriptors::VELOCITY>();
          blockLattice.get(iX,iY,iZ).computeU(u.data());
          tPartner->get(iX,iY,iZ).template setField<descriptors::VELOCITY>(u);

          // compute relaxation time, save to external field
          tau[0] = viscosity / rho * descriptors::invCs2<T,L>() + 0.5;
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PhaseFieldCouplingPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
PhaseFieldCouplingGenerator3D<T,DESCRIPTOR>::
PhaseFieldCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                              T rho_L, T rho_H, T mu_L, T mu_H, T surface_tension, T interface_thickness)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    _rho_L(rho_L), _rho_H(rho_H), _mu_L(mu_L), _mu_H(mu_H), _surface_tension(surface_tension), _interface_thickness(interface_thickness)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* PhaseFieldCouplingGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<3>* > partners) const
{
  return new PhaseFieldCouplingPostProcessor3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, _rho_L, _rho_H, _mu_L, _mu_H, _surface_tension, _interface_thickness, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* PhaseFieldCouplingGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new PhaseFieldCouplingGenerator3D<T,DESCRIPTOR>(*this);
}


//=====================================================================================
//==============  SmagorinskyBoussinesqCouplingPostProcessor3D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
SmagorinskyBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>::
SmagorinskyBoussinesqCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_, T smagoPrefactor_,
    std::vector<BlockStructureD<3>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), PrTurb(PrTurb_), smagoPrefactor(smagoPrefactor_), partners(partners_)
{
  this->getName() = "SmagorinskyBoussinesqCouplingPostProcessor3D";
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

  tauTurbADPrefactor = descriptors::invCs2<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::TAU_EFF>>() / descriptors::invCs2<T,DESCRIPTOR>() / PrTurb;
  tPartner = static_cast<BlockLattice<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::TAU_EFF>> *>(partners[0]);
}


template<typename T, typename DESCRIPTOR>
void SmagorinskyBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {

          // computation of the bousinessq force
          auto force = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::FORCE>();
          T temperatureDifference = tPartner->get(iX,iY,iZ).computeRho() - T0;
          for (unsigned iD = 0; iD < L::d; ++iD) {
            force[iD] = forcePrefactor[iD] * temperatureDifference;
          }

          // Velocity coupling
          auto u = tPartner->get(iX,iY,iZ).template getField<descriptors::VELOCITY>();

          // tau coupling
          auto tauNS = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::TAU_EFF>();
          auto tauAD = tPartner->get(iX,iY,iZ).template getFieldPointer<descriptors::TAU_EFF>();

          T rho, pi[util::TensorVal<DESCRIPTOR >::n];
          blockLattice.get(iX,iY,iZ).computeAllMomenta(rho, u.data(), pi);
          tPartner->get(iX,iY,iZ).template setField<descriptors::VELOCITY>(u);
          T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
          if (util::TensorVal<DESCRIPTOR >::n == 6) {
            PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
          }
          T PiNeqNorm    = util::sqrt(PiNeqNormSqr);
          /// Molecular realaxation time
          T tau_mol_NS = 1. / blockLattice.get(iX,iY,iZ).getDynamics()->getParameters(blockLattice).template getOrFallback<descriptors::OMEGA>(0);
          T tau_mol_AD = 1. / tPartner->get(iX,iY,iZ).getDynamics()->getParameters(*tPartner).template getOrFallback<descriptors::OMEGA>(0);
          /// Turbulent realaxation time
          T tau_turb_NS = 0.5*(util::sqrt(tau_mol_NS*tau_mol_NS + smagoPrefactor/rho*PiNeqNorm) - tau_mol_NS);
          if (tau_turb_NS != tau_turb_NS) {
            tau_turb_NS = 0.0;
          }
          /// Effective realaxation time
          tauNS[0] = tau_mol_NS+tau_turb_NS;

          T tau_turb_AD = tau_turb_NS * tauTurbADPrefactor;
          tauAD[0] = tau_mol_AD+tau_turb_AD;
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SmagorinskyBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
SmagorinskyBoussinesqCouplingGenerator3D<T,DESCRIPTOR>::
SmagorinskyBoussinesqCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_, T smagoPrefactor_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_), PrTurb(PrTurb_), smagoPrefactor(smagoPrefactor_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* SmagorinskyBoussinesqCouplingGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<3>* > partners) const
{
  return new SmagorinskyBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, gravity, T0, deltaTemp, dir, PrTurb, smagoPrefactor, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* SmagorinskyBoussinesqCouplingGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new SmagorinskyBoussinesqCouplingGenerator3D<T,DESCRIPTOR>(*this);
}



//=====================================================================================
//==============  AdvectionDiffusionParticleCouplingPostProcessor3D ===========
//=====================================================================================

template<typename T, typename DESCRIPTOR, typename ADLattice, typename FIELD_A, typename FIELD_B>
AdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice,FIELD_A,FIELD_B>::
AdvectionDiffusionParticleCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int iC_,
    std::vector<BlockStructureD<3>* > partners_,
    std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice>>> forces_)
  :  _forces(forces_),
     x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), iC(iC_),
     _partnerLattice(static_cast<BlockLattice<T,ADLattice> *>(partners_[0])),
     _cell(_partnerLattice->get(x0,y0,z0)),
     _cellXp(_partnerLattice->get(x0+1,y0,z0)),
     _cellXn(_partnerLattice->get(x0-1,y0,z0)),
     _cellYp(_partnerLattice->get(x0,y0+1,z0)),
     _cellYn(_partnerLattice->get(x0,y0-1,z0)),
     _cellZp(_partnerLattice->get(x0,y0,z0+1)),
     _cellZn(_partnerLattice->get(x0,y0,z0-1))
{
  this->getName() = "AdvectionDiffusionParticleCouplingPostProcessor3D";
}

template<typename T, typename DESCRIPTOR, typename ADLattice, typename FIELD_A, typename FIELD_B>
void AdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice,FIELD_A,FIELD_B>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  auto vel     = par ? _cell.template getField<FIELD_B>() : _cell.template getField<FIELD_A>();
  //auto vel_new = par ? _cell.template getFieldPointer<FIELD_A>() : _cell.template getFieldPointer<FIELD_B>();

  auto velXp = par ? _cellXp.template getField<FIELD_B>() : _cellXp.template getField<FIELD_A>();
  auto velXn = par ? _cellXn.template getField<FIELD_B>() : _cellXn.template getField<FIELD_A>();
  auto velYp = par ? _cellYp.template getField<FIELD_B>() : _cellYp.template getField<FIELD_A>();
  auto velYn = par ? _cellYn.template getField<FIELD_B>() : _cellYn.template getField<FIELD_A>();
  auto velZp = par ? _cellZp.template getField<FIELD_B>() : _cellZp.template getField<FIELD_A>();
  auto velZn = par ? _cellZn.template getField<FIELD_B>() : _cellZn.template getField<FIELD_A>();

  int newX0, newX1, newY0, newY1, newZ0, newZ1;

  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          int latticeR[4] = {iC, iX, iY, iZ};
          T velGrad[3] = {0.,0.,0.};
          T forceValue[3] = {0.,0.,0.};
          T velF[3] = {0.,0.,0.};

          auto nsCell = blockLattice.get(iX,iY,iZ);

          if (_forces.begin() != _forces.end()) {
            // calculating upwind Gradient
            // vel contains velocity information on ADlattice
            // velGrad contains upwind vel on ADlattice
            if (vel[0]<0.) {
              velGrad[0] = vel[0]*(velXp[0]-vel[0]);
              velGrad[1] = vel[0]*(velXp[1]-vel[1]);
              velGrad[2] = vel[0]*(velXp[2]-vel[2]);
            }
            else {
              velGrad[0] = vel[0]*(vel[0]-velXn[0]);
              velGrad[1] = vel[0]*(vel[1]-velXn[1]);
              velGrad[2] = vel[0]*(vel[2]-velXn[2]);
            }
            if (vel[1]<0.) {
              velGrad[0] += vel[1]*(velYp[0]-vel[0]);
              velGrad[1] += vel[1]*(velYp[1]-vel[1]);
              velGrad[2] += vel[1]*(velYp[2]-vel[2]);
            }
            else {
              velGrad[0] += vel[1]*(vel[0]-velYn[0]);
              velGrad[1] += vel[1]*(vel[1]-velYn[1]);
              velGrad[2] += vel[1]*(vel[2]-velYn[2]);
            }
            if (vel[2]<0.) {
              velGrad[0] += vel[2]*(velZp[0]-vel[0]);
              velGrad[1] += vel[2]*(velZp[1]-vel[1]);
              velGrad[2] += vel[2]*(velZp[2]-vel[2]);
            }
            else {
              velGrad[0] += vel[2]*(vel[0]-velZn[0]);
              velGrad[1] += vel[2]*(vel[1]-velZn[1]);
              velGrad[2] += vel[2]*(vel[2]-velZn[2]);
            }

            for (AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice>& f : _forces) {
              // writes force in forceValues, vel refers to ADlattice
              auto adCell = _partnerLattice->get(x0,y0,z0);
              f.applyForce(forceValue, &nsCell, &adCell, vel.data(), latticeR);
              if (par) {
                _cell.template setField<FIELD_B>(vel);
              }
              else {
                _cell.template setField<FIELD_A>(vel);
              }
            }

            // compute new particle velocity
            Vector<T,DESCRIPTOR::d> newVel;
            for (int i=0; i < DESCRIPTOR::d; i++) {
              newVel[i] = vel[i] + forceValue[i] - velGrad[i];
            }
            if (par) {
              _cell.template setField<FIELD_A>(newVel);
            } else {
              _cell.template setField<FIELD_B>(newVel);
            }
          }
          else {   // set particle velocity to fluid velocity
            nsCell.computeU(velF);
            Vector<T,DESCRIPTOR::d> newVel;
            for (int i = 0; i < DESCRIPTOR::d; i++) {
              newVel[i] = velF[i];
            }
            if (par) {
              _cell.template setField<FIELD_A>(newVel);
            } else {
              _cell.template setField<FIELD_B>(newVel);
            }
          }
        }
      }
    }
  }
  par = !par;
}

template<typename T, typename DESCRIPTOR, typename ADLattice, typename FIELD_A, typename FIELD_B>
void AdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice,FIELD_A,FIELD_B>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
NavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::
NavierStokesAdvectionDiffusionCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_,int z0_, int z1_, T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<3>* > partners) const
{
  return new NavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* NavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new NavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>(*this);
}

// LatticeCouplingGenerator for one-way advectionDiffusion coupling with Stokes drag

template<typename T, typename DESCRIPTOR, typename ADLattice, typename FIELD_A, typename FIELD_B>
AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice,FIELD_A,FIELD_B>::
AdvectionDiffusionParticleCouplingGenerator3D()
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0)
{ }

template<typename T, typename DESCRIPTOR, typename ADLattice, typename FIELD_A, typename FIELD_B>
PostProcessor3D<T,DESCRIPTOR>* AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice,FIELD_A,FIELD_B>::generate(
  std::vector<BlockStructureD<3>* > partners) const
{
  return new AdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,ADLattice,FIELD_A,FIELD_B>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, this->iC, partners, ADforces);
}

template<typename T, typename DESCRIPTOR, typename ADLattice, typename FIELD_A, typename FIELD_B>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice,FIELD_A,FIELD_B>::clone() const
{
  return new AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice,FIELD_A,FIELD_B>(*this);
}

template<typename T, typename DESCRIPTOR, typename ADLattice, typename FIELD_A, typename FIELD_B>
void AdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,ADLattice,FIELD_A,FIELD_B>::addForce(
  AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> &force)
{
  ADforces.push_back(force);
}


//=====================================================================================
//==============  PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D ===========
//=====================================================================================

template<typename T, typename DESCRIPTOR>
PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<BlockStructureD<3>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), partners(partners_)
{
  this->getName() = "PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D";
  // we normalize the direction of force vector
  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD]*dir[iD];
  }
  normDir = util::sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }
}

template<typename T, typename DESCRIPTOR>
void PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  typedef DESCRIPTOR L;
  enum {x,y,z};

  BlockLattice<T,descriptors::PorousAdvectionDiffusionD3Q7Descriptor> *tPartner =
    static_cast<BlockLattice<T,descriptors::PorousAdvectionDiffusionD3Q7Descriptor> *>(partners[0]);

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          //                  Velocity coupling
          auto u = tPartner->get(iX,iY,iZ).template getField<descriptors::VELOCITY>();
          blockLattice.get(iX,iY,iZ).computeU(u.data());
          tPartner->get(iX,iY,iZ).template setField<descriptors::VELOCITY>(u);

          //coupling between the temperature and navier stokes.

          auto force = blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::FORCE>();
          // this should return the interpolated solid-fluid temperature
          T temperature = tPartner->get(iX,iY,iZ).computeRho();
          T rho = blockLattice.get(iX,iY,iZ).computeRho();
          for (unsigned iD = 0; iD < L::d; ++iD) {
            force[iD] = gravity * rho * (temperature - T0) / deltaTemp * dir[iD];
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR>
PorousNavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::
PorousNavierStokesAdvectionDiffusionCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_,int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* PorousNavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<3>* > partners) const
{
  return new PorousNavierStokesAdvectionDiffusionCouplingPostProcessor3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* PorousNavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new PorousNavierStokesAdvectionDiffusionCouplingGenerator3D<T,DESCRIPTOR>(*this);
}


//=====================================================================================
//==============  MixedScaleBoussinesqCouplingPostProcessor3D ===============
//=====================================================================================

template<typename T, typename DESCRIPTOR>
MixedScaleBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>::
MixedScaleBoussinesqCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
    std::vector<BlockStructureD<3>* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
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

  tauTurbADPrefactor = descriptors::invCs2<T,descriptors::D3Q7<>>() / descriptors::invCs2<T,DESCRIPTOR>() / PrTurb;
  tPartner = static_cast<BlockLattice<T,descriptors::D3Q7<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>> *>(partners[0]);
}

template<typename T, typename DESCRIPTOR>
void MixedScaleBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{

  const T C_nu = 0.04;
  const T C_alpha = 0.5;
  // const T deltaT = 1.0;

  const T invCs2 = descriptors::invCs2<T,DESCRIPTOR>();
  const T invCs2_g = descriptors::invCs2<T,descriptors::D3Q7<>>();

  int newX0, newX1, newY0, newY1, newZ1, newZ0;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ1, newZ0 ) ) {

    auto& heatFluxCache = blockLattice.template getField<HEAT_FLUX_CACHE>()[0];
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        for (int iZ=newZ0-1; iZ<=newZ1+1; ++iZ) {
          const T temperature = tPartner->get(iX,iY,iZ).computeRho();
          heatFluxCache[blockLattice.getCellId(iX, iY, iZ)] = temperature;

          // computation of the bousinessq force
          T temperatureDifference = temperature - T0;
          blockLattice.get(iX,iY,iZ).template setField<descriptors::FORCE>(temperatureDifference*forcePrefactor);

          FieldD<T,DESCRIPTOR,descriptors::VELOCITY> u;
          blockLattice.get(iX,iY,iZ).computeU(u.data());
          tPartner->get(iX,iY,iZ).template setField<descriptors::VELOCITY>(u);
        }
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {

          auto u_ppp = tPartner->get(iX+1, iY+1, iZ+1).template getField<descriptors::VELOCITY>();
          auto u_0pp = tPartner->get(iX, iY+1, iZ+1).template getField<descriptors::VELOCITY>();
          auto u_npp = tPartner->get(iX-1, iY+1, iZ+1).template getField<descriptors::VELOCITY>();
          auto u_p0p = tPartner->get(iX+1, iY, iZ+1).template getField<descriptors::VELOCITY>();
          auto u_00p = tPartner->get(iX, iY, iZ+1).template getField<descriptors::VELOCITY>();
          auto u_n0p = tPartner->get(iX-1, iY, iZ+1).template getField<descriptors::VELOCITY>();
          auto u_pnp = tPartner->get(iX+1, iY-1, iZ+1).template getField<descriptors::VELOCITY>();
          auto u_0np = tPartner->get(iX,   iY-1, iZ+1).template getField<descriptors::VELOCITY>();
          auto u_nnp = tPartner->get(iX-1, iY-1, iZ+1).template getField<descriptors::VELOCITY>();

          auto u_pp0 = tPartner->get(iX+1, iY+1, iZ  ).template getField<descriptors::VELOCITY>();
          auto u_0p0 = tPartner->get(iX, iY+1, iZ  ).template getField<descriptors::VELOCITY>();
          auto u_np0 = tPartner->get(iX-1, iY+1, iZ  ).template getField<descriptors::VELOCITY>();
          auto u_p00 = tPartner->get(iX+1, iY, iZ  ).template getField<descriptors::VELOCITY>();
          auto u_000 = tPartner->get(iX, iY, iZ  ).template getField<descriptors::VELOCITY>();
          auto u_n00 = tPartner->get(iX-1, iY, iZ  ).template getField<descriptors::VELOCITY>();
          auto u_pn0 = tPartner->get(iX+1, iY-1, iZ  ).template getField<descriptors::VELOCITY>();
          auto u_0n0 = tPartner->get(iX,   iY-1, iZ  ).template getField<descriptors::VELOCITY>();
          auto u_nn0 = tPartner->get(iX-1, iY-1, iZ  ).template getField<descriptors::VELOCITY>();

          auto u_ppn = tPartner->get(iX+1, iY+1, iZ-1).template getField<descriptors::VELOCITY>();
          auto u_0pn = tPartner->get(iX, iY+1, iZ-1).template getField<descriptors::VELOCITY>();
          auto u_npn = tPartner->get(iX-1, iY+1, iZ-1).template getField<descriptors::VELOCITY>();
          auto u_p0n = tPartner->get(iX+1, iY, iZ-1).template getField<descriptors::VELOCITY>();
          auto u_00n = tPartner->get(iX, iY, iZ-1).template getField<descriptors::VELOCITY>();
          auto u_n0n = tPartner->get(iX-1, iY, iZ-1).template getField<descriptors::VELOCITY>();
          auto u_pnn = tPartner->get(iX+1, iY-1, iZ-1).template getField<descriptors::VELOCITY>();
          auto u_0nn = tPartner->get(iX,   iY-1, iZ-1).template getField<descriptors::VELOCITY>();
          auto u_nnn = tPartner->get(iX-1, iY-1, iZ-1).template getField<descriptors::VELOCITY>();

          const T *h_ppp = & heatFluxCache[blockLattice.getCellId(iX+1, iY+1, iZ+1)];
          const T *h_0pp = & heatFluxCache[blockLattice.getCellId(iX, iY+1, iZ+1)];
          const T *h_npp = & heatFluxCache[blockLattice.getCellId(iX-1, iY+1, iZ+1)];
          const T *h_p0p = & heatFluxCache[blockLattice.getCellId(iX+1, iY, iZ+1)];
          const T *h_pnp = & heatFluxCache[blockLattice.getCellId(iX+1, iY-1, iZ+1)];
          const T *h_00p = & heatFluxCache[blockLattice.getCellId(iX, iY, iZ+1)];
          const T *h_0np = & heatFluxCache[blockLattice.getCellId(iX, iY-1, iZ+1)];
          const T *h_n0p = & heatFluxCache[blockLattice.getCellId(iX-1, iY, iZ+1)];
          const T *h_nnp = & heatFluxCache[blockLattice.getCellId(iX-1, iY-1, iZ+1)];

          const T *h_pp0 = & heatFluxCache[blockLattice.getCellId(iX+1, iY+1, iZ  )];
          const T *h_0p0 = & heatFluxCache[blockLattice.getCellId(iX, iY+1, iZ  )];
          const T *h_np0 = & heatFluxCache[blockLattice.getCellId(iX-1, iY+1, iZ  )];
          const T *h_p00 = & heatFluxCache[blockLattice.getCellId(iX+1, iY, iZ  )];
          const T *h_pn0 = & heatFluxCache[blockLattice.getCellId(iX+1, iY-1, iZ  )];
          const T *h_000 = & heatFluxCache[blockLattice.getCellId(iX, iY, iZ  )];
          const T *h_0n0 = & heatFluxCache[blockLattice.getCellId(iX, iY-1, iZ  )];
          const T *h_n00 = & heatFluxCache[blockLattice.getCellId(iX-1, iY, iZ  )];
          const T *h_nn0 = & heatFluxCache[blockLattice.getCellId(iX-1, iY-1, iZ  )];

          const T *h_ppn = & heatFluxCache[blockLattice.getCellId(iX+1, iY+1, iZ-1)];
          const T *h_0pn = & heatFluxCache[blockLattice.getCellId(iX, iY+1, iZ-1)];
          const T *h_npn = & heatFluxCache[blockLattice.getCellId(iX-1, iY+1, iZ-1)];
          const T *h_p0n = & heatFluxCache[blockLattice.getCellId(iX+1, iY, iZ-1)];
          const T *h_pnn = & heatFluxCache[blockLattice.getCellId(iX+1, iY-1, iZ-1)];
          const T *h_00n = & heatFluxCache[blockLattice.getCellId(iX, iY, iZ-1)];
          const T *h_0nn = & heatFluxCache[blockLattice.getCellId(iX, iY-1, iZ-1)];
          const T *h_n0n = & heatFluxCache[blockLattice.getCellId(iX-1, iY, iZ-1)];
          const T *h_nnn = & heatFluxCache[blockLattice.getCellId(iX-1, iY-1, iZ-1)];

//        cout<<"h_ppn= "<<h_ppp[0]<<endl;
//        cout<<"h_0pn= "<<h_0pp[0]<<endl;
//        cout<<"h_npn= "<<h_npp[0]<<endl;
//        cout<<"h_p0n= "<<h_p0p[0]<<endl;
//        cout<<"h_00n= "<<h_00p[0]<<endl;
//        cout<<"h_n0n= "<<h_n0p[0]<<endl;
//        cout<<"h_pnn= "<<h_pnp[0]<<endl;
//        cout<<"h_0nn= "<<h_0np[0]<<endl;
//        cout<<"h_nnn= "<<h_nnp[0]<<endl;

          //Testfilter h 3D
          Vector<T,3> filtered_u;
          T filtered_h = (
                           (
                             ( h_ppp[0] + 2. * h_0pp[0] + h_npp[0]) + 2.*
                             ( h_p0p[0] + 2. * h_00p[0] + h_n0p[0]) +
                             ( h_pnp[0] + 2. * h_0np[0] + h_nnp[0])
                           ) * 0.25 * 0.25 + 2. *

                           (
                             ( h_pp0[0] + 2. * h_0p0[0] + h_np0[0]) + 2.*
                             ( h_p00[0] + 2. * h_000[0] + h_n00[0]) +
                             ( h_pn0[0] + 2. * h_0n0[0] + h_nn0[0])
                           ) * 0.25 * 0.25  +

                           (
                             ( h_ppn[0] + 2. * h_0pn[0] + h_npn[0]) + 2.*
                             ( h_p0n[0] + 2. * h_00n[0] + h_n0n[0]) +
                             ( h_pnn[0] + 2. * h_0nn[0] + h_nnn[0])
                           ) * 0.25 * 0.25
                         ) * 0.25;

          //Testfilter u 3D
          filtered_u = (
                         (
                           ( u_ppp + 2. * u_0pp + u_npp) + 2.*
                           ( u_p0p + 2. * u_00p + u_n0p) +
                           ( u_pnp + 2. * u_0np + u_nnp)
                         ) * 0.25 * 0.25 + 2. *

                         (
                           ( u_pp0 + 2. * u_0p0 + u_np0) + 2.*
                           ( u_p00 + 2. * u_000 + u_n00) +
                           ( u_pn0 + 2. * u_0n0 + u_nn0)
                         ) * 0.25 * 0.25 +

                         (
                           ( u_ppn + 2. * u_0pn + u_npn) + 2.*
                           ( u_p0n + 2. * u_00n + u_n0n) +
                           ( u_pnn + 2. * u_0nn + u_nnn)
                         ) * 0.25 * 0.25
                       ) * 0.25;

          //Ende Testfilter u
          // filtered_u[i] -= u_00[i];

//        cout << "u["<<i<<"]= "<< filtered_u[i] << std::endl;
//        cout << "h[0]= "<< filtered_h[0] << std::endl;
          const Vector<T,3> filtered_u_reduced = u_000 - filtered_u;
          const Vector<T,3> filtered_heatFlux_reduced = u_000*h_000[0] - filtered_u*filtered_h;
          T cutoffKinEnergy = filtered_u_reduced[0] * filtered_u_reduced[0]
                              + filtered_u_reduced[1] * filtered_u_reduced[1]
                              + filtered_u_reduced[2] * filtered_u_reduced[2];
          T cutoffHeatFlux = filtered_heatFlux_reduced[0] * filtered_heatFlux_reduced[0]
                             + filtered_heatFlux_reduced[1] * filtered_heatFlux_reduced[1]
                             + filtered_heatFlux_reduced[2] * filtered_heatFlux_reduced[2];
//          cout << "filtered_u_reduced= " << filtered_u_reduced << std::endl;
//          cout << "filtered_heatFlux_reduced= " << filtered_heatFlux_reduced << std::endl;
//          cout << "cutoffKinEnergy= " << cutoffKinEnergy << std::endl;
//          cout << "cutoffHeatFlux= " << cutoffHeatFlux << std::endl;

          blockLattice.get(iX,iY,iZ).template setField<descriptors::CUTOFF_KIN_ENERGY>(util::pow(0.5*cutoffKinEnergy, 0.25));
          tPartner->get(iX,iY,iZ).template setField<descriptors::CUTOFF_HEAT_FLUX>(util::pow(0.5*cutoffHeatFlux, 0.25));

//        std::cout<<"cutoffKinEnergy_14= "<<cutoffKinEnergy_14[0]<<std::endl;
//          std::cout<<"cutoffHeatFlux_14= "<<cutoffHeatFlux_14[0]<<std::endl;


          // tau coupling
          auto tauNS = blockLattice.get(iX,iY,iZ).template getField<descriptors::TAU_EFF>();
          auto tauAD = tPartner->get(iX,iY,iZ).template getField<descriptors::TAU_EFF>();


//        std::cout<<"tau_NS= "<<tauNS[0]<<std::endl;
//        std::cout<<"tau_AD= "<<tauAD[0]<<std::endl;

          /// Molecular realaxation time
          T tau_mol_NS = 1. / blockLattice.get(iX, iY, iZ).getDynamics()->getParameters(blockLattice).template getOrFallback<descriptors::OMEGA>(0);
          T tau_mol_AD = 1. / tPartner->get(iX, iY, iZ).getDynamics()->getParameters(*tPartner).template getOrFallback<descriptors::OMEGA>(0);

          const T temperature = h_000[0];

          // computation of the bousinessq force
          auto force = blockLattice.get(iX,iY,iZ).template getField<descriptors::FORCE>();
          // T temperatureDifference = temperature - T0;
          // for (unsigned iD = 0; iD < L::d; ++iD) {
          //   force[iD] = forcePrefactor[iD] * temperatureDifference;
          // }

          auto u = tPartner->get(iX,iY,iZ).template getFieldPointer<descriptors::VELOCITY>();
          T rho, pi[util::TensorVal<DESCRIPTOR>::n], j[DESCRIPTOR::d];
          // blockLattice.get(iX,iY).computeAllMomenta(rho, u, pi);
          rho = blockLattice.get(iX,iY,iZ).computeRho();
          blockLattice.get(iX,iY,iZ).computeStress(pi);

          int iPi = 0;
          const int dim = DESCRIPTOR::d;
          for (int Alpha=0; Alpha<dim; ++Alpha) {
            for (int Beta=Alpha; Beta<dim; ++Beta) {
              pi[iPi] += rho/2.*(force[Alpha]*u[Beta] + u[Alpha]*force[Beta]);
              ++iPi;
            }
          }
          const T piSqr[6] = {pi[0]*pi[0], pi[1]*pi[1], pi[2]*pi[2], pi[3]*pi[3], pi[4]*pi[4], pi[5]*pi[5]};
          const T PiNeqNormSqr = piSqr[0] + 2.0 * piSqr[1] + 2.0 * piSqr[2] + piSqr[3] + 2.0 * piSqr[4] + piSqr[5];
          const T PiNeqNorm    = util::sqrt(PiNeqNormSqr);

          tPartner->get(iX, iY, iZ).computeJ(j);

          const Vector<T,3> jNeq = {(j[0] - temperature * u[0]), (j[1] - temperature * u[1]), (j[2] - temperature * u[2])};

          const Vector<T,6> t = { 2.0 * jNeq[0] * pi[0], (jNeq[0] + jNeq[1]) * pi[1], (jNeq[0] + jNeq[2]) * pi[2], 2.0 * jNeq[1] * pi[3], (jNeq[1] + jNeq[2]) * pi[4], 2.0 * jNeq[2] * pi[5]};
          const Vector<T,6> tSqr = {t[0]*t[0], t[1]*t[1], t[2]*t[2], t[3]*t[3], t[4]*t[4], t[5]*t[5]};
          const T TnormSqr = tSqr[0] + 2.0 * tSqr[1] + 2.0 * tSqr[2] + tSqr[3] + 2.0 * tSqr[4] + tSqr[5];
          const T Tnorm    = util::sqrt(2.0 * 0.25 * TnormSqr);

          auto cutoffKinEnergy_14 = blockLattice.get(iX,iY,iZ).template getField<descriptors::CUTOFF_KIN_ENERGY>();
          auto cutoffHeatFlux_14 = tPartner->get(iX,iY,iZ).template getField<descriptors::CUTOFF_HEAT_FLUX>();

          const T tau_turb_NS = C_nu * util::sqrt(util::sqrt(2.)/2.) * invCs2 * invCs2 * util::sqrt(PiNeqNorm / rho / tauNS) * cutoffKinEnergy_14;
          const T tau_turb_AD = C_alpha * invCs2 / rho * util::sqrt(Tnorm * invCs2_g / tauNS / tauAD) * cutoffHeatFlux_14;

          /// Effective realaxation time
          blockLattice.get(iX,iY,iZ).template setField<descriptors::TAU_EFF>(tau_mol_NS + tau_turb_NS);
          tPartner->get(iX,iY,iZ).template setField<descriptors::TAU_EFF>(tau_mol_AD + tau_turb_AD);

        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void MixedScaleBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, typename DESCRIPTOR>
MixedScaleBoussinesqCouplingGenerator3D<T,DESCRIPTOR>::
MixedScaleBoussinesqCouplingGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                        T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_)
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_), PrTurb(PrTurb_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* MixedScaleBoussinesqCouplingGenerator3D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<3>* > partners) const
{
  return new MixedScaleBoussinesqCouplingPostProcessor3D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, this-> z0, this->z1, gravity, T0, deltaTemp, dir, PrTurb, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* MixedScaleBoussinesqCouplingGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new MixedScaleBoussinesqCouplingGenerator3D<T,DESCRIPTOR>(*this);
}


//=====================================================================================
//==============  VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingPostProcessor3D ===========
//=====================================================================================

template<typename T, typename DESCRIPTOR, typename POROSITY, typename ADLattice, typename FIELD_A, typename FIELD_B>
VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,POROSITY,ADLattice,FIELD_A,FIELD_B>::
VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int iC_,
    std::vector<BlockStructureD<3>* > partners_,
    std::vector<std::reference_wrapper<AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice>>> forces_)
  :  _forces(forces_),
     x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), iC(iC_),
     _partnerLattice(static_cast<BlockLattice<T,ADLattice> *>(partners_[0])),
     _cell(_partnerLattice->get(x0,y0,z0)),
     _cellXp(_partnerLattice->get(x0+1,y0,z0)),
     _cellXn(_partnerLattice->get(x0-1,y0,z0)),
     _cellYp(_partnerLattice->get(x0,y0+1,z0)),
     _cellYn(_partnerLattice->get(x0,y0-1,z0)),
     _cellZp(_partnerLattice->get(x0,y0,z0+1)),
     _cellZn(_partnerLattice->get(x0,y0,z0-1))
{ 
  this->getName() = "VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingPostProcessor3D";  
}

template<typename T, typename DESCRIPTOR, typename POROSITY, typename ADLattice, typename FIELD_A, typename FIELD_B>
void VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,POROSITY,ADLattice,FIELD_A,FIELD_B>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  auto vel     = par ? _cell.template getField<FIELD_B>() : _cell.template getField<FIELD_A>();
  auto vel_new = par ? _cell.template getFieldPointer<FIELD_A>() : _cell.template getFieldPointer<FIELD_B>();

  auto velXp = par ? _cellXp.template getFieldPointer<FIELD_B>() : _cellXp.template getFieldPointer<FIELD_A>();
  auto velXn = par ? _cellXn.template getFieldPointer<FIELD_B>() : _cellXn.template getFieldPointer<FIELD_A>();
  auto velYp = par ? _cellYp.template getFieldPointer<FIELD_B>() : _cellYp.template getFieldPointer<FIELD_A>();
  auto velYn = par ? _cellYn.template getFieldPointer<FIELD_B>() : _cellYn.template getFieldPointer<FIELD_A>();
  auto velZp = par ? _cellZp.template getFieldPointer<FIELD_B>() : _cellZp.template getFieldPointer<FIELD_A>();
  auto velZn = par ? _cellZn.template getFieldPointer<FIELD_B>() : _cellZn.template getFieldPointer<FIELD_A>();

  auto forceNS = blockLattice.get(x0,y0,z0).template getFieldPointer<descriptors::FORCE>();
  
  auto _cellNSXp = blockLattice.get(x0+1,y0,z0);
  auto _cellNSXn = blockLattice.get(x0-1,y0,z0);
  auto _cellNSYp = blockLattice.get(x0,y0+1,z0);
  auto _cellNSYn = blockLattice.get(x0,y0-1,z0);
  auto _cellNSZp = blockLattice.get(x0,y0,z0+1);
  auto _cellNSZn = blockLattice.get(x0,y0,z0-1);

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  T porosityForce[3];
  T pressure;

  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
	  auto nsCell = blockLattice.get(iX,iY,iZ);
	  T porosityNS = 1. - _cell.computeRho();
	  nsCell.template setField<POROSITY>(porosityNS);
	}
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          int latticeR[4] = {iC, iX, iY, iZ};
          T velGrad[3] = {0.,0.,0.};
          T forceValue[3] = {0.,0.,0.};
          T velF[3] = {0.,0.,0.};

          auto nsCell = blockLattice.get(iX,iY,iZ);	  
          pressure = nsCell.computeRho() / descriptors::invCs2<T,DESCRIPTOR>();
 
          if (_forces.begin() != _forces.end()) {
            // calculating upwind Gradient
            // vel contains velocity information on ADlattice
            // velGrad contains upwind vel on ADlattice
            if (vel[0]<0.) {
              velGrad[0] = vel[0]*(velXp[0]-vel[0]);
              velGrad[1] = vel[0]*(velXp[1]-vel[1]);
              velGrad[2] = vel[0]*(velXp[2]-vel[2]);
            }
            else {
              velGrad[0] = vel[0]*(vel[0]-velXn[0]);
              velGrad[1] = vel[0]*(vel[1]-velXn[1]);
              velGrad[2] = vel[0]*(vel[2]-velXn[2]);
            }
            if (vel[1]<0.) {
              velGrad[0] += vel[1]*(velYp[0]-vel[0]);
              velGrad[1] += vel[1]*(velYp[1]-vel[1]);
              velGrad[2] += vel[1]*(velYp[2]-vel[2]);
            }
            else {
              velGrad[0] += vel[1]*(vel[0]-velYn[0]);
              velGrad[1] += vel[1]*(vel[1]-velYn[1]);
              velGrad[2] += vel[1]*(vel[2]-velYn[2]);
            }
            if (vel[2]<0.) {
              velGrad[0] += vel[2]*(velZp[0]-vel[0]);
              velGrad[1] += vel[2]*(velZp[1]-vel[1]);
              velGrad[2] += vel[2]*(velZp[2]-vel[2]);
            }
            else {
              velGrad[0] += vel[2]*(vel[0]-velZn[0]);
              velGrad[1] += vel[2]*(vel[1]-velZn[1]);
              velGrad[2] += vel[2]*(vel[2]-velZn[2]);
            }

            for (AdvectionDiffusionForce3D<T, DESCRIPTOR, ADLattice>& f : _forces) {
              // writes force in forceValues, vel refers to ADlattice
              auto adCell = _partnerLattice->get(x0,y0,z0);
              f.applyForce(forceValue, &nsCell, &adCell, vel.data(), latticeR);
              if (par) {
                _cell.template setField<FIELD_B>(vel);
              } else {
                _cell.template setField<FIELD_A>(vel);
              }
            }

            // compute new particle velocity and opposite fluid force
            for (int i=0; i < DESCRIPTOR::d; i++) {
              vel_new[i] = vel[i] + forceValue[i] - velGrad[i];
	      forceNS[i] = -forceValue[i];    
            }
          }
          else {
	    
	    T porXp = _cellNSXp.template getField<POROSITY>();
            T porXn = _cellNSXn.template getField<POROSITY>();
            T porYp = _cellNSYp.template getField<POROSITY>();
            T porYn = _cellNSYn.template getField<POROSITY>();
            T porZp = _cellNSZp.template getField<POROSITY>();
            T porZn = _cellNSZn.template getField<POROSITY>();
  
            porosityForce[0] = 0.5 * pressure * (porXp - porXn);
	    porosityForce[1] = 0.5 * pressure * (porYp - porYn);
            porosityForce[2] = 0.5 * pressure * (porZp - porZn);
	    
	    for (int i = 0; i < DESCRIPTOR::d; i++) {
	      forceNS[i] += porosityForce[i];
	    }

            nsCell.computeU(velF);

            for (int i = 0; i < DESCRIPTOR::d; i++) {   // set particle velocity to fluid velocity
              vel_new[i] = velF[i];
            }
          }
        }
      }
    }
  }
  par = !par;
}

template<typename T, typename DESCRIPTOR, typename POROSITY, typename ADLattice, typename FIELD_A, typename FIELD_B>
void VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,POROSITY,ADLattice,FIELD_A,FIELD_B>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

// LatticeCouplingGenerator for one-way VANS advectionDiffusion coupling with Stokes drag

template<typename T, typename DESCRIPTOR, typename POROSITY, typename ADLattice, typename FIELD_A, typename FIELD_B>
VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,POROSITY,ADLattice,FIELD_A,FIELD_B>::
VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingGenerator3D()
  : LatticeCouplingGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0)
{ }

template<typename T, typename DESCRIPTOR, typename POROSITY, typename ADLattice, typename FIELD_A, typename FIELD_B>
PostProcessor3D<T,DESCRIPTOR>* VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,POROSITY,ADLattice,FIELD_A,FIELD_B>::generate(
  std::vector<BlockStructureD<3>* > partners) const
{
  return new VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingPostProcessor3D<T,DESCRIPTOR,POROSITY,ADLattice,FIELD_A,FIELD_B>(
																	  this->x0,this->x1,this->y0,this->y1,this->z0,this->z1, this->iC, partners, ADforces);
}

template<typename T, typename DESCRIPTOR, typename POROSITY, typename ADLattice, typename FIELD_A, typename FIELD_B>
LatticeCouplingGenerator3D<T,DESCRIPTOR>* VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,POROSITY,ADLattice,FIELD_A,FIELD_B>::clone() const
{
  return new VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,POROSITY,ADLattice,FIELD_A,FIELD_B>(*this);
}

template<typename T, typename DESCRIPTOR, typename POROSITY, typename ADLattice, typename FIELD_A, typename FIELD_B>
void VolumeAveragedNavierStokesAdvectionDiffusionParticleCouplingGenerator3D<T,DESCRIPTOR,POROSITY,ADLattice,FIELD_A,FIELD_B>::addForce(
  AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> &force)
{
  ADforces.push_back(force);
}

}  // namespace olb

#endif
