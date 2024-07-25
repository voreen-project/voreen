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

#ifndef NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_H
#define NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_H

#include "core/blockStructure.h"
#include "core/postProcessing.h"
#include "utilities/omath.h"


namespace olb {


/**
* Class for the coupling between a Navier-Stokes (NS) lattice and an
* Advection-Diffusion (AD) lattice.
*/

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
// This coupling must be necessarily be put on the Navier-Stokes lattice!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

//======================================================================
// ======== AD coupling with Boussinesq bouancy 2D ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class NavierStokesAdvectionDiffusionCouplingPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  NavierStokesAdvectionDiffusionCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
      std::vector<BlockStructureD<2>* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY>> *tPartner;
  T forcePrefactor[L::d];

  std::vector<BlockStructureD<2>*> partners;
};

template<typename T, typename DESCRIPTOR>
class NavierStokesAdvectionDiffusionCouplingGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  NavierStokesAdvectionDiffusionCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_);
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<BlockStructureD<2>* > partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
};


//======================================================================
// ======== Total enthalpy coupling with Boussinesq bouancy 2D and phase change====================//
//======================================================================
template<typename T, typename DESCRIPTOR, typename DYNAMICS>
class TotalEnthalpyPhaseChangeCouplingPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  TotalEnthalpyPhaseChangeCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
      std::vector<BlockStructureD<2>* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TEMPERATURE>> *tPartner;
  T forcePrefactor[L::d];

  std::vector<BlockStructureD<2>*> partners;
};

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class TotalEnthalpyPhaseChangeCouplingGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  TotalEnthalpyPhaseChangeCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_);
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<BlockStructureD<2>* > partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
};


//======================================================================
// ======== Phase field coupling without bouancy 2D ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class PhaseFieldCouplingPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  PhaseFieldCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
                                    T rho_L, T rho_H, T mu_L, T mu_H, T surface_tension, T interface_thickness,
                                    std::vector<BlockStructureD<2>* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_) override;
private:
  using L = DESCRIPTOR;
  using PHI_CACHE = descriptors::FIELD_BASE<1, 0, 0>;

  int x0, x1, y0, y1;

  T _rho_L, _rho_H, _delta_rho;
  T _mu_L, _mu_H;
  T _surface_tension, _interface_thickness;
  T _beta, _kappa;

  BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::INTERPHASE_NORMAL>> *tPartner;

  std::vector<BlockStructureD<2>*> partners;
};

template<typename T, typename DESCRIPTOR>
class PhaseFieldCouplingGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  PhaseFieldCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
                                T rho_L, T rho_H, T mu_L, T mu_H, T surface_tension, T interface_thickness);
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<BlockStructureD<2>* > partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;

private:
  T _rho_L, _rho_H, _delta_rho;
  T _mu_L, _mu_H;
  T _surface_tension, _interface_thickness;
};


//======================================================================
// ======== AD coupling with Boussinesq bouancy for Smagorinsky-LES ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class SmagorinskyBoussinesqCouplingPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  SmagorinskyBoussinesqCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_, T smagoPrefactor_,
      std::vector<BlockStructureD<2>* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_) override;
private:
  typedef DESCRIPTOR L;
  int x0, x1, y0, y1;
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
  BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF>> *tPartner;
  T forcePrefactor[L::d];
  T tauTurbADPrefactor;
  T smagoPrefactor;

  std::vector<BlockStructureD<2>*> partners;
};

template<typename T, typename DESCRIPTOR>
class SmagorinskyBoussinesqCouplingGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  SmagorinskyBoussinesqCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_, T smagoPrefactor_);
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<BlockStructureD<2>* > partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp;
  std::vector<T> dir;
  T PrTurb;
  T smagoPrefactor;
};

//======================================================================
// ======== AD coupling with Boussinesq bouancy for Mixed Scale-LES ====================//
//======================================================================
template<typename T, typename DESCRIPTOR>
class MixedScaleBoussinesqCouplingPostProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  MixedScaleBoussinesqCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
      T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_,
      std::vector<BlockStructureD<2>* > partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_) override;
private:
  typedef DESCRIPTOR L;
  using HEAT_FLUX_CACHE = descriptors::FIELD_BASE<1, 0, 0>;
  int x0, x1, y0, y1;
  T gravity, T0, deltaTemp, PrTurb;
  std::vector<T> dir;
  BlockLattice<T,descriptors::D2Q5<descriptors::VELOCITY,descriptors::TAU_EFF,descriptors::CUTOFF_HEAT_FLUX>> *tPartner;
  Vector<T, L::d> forcePrefactor;
  T tauTurbADPrefactor;

  std::vector<BlockStructureD<2>*> partners;
};

template<typename T, typename DESCRIPTOR>
class MixedScaleBoussinesqCouplingGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  MixedScaleBoussinesqCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
                                          T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_, T PrTurb_);
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<BlockStructureD<2>* > partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;

private:
  T gravity, T0, deltaTemp, PrTurb;
  std::vector<T> dir;
};

}

#endif
