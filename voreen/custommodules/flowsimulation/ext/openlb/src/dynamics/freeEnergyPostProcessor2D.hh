/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Robin Trunk, Sam Avis
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

#ifndef FREE_ENERGY_POST_PROCESSOR_2D_HH
#define FREE_ENERGY_POST_PROCESSOR_2D_HH

#include "freeEnergyPostProcessor2D.h"

namespace olb {

////////  FreeEnergyChemicalPotentialCoupling2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialCoupling2D <T,DESCRIPTOR>::FreeEnergyChemicalPotentialCoupling2D (
  int x0_, int x1_, int y0_, int y1_, T alpha_, T kappa1_, T kappa2_, T kappa3_,
  std::vector<BlockStructureD<2>*> partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), alpha(alpha_), kappa1(kappa1_),
     kappa2(kappa2_), kappa3(kappa3_), partners(partners_)
{
  this->getName() = "FreeEnergyChemicalPotentialCoupling2D";
}

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialCoupling2D <T,DESCRIPTOR>::FreeEnergyChemicalPotentialCoupling2D (
  T alpha_, T kappa1_, T kappa2_, T kappa3_, std::vector<BlockStructureD<2>*> partners_)
  :  x0(0), x1(0), y0(0), y1(0), alpha(alpha_), kappa1(kappa1_), kappa2(kappa2_),
     kappa3(kappa3_), partners(partners_)
{
  this->getName() = "FreeEnergyChemicalPotentialCoupling2D";
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyChemicalPotentialCoupling2D<T,DESCRIPTOR>::processSubDomain (
  BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_ )
{
  // If partners.size() == 1: two fluid components
  // If partners.size() == 2: three fluid components
  BlockLattice<T,DESCRIPTOR> *partnerLattice1 = static_cast<BlockLattice<T,DESCRIPTOR> *>(partners[0]);
  BlockLattice<T,DESCRIPTOR> *partnerLattice2 = 0;
  if (partners.size() > 1) {
    partnerLattice2 = static_cast<BlockLattice<T,DESCRIPTOR> *>(partners[1]);
  }

  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( x0, x1, y0, y1,
                         x0_, x1_, y0_, y1_,
                         newX0, newX1, newY0, newY1 ) ) {

    auto& rhoField = blockLattice.template getField<RHO_CACHE>();

    for (int iX=newX0-1; iX<=newX1+1; ++iX)
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        rhoField[0][blockLattice.getCellId(iX, iY)] = blockLattice.get(iX,iY).computeRho();
      }
    for (int iX=newX0-1; iX<=newX1+1; ++iX)
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        rhoField[1][blockLattice.getCellId(iX, iY)] = partnerLattice1->get(iX,iY).computeRho();
      }
    if (partners.size() > 1) {
      for (int iX=newX0-1; iX<=newX1+1; ++iX)
        for (int iY=newY0-1; iY<=newY1+1; ++iY) {
          rhoField[2][blockLattice.getCellId(iX, iY)] = partnerLattice2->get(iX,iY).computeRho();
        }
    }

    // calculate chemical potential
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        T densitySum = rhoField[0][blockLattice.getCellId(iX, iY)]
                       + rhoField[1][blockLattice.getCellId(iX, iY)];
        T densityDifference = rhoField[0][blockLattice.getCellId(iX, iY)]
                              - rhoField[1][blockLattice.getCellId(iX, iY)];
        if (partners.size() > 1) {
          densitySum -= rhoField[2][blockLattice.getCellId(iX, iY)];
          densityDifference -= rhoField[2][blockLattice.getCellId(iX, iY)];
        }
        T term1 = 0.125 * kappa1 * (densitySum)
                  * (densitySum-1.) * (densitySum-2.);
        T term2 = 0.125 * kappa2 * (densityDifference)
                  * (densityDifference-1.) * (densityDifference-2.);
        T term3 = 0.;
        if (partners.size() > 1) {
          T rho3 = rhoField[2][blockLattice.getCellId(iX, iY)];
          term3 = kappa3 * rho3 * (rho3 - 1.) * (2.*rho3 - 1.);
        }

        T laplaceRho1 = 0.25 * (
                          rhoField[0][blockLattice.getCellId(iX-1, iY-1)]
                          + 2. * rhoField[0][blockLattice.getCellId(iX, iY-1)]
                          +      rhoField[0][blockLattice.getCellId(iX+1, iY-1)]
                          + 2. * rhoField[0][blockLattice.getCellId(iX-1, iY)]
                          -12. * rhoField[0][blockLattice.getCellId(iX, iY)]
                          + 2. * rhoField[0][blockLattice.getCellId(iX+1, iY)]
                          +      rhoField[0][blockLattice.getCellId(iX-1, iY+1)]
                          + 2. * rhoField[0][blockLattice.getCellId(iX, iY+1)]
                          +      rhoField[0][blockLattice.getCellId(iX+1, iY+1)]
                        );
        T laplaceRho2 = 0.25 * (
                          rhoField[1][blockLattice.getCellId(iX-1, iY-1)]
                          + 2. * rhoField[1][blockLattice.getCellId(iX, iY-1)]
                          +      rhoField[1][blockLattice.getCellId(iX+1, iY-1)]
                          + 2. * rhoField[1][blockLattice.getCellId(iX-1, iY)]
                          -12. * rhoField[1][blockLattice.getCellId(iX, iY)]
                          + 2. * rhoField[1][blockLattice.getCellId(iX+1, iY)]
                          +      rhoField[1][blockLattice.getCellId(iX-1, iY+1)]
                          + 2. * rhoField[1][blockLattice.getCellId(iX, iY+1)]
                          +      rhoField[1][blockLattice.getCellId(iX+1, iY+1)]
                        );
        T laplaceRho3 = 0.;
        if (partners.size() > 1) {
          laplaceRho3 = 0.25 * (
                          rhoField[2][blockLattice.getCellId(iX-1, iY-1)]
                          + 2. * rhoField[2][blockLattice.getCellId(iX, iY-1)]
                          +      rhoField[2][blockLattice.getCellId(iX+1, iY-1)]
                          + 2. * rhoField[2][blockLattice.getCellId(iX-1, iY)]
                          -12. * rhoField[2][blockLattice.getCellId(iX, iY)]
                          + 2. * rhoField[2][blockLattice.getCellId(iX+1, iY)]
                          +      rhoField[2][blockLattice.getCellId(iX-1, iY+1)]
                          + 2. * rhoField[2][blockLattice.getCellId(iX, iY+1)]
                          +      rhoField[2][blockLattice.getCellId(iX+1, iY+1)]
                        );
        }

        // setting chemical potential to the respective lattices
        blockLattice.get(iX, iY).template setField<descriptors::CHEM_POTENTIAL>(term1 + term2
            + 0.25*alpha*alpha*( (kappa2 - kappa1) * laplaceRho2
                                 +(kappa2 + kappa1) * (laplaceRho3 - laplaceRho1) ));
        partnerLattice1->get(iX, iY).template setField<descriptors::CHEM_POTENTIAL>(term1 - term2
            + 0.25*alpha*alpha*( (kappa2 - kappa1) * (laplaceRho1 - laplaceRho3)
                                 -(kappa2 + kappa1) * laplaceRho2 ));
        if (partners.size() > 1) {
          partnerLattice2->get(iX, iY).template setField<descriptors::CHEM_POTENTIAL>(- term1 - term2 + term3
              + 0.25*alpha*alpha*( (kappa2 + kappa1) * laplaceRho1
                                   -(kappa2 - kappa1) * laplaceRho2
                                   -(kappa2 + kappa1 + 4.*kappa3) * laplaceRho3 ));
        }
      }
    }

  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyChemicalPotentialCoupling2D<T,DESCRIPTOR>::process (
  BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


////////  FreeEnergyForceCoupling2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyForceCoupling2D <T,DESCRIPTOR>::FreeEnergyForceCoupling2D (
  int x0_, int x1_, int y0_, int y1_,
  std::vector<BlockStructureD<2>*> partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), partners(partners_)
{
  this->getName() = "FreeEnergyForceCoupling2D";
}

template<typename T, typename DESCRIPTOR>
FreeEnergyForceCoupling2D <T,DESCRIPTOR>::FreeEnergyForceCoupling2D (
  std::vector<BlockStructureD<2>*> partners_)
  :  x0(0), x1(0), y0(0), y1(0), partners(partners_)
{
  this->getName() = "FreeEnergyForceCoupling2D";
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyForceCoupling2D<T,DESCRIPTOR>::processSubDomain (
  BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_ )
{
  // If partners.size() == 1: two fluid components
  // If partners.size() == 2: three fluid components
  BlockLattice<T,DESCRIPTOR> *partnerLattice1 = static_cast<BlockLattice<T,DESCRIPTOR> *>(partners[0]);
  BlockLattice<T,DESCRIPTOR> *partnerLattice2 = nullptr;
  if (partners.size() > 1) {
    partnerLattice2 = static_cast<BlockLattice<T,DESCRIPTOR> *>(partners[1]);
  }

  int newX0, newX1, newY0, newY1;
  if ( util::intersect(x0, x1, y0, y1,
                       x0_, x1_, y0_, y1_,
                       newX0, newX1, newY0, newY1) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        T phi = blockLattice.get(iX,iY).computeRho();
        T rho = partnerLattice1->get(iX,iY).computeRho();
        T gradMuPhiX = 1./12. * ( -blockLattice.get(iX-1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  - 4.* blockLattice.get(iX-1,iY  ).template getField<descriptors::CHEM_POTENTIAL>()
                                  -  blockLattice.get(iX-1,iY+1).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  blockLattice.get(iX+1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  + 4.* blockLattice.get(iX+1,iY  ).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  blockLattice.get(iX+1,iY+1).template getField<descriptors::CHEM_POTENTIAL>() );
        T gradMuPhiY = 1./12. * ( -blockLattice.get(iX-1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  - 4.* blockLattice.get(iX,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  -  blockLattice.get(iX+1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  blockLattice.get(iX-1,iY+1).template getField<descriptors::CHEM_POTENTIAL>()
                                  + 4.* blockLattice.get(iX,iY+1).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  blockLattice.get(iX+1,iY+1).template getField<descriptors::CHEM_POTENTIAL>() );
        T gradMuRhoX = 1./12. * ( -partnerLattice1->get(iX-1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  - 4.* partnerLattice1->get(iX-1,iY  ).template getField<descriptors::CHEM_POTENTIAL>()
                                  -  partnerLattice1->get(iX-1,iY+1).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  partnerLattice1->get(iX+1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  + 4.* partnerLattice1->get(iX+1,iY  ).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  partnerLattice1->get(iX+1,iY+1).template getField<descriptors::CHEM_POTENTIAL>() );
        T gradMuRhoY = 1./12. * ( -partnerLattice1->get(iX-1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  - 4.* partnerLattice1->get(iX,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  -  partnerLattice1->get(iX+1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  partnerLattice1->get(iX-1,iY+1).template getField<descriptors::CHEM_POTENTIAL>()
                                  + 4.* partnerLattice1->get(iX,iY+1).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  partnerLattice1->get(iX+1,iY+1).template getField<descriptors::CHEM_POTENTIAL>() );
        T psi = 0.;
        T gradMuPsiX = 0.;
        T gradMuPsiY = 0.;
        if (partners.size() > 1) {
          psi = partnerLattice2->get(iX,iY).computeRho();
          gradMuPsiX = 1./12. * ( -partnerLattice2->get(iX-1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  - 4.* partnerLattice2->get(iX-1,iY  ).template getField<descriptors::CHEM_POTENTIAL>()
                                  -  partnerLattice2->get(iX-1,iY+1).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  partnerLattice2->get(iX+1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  + 4.* partnerLattice2->get(iX+1,iY  ).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  partnerLattice2->get(iX+1,iY+1).template getField<descriptors::CHEM_POTENTIAL>() );
          gradMuPsiY = 1./12. * ( -partnerLattice2->get(iX-1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  - 4.* partnerLattice2->get(iX,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  -  partnerLattice2->get(iX+1,iY-1).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  partnerLattice2->get(iX-1,iY+1).template getField<descriptors::CHEM_POTENTIAL>()
                                  + 4.* partnerLattice2->get(iX,iY+1).template getField<descriptors::CHEM_POTENTIAL>()
                                  +  partnerLattice2->get(iX+1,iY+1).template getField<descriptors::CHEM_POTENTIAL>() );
        }

        auto partnerCell = partnerLattice1->get(iX, iY);
        partnerCell.template setField<descriptors::FORCE>({
          -rho*gradMuRhoX - phi*gradMuPhiX - psi*gradMuPsiX,
            -rho*gradMuRhoY - phi*gradMuPhiY - psi*gradMuPsiY
          });
        T u[2] { };
        partnerCell.computeU(u);
        blockLattice.get(iX, iY).template setField<descriptors::FORCE>(u);
        if (partners.size() > 1) {
          partnerLattice2->get(iX, iY).template setField<descriptors::FORCE>(u);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyForceCoupling2D<T,DESCRIPTOR>::process(
  BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


////////  FreeEnergyInletOutletCoupling2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyInletOutletCoupling2D <T,DESCRIPTOR>::FreeEnergyInletOutletCoupling2D (
  int x0_, int x1_, int y0_, int y1_, std::vector<BlockStructureD<2>*> partners_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), partners(partners_)
{
  this->getName() = "FreeEnergyInletOutletCoupling2D";
}

template<typename T, typename DESCRIPTOR>
FreeEnergyInletOutletCoupling2D <T,DESCRIPTOR>::FreeEnergyInletOutletCoupling2D (
  std::vector<BlockStructureD<2>*> partners_)
  : x0(0), x1(0), y0(0), y1(0), partners(partners_)
{
  this->getName() = "FreeEnergyInletOutletCoupling2D";
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyInletOutletCoupling2D<T,DESCRIPTOR>::processSubDomain (
  BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_ )
{
  // If partners.size() == 1: two fluid components
  // If partners.size() == 2: three fluid components
  BlockLattice<T,DESCRIPTOR> *partnerLattice1 = static_cast<BlockLattice<T,DESCRIPTOR> *>(partners[0]);
  BlockLattice<T,DESCRIPTOR> *partnerLattice2 = 0;
  if (partners.size() > 1) {
    partnerLattice2 = static_cast<BlockLattice<T,DESCRIPTOR> *>(partners[1]);
  }

  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( x0, x1, y0, y1,
                         x0_, x1_, y0_, y1_,
                         newX0, newX1, newY0, newY1 ) ) {
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        T u[DESCRIPTOR::d];
        partnerLattice1->get(iX,iY).computeU(u);
        blockLattice.get(iX,iY).defineU(u);
        if (partners.size() > 1) {
          partnerLattice2->get(iX,iY).defineU(u);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyInletOutletCoupling2D<T,DESCRIPTOR>::process(
  BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


////////  FreeEnergyDensityOutletCoupling2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyDensityOutletCoupling2D <T,DESCRIPTOR>::FreeEnergyDensityOutletCoupling2D (
  int x0_, int x1_, int y0_, int y1_, T rho_,
  std::vector<BlockStructureD<2>*> partners_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), rho(rho_), partners(partners_)
{
  this->getName() = "FreeEnergyDensityOutletCoupling2D";
}

template<typename T, typename DESCRIPTOR>
FreeEnergyDensityOutletCoupling2D <T,DESCRIPTOR>::FreeEnergyDensityOutletCoupling2D (
  T rho_, std::vector<BlockStructureD<2>*> partners_)
  : x0(0), x1(0), y0(0), y1(0), rho(rho_), partners(partners_)
{
  this->getName() = "FreeEnergyDensityOutletCoupling2D";
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyDensityOutletCoupling2D<T,DESCRIPTOR>::processSubDomain (
  BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_ )
{
  // If partners.size() == 1: two fluid components
  // If partners.size() == 2: three fluid components
  BlockLattice<T,DESCRIPTOR> *partnerLattice1 = static_cast<BlockLattice<T,DESCRIPTOR> *>(partners[0]);
  BlockLattice<T,DESCRIPTOR> *partnerLattice2 = nullptr;
  if (partners.size() > 1) {
    partnerLattice2 = static_cast<BlockLattice<T,DESCRIPTOR> *>(partners[1]);
  }

  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( x0, x1, y0, y1,
                         x0_, x1_, y0_, y1_,
                         newX0, newX1, newY0, newY1 ) ) {
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        T rho0, phi, psi;
        rho0 = blockLattice.get(iX,iY).computeRho();
        phi = partnerLattice1->get(iX,iY).computeRho();
        blockLattice.get(iX,iY).defineRho(rho);
        partnerLattice1->get(iX,iY).defineRho(phi * rho / rho0);
        if (partners.size() > 1) {
          psi = partnerLattice2->get(iX,iY).computeRho();
          partnerLattice2->get(iX,iY).defineRho(psi * rho / rho0);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyDensityOutletCoupling2D<T,DESCRIPTOR>::process(
  BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


////////  FreeEnergyChemicalPotentialGenerator2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialGenerator2D<T,DESCRIPTOR>::FreeEnergyChemicalPotentialGenerator2D (
  int x0_, int x1_, int y0_, int y1_, T alpha_, T kappa1_, T kappa2_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), alpha(alpha_),
    kappa1(kappa1_), kappa2(kappa2_), kappa3(0)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialGenerator2D<T,DESCRIPTOR>::FreeEnergyChemicalPotentialGenerator2D (
  T alpha_, T kappa1_, T kappa2_ )
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(0, 0, 0, 0), alpha(alpha_),
    kappa1(kappa1_), kappa2(kappa2_), kappa3(0)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialGenerator2D<T,DESCRIPTOR>::FreeEnergyChemicalPotentialGenerator2D (
  int x0_, int x1_, int y0_, int y1_, T alpha_, T kappa1_, T kappa2_, T kappa3_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), alpha(alpha_),
    kappa1(kappa1_), kappa2(kappa2_), kappa3(kappa3_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyChemicalPotentialGenerator2D<T,DESCRIPTOR>::FreeEnergyChemicalPotentialGenerator2D (
  T alpha_, T kappa1_, T kappa2_, T kappa3_ )
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(0, 0, 0, 0), alpha(alpha_),
    kappa1(kappa1_), kappa2(kappa2_), kappa3(kappa3_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* FreeEnergyChemicalPotentialGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<2>*> partners) const
{
  return new FreeEnergyChemicalPotentialCoupling2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, alpha, kappa1, kappa2, kappa3, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>*
FreeEnergyChemicalPotentialGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyChemicalPotentialGenerator2D<T,DESCRIPTOR>(*this);
}

////////  FreeEnergyForceGenerator2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyForceGenerator2D<T,DESCRIPTOR>::FreeEnergyForceGenerator2D (
  int x0_, int x1_, int y0_, int y1_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyForceGenerator2D<T,DESCRIPTOR>::FreeEnergyForceGenerator2D ( )
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(0, 0, 0, 0)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* FreeEnergyForceGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<2>*> partners) const
{
  return new FreeEnergyForceCoupling2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* FreeEnergyForceGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyForceGenerator2D<T,DESCRIPTOR>(*this);
}

////////  FreeEnergyInletOutletGenerator2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyInletOutletGenerator2D<T,DESCRIPTOR>::FreeEnergyInletOutletGenerator2D (
  int x0_, int x1_, int y0_, int y1_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyInletOutletGenerator2D<T,DESCRIPTOR>::FreeEnergyInletOutletGenerator2D ( )
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(0, 0, 0, 0)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* FreeEnergyInletOutletGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<2>*> partners) const
{
  return new FreeEnergyInletOutletCoupling2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* FreeEnergyInletOutletGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyInletOutletGenerator2D<T,DESCRIPTOR>(*this);
}

////////  FreeEnergyDensityOutletGenerator2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyDensityOutletGenerator2D<T,DESCRIPTOR>::FreeEnergyDensityOutletGenerator2D (
  int x0_, int x1_, int y0_, int y1_, T rho_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), rho(rho_)
{ }

template<typename T, typename DESCRIPTOR>
FreeEnergyDensityOutletGenerator2D<T,DESCRIPTOR>::FreeEnergyDensityOutletGenerator2D (
  T rho_)
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(0, 0, 0, 0), rho(rho_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* FreeEnergyDensityOutletGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<2>*> partners) const
{
  return new FreeEnergyDensityOutletCoupling2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1, rho, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* FreeEnergyDensityOutletGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyDensityOutletGenerator2D<T,DESCRIPTOR>(*this);
}


}  // namespace olb

#endif
