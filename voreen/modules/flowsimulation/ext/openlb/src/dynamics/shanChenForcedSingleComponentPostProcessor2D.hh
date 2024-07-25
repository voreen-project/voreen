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

#ifndef SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_2D_HH
#define SHAN_CHEN_FORCED_SINGLE_COMPONENT_POST_PROCESSOR_2D_HH

#include "shanChenForcedSingleComponentPostProcessor2D.h"
#include "interactionPotential.h"
#include "core/util.h"
#include "core/finiteDifference2D.h"

namespace olb {

////////  ShanChenForcedSingleComponentPostProcessor2D ///////////////////////////////////


template<typename T, typename DESCRIPTOR>
ShanChenForcedSingleComponentPostProcessor2D <T,DESCRIPTOR>::
ShanChenForcedSingleComponentPostProcessor2D(int x0_, int x1_, int y0_, int y1_, T G_,
    std::vector<T> rho0_,
    AnalyticalF<1,T,T>& iP_,
    std::vector<BlockStructureD<2>*> partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_), G(G_), rho0(rho0_), interactionPotential(iP_), partners(partners_)
{
  this->getName() = "ShanChenForcedSingleComponentPostProcessor2D";
}

template<typename T, typename DESCRIPTOR>
ShanChenForcedSingleComponentPostProcessor2D <T,DESCRIPTOR>::
ShanChenForcedSingleComponentPostProcessor2D(T G_,
    std::vector<T> rho0_,
    AnalyticalF<1,T,T>& iP_,
    std::vector<BlockStructureD<2>*> partners_)
  :  x0(0), x1(0), y0(0), y1(0), G(G_), rho0(rho0_), interactionPotential(iP_), partners(partners_)
{
  this->getName() = "ShanChenForcedSingleComponentPostProcessor2D";
}

template<typename T, typename DESCRIPTOR>
void ShanChenForcedSingleComponentPostProcessor2D<T,DESCRIPTOR>::
processSubDomain( BlockLattice<T,DESCRIPTOR>& blockLattice,
                  int x0_, int x1_, int y0_, int y1_ )
{
  typedef DESCRIPTOR L;

  int newX0, newX1, newY0, newY1;
  if ( util::intersect ( x0, x1, y0, y1,
                         x0_, x1_, y0_, y1_,
                         newX0, newX1, newY0, newY1 ) ) {

    auto& rhoField = blockLattice.template getField<RHO_CACHE>()[0];

    // Compute density and velocity on every site of first lattice, and store result
    //   in external scalars; envelope cells are included, because they are needed
    //   to compute the interaction potential in what follows.
    for (int iX=newX0-1; iX<=newX1+1; ++iX) {
      for (int iY=newY0-1; iY<=newY1+1; ++iY) {
        Cell<T,DESCRIPTOR> cell = blockLattice.get(iX,iY);
        rhoField[cell.getCellId()] = cell.computeRho()*rho0[0];
      }
    }

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        Cell<T,DESCRIPTOR> blockCell = blockLattice.get(iX,iY);

        auto j = blockCell.template getField<descriptors::VELOCITY>();
        lbm<DESCRIPTOR>::computeJ(blockCell,j);
        blockCell.template setField<descriptors::VELOCITY>(j);

        T blockOmega   = blockCell.getDynamics()->getParameters(blockLattice).template getOrFallback<descriptors::OMEGA>(0);

        // Computation of the common velocity, shared among the two populations
        T rhoTot = rhoField[blockCell.getCellId()]*blockOmega;

        Vector<T, 2> uTot;
        auto blockU = blockCell.template getField<descriptors::VELOCITY>();     // contains precomputed value rho*u
        uTot = (blockU*rho0[0]*blockOmega) / rhoTot;

        // Computation of the interaction potential
        Vector<T, 2> rhoBlockContribution;
        T psi;
        interactionPotential(&psi, &rhoField[blockCell.getCellId()]);
        for (int iPop = 0; iPop < L::q; ++iPop) {
          int nextX = iX + descriptors::c<L>(iPop,0);
          int nextY = iY + descriptors::c<L>(iPop,1);
          T blockRho;
          interactionPotential(&blockRho, &rhoField[blockLattice.getCellId(nextX, nextY)]);//rho0[0];
          rhoBlockContribution   += psi * blockRho * descriptors::c<L>(iPop)* descriptors::t<T,L>(iPop);
        }

        // Computation and storage of the final velocity, consisting
        //   of u and the momentum difference due to interaction
        //   potential plus external force
        auto externalBlockForce   = blockCell.template getField<descriptors::EXTERNAL_FORCE>();

        blockCell.template setField<descriptors::VELOCITY>(uTot);
        blockCell.template setField<descriptors::FORCE>(externalBlockForce
            - G*rhoBlockContribution/rhoField[blockCell.getCellId()]);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ShanChenForcedSingleComponentPostProcessor2D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


/// LatticeCouplingGenerator for NS coupling

template<typename T, typename DESCRIPTOR>
ShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>::ShanChenForcedSingleComponentGenerator2D (
  int x0_, int x1_, int y0_, int y1_, T G_, std::vector<T> rho0_, AnalyticalF<1,T,T>& iP_ )
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), G(G_), rho0(rho0_), interactionPotential(iP_)
{ }

template<typename T, typename DESCRIPTOR>
ShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>::ShanChenForcedSingleComponentGenerator2D (
  T G_, std::vector<T> rho0_, AnalyticalF<1,T,T>& iP_ )
  : LatticeCouplingGenerator2D<T,DESCRIPTOR>(0, 0, 0, 0), G(G_), rho0(rho0_), interactionPotential(iP_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* ShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>::generate (
  std::vector<BlockStructureD<2>*> partners) const
{
  return new ShanChenForcedSingleComponentPostProcessor2D<T,DESCRIPTOR>(
           this->x0,this->x1,this->y0,this->y1,G, rho0, interactionPotential, partners);
}

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>* ShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new ShanChenForcedSingleComponentGenerator2D<T,DESCRIPTOR>(*this);
}



}  // namespace olb

#endif
