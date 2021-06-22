/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 the OpenLB project
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

#ifndef DEFINE_RHO_2D_HH
#define DEFINE_RHO_2D_HH

#include <vector>
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"
#include "geometry/superGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "core/superLattice2D.h"
#include "core/blockLatticeStructure2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometry2D.h"

//defineRho for offLatticeBoundaryConditions
namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineRho(SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry, int material,
               AnalyticalF2D<T,T>& rho,
               std::vector<int> bulkMaterials)
{
  defineRho<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material),
                          superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
                          rho);
}


template<typename T, typename DESCRIPTOR>
void defineRho(SuperLattice2D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
               FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
               AnalyticalF2D<T,T>&                rho)
{

  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    defineRho<T,DESCRIPTOR>(sLattice.getExtendedBlockIndicator(iCloc), indicator->getExtendedBlockIndicatorF(iCloc),
                            bulkIndicator->getExtendedBlockIndicatorF(iCloc),
                            rho);
  }
}


////////// BlockLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineRho(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator, BlockIndicatorF2D<T>& bulkIndicator, AnalyticalF2D<T,T>& rho)
{

  if ( indicator.isEmpty() ) {
    return;
  }

  const Vector<int,2> min = indicator.getMin();
  const Vector<int,2> max = indicator.getMax();

  for (int iX = min[0]; iX <= max[0]; ++iX) {
    for (int iY = min[1]; iY <= max[1]; ++iY) {
      if (indicator(iX, iY)) {
        for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
          // Get direction
          int iXn = iX + descriptors::c<DESCRIPTOR >(iPop,0);
          int iYn = iY + descriptors::c<DESCRIPTOR >(iPop,1);
          if (bulkIndicator(iXn, iYn)) {
            T intersection[] = { T(), T() }; // coord. of intersection
            int opp = util::opposite<DESCRIPTOR >(iPop);
            //getBoundaryIntersection is defined in defineU2D.h/hh
            if (getBoundaryIntersection<T,DESCRIPTOR>(iX, iY, opp, intersection) ) {
              T rhoLocal[]= {T(1)};
              rho(rhoLocal,intersection);
              defineRho<T,DESCRIPTOR>(iX, iY, opp, rhoLocal[0]);
            }
          }
        }
      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
void defineRho(BlockLatticeStructure2D<T,DESCRIPTOR>& block, int iX, int iY, int iPop, const T rho)
{
  bool _output = false;
  OstreamManager clout(std::cout, "defineRho");
  block.getDynamics(iX, iY)->defineRho(iPop, rho);
  if (_output) {
    clout << "defineRho(" << iX << ", " << iY << " )" << std::endl;
  }
}


}//namespace olb


#endif
