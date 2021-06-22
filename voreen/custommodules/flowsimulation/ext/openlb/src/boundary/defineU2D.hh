/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

#ifndef DEFINE_U_2D_HH
#define DEFINE_U_2D_HH

#include "defineU2D.h"

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry, int material,
                    AnalyticalF2D<T,T>& u,
                    std::vector<int> bulkMaterials)
{
  defineUBouzidi<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material),
                               superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
                               u);
}

template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice2D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                    FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                    AnalyticalF2D<T,T>& u)
{

  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    defineUBouzidi<T,DESCRIPTOR>(sLattice.getExtendedBlockIndicator(iCloc), indicator->getExtendedBlockIndicatorF(iCloc),
                                 bulkIndicator->getExtendedBlockIndicatorF(iCloc),
                                 u);
  }
}

////////// BlockLattice Domain  /////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator, BlockIndicatorF2D<T>& bulkIndicator, AnalyticalF2D<T,T>& u)
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
            if (getBoundaryIntersection<T,DESCRIPTOR>(iX, iY, opp, intersection) ) {
              T vel[]= {T(),T()};
              u(vel,intersection);
              defineUBouzidi<T,DESCRIPTOR>(iX, iY, opp, vel);
            }
          }
        }
      }
    }
  }

}


template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLatticeStructure2D<T,DESCRIPTOR>& block, int iX, int iY, int iPop, const T u[DESCRIPTOR::d])
{
  bool _output = false;
  OstreamManager clout(std::cout, "defineUBouzidi");
  block.getDynamics(iX, iY)->defineU(iPop, u);
  if (_output) {
    clout << "defineUBouzidi(" << iX << ", " << iY << " )" << std::endl;
  }
}



template<typename T, typename DESCRIPTOR>
bool getBoundaryIntersection(BlockLatticeStructure2D<T,DESCRIPTOR>& block, int iX, int iY, int iPop, T point[DESCRIPTOR::d])
{
  return block.getDynamics(iX, iY)->getBoundaryIntersection(iPop, point);
}

template<typename T, typename DESCRIPTOR>
void setBoundaryIntersection(BlockLatticeStructure2D<T,DESCRIPTOR>& block, int iX, int iY, int iPop, T distance)
{
  bool _output = false;
  OstreamManager clout(std::cout, "setBoundaryIntersection");
  block.getDynamics(iX, iY)->setBoundaryIntersection(iPop, distance);
  if (_output) {
    clout << "setBoundaryIntersection(" << iX << ", " << iY << " )" << std::endl;
  }
}


}//namespace olb


#endif
