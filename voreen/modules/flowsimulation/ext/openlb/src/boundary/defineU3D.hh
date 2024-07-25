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

#ifndef DEFINE_U_3D_HH
#define DEFINE_U_3D_HH

#include "defineU3D.h"

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material,
                    AnalyticalF3D<T,T>& u, std::vector<int> bulkMaterials)
{
  defineUBouzidi<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material), u, bulkMaterials);
}


template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice<T,DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                    AnalyticalF3D<T,T>& u, std::vector<int> bulkMaterials)
{
  defineUBouzidi<T,DESCRIPTOR>(sLattice, std::forward<decltype(boundaryIndicator)>(boundaryIndicator),
                               boundaryIndicator->getSuperGeometry().getMaterialIndicator(std::move(bulkMaterials)),
                               u);
}

template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice<T,DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                    FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator, AnalyticalF3D<T,T>& u)
{

  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    defineUBouzidi<T,DESCRIPTOR>(sLattice.getBlock(iCloc),
                                 boundaryIndicator->getBlockIndicatorF(iCloc),
                                 bulkIndicator->getBlockIndicatorF(iCloc),u);
  }

}


////////// BlockLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, BlockIndicatorF3D<T>& bulkIndicator, AnalyticalF3D<T,T>& u)
{
  _block.forSpatialLocations([&](auto iX, auto iY, auto iZ) {
    if (indicator(iX,iY,iZ)) {
      for (int q = 1; q < DESCRIPTOR::q ; ++q) {
        // Get direction
        const int iXn = iX + descriptors::c<DESCRIPTOR>(q,0);
        const int iYn = iY + descriptors::c<DESCRIPTOR>(q,1);
        const int iZn = iZ + descriptors::c<DESCRIPTOR>(q,2);
        if (_block.isInside({iXn,iYn,iZn}) && bulkIndicator(iXn,iYn,iZn)) {
          T intersection[3] = { };
          const int opp = util::opposite<DESCRIPTOR>(q);
          if (getBoundaryIntersection<T,DESCRIPTOR>(_block, iX, iY, iZ, opp, intersection)) {
            T vel[3]= { };
            u(vel, intersection);
            defineUBouzidi<T,DESCRIPTOR>(_block, iX, iY, iZ, opp, vel);
          }
        }
      }
    }
  });
}



template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLattice<T,DESCRIPTOR>& _block, int iX, int iY, int iZ, int iPop, const T u[DESCRIPTOR::d])
{
  OstreamManager clout(std::cout, "defineUBouzidi");
  bool _output = false;
  static_cast<legacy::OffDynamics<T,DESCRIPTOR>*>(_block.getDynamics(iX, iY, iZ))->defineU(iPop, u);
  if (_output) {
    clout << "defineUBouzidi(" << iX << ", " << iY << ", " << iZ << " )" << std::endl;
  }

}

template<typename T, typename DESCRIPTOR>
bool getBoundaryIntersection(BlockLattice<T,DESCRIPTOR>& _block, int iX, int iY, int iZ, int iPop, T point[DESCRIPTOR::d])
{
  return static_cast<legacy::OffDynamics<T,DESCRIPTOR>*>(_block.getDynamics(iX, iY, iZ))->getBoundaryIntersection(iPop, point);
}

template<typename T, typename DESCRIPTOR>
void setBoundaryIntersection(BlockLattice<T,DESCRIPTOR>& _block, int iX, int iY, int iZ, int iPop, T distance)
{
  bool _output = false;
  OstreamManager clout(std::cout, "setBoundaryIntersection");
  _block.getDynamics(iX, iY, iZ)->setBoundaryIntersection(iPop, distance);
  if (_output) {
    clout << "setBoundaryIntersection(" << iX << ", " << iY << ", " << iZ << " )" << std::endl;
  }
}



}//namespace olb


#endif
