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
void defineUBouzidi(SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry, int material,
                    AnalyticalF3D<T,T>& u, std::vector<int> bulkMaterials)
{
  defineUBouzidi<T,DESCRIPTOR>(sLattice, superGeometry.getMaterialIndicator(material), u, bulkMaterials);
}


template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice3D<T,DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                    AnalyticalF3D<T,T>& u, std::vector<int> bulkMaterials)
{
  defineUBouzidi<T,DESCRIPTOR>(sLattice, std::forward<decltype(boundaryIndicator)>(boundaryIndicator),
                               boundaryIndicator->getSuperGeometry().getMaterialIndicator(std::move(bulkMaterials)),
                               u);
}

template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice3D<T,DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                    FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator, AnalyticalF3D<T,T>& u)
{

  for (int iCloc = 0; iCloc < sLattice.getLoadBalancer().size(); ++iCloc) {
    defineUBouzidi<T,DESCRIPTOR>(sLattice.getExtendedBlockLattice(iCloc),boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
                                 bulkIndicator->getExtendedBlockIndicatorF(iCloc),u);
  }

}


////////// BlockLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, BlockIndicatorF3D<T>& bulkIndicator, AnalyticalF3D<T,T>& u)
{
  if ( !indicator.isEmpty() ) {
    const Vector<int,3> min = indicator.getMin();
    const Vector<int,3> max = indicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        for (int iZ = min[2]; iZ <= max[2]; ++iZ) {
          if (indicator(iX, iY, iZ)) {
            for (int q = 1; q < DESCRIPTOR::q ; ++q) {
              // Get direction
              const int iXn = iX + descriptors::c<DESCRIPTOR>(q,0);
              const int iYn = iY + descriptors::c<DESCRIPTOR>(q,1);
              const int iZn = iZ + descriptors::c<DESCRIPTOR>(q,2);
              if (bulkIndicator.getBlockGeometryStructure().isInside(iXn,iYn,iZn) &&
                  bulkIndicator(iXn,iYn,iZn)) {
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
        }
      }
    }
  }

}



template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, int iX, int iY, int iZ, int iPop, const T u[DESCRIPTOR::d])
{
  OstreamManager clout(std::cout, "defineUBouzidi");
  bool _output = false;
  _block.getDynamics(iX, iY, iZ)->defineU(iPop, u);
  if (_output) {
    clout << "defineUBouzidi(" << iX << ", " << iY << ", " << iZ << " )" << std::endl;
  }

}

template<typename T, typename DESCRIPTOR>
bool getBoundaryIntersection(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, int iX, int iY, int iZ, int iPop, T point[DESCRIPTOR::d])
{
  return _block.getDynamics(iX, iY, iZ)->getBoundaryIntersection(iPop, point);
}

template<typename T, typename DESCRIPTOR>
void setBoundaryIntersection(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, int iX, int iY, int iZ, int iPop, T distance)
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
