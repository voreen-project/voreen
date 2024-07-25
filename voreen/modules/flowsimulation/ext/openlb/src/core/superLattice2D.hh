/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause
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
 * The description of a 2D super lattice -- generic implementation.
 */


#ifndef SUPER_LATTICE_2D_HH
#define SUPER_LATTICE_2D_HH

#include <limits>
#include <numeric>

#include "communication/mpiManager.h"
#include "cell.h"
#include "geometry/cuboidGeometry2D.h"
#include "geometry/superGeometry.h"
#include "communication/loadBalancer.h"
#include "superLattice2D.h"
#include "io/base64.h"
#include "functors/lattice/superBaseF2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "io/serializerIO.h"

namespace olb {

//TODO: 200116 preliminary version
template<typename T, typename DESCRIPTOR>
void setSuperExternalPSMParticleField( SuperGeometry<T,2>& sGeometry, int material, AnalyticalF2D<T,T>& velocity, 
                                    T size,
                                    SuperLatticeF2D<T,DESCRIPTOR>& epsilon,
                                    SuperLattice<T, DESCRIPTOR>& sLattice )
{
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator = sGeometry.getMaterialIndicator(material);
  const int overlap = indicator->getSuperGeometry().getOverlap();
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) { 
    int globIC = sLattice.getLoadBalancer().glob(iC);
    BlockIndicatorF2D<T>& blockIndicator = indicator->getBlockIndicatorF(iC);
    setBlockExternalPSMParticleField( sGeometry.getBlockGeometry(iC), velocity, size, epsilon,
    sLattice.getBlock(iC), blockIndicator, globIC, overlap);
  }
}



//Set Zeta-Field (Geng2019)
/* TODO: Change that it uses the new particle interface
template<typename T, typename DESCRIPTOR>
void setSuperZetaParticleField( SuperGeometry<T,2>& sGeometry, AnalyticalF<2,T,T>& velocity,
                                SmoothIndicatorF2D<T,T,true>& sIndicator,
                                SuperLattice<T, DESCRIPTOR>& sLattice )
{
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setBlockZetaParticleField( sGeometry.getBlockGeometry(iC), velocity, sIndicator,
                               sLattice.getBlock(iC) );
  }
}
*/

} // namespace olb

#endif
