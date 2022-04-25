/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007-2014 Mathias J. Krause
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
 * The description of a 2D super lattice -- header file.
 */

#ifndef SUPER_LATTICE_2D_H
#define SUPER_LATTICE_2D_H

#include <memory>
#include <vector>
#include <map>

#include "superLattice.hh"
#include "cellD.h"
#include "communication/superCommunicator.h"
#include "postProcessing.h"
#include "serializer.h"
#include "communication/superStructure.h"
#include "utilities/functorPtr.h"
#include "functors/analytical/analyticalF.h"

#include "core/olbDebug.h"

// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class CuboidGeometry2D;
template<typename T, typename DESCRIPTOR> class SuperLattice;
template<typename T> class LoadBalancer;
template<typename T, unsigned D> class SuperGeometry;
template<typename T, typename DESCRIPTOR> class SuperLatticeF2D;
template<typename T> class SuperStructure2D;
template<typename T> class SuperIndicatorF2D;

//TODO: 200116 preliminary version
template<typename T, typename DESCRIPTOR>
void setSuperExternalPSMParticleField( SuperGeometry<T,2>& sGeometry, int material, AnalyticalF2D<T,T>& velocity,
                                    T size,
                                    SuperLatticeF2D<T,DESCRIPTOR>& epsilon,
                                    SuperLattice<T, DESCRIPTOR>& sLattice );

//Geng2019
/* TODO: Change that it uses the new particle interface
template<typename T, typename DESCRIPTOR>
void setSuperZetaParticleField( SuperGeometry<T,2>& sGeometry, AnalyticalF<2,T,T>& velocity,
                                SmoothIndicatorF2D<T,T,true>& sIndicator,
                                SuperLattice<T, DESCRIPTOR>& sLattice );
*/

} // namespace olb

#endif
