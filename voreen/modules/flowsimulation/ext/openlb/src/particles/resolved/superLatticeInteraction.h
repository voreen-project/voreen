/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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


#ifndef SUPER_LATTICE_INTERACTION_H
#define SUPER_LATTICE_INTERACTION_H


// All OpenLB code is contained in this namespace.
namespace olb {

namespace particles {


/// Set particle field
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void setSuperParticleField( SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                            AnalyticalF<DESCRIPTOR::d,T,T>& velocity,
                            SuperLattice<T, DESCRIPTOR>& sLattice,
                            Particle<T,PARTICLETYPE>& particle );

} // namespace particles

} // namespace olb

#endif
