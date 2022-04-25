/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
 *                2008-2020 Mathias Krause
 *                2020      Adrian Kummerlaender
 *                2021      Nicolas Hafen
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

#ifndef BLOCK_LATTICE_INTERACTION_H
#define BLOCK_LATTICE_INTERACTION_H

namespace olb {

namespace particles {

//Forward declaration
template<typename T, typename PARTICLETYPE> class ParticleSystem;
template<typename T, typename PARTICLETYPE> class Particle;


//Get block particle intersection
template <typename T, unsigned D>
bool getBlockParticleIntersection( BlockGeometry<T,D>& blockGeometry,
                                   T invDeltaX,
                                   LatticeR<D>& start, LatticeR<D>& end,
                                   Vector<T,D> position, T circumRadius );

//Check whether smoothIndicator is out of geometry
template<typename T, unsigned D>
void checkSmoothIndicatorOutOfGeometry( bool& outOfGeometry, Vector<T,D>& ghostPos,
                                        const PhysR<T,D>& cellMin,
                                        const PhysR<T,D>& cellMax,
                                        const Vector<T,D>& position, T circumRadius,
                                        const Vector<bool,D>& periodic);

//Iterate over spacial locations in block particle intersection with lambda expression f
template <typename T, typename DESCRIPTOR, typename F>
void forSpatialLocationsInBlockParticleIntersection( BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
    BlockLattice<T,DESCRIPTOR>& blockLattice,
    Vector<T,DESCRIPTOR::d> position,
    T circumRadius, F f );

//Set block particle field
template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void setBlockParticleField( BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                            AnalyticalF<DESCRIPTOR::d,T,T>& eccentricVelocity,
                            BlockLattice<T,DESCRIPTOR>& blockLattice,
                            Particle<T,PARTICLETYPE>& particle );

} //namespace particles

} //namespace olb

#endif
