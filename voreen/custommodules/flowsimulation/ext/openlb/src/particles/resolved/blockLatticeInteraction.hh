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

#ifndef BLOCK_LATTICE_INTERACTION_HH
#define BLOCK_LATTICE_INTERACTION_HH


namespace olb {

namespace particles {

template <typename T, unsigned D>
bool getBlockParticleIntersection( BlockGeometry<T,D>& blockGeometry,
                                   T invDeltaX,
                                   LatticeR<D>& start, LatticeR<D>& end,
                                   PhysR<T,D> position, T circumRadius )
{

  /// Setting block bounds excluding(!) padding
  LatticeR<D> blockMin( 0. );
  LatticeR<D> blockMax( blockGeometry.getExtent()-1 );

  /// Calculate relative start/end in block domain
  PhysR<T,D> particleMin( position - circumRadius );
  PhysR<T,D> relStart( particleMin - blockGeometry.getOrigin() );
  PhysR<T,D> particleMax( position + circumRadius );
  PhysR<T,D> relEnd( particleMax - blockGeometry.getOrigin() );

  /// Set latticeR start/end and check validity
  bool validRange = true;
  for (unsigned iDim=0; iDim<D; ++iDim) {
    start[iDim] = util::max( util::floor( invDeltaX * relStart[iDim]), blockMin[iDim] );
    end[iDim] = util::min( util::ceil( invDeltaX * relEnd[iDim] ), blockMax[iDim] );
    T intersectionRange = end[iDim] - start[iDim];
    validRange = validRange && (intersectionRange >= 0.);
  }

  //Interpret valid range
  if (!validRange) {
    return false;
  }
  else {
    return true;
  }
}



template<typename T, unsigned D>
void checkSmoothIndicatorOutOfGeometry( bool& outOfGeometry, Vector<T,D>& ghostPos,
                                        const PhysR<T,D>& cellMin,
                                        const PhysR<T,D>& cellMax,
                                        const Vector<T,D>& position, T circumRadius,
                                        const Vector<bool,D>& periodic)
{
  Vector<bool,D>dir = Vector<bool,D> (false);
  outOfGeometry = false;
  for (unsigned i=0; i<D; i++) {
    T const particleMin = position[i] - circumRadius;
    T const particleMax = position[i] + circumRadius;
    if (particleMin < cellMin[i] && periodic[i]) {
      outOfGeometry = true;
      dir[i] = true;
      ghostPos[i] = cellMax[i] - (cellMin[i] - position[i]);
    }
    else if (particleMax > cellMax[i] && periodic[i]) {
      outOfGeometry = true;
      dir[i] = true;
      ghostPos[i] = cellMin[i] + (position[i] - cellMax[i]);
    }
    if (!dir[i]) {
      ghostPos[i] = position[i];
    }
  }
}

//Interation for block particle intersection for provided lambda function F
template <typename T, typename DESCRIPTOR, typename F>
void forSpatialLocationsInBlockParticleIntersection( BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
    BlockLattice<T,DESCRIPTOR>& blockLattice,
    int padding,
    Vector<T,DESCRIPTOR::d> position,
    T circumRadius, F f)
{
  constexpr unsigned D = DESCRIPTOR::d;
  LatticeR<D> start;
  LatticeR<D> end;
  if (getBlockParticleIntersection( blockGeometry, T{1}/blockGeometry.getDeltaR(),
                                    start, end, position, circumRadius )) {
    //Extend range by padding if desired
    start -= LatticeR<D>(padding);
    end += LatticeR<D>(padding);
    //For all cells in subdomain defined by start/end
    blockLattice.forSpatialLocations(start, end, [&](LatticeR<D> latticeR) {
      f(latticeR);
    });

  }
}



template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void setBlockParticleField( BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                            AnalyticalF<DESCRIPTOR::d,T,T>& eccentricVelocity,
                            BlockLattice<T,DESCRIPTOR>& blockLattice,
                            Particle<T,PARTICLETYPE>& particle )
{
  constexpr unsigned D = DESCRIPTOR::d;

  using namespace descriptors;
  auto circumRadius = particle.template getField<SURFACE,SINDICATOR>()->getCircumRadius();
  auto position = particle.template getField<GENERAL,POSITION>();

  //For all cells in block particle intersection
  int padding = blockGeometry.getPadding();
  forSpatialLocationsInBlockParticleIntersection( blockGeometry, blockLattice, padding,
      position, circumRadius,
  [&](LatticeR<D> latticeR) {
    //Get phys position
    T physR[D] = { };
    blockGeometry.getPhysR(physR, latticeR);
    //Retrieve porosity
    FieldD<T,DESCRIPTOR,descriptors::POROSITY> porosity { };
    particles::resolved::evalPorosity(porosity.data(), physR, particle);
    //For cells containing particle bits
    if (!util::nearZero(porosity[0])) {
      //Retrieve local velocity
      FieldD<T,DESCRIPTOR,descriptors::VELOCITY> velocity { };
      eccentricVelocity( velocity.data(), physR );
      //Apply weighting by porosity
      for (int iDim=0; iDim!=D; ++iDim) {
        velocity[iDim]*=porosity[0];
      }
      //Calculate solid volume fraction
      T solidVolumeFraction = 1. - porosity[0];
      //Write to fields
      {
        auto cell = blockLattice.get(latticeR);
        cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()   += velocity;
        cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>() += porosity;
        cell.template getFieldPointer<descriptors::POROSITY>() *= solidVolumeFraction;
      }
    }
  });
}

} //namespace particles

}  // namespace olb

#endif
