/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Jan E. Marquardt, Mathias J. Krause
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

#ifndef PARTICLE_CREATOR_HELPER_FUNCTIONS_H
#define PARTICLE_CREATOR_HELPER_FUNCTIONS_H

namespace olb {

namespace particles {

namespace creators {

template<typename PARTICLETYPE, bool ROTATION_IS_OPTIONAL=false>
void checkForErrors()
{
  using namespace descriptors;

  static_assert(ROTATION_IS_OPTIONAL || PARTICLETYPE::template providesNested<SURFACE,ROT_MATRIX>(),
                "A rotation matrix is necessary but not provided.");
}

/// Set new particle for existing surface
template<typename T, typename PARTICLETYPE, bool ROTATION_IS_OPTIONAL=false>
void setResolvedObject(
  ParticleSystem<T,PARTICLETYPE>& particleSystem,
  std::size_t idxParticle, std::size_t idxSurface,
  const Vector<T,PARTICLETYPE::d>& position,
  T density,
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angle,
  const Vector<T,PARTICLETYPE::d>& velocity)
{
  particles::creators::checkForErrors<PARTICLETYPE,ROTATION_IS_OPTIONAL>();

  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  typedef SmoothIndicatorF<T,T,D,true> SIndicatorBaseType;

  //Retrieve smart pointer from particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::shared_ptr<SIndicatorBaseType>>>();
  auto sIndicatorSPtr = vectorOfIndicators.at( idxSurface );

  //Initialize fields (Seems to be necessary for clang-1000.10.44.4 but not for gcc)
  dynamics::initializeParticle<T,PARTICLETYPE>(particleSystem.get().data(), idxParticle);

  //Calculate moment of inertia and mass
  Vector<T,utilities::dimensions::convert<D>::rotation> momentOfInertia;
  T mass;
  if constexpr(D==3) {
    for (unsigned iD=0; iD<D; ++iD) {
      momentOfInertia[iD] = sIndicatorSPtr->calcMofiAndMass(density)[iD];
    }
    mass = sIndicatorSPtr->calcMofiAndMass(density)[D];
  }
  else {
    momentOfInertia[0] = sIndicatorSPtr->calcMofiAndMass(density)[0];
    mass = sIndicatorSPtr->calcMofiAndMass(density)[1];
  }

  //Set values
  auto particle = particleSystem.get( idxParticle );
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<SURFACE,SINDICATOR>( sIndicatorSPtr.get() );
  particle.template setField<PHYSPROPERTIES,MASS>( mass );
  particle.template setField<MOBILITY,VELOCITY>( velocity );
  particle.template setField<PHYSPROPERTIES,MOFI>( utilities::dimensions::convert<D>::serialize_rotation(momentOfInertia) );

  angle *= (M_PI/180.);
  particle.template setField<SURFACE,ANGLE>( utilities::dimensions::convert<D>::serialize_rotation(angle) );
  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ROT_MATRIX>() ) {
    const Vector<T,utilities::dimensions::convert<D>::matrix> rotationMatrix = dynamics::rotation_matrix<D,T>::calculate(
          utilities::dimensions::convert<D>::serialize_rotation(angle) );
    particle.template setField<SURFACE,ROT_MATRIX>( rotationMatrix );
  }
}

/// Add resolved object as new particle with existing surface
template<typename T, typename PARTICLETYPE>
void addResolvedObject(
  ParticleSystem<T,PARTICLETYPE>& particleSystem, std::size_t idxSurface,
  const Vector<T,PARTICLETYPE::d>& position, T density=0.,
  const Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>& angle = Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> (0.),
  const Vector<T,PARTICLETYPE::d>& velocity = Vector<T,PARTICLETYPE::d> (0.))
{
  //Retrieve new index
  std::size_t idxParticle = particleSystem.size();

  //Initialize particle address
  particleSystem.extend();

  /// Set resolved object 3D at given index
  setResolvedObject( particleSystem, idxParticle, idxSurface, position, density, angle, velocity );
}

} //namespace creators

} //namespace particles

} //namespace olb

#endif