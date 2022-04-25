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

#ifndef PARTICLE_DYNAMICS_UTILITIES_H
#define PARTICLE_DYNAMICS_UTILITIES_H

#include <cassert>

namespace olb {

namespace particles {

//Forward declaration of particle
template <typename T, typename DESCRIPTOR> class Particle;

namespace dynamics {

/// Common nuclear operations provided for convenience

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getPosition( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned D = PARTICLETYPE::d;
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>() );
  return position;
}

template<typename T, typename PARTICLETYPE>
Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> getAngle( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned Drot = utilities::dimensions::convert<PARTICLETYPE::d>::rotation;
  Vector<T,Drot> angle( particle.template getField<SURFACE,ANGLE>() );
  return angle;
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getVelocity( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned D = PARTICLETYPE::d;
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() );
  return velocity;
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getForce( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned D = PARTICLETYPE::d;
  Vector<T,D> force( particle.template getField<FORCING,FORCE>() );
  return force;
}

template<typename T, typename PARTICLETYPE>
Vector<T,2> getAdhesion( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  Vector<T,2> adhesion( particle.template getField<FORCING,ADHESION>() );
  return adhesion;
}

template<typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> getAcceleration( Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned D = PARTICLETYPE::d;
  Vector<T,D> acceleration;
  Vector<T,D> force( particle.template getField<FORCING,FORCE>() );
  T mass = particle.template getField<PHYSPROPERTIES,MASS>();
  for (unsigned iDim=0; iDim<D; ++iDim) {
    acceleration[iDim] = force[iDim] / mass;
  }
  return acceleration;
}

template<typename T, typename PARTICLETYPE>
Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> getAngAcceleration(
  Particle<T,PARTICLETYPE> particle )
{
  using namespace descriptors;
  const unsigned Drot = utilities::dimensions::convert<PARTICLETYPE::d>::rotation;
  Vector<T,Drot> angularAcceleration;
  Vector<T,Drot> torque( particle.template getField<FORCING,TORQUE>() );
  Vector<T,Drot> momentOfInertia( particle.template getField<PHYSPROPERTIES,MOFI>() );
  for (unsigned iRot=0; iRot<Drot; ++iRot) {
    angularAcceleration[iRot] = torque[iRot] / momentOfInertia[iRot];
  }
  return angularAcceleration;
}

template<typename T, typename PARTICLETYPE>
T getRadius( Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;
  T radius;
  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,SINDICATOR>() ) {
    radius = particle.template getField<SURFACE,SINDICATOR>()->getCircumRadius();
  }
  else if constexpr ( PARTICLETYPE::template providesNested<PHYSPROPERTIES,RADIUS>() ) {
    radius = particle.template getField<PHYSPROPERTIES,RADIUS>();
  }
  else {
    radius = 0.;
  }
  return radius;
}

template<typename T, typename PARTICLETYPE>
void updateRotationMatrix( Particle<T,PARTICLETYPE>& particle )
{
  using namespace descriptors;
  auto rotationMatrix = particles::dynamics::rotation_matrix<PARTICLETYPE::d,T>::calculate(
                          particle.template getField<SURFACE,ANGLE>() );
  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ROT_MATRIX>() ) {
    particle.template setField<SURFACE,ROT_MATRIX>( rotationMatrix );
  } else {
    std::cerr << "ERROR: The field ROT_MATRIX must be provided." << std::endl;
  }
}

template<typename PARTICLETYPE>
constexpr bool providesAngle()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<SURFACE,ANGLE>();
}

template<typename PARTICLETYPE>
constexpr bool providesActive()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<NUMERICPROPERTIES,ACTIVE>();
}

template<typename PARTICLETYPE>
constexpr bool providesRotationMatrix()
{
  using namespace descriptors;
  return PARTICLETYPE::template providesNested<SURFACE,ROT_MATRIX>();
}

/// Helper functions

template<typename T, typename PARTICLETYPE, typename F>
void doAtCubicBoundPenetration(
  Particle<T,PARTICLETYPE>& particle,
  Vector<T,PARTICLETYPE::d> domainMin,
  Vector<T,PARTICLETYPE::d> domainMax,
  F boundTreatment )
{
  auto radius = getRadius<T,PARTICLETYPE>( particle );
  auto position = getPosition<T,PARTICLETYPE>( particle );
  for (int iDim=0; iDim<PARTICLETYPE::d; ++iDim) {
    T distToLowerBound = (position[iDim] - radius) - domainMin[iDim];
    T distToUpperBound = domainMax[iDim] - (position[iDim] + radius);
    Vector<T,PARTICLETYPE::d> normal(0.);
    if ( distToLowerBound <= 0. ) {
      normal[iDim] = 1.;
      boundTreatment(iDim, normal, distToLowerBound);
    }
    if ( distToUpperBound <= 0. ) {
      normal[iDim] = -1.;
      boundTreatment(iDim, normal, distToUpperBound);
    }
  }
}

template<typename T, typename PARTICLETYPE>
void resetDirection( Particle<T,PARTICLETYPE>& particle,
                     Vector<T,PARTICLETYPE::d> positionPre, int iDir )
{
  using namespace descriptors;
  const unsigned D = PARTICLETYPE::d;
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>() );
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() );
  position[iDir] = positionPre[iDir];
  velocity[iDir] = 0.;
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<MOBILITY,VELOCITY>( velocity );
}
 
template<typename T, typename PARTICLETYPE>
void resetMovement( Particle<T,PARTICLETYPE>& particle,
                    Vector<T,PARTICLETYPE::d> positionPre,
                    Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> anglePre
                      = Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>(0.) )
{
  using namespace descriptors;
  const unsigned D = PARTICLETYPE::d;
  const unsigned Drot = utilities::dimensions::convert<PARTICLETYPE::d>::rotation;
  particle.template setField<GENERAL,POSITION>( positionPre );
  particle.template setField<MOBILITY,VELOCITY>( Vector<T,D>(0.) );
  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ANGLE>() ) {
    particle.template setField<SURFACE,ANGLE>( anglePre );
    particle.template setField<MOBILITY,ANG_VELOCITY>( Vector<T,Drot>(0.) );
  }
} 

template<typename T, typename PARTICLETYPE>
struct ParticleDynamicStateNoAngle{
  Vector<T,PARTICLETYPE::d> position;
  static const bool hasAngle = false;

  ParticleDynamicStateNoAngle<T, PARTICLETYPE>()
  {}

  ParticleDynamicStateNoAngle<T, PARTICLETYPE>(
    Vector<T,PARTICLETYPE::d> position )
    : position(position)
  {}
};

template<typename T, typename PARTICLETYPE>
struct ParticleDynamicStateAngle{
  Vector<T,PARTICLETYPE::d> position;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angle;  
  static const bool hasAngle = true;

  ParticleDynamicStateAngle<T,PARTICLETYPE>()
  {}

  ParticleDynamicStateAngle<T,PARTICLETYPE>(
    Vector<T,PARTICLETYPE::d> position )
    : position(position)
  {}

  ParticleDynamicStateAngle<T,PARTICLETYPE>(
    Vector<T,PARTICLETYPE::d> position,
    Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angle )
    : position(position), angle(angle)
  {}
};

template <typename T, typename PARTICLETYPE>
using DynState = std::conditional_t<
  PARTICLETYPE::template providesNested<descriptors::SURFACE,descriptors::ANGLE>(),
  ParticleDynamicStateAngle<T,PARTICLETYPE>,
  ParticleDynamicStateNoAngle<T,PARTICLETYPE>
>;

} //namespace dynamics
  
} //namespace particles

} //namespace olb

#endif
