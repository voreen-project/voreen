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


/* This file contains functions used for the calculation of particle related dynamics.
 *
*/

#ifndef PARTICLE_MOTION_FUNCTIONS_H
#define PARTICLE_MOTION_FUNCTIONS_H


namespace olb {

namespace particles {

namespace dynamics {

/////////// Basic Newtons Movement Functions /////////

// Velocity verlet integration according to
// Verlet1967 (https://doi.org/10.1103/PhysRev.159.98)
// Swope1982 (https://doi.org/10.1063/1.442716)
template<typename T, typename PARTICLETYPE>
void velocityVerletIntegration( Particle<T,PARTICLETYPE>& particle, T delTime,
                                Vector<T,PARTICLETYPE::d> acceleration,
                                Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angularAcceleration
                                = Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>(0.) )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  static_assert(PARTICLETYPE::template providesNested<GENERAL,POSITION>(), "Field POSITION has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,VELOCITY>(), "Field VELOCITY has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,ACCELERATION_STRD>(), "Field ACCELERATION_STRD has to be provided");
  //Calculate cartesian directions
  T delTime2 = delTime*delTime;
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>()
                        + particle.template getField<MOBILITY,VELOCITY>() * delTime
                        + (0.5 * particle.template getField<MOBILITY,ACCELERATION_STRD>() * delTime2) );
  Vector<T,D> avgAcceleration( .5*(particle.template getField<MOBILITY,ACCELERATION_STRD>() + acceleration) );
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() + avgAcceleration * delTime );
  //Update values
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<MOBILITY,VELOCITY>( velocity );
  particle.template setField<MOBILITY,ACCELERATION_STRD>( acceleration );

  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ANGLE>() ) {
    const unsigned Drot = utilities::dimensions::convert<D>::rotation;
    //Calculate rotation
    Vector<T,Drot> angle( particle.template getField<SURFACE,ANGLE>()
                          + particle.template getField<MOBILITY,ANG_VELOCITY>() * delTime
                          + (0.5 * particle.template getField<MOBILITY,ANG_ACC_STRD>() * delTime2) );
    Vector<T,Drot> avgAngularAcceleration( .5*(particle.template getField<MOBILITY,ANG_ACC_STRD>() + angularAcceleration) );
    Vector<T,Drot> angularVelocity( particle.template getField<MOBILITY,ANG_VELOCITY>() + avgAngularAcceleration * delTime );
    //Treat full rotation
    for (unsigned iRot=0; iRot<Drot; ++iRot) {
      angle[iRot] = std::fmod( angle[iRot], 2.*M_PI );
    }
    //Update values
    particle.template setField<SURFACE,ANGLE>(
      utilities::dimensions::convert<D>::serialize_rotation(angle) );
    particle.template setField<MOBILITY,ANG_VELOCITY>(
      utilities::dimensions::convert<D>::serialize_rotation(angularVelocity) );
    particle.template setField<MOBILITY,ANG_ACC_STRD>(
      utilities::dimensions::convert<D>::serialize_rotation(angularAcceleration) );

  }
}


/// Euler integration
template<typename T, typename PARTICLETYPE>
void eulerIntegration( Particle<T,PARTICLETYPE>& particle, T delTime,
                       Vector<T,PARTICLETYPE::d> acceleration,
                       Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angularAcceleration
                       = Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>(0.) )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  static_assert(PARTICLETYPE::template providesNested<GENERAL,POSITION>(), "Field POSITION has to be provided");
  static_assert(PARTICLETYPE::template providesNested<MOBILITY,VELOCITY>(), "Field VELOCITY has to be provided");
  //Calculate cartesian directions
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() + acceleration * delTime );
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>() + velocity * delTime );
  //Update values
  particle.template setField<GENERAL,POSITION>( position );
  particle.template setField<MOBILITY,VELOCITY>( velocity );

  if constexpr ( PARTICLETYPE::template providesNested<SURFACE,ANGLE>() ) {
    const unsigned Drot = utilities::dimensions::convert<D>::rotation;
    //Calculate rotational directions
    Vector<T,Drot> angularVelocity( particle.getMobility().getAngularVelocity()
                                    + angularAcceleration * delTime );
    Vector<T,Drot> angle( particle.getSurface().getAngle() + angularVelocity * delTime );
    //Update values
    particle.template setField<SURFACE,ANGLE>( angle );
    particle.template setField<MOBILITY,ANG_VELOCITY>( angularVelocity );
  }
}


} //namespace dynamics

} //namespace particles

} //namespace olb


#endif
