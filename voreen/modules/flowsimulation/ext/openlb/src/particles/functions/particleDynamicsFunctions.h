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

#ifndef PARTICLE_DYNAMICS_FUNCTIONS_H
#define PARTICLE_DYNAMICS_FUNCTIONS_H


#include <cassert>

namespace olb {

namespace particles {

namespace dynamics {

////////////// Dimension sensitive Functions ////////////

//Calculate rotation matrix
template<unsigned D, typename T>
struct rotation_matrix;

template<typename T>
struct rotation_matrix<2,T> {
  static constexpr Vector<T,4> calculate( Vector<T,1> angle )
  {
    Vector<T,4> rotationMatrix;

    T cos = util::cos(angle[0]);
    T sin = util::sin(angle[0]);

    rotationMatrix[0] = cos;
    rotationMatrix[1] = sin;
    rotationMatrix[2] = -sin;
    rotationMatrix[3] = cos;

    return rotationMatrix;
  }

  static constexpr Vector<T,4> calculateInverse( Vector<T,1> angle )
  {
    return calculate(-1 * angle);
  }
};

template<typename T>
struct rotation_matrix<3,T> {
  static constexpr Vector<T,9> calculate( Vector<T,3> angle )
  {
    Vector<T,9> rotationMatrix;

    T cos0 = util::cos(angle[0]);
    T cos1 = util::cos(angle[1]);
    T cos2 = util::cos(angle[2]);
    T sin0 = util::sin(angle[0]);
    T sin1 = util::sin(angle[1]);
    T sin2 = util::sin(angle[2]);

    rotationMatrix[0] = cos1*cos2;
    rotationMatrix[1] = sin0*sin1*cos2 - cos0*sin2;
    rotationMatrix[2] = cos0*sin1*cos2 + sin0*sin2;
    rotationMatrix[3] = cos1*sin2;
    rotationMatrix[4] = sin0*sin1*sin2 + cos0*cos2;
    rotationMatrix[5] = cos0*sin1*sin2 - sin0*cos2;
    rotationMatrix[6] = -sin1;
    rotationMatrix[7] = sin0*cos1;
    rotationMatrix[8] = cos0*cos1;

    return rotationMatrix;
  }

  static constexpr Vector<T,9> calculateInverse( Vector<T,3> angle )
  {
    return calculate(-1 * angle);
  }
};

//Calculate torque from force and lever
template<unsigned D, typename T>
struct torque_from_force;

template<typename T>
struct torque_from_force<2,T> {
  static constexpr Vector<T,1> calculate(
    Vector<T,2> force,
    PhysR<T,2> lever )
  {
    return Vector<T,1>( lever[0]*force[1]-lever[1]*force[0] );
  }
};

template<typename T>
struct torque_from_force<3,T> {
  static constexpr Vector<T,3> calculate(
    Vector<T,3> force,
    PhysR<T,3> lever )
  {
    return Vector<T,3>( lever[1]*force[2]-lever[2]*force[1],
                        lever[2]*force[0]-lever[0]*force[2],
                        lever[0]*force[1]-lever[1]*force[0]
                      );
  }
};

//Calculate local velocity
template <typename T, typename PARTICLETYPE>
Vector<T,PARTICLETYPE::d> calculateLocalVelocity(Particle<T,PARTICLETYPE>& particle, const PhysR<T,PARTICLETYPE::d>& input)
{
  using namespace descriptors;

  auto velocity = particle.template getField<MOBILITY,VELOCITY>();
  auto position = particle.template getField<GENERAL,POSITION>();
  auto angVelocity = Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>(
    particle.template getField<MOBILITY,ANG_VELOCITY>() );

  return util::calculateLocalVelocity(position, velocity, angVelocity, input);
}


////////////// Force and Utility Functions ////////////

/// Unserialize force field provieded by force integration functor (e.g. momentumExchange)
template<typename T, typename PARTICLETYPE>
void unserializeForceTorqueVoxels( Vector<T,PARTICLETYPE::d>& force,
                                   Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation>& torque,
                                   T serializedFunctorForceField[], int iP )
{
  constexpr unsigned D = PARTICLETYPE::d;
  const unsigned Drot = utilities::dimensions::convert<D>::rotation;
  const int serialSize = D+Drot+1;

  const int idxForce = iP*(serialSize);
  const int idxTorque = idxForce+D;
  const int idxTorqueEnd = idxTorque+Drot;

  //Get force
  force = std::vector<T>(  serializedFunctorForceField+idxForce,
                           serializedFunctorForceField+idxTorque);
  //Get torque
  torque = std::vector<T>( serializedFunctorForceField+idxTorque,
                           serializedFunctorForceField+idxTorqueEnd );
}


/// Unserialize force field provieded by force integration functor (e.g. stokesDragForce)
template<typename T, typename PARTICLETYPE>
void unserializeForce( Vector<T,PARTICLETYPE::d>& force,
                       T serializedFunctorForceField[], int iP )
{
  const unsigned D = PARTICLETYPE::d;
  const int serialSize = D;

  const int idxForce = iP*(serialSize);
  const int idxForceEnd = idxForce+D;

  //Get force
  force = std::vector<T>(  serializedFunctorForceField+idxForce,
                           serializedFunctorForceField+idxForceEnd);
}

/// Apply boundary force provided by force functor to the particle center as torque and force
template<typename T, typename PARTICLETYPE>
void applyParticleForce( SuperF<PARTICLETYPE::d,T,T>& forceF, ParticleSystem<T,PARTICLETYPE>& particleSystem, std::size_t iP0=0 )
{
  constexpr unsigned D = PARTICLETYPE::d;
  const unsigned Drot = utilities::dimensions::convert<D>::rotation;

  //Retrieve boundary force field
  int input[1];
  T serializedFunctorForceField[forceF.getTargetDim()];
  forceF(serializedFunctorForceField, input);

  //Loop over particles and apply individual force and torque contribution
  for (std::size_t iP=iP0; iP<particleSystem.size(); iP++) {
    auto particle = particleSystem.get(iP);

    Vector<T,D> force;
    std::size_t iPeval = iP-iP0; //Shift, if iP!=0
    if constexpr ( PARTICLETYPE::template providesNested<descriptors::FORCING,descriptors::TORQUE>() ) {
      Vector<T,Drot> torque;
      unserializeForceTorqueVoxels<T,PARTICLETYPE>( force, torque, serializedFunctorForceField, iPeval );
      particle.template setField<descriptors::FORCING,descriptors::TORQUE>(
        utilities::dimensions::convert<D>::serialize_rotation(torque) );
    }
    else {
      unserializeForce<T,PARTICLETYPE>( force, serializedFunctorForceField, iPeval );
    }
    particle.template setField<descriptors::FORCING,descriptors::FORCE>( force );
  }
}


/// Initialize all fields in particle (necessary for clang)
template<typename T, typename PARTICLETYPE>
void initializeParticle( DynamicFieldGroupsD<T, typename PARTICLETYPE::fieldList>& dynamicFieldGroups, std::size_t iP )
{
  //Define init lambda expression
  typedef std::function<bool(const std::type_info&,int,std::string)> FunctionType;
  FunctionType initFunction = [](const std::type_info& typeInfo, int fieldSize, std::string fieldContentStr) {
    return true; //resetField=true
  };
  //Call recursive field traversal function with lambda expression
  descriptors::access_field_content<FunctionType,T,PARTICLETYPE,typename PARTICLETYPE::fieldList>::fieldsL2(
    initFunction, dynamicFieldGroups, iP );
}

//Apply stokes drag force
template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void applyStokesForce( ParticleSystem<T,PARTICLETYPE>& particleSystem,
                       SuperLattice<T, DESCRIPTOR>& sLattice,
                       SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                       UnitConverter<T,DESCRIPTOR> const& converter,
                       SuperLatticeInterpPhysVelocity<T, DESCRIPTOR>& interpPhysVelF )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;

  //Calculate general constants
  T dTinv = 1./converter.getPhysDeltaT();
  T C1 = 6. * M_PI * converter.getPhysViscosity()
         * converter.getPhysDensity() * converter.getConversionFactorTime();

  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);

    //Retrieve particle quantities
    Vector<T,D> position = particle.template getField<GENERAL,POSITION>();
    Vector<T,D> velocity = particle.template getField<MOBILITY,VELOCITY>();
    T radius = particle.template getField<PHYSPROPERTIES,RADIUS>();
    T mass = particle.template getField<PHYSPROPERTIES,MASS>();
    T* positionArray = position.data();

    //Calculate particle coefficiants
    T c = C1 * radius * 1./mass;
    T C2 = 1. / (1. + c);

    //Create fluid velocity array
    T fluidVelArray[D] = {0.};

    //Loop over iCs
    int maxC = sLattice.getLoadBalancer().size();
    for (int iC = 0; iC < maxC; iC++) {
      int globiC = sLattice.getLoadBalancer().glob(iC);

      //Retrieve block bounds
      LatticeR<D> extend = sGeometry.getBlockGeometry(iC).getExtent();
      PhysR<T,D> physExtend = converter.getPhysDeltaX()*extend;
      PhysR<T,D> origin = sGeometry.getBlockGeometry(iC).getOrigin();
      PhysR<T,D> end = origin + physExtend;

      //Check whether inside cuboid
      bool inside = true;
      for (unsigned iDim=0; iDim<D; ++iDim) {
        inside = inside
                 && (positionArray[iDim] > origin[iDim])
                 && (positionArray[iDim] < end[iDim]);
      }
      if (inside) {
        interpPhysVelF(fluidVelArray, positionArray, globiC);
      }
    }

    //Synchronize, and calculate stokes force
    Vector<T,3> force;
    for (int iDim = 0; iDim < PARTICLETYPE::d; ++iDim) {
      //Communicate
#ifdef PARALLEL_MODE_MPI
      singleton::mpi().reduceAndBcast(fluidVelArray[iDim], MPI_SUM);
#endif
      //Calculate force
      force[iDim] = mass * dTinv
                    * ((c * fluidVelArray[iDim] + velocity[iDim]) * C2 - velocity[iDim]);
    }
    particle.template setField<FORCING,FORCE>(force);
  }
}

} //namespace dynamics

} //namespace particles

} //namespace olb


#endif
