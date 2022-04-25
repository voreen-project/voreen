/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Robin Trunk, Mathias J. Kraus
 *                2021 Nicolas Hafen
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

#ifndef LATTICE_MOMENTUM_EXCHANGE_FORCE_HH
#define LATTICE_MOMENTUM_EXCHANGE_FORCE_HH

#include "particles/resolved/blockLatticeInteraction.hh"

namespace olb {


template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
SuperLatticeMomentumExchangeForce<T,DESCRIPTOR,PARTICLETYPE>::SuperLatticeMomentumExchangeForce(
  SuperLattice<T,DESCRIPTOR>& sLattice,
  SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter,
  Vector<bool,DESCRIPTOR::d> periodic, std::size_t iP0 )
  : SuperLatticePhysF<T,DESCRIPTOR>(sLattice,converter,
                                    (DESCRIPTOR::d+utilities::dimensions::convert<DESCRIPTOR::d>::rotation+1)*(particleSystem.size()-iP0))
{
  this->getName() = "physMomentumExchangeForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  PhysR<T,DESCRIPTOR::d> min = superGeometry.getStatistics().getMinPhysR( 1 ); //TODO: disable hardcoded 1
  PhysR<T,DESCRIPTOR::d> max = superGeometry.getStatistics().getMaxPhysR( 1 ); //TODO: disable hardcoded 1
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeMomentumExchangeForce<T,DESCRIPTOR, PARTICLETYPE>(
                                  this->_sLattice.getBlock(iC),
                                  superGeometry.getBlockGeometry(iC),
                                  particleSystem, converter, min, max, periodic, iP0));
  }
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool SuperLatticeMomentumExchangeForce<T,DESCRIPTOR,PARTICLETYPE>::operator() (T output[],
    const int input[])
{
  for (int iS=0; iS<this->getTargetDim(); iS++) {
    output[iS] = 0.;
  }
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    int globiC = this->_sLattice.getLoadBalancer().glob(iC);
    if ( this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank() ) {
      this->getBlockF(iC)(output,&input[1]);
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int iS = 0; iS < this->getTargetDim(); ++iS) {
    singleton::mpi().reduceAndBcast(output[iS], MPI_SUM);
  }
#endif
  return true;

}


template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
BlockLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>::BlockLatticeMomentumExchangeForce(
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter,
  PhysR<T,DESCRIPTOR::d>cellMin, PhysR<T,DESCRIPTOR::d> cellMax,
  Vector<bool,DESCRIPTOR::d> periodic, std::size_t iP0 )
  : BlockLatticePhysF<T,DESCRIPTOR>(blockLattice, converter, particleSystem.size()),
    _blockGeometry(blockGeometry), _blockLattice(blockLattice),
    _particleSystem(particleSystem),
    _cellMin(cellMin), _cellMax(cellMax), _periodic(periodic), _iP0(iP0)
{
  this->getName() = "physMomentumExchangeForce";
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void BlockLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>::evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP)
{
  constexpr unsigned D = DESCRIPTOR::d;
  const unsigned Drot = utilities::dimensions::convert<D>::rotation;
  const int serialSize = D+Drot+1;

  using namespace descriptors;
  auto position = particle.template getField<GENERAL,POSITION>();
  auto circumRadius = particle.template getField<SURFACE,SINDICATOR>()->getCircumRadius();

  //For all cells in block particle intersection
  int numVoxels = 0;
  int padding = 0;
  std::size_t iPeval = iP-_iP0; //Shifted output index, if iP!=0
  particles::forSpatialLocationsInBlockParticleIntersection( _blockGeometry, _blockLattice, padding,
      position, circumRadius,
  [&](LatticeR<D> latticeRinner) {

    T tmpForce[D] = {0.};
    PhysR<T,DESCRIPTOR::d> lever(0.);
    if ( particles::resolved::momentumExchangeAtSurfaceLocation(tmpForce,lever,latticeRinner,
         this->_blockGeometry, this->_blockLattice, this->_converter, particle) ){

      // count cells of considered object
      numVoxels++;

      //Calculate torque
      auto torque = particles::dynamics::torque_from_force<D,T>::calculate( Vector<T,D>(tmpForce), lever );

      //Add force and torque to output
      for (unsigned iDim=0; iDim<D; iDim++) {
        output[iDim+iPeval*serialSize] += tmpForce[iDim];
      }
      for (unsigned iRot=0; iRot<Drot; iRot++) {
        output[(iRot+D)+iPeval*serialSize] += torque[iRot];
      }
    }
  });

  output[D+Drot+iPeval*serialSize] = numVoxels;
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool BlockLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE>::operator()(T output[], const int input[])
{
  constexpr unsigned D = DESCRIPTOR::d;

  using namespace descriptors;

  // iterate over all particles in _particleSystem
  for (std::size_t iP=_iP0; iP!=_particleSystem.size(); iP++) {

    auto particle = _particleSystem.get(iP);
    auto circumRadius = particle.template getField<SURFACE,SINDICATOR>()->getCircumRadius();
    auto position = particle.template getField<GENERAL,POSITION>();

    //Evaluate periodicity
    bool anyPeriodic = false;
    for (unsigned iDim=0; iDim<D; iDim++) {
      anyPeriodic = anyPeriodic || _periodic[iDim];
    }

    //If any direction requires periodic treatment
    if (anyPeriodic) {
      bool outOfGeometry = false;
      PhysR<T,D> ghostPos = PhysR<T,D> (0.);
      particles::checkSmoothIndicatorOutOfGeometry(outOfGeometry, ghostPos,
        _cellMin, _cellMax, position, circumRadius, _periodic);
      PhysR<T,D>particlePosition = position;
      // Sets the particle to the other domainside if it leaves the domain
      bool outOfDomain = false;
      for (unsigned iDim=0; iDim<D; iDim++) {
        outOfDomain = outOfDomain || particlePosition[iDim] <_cellMin[iDim] || particlePosition[iDim]>_cellMax[iDim];
      }

      if ( outOfDomain) {
        particle.template setField<GENERAL,POSITION>(ghostPos);
        ghostPos = particlePosition;
        particlePosition = particle.template getField<GENERAL,POSITION>();
      }

      if (!outOfGeometry) {
        evaluate(output, particle, iP);
      }
      else {
        // Calculate force for the ghost particle
        particle.template setField<GENERAL,POSITION>(ghostPos);
        evaluate(output, particle, iP);
        // Calculate force of actual particle
        particle.template setField<GENERAL,POSITION>(particlePosition);
        evaluate(output, particle, iP);
      }
    }
    else {
      evaluate(output, particle, iP);
    }
  }
  return true;
}




template<typename T,typename DESCRIPTOR,typename PARTICLETYPE>
SuperLatticeMomentumExchangeForceLocal<T,DESCRIPTOR,PARTICLETYPE>::SuperLatticeMomentumExchangeForceLocal(
  SuperLattice<T,DESCRIPTOR>& sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem )
  : SuperLatticePhysF<T,DESCRIPTOR>(sLattice, converter, DESCRIPTOR::d)
{
  this->getName() = "localMomentumExchange";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeMomentumExchangeForceLocal<T,DESCRIPTOR,PARTICLETYPE>(
      this->_sLattice.getBlock(iC),
      superGeometry.getBlockGeometry(iC),
      particleSystem,
      this->_converter));
  }
}


template <typename T, typename DESCRIPTOR,typename PARTICLETYPE>
BlockLatticeMomentumExchangeForceLocal<T,DESCRIPTOR,PARTICLETYPE>::BlockLatticeMomentumExchangeForceLocal(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF<T,DESCRIPTOR>(blockLattice,converter,DESCRIPTOR::d),
    _blockLattice(blockLattice),
    _blockGeometry(blockGeometry),
    _particleSystem(particleSystem)
{
  this->getName() = "localMomentumExchange";
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool BlockLatticeMomentumExchangeForceLocal<T,DESCRIPTOR,PARTICLETYPE>::operator() (T output[], const int input[])
{
  using namespace descriptors;

  for (unsigned iDim=0; iDim<DESCRIPTOR::d; iDim++) {
    output[iDim] = 0.;
  }

  // iterate over all particles in _particleSystem
  int _iP0 = 0;
  for (std::size_t iP=_iP0; iP!=_particleSystem.size(); iP++) {
    auto particle = _particleSystem.get(iP);
    auto circumRadius = particle.template getField<SURFACE,SINDICATOR>()->getCircumRadius();
    auto position = particle.template getField<GENERAL,POSITION>();

    LatticeR<DESCRIPTOR::d> start, end;
    if ( particles::getBlockParticleIntersection(this->_blockGeometry, 1./this->_converter.getPhysDeltaX(), start, end,
                                      position, circumRadius) ){

      PhysR<T,DESCRIPTOR::d> lever(0.); //Dummy here
      particles::resolved::momentumExchangeAtSurfaceLocation( output,lever,input,
        this->_blockGeometry,this->_blockLattice,this->_converter,particle);
    } 
  }
  return true;
}



}
#endif
