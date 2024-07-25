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

#ifndef LATTICE_STOKES_DRAG_FORCE_HH
#define LATTICE_STOKES_DRAG_FORCE_HH

namespace olb {


template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
SuperLatticeStokesDragForce<T,DESCRIPTOR,PARTICLETYPE>::SuperLatticeStokesDragForce(
  SuperLattice<T,DESCRIPTOR>& sLattice,
  SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter,
  Vector<bool,DESCRIPTOR::d> periodic,
  std::size_t iP0 )
  : SuperLatticePhysF<T,DESCRIPTOR>(sLattice,converter,
      (DESCRIPTOR::d)*(particleSystem.size()-iP0))
{
  this->getName() = "physStokesDragForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  PhysR<T,DESCRIPTOR::d> min = superGeometry.getStatistics().getMinPhysR( 1 ); //TODO: disable hardcoded 1
  PhysR<T,DESCRIPTOR::d> max = superGeometry.getStatistics().getMaxPhysR( 1 ); //TODO: disable hardcoded 1
  for (int iC = 0; iC < maxC; iC++) {
    int globiC = this->_sLattice.getLoadBalancer().glob(iC);
    this->_blockF.emplace_back( new BlockLatticeStokesDragForce<T,DESCRIPTOR, PARTICLETYPE>(
      this->_sLattice.getBlock(iC),
      superGeometry.getBlockGeometry(iC),
      &this->_sLattice.getCuboidGeometry().get(globiC),
      particleSystem, converter, min, max, periodic, iP0));
  }
}

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool SuperLatticeStokesDragForce<T,DESCRIPTOR,PARTICLETYPE>::operator() (T output[],
    const int input[])
{
  for (std::size_t iS=0; iS<this->getTargetDim(); ++iS) {
    output[iS] = 0.;
  }
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    int globiC = this->_sLattice.getLoadBalancer().glob(iC);
    if ( this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank() ) {
      this->getBlockF(iC)(output,&input[1]);
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (std::size_t iS = 0; iS < this->getTargetDim(); ++iS) {
    singleton::mpi().reduceAndBcast(output[iS], MPI_SUM);
  }
#endif
  return true;

}




template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
BlockLatticeStokesDragForce<T, DESCRIPTOR, PARTICLETYPE>::BlockLatticeStokesDragForce(
  BlockLattice<T, DESCRIPTOR>& blockLattice,
  BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
  Cuboid<T,DESCRIPTOR::d>* cuboid,
  particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
  const UnitConverter<T,DESCRIPTOR>& converter,
  PhysR<T,DESCRIPTOR::d>cellMin, PhysR<T,DESCRIPTOR::d> cellMax,
  Vector<bool,DESCRIPTOR::d> periodic,
  std::size_t iP0 )
  : BlockLatticePhysF<T,DESCRIPTOR>(blockLattice, converter, particleSystem.size()),
    _blockGeometry(blockGeometry), _blockLattice(blockLattice), _cuboid(cuboid),
    _particleSystem(particleSystem),
    _cellMin(cellMin), _cellMax(cellMax), _periodic(periodic), _iP0(iP0)
{
  this->getName() = "physStokesDragForce";
}

template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
void BlockLatticeStokesDragForce<T, DESCRIPTOR, PARTICLETYPE>::evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP)
{
  //TODO: adapt serialization size etc!
  const unsigned D = DESCRIPTOR::d;
  const int serialSize = D;

  using namespace descriptors;

  //Retrieve particle quantities
  auto position = particle.template getField<GENERAL,POSITION>();
  auto radius = particle.template getField<PHYSPROPERTIES,RADIUS>();
  auto velocity = particle.template getField<MOBILITY,VELOCITY>();
  T mass = particle.template getField<PHYSPROPERTIES,MASS>();
  T* positionArray = position.data();

  //TODO: check, whether creation can be avoided by using a second _blockF list
  BlockLatticeInterpPhysVelocity<T,DESCRIPTOR> blockInterpPhysVelF( _blockLattice, this->_converter, _cuboid); 

  //Calculate general constants
  T dTinv = 1./this->_converter.getPhysDeltaT();
  T C1 = 6. * M_PI * this->_converter.getPhysViscosity()
       * this->_converter.getPhysDensity() * this->_converter.getConversionFactorTime();

  //Calculate particle coefficiants
  T c = C1 * radius * 1./mass;
  T C2 = 1. / (1. + c);

  //Retrieve block bounds    
  LatticeR<D> extend = _blockGeometry.getExtent();
  PhysR<T,D> physExtend = this->_converter.getPhysDeltaX()*extend;
  PhysR<T,D> origin = _blockGeometry.getOrigin();
  PhysR<T,D> end = origin + physExtend;

  //Check whether inside cuboid
  bool inside = true;
  for (int iDim=0; iDim<D; ++iDim){
    inside = inside 
          && (positionArray[iDim] > origin[iDim])
          && (positionArray[iDim] < end[iDim]);
  }
  if (inside){

    //Retrieve interpolated velocity at position
    T fluidVelArray[D] = {0.};
    blockInterpPhysVelF(fluidVelArray, positionArray);

    //Calculate stokes force
    Vector<T,D> force(0.);
    for (int iDim = 0; iDim < PARTICLETYPE::d; ++iDim) {
      force[iDim] = mass * dTinv
        * ((c * fluidVelArray[iDim] + velocity[iDim]) * C2 - velocity[iDim]);
    }
   
    //Add force and torque to output
    std::size_t iPeval = iP-_iP0; //Shifted output index, if iP!=0
    for (unsigned iDim=0; iDim<D; iDim++) {
      output[iDim+iPeval*serialSize] += force[iDim];
    }

  }
}



template<typename T, typename DESCRIPTOR, typename PARTICLETYPE>
bool BlockLatticeStokesDragForce<T, DESCRIPTOR, PARTICLETYPE>::operator()(T output[], const int input[])
{
  using namespace descriptors;
  // iterate over all particles in _indicator //TODO: add periodic treatment analogous to momentum exchange
  for (std::size_t iP=_iP0; iP!=_particleSystem.size(); iP++) {
    auto particle = _particleSystem.get(iP);
    if constexpr( PARTICLETYPE::template providesNested<NUMERICPROPERTIES,ACTIVE>() ) {
      if (particle.template getField<NUMERICPROPERTIES,ACTIVE>()){
        evaluate(output, particle, iP);
      }
    } else {
      evaluate(output, particle, iP);
    }
  }
  return true;
}



}
#endif
