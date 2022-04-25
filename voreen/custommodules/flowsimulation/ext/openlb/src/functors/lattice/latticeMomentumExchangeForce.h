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

#ifndef LATTICE_MOMENTUM_EXCHANGE_FORCE_H
#define LATTICE_MOMENTUM_EXCHANGE_FORCE_H



namespace olb {

//Forward declaration
namespace particles{
template<typename T, typename PARTICLETYPE> class ParticleSystem;
template<typename T, typename PARTICLETYPE> class Particle;
}

//TODO: adapt description
/** Functor that returns forces acting on a particle surface, returns data in output for every particle in a row(described are return values for the first particle).
 * \return output[0]-output[2] translational force - physical units
 * \return output[3]-output[5] torque - physical units
 * \return output[7] number of voxels
 */
template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
class SuperLatticeMomentumExchangeForce final : public SuperLatticePhysF<T,DESCRIPTOR> {
public:
  SuperLatticeMomentumExchangeForce( SuperLattice<T,DESCRIPTOR>& sLattice,
                                     SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
                                     particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
                                     const UnitConverter<T,DESCRIPTOR>& converter,
                                     Vector<bool,DESCRIPTOR::d> periodic = Vector<bool,DESCRIPTOR::d> (false),
                                     std::size_t iP0=0 );
  bool operator() (T output[], const int input[]) override;
};

//TODO: adapt description
/** Functor that returns forces acting on a particle surface, returns data in output for every particle in a row(described are return values for the first particle).
 * \return output[0]-output[2] translational force - physical units
 * \return output[3]-output[5] torque - physical units
 * \return output[7] number of voxels
 */

template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
class BlockLatticeMomentumExchangeForce final : public BlockLatticePhysF<T,DESCRIPTOR> {
private:
  BlockGeometry<T,DESCRIPTOR::d>& _blockGeometry;
  BlockLattice<T,DESCRIPTOR>& _blockLattice;
  particles::ParticleSystem<T,PARTICLETYPE>& _particleSystem;
  PhysR<T,DESCRIPTOR::d> _cellMin;
  PhysR<T,DESCRIPTOR::d> _cellMax;
  Vector<bool,DESCRIPTOR::d> _periodic;
  std::size_t _iP0;
public:
  BlockLatticeMomentumExchangeForce( BlockLattice<T,DESCRIPTOR>& blockLattice,
                                     BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                                     particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
                                     const UnitConverter<T,DESCRIPTOR>& converter,
                                     PhysR<T,DESCRIPTOR::d> cellMin = PhysR<T,DESCRIPTOR::d> (0.),
                                     PhysR<T,DESCRIPTOR::d> cellMax = PhysR<T,DESCRIPTOR::d> (0.),
                                     Vector<bool,DESCRIPTOR::d> periodic = Vector<bool,DESCRIPTOR::d> (false),
                                     std::size_t iP0=0 );
  void evaluate(T output[], particles::Particle<T,PARTICLETYPE>& particle, int iP);
  bool operator() (T output[], const int input[]) override;
};





///The following are functors that work in the traditional (output[], input[]) sense,
///They can therefore be used e.g. in the vtk writer as well

/// functor to get pointwise momentum exchange on local lattice
template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
class SuperLatticeMomentumExchangeForceLocal final : public SuperLatticePhysF<T,DESCRIPTOR> {
public:
  SuperLatticeMomentumExchangeForceLocal( SuperLattice<T,DESCRIPTOR>& sLattice,
                                          const UnitConverter<T,DESCRIPTOR>& converter,
                                          SuperGeometry<T,DESCRIPTOR::d>& superGeometry,
                                          particles::ParticleSystem<T,PARTICLETYPE>& particleSystem );
};


template <typename T, typename DESCRIPTOR, typename PARTICLETYPE>
class BlockLatticeMomentumExchangeForceLocal final : public BlockLatticePhysF<T,DESCRIPTOR> {
private:
  BlockLattice<T,DESCRIPTOR>& _blockLattice;
  BlockGeometry<T,DESCRIPTOR::d>& _blockGeometry;
  particles::ParticleSystem<T,PARTICLETYPE>& _particleSystem;
public:
  BlockLatticeMomentumExchangeForceLocal(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                         BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                                         particles::ParticleSystem<T,PARTICLETYPE>& particleSystem,
                                         const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};




}
#endif
