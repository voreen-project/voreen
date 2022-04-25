/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Particle Reynold number classes for Lagrangian two-way coupling methods -- generic implementation.
 */

#ifndef LB_REYNOLDS_NUMBERS_3D_HH
#define LB_REYNOLDS_NUMBERS_3D_HH

namespace olb {

////////////////////// Class ParticleReynoldsNumberBase ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
ParticleReynoldsNumberBase<T,Lattice,Particle>::ParticleReynoldsNumberBase(UnitConverter<T, Lattice>& converter)
  : _converter(converter)
{}


////////////////////// Class NewtonianParticleReynoldsNumber ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NewtonianParticleReynoldsNumber<T,Lattice,Particle>::NewtonianParticleReynoldsNumber(UnitConverter<T, Lattice>& converter)
  : ParticleReynoldsNumberBase<T,Lattice,Particle>(converter)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
T NewtonianParticleReynoldsNumber<T,Lattice,Particle>::operator() ( Particle<T>* p, T magU, int globicFull[])
{
  T ReP = 2. * p->getRad() * magU / this->_converter.getPhysViscosity();
  return ReP > this->_RePmin ? ReP : this->_RePmin;
}


////////////////////// Class PowerLawParticleReynoldsNumber ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
PowerLawParticleReynoldsNumber<T,Lattice,Particle>::PowerLawParticleReynoldsNumber (
  UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice )
  : ParticleReynoldsNumberBase<T,Lattice,Particle>(converter),
    _sLattice(sLattice)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
T PowerLawParticleReynoldsNumber<T,Lattice,Particle>::operator() ( Particle<T>* p, T magU, int globicFull[])
{
  // loc() indicates the local cuboid number locIC of the actual processing thread,
  // for given global cuboid number iC
  // this is to get appropriate particle system on locIC
  int locIC = _sLattice.getLoadBalancer().loc(globicFull[0]);

  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = { globicFull[1], globicFull[2], globicFull[3] };

  // getting the power-law relaxation frequency form the dynamics's external field
  T omega = _sLattice.getBlock(locIC).get (
              latticeRoundedPosP[0],
              latticeRoundedPosP[1],
              latticeRoundedPosP[2] ).template getField<descriptors::OMEGA>();

  // physical viscosity from relaxation time
  T nu = this->_converter.getPhysViscosity (
           (1./omega - 0.5) / descriptors::invCs2<T,Lattice>() );

  T ReP = 2. * p->getRad() * magU / nu;
  return ReP > this->_RePmin ? ReP : this->_RePmin;
}


}

#endif
