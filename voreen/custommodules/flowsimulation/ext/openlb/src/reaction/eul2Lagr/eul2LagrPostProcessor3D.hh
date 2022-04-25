/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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

/** \file
 * Postprocessor to store the particle density, converted to Euler
 *  -- generic implementation
 */
#ifndef EUL2LAGR_POST_PROCESSOR_HH
#define EUL2LAGR_POST_PROCESSOR_HH

namespace olb {


//////// Eul2LagrOperator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, template<typename U> class PARTICLETYPE>
Eul2LagrOperator3D<T,DESCRIPTOR,PARTICLETYPE>::Eul2LagrOperator3D(ParticleSystem3D<T,PARTICLETYPE>& pSystem, SuperGeometry<T,3>& superGeometry)
  : _superGeometry(superGeometry),
    _pSystem(pSystem)
{}

template<typename T, typename DESCRIPTOR, template<typename U> class PARTICLETYPE>
bool Eul2LagrOperator3D<T,DESCRIPTOR,PARTICLETYPE>::operator()(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
#ifdef TEST_L
  OstreamManager clout( std::cout,"eul2LagrOperator::operator()" );
#endif
  // Resetting the external field
  for (int iX=0; iX<blockLattice.getNx(); ++iX) {
    for (int iY=0; iY<blockLattice.getNy(); ++iY) {
      for (int iZ=0; iZ<blockLattice.getNz(); ++iZ) {
        blockLattice.get(iX,iY,iZ).template setField<descriptors::EUL2LAGR>(0.);
      }
    }
  }

  // Populating each voxel
  int overlap = _superGeometry.getOverlap();
  for (auto p : _pSystem.getParticlesPointer()) {
    if (p->getActive()) {
      T physPosP[] { p->getPos()[0], p->getPos()[1], p->getPos()[2] };
      int latticeRoundedPosP[] {0, 0, 0, 0};
      _superGeometry.getCuboidGeometry().getLatticeR (latticeRoundedPosP, physPosP);
      latticeRoundedPosP[1] += overlap;
      latticeRoundedPosP[2] += overlap;
      latticeRoundedPosP[3] += overlap;
      //Populating the lattice
      T eul2LagrRho = blockLattice.get(&latticeRoundedPosP[1]).template getField<descriptors::EUL2LAGR>();
      eul2LagrRho += 1.;
      blockLattice.get(&latticeRoundedPosP[1]).template setField<descriptors::EUL2LAGR>(eul2LagrRho);
    }
  }

  return true;
}


////////  Eul2LagrPostProcessor3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
Eul2LagrPostProcessor3D <T,DESCRIPTOR>::
Eul2LagrPostProcessor3D ( int x0, int x1, int y0, int y1, int z0, int z1,
                          std::shared_ptr<Eul2LagrOperatorBase3D<T,DESCRIPTOR>> eul2LagrOperator )
  :  _x0(x0), _x1(x1), _y0(y0), _y1(y1), _z0(z0), _z1(z1), _eul2LagrOperator(eul2LagrOperator)
{
  this->getName() = "Eul2LagrPostProcessor3D";
}

template<typename T, typename DESCRIPTOR>
Eul2LagrPostProcessor3D <T,DESCRIPTOR>::
Eul2LagrPostProcessor3D (std::shared_ptr<Eul2LagrOperatorBase3D<T,DESCRIPTOR>> eul2LagrOperator)
  :  _x0(0), _x1(0), _y0(0), _y1(0), _z0(0), _z1(0), _eul2LagrOperator(eul2LagrOperator)
{
  this->getName() = "Eul2LagrPostProcessor3D";
}

template<typename T, typename DESCRIPTOR>
void Eul2LagrPostProcessor3D<T,DESCRIPTOR>::
processSubDomain( BlockLattice<T,DESCRIPTOR>& blockLattice,
                  int x0, int x1, int y0, int y1, int z0, int z1 )
{
  _eul2LagrOperator->operator()(blockLattice);
}

template<typename T, typename DESCRIPTOR>
void Eul2LagrPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, _x0, _x1, _y0, _y1, _z0, _z1);
}


//////// Eul2LagrPostProcessorGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, template<typename U> class PARTICLETYPE>
Eul2LagrPostProcessorGenerator3D<T,DESCRIPTOR,PARTICLETYPE>::Eul2LagrPostProcessorGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  SuperParticleSystem3D<T,PARTICLETYPE>& spSys, SuperGeometry<T,3>& superGeometry )
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_)
{
  _eul2LagrOperator = std::make_shared<Eul2LagrOperator3D<T,DESCRIPTOR,PARTICLETYPE>>(spSys[singleton::mpi().getRank()], superGeometry);
}

template<typename T, typename DESCRIPTOR, template<typename U> class PARTICLETYPE>
Eul2LagrPostProcessorGenerator3D<T,DESCRIPTOR,PARTICLETYPE>::Eul2LagrPostProcessorGenerator3D (
  SuperParticleSystem3D<T,PARTICLETYPE>& spSys, SuperGeometry<T,3>& superGeometry )
  : PostProcessorGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0)
{
  _eul2LagrOperator = std::make_shared<Eul2LagrOperator3D<T,DESCRIPTOR,PARTICLETYPE>>(spSys[singleton::mpi().getRank()], superGeometry);
}

template<typename T, typename DESCRIPTOR, template<typename U> class PARTICLETYPE>
PostProcessor3D<T,DESCRIPTOR>* Eul2LagrPostProcessorGenerator3D<T,DESCRIPTOR,PARTICLETYPE>::generate() const
{
  return new Eul2LagrPostProcessor3D<T,DESCRIPTOR>(this->x0,this->x1,this->y0,this->y1,this->z0,this->z1,
         this->_eul2LagrOperator);
}

template<typename T, typename DESCRIPTOR, template<typename U> class PARTICLETYPE>
PostProcessorGenerator3D<T,DESCRIPTOR>* Eul2LagrPostProcessorGenerator3D<T,DESCRIPTOR,PARTICLETYPE>::clone() const
{
  return new Eul2LagrPostProcessorGenerator3D<T,DESCRIPTOR,PARTICLETYPE>(*this);
}


////////  Eul2LagrNormDistrPostProcessor3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
Eul2LagrNormDistrPostProcessor3D <T,DESCRIPTOR>::
Eul2LagrNormDistrPostProcessor3D(int x0, int x1, int y0, int y1, int z0, int z1, T mean, T stdDev, SuperGeometry<T,3>& superGeometry)
  :  _x0(x0), _x1(x1), _y0(y0), _y1(y1), _z0(z0), _z1(z1), _randomNormal(mean, stdDev), _superGeometry(superGeometry)
{
  this->getName() = "Eul2LagrNormDistrPostProcessor3D";
}

template<typename T, typename DESCRIPTOR>
Eul2LagrNormDistrPostProcessor3D <T,DESCRIPTOR>::
Eul2LagrNormDistrPostProcessor3D (T mean, T stdDev, SuperGeometry<T,3>& superGeometry)
  :  _x0(0), _x1(0), _y0(0), _y1(0), _z0(0), _z1(0), _randomNormal(mean, stdDev), _superGeometry(superGeometry)
{
  this->getName() = "Eul2LagrNormDistrPostProcessor3D";
}

template<typename T, typename DESCRIPTOR>
void Eul2LagrNormDistrPostProcessor3D<T,DESCRIPTOR>::
processSubDomain( BlockLattice<T,DESCRIPTOR>& blockLattice,
                  int x0, int x1, int y0, int y1, int z0, int z1 )
{
  std::cout << "ciao" << std::endl;
  T output[1];
  T input[3];
  _randomNormal(output, input);

  T physPos[] {output[0], 0., 0.};
  int latticeRoundedPos[] {0, 0, 0, 0};
  _superGeometry.getCuboidGeometry().getLatticeR (latticeRoundedPos, physPos);
  int overlap = _superGeometry.getOverlap();
  latticeRoundedPos[1] += overlap;
  latticeRoundedPos[2] += overlap;
  latticeRoundedPos[3] += overlap;
  //Populating the lattice
  T eul2LagrRho = blockLattice.get(&latticeRoundedPos[1]).template getField<descriptors::EUL2LAGR>();
  eul2LagrRho += 1.;
  blockLattice.get(&latticeRoundedPos[1]).template setField<descriptors::EUL2LAGR>(eul2LagrRho);
}

template<typename T, typename DESCRIPTOR>
void Eul2LagrNormDistrPostProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, _x0, _x1, _y0, _y1, _z0, _z1);
}


//////// Eul2LagrNormDistrPostProcessorGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR>
Eul2LagrNormDistrPostProcessorGenerator3D<T,DESCRIPTOR>::Eul2LagrNormDistrPostProcessorGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T mean, T stdDev, SuperGeometry<T,3>& superGeometry )
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), _mean(mean), _stdDev(stdDev), _superGeometry(superGeometry)
{ }

template<typename T, typename DESCRIPTOR>
Eul2LagrNormDistrPostProcessorGenerator3D<T,DESCRIPTOR>::Eul2LagrNormDistrPostProcessorGenerator3D(T mean, T stdDev, SuperGeometry<T,3>& superGeometry)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(0, 0, 0, 0, 0, 0), _mean(mean), _stdDev(stdDev), _superGeometry(superGeometry)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* Eul2LagrNormDistrPostProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new Eul2LagrNormDistrPostProcessor3D<T,DESCRIPTOR>(this->x0,this->x1,this->y0,this->y1,this->z0,this->z1,
         this->_mean,this->_stdDev,this->_superGeometry);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>* Eul2LagrNormDistrPostProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new Eul2LagrNormDistrPostProcessorGenerator3D<T,DESCRIPTOR>(*this);
}


} // olb

#endif
