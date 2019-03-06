/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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

#ifndef SUPER_LATTICE_LOCAL_F_2D_HH
#define SUPER_LATTICE_LOCAL_F_2D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<limits>

#include "superLatticeLocalF2D.h"
#include "blockLatticeLocalF2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry2D.h"

namespace olb {

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticeDissipation2D<T,DESCRIPTOR>::SuperLatticeDissipation2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1), _converter(converter)
{
  this->getName() = "dissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeDissipation2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticeDissipation2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{

  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  ////  int lociz= input[3];
  //  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
  //    // local coords are given, fetch local cell and compute value(s)
  //    T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  //    int overlap = this->sLattice.getOverlap();
  //    this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap/*, lociz+overlap*/).computeAllMomenta(rho, uTemp, pi);

  //    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  //    if (util::TensorVal<DESCRIPTOR<T> >::n == 6)
  //      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];

  //    T nuLattice = converter.getLatticeNu();
  //    T omega = converter.getOmega();
  //    T finalResult = PiNeqNormSqr*nuLattice*pow(omega*DESCRIPTOR<T>::invCs2,2)/rho/2.;

  //    return std::vector<T>(1,finalResult);
  //  } else {
  //    return std::vector<T>(); // empty vector
  //  }
  return false;
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysDissipation2D<T,DESCRIPTOR>::SuperLatticePhysDissipation2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physDissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysDissipation2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysDissipation2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticeDensity2D<T,DESCRIPTOR>::SuperLatticeDensity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "density";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeDensity2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)) );
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticeDensity2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticeVelocity2D<T,DESCRIPTOR>::SuperLatticeVelocity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 2)
{
  this->getName() = "velocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVelocity2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticeVelocity2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysStrainRate2D<T,DESCRIPTOR>::SuperLatticePhysStrainRate2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 4)
{
  this->getName() = "physStrainRate";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysStrainRate2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysStrainRate2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{

  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysWallShearStress2D<T, DESCRIPTOR>::SuperLatticePhysWallShearStress2D(
  SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physWallShearStress";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysWallShearStress2D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC), _superGeometry.getBlockGeometry(iC), _material, this->_converter));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysWallShearStress2D<T, DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticeGeometry2D<T,DESCRIPTOR>::SuperLatticeGeometry2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1), _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "geometry";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new  BlockLatticeGeometry2D<T,DESCRIPTOR>(
                                  this->_sLattice.getBlockLattice(iC),
                                  this->_superGeometry.getBlockGeometry(iC),
                                  _material) );
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticeGeometry2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticeRank2D<T,DESCRIPTOR>::SuperLatticeRank2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "rank";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeRank2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)) );
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticeRank2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
    return true;
  } else {
    return false;
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticeCuboid2D<T,DESCRIPTOR>::SuperLatticeCuboid2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "cuboid";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeCuboid2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                                this->_sLattice.getLoadBalancer().glob(iC)) );
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticeCuboid2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
    return true;
  } else {
    return false;
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysPressure2D<T,DESCRIPTOR>::SuperLatticePhysPressure2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physPressure";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysPressure2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_converter));
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysPressure2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  } else {
    return false;
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysVelocity2D<T,DESCRIPTOR>::SuperLatticePhysVelocity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "physVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysVelocity2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                               this->_converter));
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysVelocity2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template <typename T,template <typename U> class DESCRIPTOR>
SuperLatticePhysExternalPorosity2D<T,DESCRIPTOR>::SuperLatticePhysExternalPorosity2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,1)
{
  this->getName() = "ExtPorosityField";
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysExternalPorosity2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  int globIC = input[0];
  if (this->_sLattice.getLoadBalancer().rank(globIC)
      == singleton::mpi().getRank()) {
    int inputLocal[2] = { };

    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    BlockLatticePhysExternalPorosity2D<T,DESCRIPTOR> blockLatticeF(
      this->_sLattice.getExtendedBlockLattice(
        this->_sLattice.getLoadBalancer().loc(globIC)),
      this->_converter);
    blockLatticeF(output, inputLocal);
    return true;
  } else {
    return false;
  }
  return false;
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysExternalVelocity2D<T,DESCRIPTOR>::SuperLatticePhysExternalVelocity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "ExtVelocityField";
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysExternalVelocity2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  int globIC = input[0];
  if (this->_sLattice.getLoadBalancer().rank(globIC)
      == singleton::mpi().getRank()) {
    int inputLocal[2] = { };

    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    BlockLatticePhysExternalVelocity2D<T,DESCRIPTOR> blockLatticeF(
      this->_sLattice.getExtendedBlockLattice(
        this->_sLattice.getLoadBalancer().loc(globIC)),
      this->_converter);
    blockLatticeF(output, inputLocal);
    return true;
  } else {
    return false;
  }
  return false;
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysExternalParticleVelocity2D<T,DESCRIPTOR>::SuperLatticePhysExternalParticleVelocity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "ExtPartVelField";
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysExternalParticleVelocity2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  int globIC = input[0];
  if (this->_sLattice.getLoadBalancer().rank(globIC)
      == singleton::mpi().getRank()) {
    int inputLocal[2] = { };

    int overlap = this->_sLattice.getOverlap();
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    BlockLatticePhysExternalParticleVelocity2D<T,DESCRIPTOR> blockLatticeF(
      this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)),
      this->_converter);
    blockLatticeF(output, inputLocal);
    return true;
  } else {
    return false;
  }
  return false;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysBoundaryForce2D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce2D(
  SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physBoundaryForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysBoundaryForce2D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC), _superGeometry.getBlockGeometry(iC), _material, this->_converter));
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysBoundaryForce2D<T, DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysBoundaryForceIndicator2D<T,DESCRIPTOR>::SuperLatticePhysBoundaryForceIndicator2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  ParticleIndicatorF2D<T,T>& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2),
    _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "physBoundaryForceIndicator";
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysBoundaryForceIndicator2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  Cuboid2D<T>* cub = &_superGeometry.getCuboidGeometry().get(input[0]);

  T posIn[3];
  T physR[2];
  T inside[1];

  cub->getPhysR(posIn, input[1], input[2]);
  output[0] = T();
  output[1] = T();
  output[2] = T();

  if (this->_sLattice.getLoadBalancer().rank(input[0])
      == singleton::mpi().getRank()) {
    physR[0] = posIn[0];
    physR[1] = posIn[1];
    _indicator(inside, physR);
    if (inside[0] > std::numeric_limits<T>::epsilon()) {
      int overlap = this->_sLattice.getOverlap();
      BlockLatticeStructure2D<T,DESCRIPTOR>* bLat = &this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(input[0]));

      for (int iPop = 1; iPop < DESCRIPTOR<T>::q; ++iPop) {
        const int* c = DESCRIPTOR<T>::c[iPop];
        /* cub->getPhysR(physR, input[1] + c[0], input[2] + c[1]);
        _indicator(inside, physR);
        if (inside[0] < 1) { */
        T f = bLat->get(input[1]+overlap + c[0], input[2]+overlap + c[1])[iPop];
        f +=  bLat->get(input[1]+overlap, input[2]+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
        output[0] -= c[0]*f;
        output[1] -= c[1]*f;
        /*}*/
      }
      output[0] = this->_converter.getPhysForce(output[0]);
      output[1] = this->_converter.getPhysForce(output[1]);
      output[2] = (posIn[0] - _indicator.getPos()[0]) * output[1]- (posIn[1] - _indicator.getPos()[1]) * output[0];
    }
    return true;
  } else {
    return false;
  }
}

/*****************************NEW*****************************/

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysVolumeForceIndicator2D<T,DESCRIPTOR>::SuperLatticePhysVolumeForceIndicator2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  SmoothIndicatorF2D<T,T>& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "physBoundaryForceIndicator";
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysVolumeForceIndicator2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  int globIC = input[0];
  int locix = input[1];
  int lociy = input[2];

  std::vector<T> posIn(3, T());
  T physR[2];
  T inside[1];
  posIn = _superGeometry.getPhysR(globIC, locix, lociy);
  output[0] = T();
  output[1] = T();
  output[2] = T();

  T* vel;
  T* d;
  T u[2];
  T rho = T();
  T rho1 = T();
  if (this->_sLattice.getLoadBalancer().rank(globIC)
      == singleton::mpi().getRank()) {
    physR[0] = posIn[0];
    physR[1] = posIn[1];
    _indicator(inside, physR);
    //    if (inside[0] != 0) {
    T u_new[2] = { 0., 0. };
    int overlap = this->_sLattice.getOverlap();

    for (int iPop = 0; iPop < DESCRIPTOR<T>::q; ++iPop) {
      const int* c = DESCRIPTOR<T>::c[iPop];
      std::vector<int> input2(input, input + 3);
      input2[1] = input[1] + c[0];
      input2[2] = input[2] + c[1];

      T f = this->_sLattice.getExtendedBlockLattice(
              this->_sLattice.getLoadBalancer().loc(globIC)).get(
              locix + overlap + c[0], lociy + overlap + c[1])[iPop];
      u_new[0] += f * c[0];
      u_new[1] += f * c[1];
      rho1 += f;
    }
    rho1 += 1.;
    u_new[0] /= rho1;
    u_new[1] /= rho1;

    this->_sLattice.getExtendedBlockLattice(
      this->_sLattice.getLoadBalancer().loc(globIC)).get(locix + overlap,
          lociy + overlap)
    .computeRhoU(rho, u);



    //        if (inside[0] == 0) {
    d = this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).
        get(locix+overlap, lociy+overlap).getExternal(0);
    vel = this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).
          get(locix+overlap, lociy+overlap).getExternal(1);
    //          this->_sLattice.getExtendedBlockLattice(this->_sLattice.getLoadBalancer().loc(globIC)).
    //              get(locix+overlap, lociy+overlap).computeRhoU(rho, u);
    ////          if (locix == 16 && lociy == 18)
    std::cout << (u_new[0] - vel[0]) << " " << (u_new[1] - vel[1]) << " " << rho1 << " " <<  rho << " " << rho1-rho <<  " " << *d << std::endl;

    //      T factor = this->_converter.getCharRho(); // *(1-*d);
    output[0] = rho * (u_new[0] - vel[0]);  //  + vel[0];
    output[1] = rho * (u_new[1] - vel[1]);  //  + vel[1];
    //        }
    //    }
    output[0] = this->_converter.getPhysForce(output[0]);
    output[1] = this->_converter.getPhysForce(output[1]);
    output[2] = (posIn[0] - _indicator.getCenter()[0]) * -output[1]
                - (posIn[1] - _indicator.getCenter()[1]) * -output[0];

    std::cout << std::numeric_limits<T>::epsilon() << std::endl;

    return true;
    //    } else {
    //      return true;
    //    }
  } else {
    return false;
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::SuperLatticePhysCorrBoundaryForce2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physCorrBoundaryForce";
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  std::vector<T> force(3, T());
  //  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
  //    int globX = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosX() + locix;
  //    int globY = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosY() + lociy;
  //    int globZ = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosZ() + lociz;
  //    if(superGeometry.getMaterial(globX, globY, globZ) == material) {
  //      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
  //        // Get direction
  //        const int* c = DESCRIPTOR<T>::c[iPop];
  //        // Get next cell located in the current direction
  //        // Check if the next cell is a fluid node
  //        if (superGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
  //          int overlap = this->sLattice.getOverlap();
  //          // Get f_q of next fluid cell where l = opposite(q)
  //          T f = this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
  //          // Get f_l of the boundary cell
  //          // Add f_q and f_opp
  //          f += this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
  //          // Update force
  //          force[0] -= c[0]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //          force[1] -= c[1]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //          force[2] -= c[2]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //        }
  //      }
  //      force[0]=this->converter.physForce(force[0]);
  //      force[1]=this->converter.physForce(force[1]);
  //      force[2]=this->converter.physForce(force[2]);
  //      return force;
  //    }
  //    else {
  //      return force;
  //    }
  //  }
  //  else {
  //    return force;
  //  }
  return false;
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePorosity2D<T,DESCRIPTOR>::SuperLatticePorosity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "porosity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePorosity2D<T,DESCRIPTOR>(
                                 this->_sLattice.getBlockLattice(iC),
                                 this->_superGeometry.getBlockGeometry(iC),
                                 _material,
                                 converter));
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePorosity2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
  // // new version
  // // how to deal with value??

  //  int globIC = input[0];
  //  std::vector<int> inputLocal(3,T());
  //  T* value = new T[1];
  //  int overlap = this->sLattice.getOverlap();
  //  inputLocal[0] = input[1] + overlap;
  //  inputLocal[1] = input[2] + overlap;
  //  inputLocal[2] = input[3] + overlap;
  //  BlockLatticePorosity2D<T,DESCRIPTOR> blockLatticeF( this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)) );
  //  std::vector<T> result(1,value[0]);
  //  delete value;
  //  return result;

  // // old version
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  T* value = new T[1];
  //  int overlap = this->sLattice.getOverlap();
  //  this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeExternalField(0,1,value);
  //  std::vector<T> result(1,value[0]);
  //  delete value;
  //  return result;
//  return false;
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysPermeability2D<T,DESCRIPTOR>::SuperLatticePhysPermeability2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "permeability";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePhysPermeability2D<T,DESCRIPTOR>(
                                  this->_sLattice.getBlockLattice(iC),
                                  this->_superGeometry.getBlockGeometry(iC),
                                  _material,
                                  this->getConverter() ) );
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysPermeability2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  ////  int lociz= input[3];

  //  T* value = new T[1];
  //  int overlap = this->sLattice.getOverlap();
  //  this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap/*, lociz+overlap*/).computeExternalField(0,1,value);
  //  std::vector<T> result(1,this->converter.physPermeability(value[0]));
  //  delete value;
  //  if (!(result[0]<42)&&!(result[0]>42)&&!(result[0]==42)) result[0] = 999999;
  //  if (isinf(result[0])) result[0] = 1e6;
  //  return result;
  return false;
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticePhysDarcyForce2D<T,DESCRIPTOR>::SuperLatticePhysDarcyForce2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "alphaU";
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticePhysDarcyForce2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  SuperLatticePhysPermeability2D<T,DESCRIPTOR> permeability(this->sLattice,this->superGeometry,this->material,this->converter);
  //  SuperLatticeVelocity2D<T,DESCRIPTOR> velocity(this->sLattice);

  //  T nu = this->converter.getCharNu();
  //  T K = permeability(input)[0];
  //  std::vector<T> u = velocity(input);

  //  std::vector<T> result(2,T());
  //  result[0] = -nu/K*u[0];
  //  result[1] = -nu/K*u[1];
  ////  result[2] = -nu/K*u[2];

  //  return result;
  return false;
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperLatticeAverage2D<T,DESCRIPTOR>::SuperLatticeAverage2D(
  SuperLatticeF2D<T,DESCRIPTOR>& f, SuperGeometry2D<T>& superGeometry,
  const int material, T radius)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice(), f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material), _radius(radius)
{
  this->getName() = "Average(" + _f.getName() + ")";
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperLatticeAverage2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  //  CuboidGeometry2D<T>& cGeometry = f.getSuperLattice().getCuboidGeometry();
  //  LoadBalancer<T>& load = f.getSuperLattice().getLoadBalancer();

  //  //create boolean indicator functor isInSphere
  //  std::vector<T> center(3,0);
  //  center[0] = (int)cGeometry.get(input[0]).get_globPosX() + input[1];
  //  center[1] = (int)cGeometry.get(input[0]).get_globPosY() + input[2];
  //  center[2] = (int)cGeometry.get(input[0]).get_globPosZ() + input[3];
  //  SphereAnalyticalF2D<bool,T> isInSphere(center,radius);

  //  // iterate over all cuboids & points and test for material && isInSphere
  //  std::vector<T> tmp( this->_n, T() );
  //  int numVoxels(0);
  //  if (this->superGeometry.getMaterial(center[0],center[1],center[2]) == material) {
  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      for (int iX=0; iX<nX; iX++) {
  //        for (int iY=0; iY<nY; iY++) {
  //          for (int iZ=0; iZ<nZ; iZ++) {
  //            std::vector<T> glob(3,0);
  //            glob[0] = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            glob[1] = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            glob[2] = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->superGeometry.getMaterial(glob[0],glob[1],glob[2]) == material
  //                && isInSphere(glob)[0]==true) {
  //              for (unsigned iD=0; iD<f(load.glob(0),0,0,0).size(); iD++) {
  //                tmp[iD]+=f(load.glob(iC),iX,iY,iZ)[iD];
  //              }
  //              numVoxels++;
  //            }
  //          }
  //        }
  //      }
  //    }

  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
  //#endif
  //    for (int iD=0; iD<f.getTargetDim(); iD++) {
  //#ifdef PARALLEL_MODE_MPI
  //      singleton::mpi().reduceAndBcast(tmp[iD], MPI_SUM);
  //#endif
  //      if (numVoxels>0) {
  //        tmp[iD] /= numVoxels;
  //      }
  //    }
  //  }
  //  return tmp;
  return false;
}

template<typename T,template<typename U> class DESCRIPTOR>
SuperEuklidNorm2D<T,DESCRIPTOR>::SuperEuklidNorm2D(
  SuperLatticeF2D<T,DESCRIPTOR>& f)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice(), 1), _f(f)
{
  this->getName() = "l2(" + _f.getName() + ")";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockEuklidNorm2D<T,DESCRIPTOR>(f.getBlockF(iC)));
  }
}

template<typename T,template<typename U> class DESCRIPTOR>
bool SuperEuklidNorm2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
SuperLatticePhysTemperature2D<T,DESCRIPTOR,ThermalDESCRIPTOR>::SuperLatticePhysTemperature2D(
  SuperLattice2D<T,ThermalDESCRIPTOR>& sLattice, ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR> const& converter)
  : SuperLatticeThermalPhysF2D<T,DESCRIPTOR,ThermalDESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physTemperature";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysTemperature2D<T,DESCRIPTOR,ThermalDESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_converter));
  }
}

template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
bool SuperLatticePhysTemperature2D<T,DESCRIPTOR,ThermalDESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  } else {
    return false;
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticeMomentumExchangeForce2D<T,DESCRIPTOR>::SuperLatticeMomentumExchangeForce2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 std::vector<ParticleIndicatorF2D<T,T>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,4*indicator.size()),
    _superGeometry(superGeometry), _vectorOfIndicator(indicator)
{
  this->getName() = "physMomentumExchangeForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeMomentumExchangeForce2D<T,DESCRIPTOR>(this->_sLattice.getExtendedBlockLattice(iC), this->_superGeometry.getBlockGeometry(iC), indicator, converter, sLattice.getOverlap()));
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticeMomentumExchangeForce2D<T,DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  for (int i=0; i<this->getTargetDim(); i++) {
    output[i] = 0.;
  }
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    int globiC = this->_sLattice.getLoadBalancer().glob(iC);
    if ( this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank() ) {
      this->getBlockF(iC)(output,&input[1]);
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim(); ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
#endif
  return true;

}

}  // end namespace olb

#endif
