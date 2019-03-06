/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_LOCAL_F_3D_HH
#define BLOCK_LATTICE_LOCAL_F_3D_HH


#include<cmath>

#include "blockLatticeLocalF3D.h"
#include "blockBaseF3D.h"
#include "core/blockLatticeStructure3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeFpop3D<T, DESCRIPTOR>::BlockLatticeFpop3D(BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, DESCRIPTOR<T>::q)
{
  this->getName() = "fPop";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeFpop3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int iPop = 0; iPop < DESCRIPTOR<T>::q; ++iPop) {
    output[iPop] =
      this->_blockLattice.get(input[0], input[1], input[2])[iPop]
      + DESCRIPTOR<T>::t[iPop];
  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeDissipation3D<T, DESCRIPTOR>::BlockLatticeDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1), _converter(converter)
{
  this->getName() = "dissipation";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T nuLattice = _converter.getLatticeViscosity();
  T omega = 1. / _converter.getLatticeRelaxationTime();
  output[0] = PiNeqNormSqr * nuLattice
              * pow(omega * DESCRIPTOR<T>::invCs2, 2) / rho / 2.;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysDissipation3D<T, DESCRIPTOR>::BlockLatticePhysDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _overlap(overlap),
    _converter(converter)
{
  this->getName() = "physDissipation";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap, input[2]+_overlap
  ).computeAllMomenta(rho, uTemp, pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T nuLattice = this->_converter.getLatticeViscosity();
  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();
  output[0] = PiNeqNormSqr * nuLattice
              * pow(omega * DESCRIPTOR<T>::invCs2 / rho, 2) / 2.
              * this->_converter.getPhysViscosity() / nuLattice / dt / dt;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::BlockLatticeEffevtiveDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
  LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter), _smagoConst(smagoConst), _LESdynamics(LESdynamics)
{
  this->getName() = "EffevtiveDissipation";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T omegaEff = _LESdynamics.getEffectiveOmega(this->_blockLattice.get(input[0], input[1], input[2]));
  T nuEff = ((1./omegaEff)-0.5)/DESCRIPTOR<T>::invCs2;

  output[0] = PiNeqNormSqr * (nuEff)
              * pow(omegaEff * DESCRIPTOR<T>::invCs2 / rho, 2)  / 2.;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::BlockLatticePhysEffevtiveDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
  LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter), _smagoConst(smagoConst), _LESdynamics(LESdynamics)
{
  this->getName() = "physEffevtiveDissipation";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR<T> >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T dt = 1./ _converter.getConversionFactorTime();
  T omegaEff = _LESdynamics.getEffectiveOmega(this->_blockLattice.get(input[0], input[1], input[2]));
  T nuEff = ((1./omegaEff)-0.5)/DESCRIPTOR<T>::invCs2;  // BGK shear viscosity definition

  output[0] = PiNeqNormSqr * nuEff
              * pow(omegaEff * DESCRIPTOR<T>::invCs2 / rho, 2) / 2.
              * _converter.getPhysViscosity() / _converter.getLatticeViscosity() / dt / dt;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeDensity3D<T, DESCRIPTOR>::BlockLatticeDensity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "density";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeDensity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = this->_blockLattice.get(input[0], input[1], input[2]).computeRho();
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeVelocity3D<T, DESCRIPTOR>::BlockLatticeVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3)
{
  this->getName() = "velocity";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho;
  this->_blockLattice.get(input[0], input[1], input[2]).computeRhoU(rho, output);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeExternalVelocity3D<T, DESCRIPTOR>::BlockLatticeExternalVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3)
{
  this->getName() = "externalVelocity";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeExternalVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T* ExtVel = this->_blockLattice.get(input[0], input[1], input[2]).getExternal(DESCRIPTOR<T>::ExternalField::velocityIsAt);
  for (int iVel=0; iVel<DESCRIPTOR<T>::d; ++iVel) {
    output[iVel] = ExtVel[iVel];
  }
  //ExtVel = NULL; // Necessary to desallouer la memoire
  //delete ExtVel;
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeStrainRate3D<T, DESCRIPTOR>::BlockLatticeStrainRate3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 9)
{
  this->getName() = "strainRate";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeStrainRate3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T omega = 1. / this->_converter.getLatticeRelaxationTime();

  output[0] = -pi[0] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[1] = -pi[1] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[2] = -pi[2] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[3] = -pi[1] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[4] = -pi[3] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[5] = -pi[4] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[6] = -pi[2] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[7] = -pi[4] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;
  output[8] = -pi[5] * omega * DESCRIPTOR<T>::invCs2 / rho / 2.;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysStrainRate3D<T, DESCRIPTOR>::BlockLatticePhysStrainRate3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 9),
    _overlap(overlap)
{
  this->getName() = "strainRate";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysStrainRate3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR<T>::d], pi[util::TensorVal<DESCRIPTOR<T> >::n];
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap, input[2]+_overlap
  ).computeAllMomenta(rho, uTemp, pi);

  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();

  output[0] = -pi[0] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[1] = -pi[1] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[2] = -pi[2] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[3] = -pi[1] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[4] = -pi[3] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[5] = -pi[4] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[6] = -pi[2] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[7] = -pi[4] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;
  output[8] = -pi[5] * omega * DESCRIPTOR<T>::invCs2 / rho / 2. / dt;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeGeometry3D<T, DESCRIPTOR>::BlockLatticeGeometry3D(BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
    BlockGeometryStructure3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "geometry";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeGeometry3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = _blockGeometry.getMaterial(input[0], input[1], input[2]);

  if (_material != -1) {
    if ( util::nearZero(_material-output[0]) ) {
      output[0] = 1.;
      return true;
    } else {
      output[0] = 0.;
      return true;
    }
  }
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeRank3D<T,DESCRIPTOR>::BlockLatticeRank3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "rank";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeRank3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = singleton::mpi().getRank() + 1;
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeCuboid3D<T,DESCRIPTOR>::BlockLatticeCuboid3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, int iC)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _iC(iC)
{
  this->getName() = "cuboid";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeCuboid3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = _iC + 1;
  return true;
}



template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysPressure3D<T, DESCRIPTOR>::BlockLatticePhysPressure3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1),
    _overlap(overlap)
{
  this->getName() = "physPressure";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysPressure3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  // lattice pressure = c_s^2 ( rho -1 )
  T latticePressure = ( this->_blockLattice.get(input[0]+_overlap, input[1]+_overlap, input[2]+_overlap).computeRho() - 1.0) / DESCRIPTOR<T>::invCs2;
  output[0] = this->_converter.getPhysPressure(latticePressure);

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysVelocity3D<T, DESCRIPTOR>::BlockLatticePhysVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter,
  bool print)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _overlap(overlap),
    _print(print)
{
  this->getName() = "physVelocity";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (_print) {
    std::cout << input[0] << " " << input[1] << " " << input[2] << " | "
              << singleton::mpi().getRank() << std::endl;
  }

  T rho;
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap, input[2]+_overlap).computeRhoU(rho, output);
  output[0] = this->_converter.getPhysVelocity(output[0]);
  output[1] = this->_converter.getPhysVelocity(output[1]);
  output[2] = this->_converter.getPhysVelocity(output[2]);

  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysExternalPorosity3D<T,DESCRIPTOR>::BlockLatticePhysExternalPorosity3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,2),
    _overlap(overlap)
{
  this->getName() = "ExtPorosityField";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysExternalPorosity3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap, input[2]+_overlap
  ).computeExternalField(DESCRIPTOR<T>::ExternalField::porosityIsAt, 1, output);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysExternalVelocity3D<T, DESCRIPTOR>::BlockLatticePhysExternalVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3)
{
  this->getName() = "physVelExtField";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysExternalVelocity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    DESCRIPTOR<T>::ExternalField::velocityBeginsAt,
    DESCRIPTOR<T>::ExternalField::sizeOfVelocity, output);
  output[0] = this->_converter.getPhysVelocity(output[0]);
  output[1] = this->_converter.getPhysVelocity(output[1]);
  output[2] = this->_converter.getPhysVelocity(output[2]);
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysExternalParticleVelocity3D<T,DESCRIPTOR>::BlockLatticePhysExternalParticleVelocity3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtParticleVelocityField";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysExternalParticleVelocity3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T* foo = this->_blockLattice.get(input[0],input[1],input[2]).getExternal(0);

  if (foo[4] > std::numeric_limits<T>::epsilon()) {
    output[0]=this->_converter.getPhysVelocity(foo[1]/foo[4]);
    output[1]=this->_converter.getPhysVelocity(foo[2]/foo[4]);
    output[2]=this->_converter.getPhysVelocity(foo[3]/foo[4]);
    return true;
  }
  output[0]=this->_converter.getPhysVelocity(foo[1]);
  output[1]=this->_converter.getPhysVelocity(foo[2]);
  output[2]=this->_converter.getPhysVelocity(foo[3]);
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeExternal3D<T, DESCRIPTOR>::BlockLatticeExternal3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, int start, int size)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, size), _start(start), _size(size)
{
  this->getName() = "extField";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeExternal3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    _start, _size, output);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysExternal3D<T, DESCRIPTOR>::BlockLatticePhysExternal3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  T convFactorToPhysUnits, int offset, int size)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
    _convFactorToPhysUnits(convFactorToPhysUnits),
    _offset(offset), _size(size)
{
  this->getName() = "physExtField";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysExternal3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    _offset, _size, output);
  output[0] = _convFactorToPhysUnits*output[0];
  output[1] = _convFactorToPhysUnits*output[1];
  output[2] = _convFactorToPhysUnits*output[2];
  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysBoundaryForce3D<T,DESCRIPTOR>::BlockLatticePhysBoundaryForce3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>& blockGeometry,
  int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physBoundaryForce";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysBoundaryForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();
  output[1] = T();
  output[2] = T();

  if (this->_blockGeometry.get(input[0],input[1],input[2]) == _material) {
    for (int iPop = 1; iPop < DESCRIPTOR<T>::q; ++iPop) {
      // Get direction
      const int* c = DESCRIPTOR<T>::c[iPop];
      // Get next cell located in the current direction
      // Check if the next cell is a fluid node
      if (_blockGeometry.get(input[0] + c[0], input[1] + c[1], input[2] + c[2]) == 1) {
        // Get f_q of next fluid cell where l = opposite(q)
        T f = this->_blockLattice.get(input[0] + c[0], input[1] + c[1], input[2] + c[2])[iPop];
        // Get f_l of the boundary cell
        // Add f_q and f_opp
        f += this->_blockLattice.get(input[0], input[1], input[2])[util::opposite<DESCRIPTOR<T> >(iPop)];
        // Update force
        output[0] -= c[0] * f;
        output[1] -= c[1] * f;
        output[2] -= c[2] * f;
      }
    }
    output[0] = this->_converter.getPhysForce(output[0]);
    output[1] = this->_converter.getPhysForce(output[1]);
    output[2] = this->_converter.getPhysForce(output[2]);
    return true;
  } else {
    return true;
  }
}

template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysWallShearStress3D<T,DESCRIPTOR>::BlockLatticePhysWallShearStress3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>& blockGeometry,
  int material, const UnitConverter<T,DESCRIPTOR>& converter, IndicatorF3D<T>& indicator)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,1),
    _blockGeometry(blockGeometry),_blockLattice(blockLattice), _material(material), _converter(converter)
{
  this->getName() = "PhysWallShearStress";
  T scaling = converter.getConversionFactorLength() * 0.1;
  T omega = 1. / _converter.getLatticeRelaxationTime();
  T dt = _converter.getConversionFactorTime();
  _physFactor = -omega * DESCRIPTOR<T>::invCs2 / dt * _converter.getPhysDensity() * _converter.getPhysViscosity();
  std::vector<int> discreteNormalOutwards(4, 0);
  for (int iX = 0; iX < _blockGeometry.getNx(); iX++) {
    _discreteNormal.resize(_blockGeometry.getNx());
    _normal.resize(_blockGeometry.getNx());
    for (int iY = 0; iY < _blockGeometry.getNy(); iY++) {
      _discreteNormal[iX].resize(_blockGeometry.getNy());
      _normal[iX].resize(_blockGeometry.getNy());
      for (int iZ = 0; iZ < _blockGeometry.getNz(); iZ++) {
        _discreteNormal[iX][iY].resize(_blockGeometry.getNz());
        _normal[iX][iY].resize(_blockGeometry.getNz());
        if (_blockGeometry.get(iX, iY, iZ) == _material) {
          discreteNormalOutwards = _blockGeometry.getStatistics().getType(iX, iY, iZ);
          _discreteNormal[iX][iY][iZ].resize(3);
          _normal[iX][iY][iZ].resize(3);
          // _discreteNormal pointing inside the fluid domain
          _discreteNormal[iX][iY][iZ][0] = -discreteNormalOutwards[1];
          _discreteNormal[iX][iY][iZ][1] = -discreteNormalOutwards[2];
          _discreteNormal[iX][iY][iZ][2] = -discreteNormalOutwards[3];

          T physR[3];
          _blockGeometry.getPhysR(physR,iX, iY, iZ);
          Vector<T,3> origin(physR[0],physR[1],physR[2]);
          T distance = 0.;
          Vector<T,3> direction(0.,0.,0.);
          Vector<T,3> normal(0.,0.,0.);
          int smallestDistance_i = 0;
          T smallestDistance = 0.;
          bool firstHit = true;
          origin[0] = physR[0];
          origin[1] = physR[1];
          origin[2] = physR[2];
          int discreteDirection[6][3];
          for (int i=0; i < 6; i++) {
            for (int j=0; j < 3; j++) {
              discreteDirection[i][j] = 0;
            }
          }
          discreteDirection[0][0] = 1;
          discreteDirection[1][0] = -1;
          discreteDirection[2][1] = 1;
          discreteDirection[3][1] = -1;
          discreteDirection[4][2] = 1;
          discreteDirection[5][2] = -1;
          for (int i=0; i < 6; i++) {
            direction[0] = discreteDirection[i][0] * scaling;
            direction[1] = discreteDirection[i][1] * scaling;
            direction[2] = discreteDirection[i][2] * scaling;
            if (indicator.distance(distance, origin, direction)) {
              if (firstHit) {
                smallestDistance = distance;
                smallestDistance_i = i;
                firstHit = false;
              } else {
                if (distance < smallestDistance) {
                  smallestDistance = distance;
                  smallestDistance_i = i;
                }
              }
            }
          }
          direction[0] = discreteDirection[smallestDistance_i][0] * scaling;
          direction[1] = discreteDirection[smallestDistance_i][1] * scaling;
          direction[2] = discreteDirection[smallestDistance_i][2] * scaling;
          Vector<T,3> direction2(direction);
          Vector<T,3> direction3(direction);
          if (smallestDistance_i == 0 || smallestDistance_i == 1 ) {
            direction2[1] = direction2[0] * scaling;
            direction3[2] = direction3[0] * scaling;
          } else if (smallestDistance_i == 2 || smallestDistance_i == 3 ) {
            direction2[0] = direction2[1] * scaling;
            direction3[2] = direction3[1] * scaling;
          } else {
            direction2[0] = direction2[2] * scaling;
            direction3[1] = direction3[2] * scaling;
          }

          Vector<S,3> directionN = direction*(1/const_cast<Vector<S,3>&> (direction).norm());
          Vector<S,3> POS(origin + smallestDistance*directionN); //Point on Surface

          indicator.distance(distance, origin, direction2);
          Vector<S,3> direction2N = direction2*(1/const_cast<Vector<S,3>&> (direction2).norm());
          Vector<S,3> POS2(origin + distance*direction2N); //Point on Surface

          indicator.distance(distance, origin, direction3);

          Vector<S,3> direction3N = direction3*(1/const_cast<Vector<S,3>&> (direction3).norm());
          Vector<S,3> POS3(origin + distance*direction3N); //Point on Surface

          Vector<S,3> vec1 (POS - POS2);
          Vector<S,3> vec2 (POS - POS3);

          normal[0] = -(vec1[1]*vec2[2] - vec1[2]*vec2[1]);
          normal[1] = -(vec1[2]*vec2[0] - vec1[0]*vec2[2]);
          normal[2] = -(vec1[0]*vec2[1] - vec1[1]*vec2[0]);

          T normalMagnitude = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
          normal[0] /= normalMagnitude;
          normal[1] /= normalMagnitude;
          normal[2] /= normalMagnitude;
          _normal[iX][iY][iZ][0] = normal[0];
          _normal[iX][iY][iZ][1] = normal[1];
          _normal[iX][iY][iZ][2] = normal[2];
        }
      }
    }
  }
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysWallShearStress3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();

  if (_blockGeometry.get(input[0],input[1],input[2]) == _material) {
    T traction[3];
    T pi[6];
    T rho = _blockLattice.get(input[0] + _discreteNormal[input[0]][input[1]][input[2]][0],
                              input[1] + _discreteNormal[input[0]][input[1]][input[2]][1],
                              input[2] + _discreteNormal[input[0]][input[1]][input[2]][2]).computeRho();
    _blockLattice.get(input[0] +   _discreteNormal[input[0]][input[1]][input[2]][0],
                      input[1] +   _discreteNormal[input[0]][input[1]][input[2]][1],
                      input[2] +   _discreteNormal[input[0]][input[1]][input[2]][2]).computeStress(pi);

    traction[0] = pi[0]*_physFactor/rho*_normal[input[0]][input[1]][input[2]][0] +
                  pi[1]*_physFactor/rho*_normal[input[0]][input[1]][input[2]][1] +
                  pi[2]*_physFactor/rho*_normal[input[0]][input[1]][input[2]][2];
    traction[1] = pi[1]*_physFactor/rho*_normal[input[0]][input[1]][input[2]][0] +
                  pi[3]*_physFactor/rho*_normal[input[0]][input[1]][input[2]][1] +
                  pi[4]*_physFactor/rho*_normal[input[0]][input[1]][input[2]][2];
    traction[2] = pi[2]*_physFactor/rho*_normal[input[0]][input[1]][input[2]][0] +
                  pi[4]*_physFactor/rho*_normal[input[0]][input[1]][input[2]][1] +
                  pi[5]*_physFactor/rho*_normal[input[0]][input[1]][input[2]][2];

    T traction_normal_SP;
    T tractionNormalComponent[3];
    // scalar product of traction and normal vector
    traction_normal_SP = traction[0] * _normal[input[0]][input[1]][input[2]][0] +
                         traction[1] * _normal[input[0]][input[1]][input[2]][1] +
                         traction[2] * _normal[input[0]][input[1]][input[2]][2];
    tractionNormalComponent[0] = traction_normal_SP * _normal[input[0]][input[1]][input[2]][0];
    tractionNormalComponent[1] = traction_normal_SP * _normal[input[0]][input[1]][input[2]][1];
    tractionNormalComponent[2] = traction_normal_SP * _normal[input[0]][input[1]][input[2]][2];

    T WSS[3];
    WSS[0] = traction[0] - tractionNormalComponent[0];
    WSS[1] = traction[1] - tractionNormalComponent[1];
    WSS[2] = traction[2] - tractionNormalComponent[2];
    // magnitude of the wall shear stress vector
    output[0] = sqrt(WSS[0]*WSS[0] + WSS[1]*WSS[1] + WSS[2]*WSS[2]);
    return true;
  } else {
    return true;
  }
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::BlockLatticePhysCorrBoundaryForce3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "physCorrBoundaryForce";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  std::vector<T> force(3, T());
  //  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
  //    int globX = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosX() + locix;
  //    int globY = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosY() + lociy;
  //    int globZ = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosZ() + lociz;
  //    if(BlockGeometry.getMaterial(globX, globY, globZ) == material) {
  //      for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
  //        // Get direction
  //        const int* c = DESCRIPTOR<T>::c[iPop];
  //        // Get next cell located in the current direction
  //        // Check if the next cell is a fluid node
  //        if (BlockGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
  //          int overlap = this->_blockLattice.getOverlap();
  //          // Get f_q of next fluid cell where l = opposite(q)
  //          T f = this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
  //          // Get f_l of the boundary cell
  //          // Add f_q and f_opp
  //          f += this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
  //          // Update force
  //          force[0] -= c[0]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //          force[1] -= c[1]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //          force[2] -= c[2]*(f-2.*DESCRIPTOR<T>::t[iPop]);
  //        }
  //      }
  //      force[0]=this->_converter.physForce(force[0]);
  //      force[1]=this->_converter.physForce(force[1]);
  //      force[2]=this->_converter.physForce(force[2]);
  //      return force;
  //    }
  //    else {
  //      return force;
  //    }
  //  }
  //  else {
  //    return force;
  //  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeExternalField3D<T, DESCRIPTOR>::BlockLatticeExternalField3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, int beginsAt, int sizeOf)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, sizeOf),
    _beginsAt(beginsAt), _sizeOf(sizeOf)
{
  this->getName() = "externalField";
}

// under construction
template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeExternalField3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(_beginsAt, _sizeOf, output);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePorosity3D<T, DESCRIPTOR>::BlockLatticePorosity3D(BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "porosity";
}

// under construction
template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePorosity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
#ifndef excludeDualDynamics
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    DESCRIPTOR<T>::ExternalField::porosityIsAt, DESCRIPTOR<T>::ExternalField::sizeOfPorosity, output);
#else
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    DESCRIPTOR<T>::ExternalField::porosityIsAt, 1, output);

#endif
  return true;
}


template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysPermeability3D<T, DESCRIPTOR>::BlockLatticePhysPermeability3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1)
{
  this->getName() = "permeability";
}

//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticePhysPermeability3D<T, DESCRIPTOR>::BlockLatticePhysPermeability3D(
//  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
//  BlockGeometry3D<T>& blockGeometry, int material,
//  const UnitConverter<T>& converter)
//  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1),
//    _blockGeometry(blockGeometry),
//    _material(material)
//{
//  this->getName() = "permeability";
//}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysPermeability3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T value;
  // get porosity
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    DESCRIPTOR<T>::ExternalField::porosityIsAt, DESCRIPTOR<T>::ExternalField::sizeOfPorosity, &value);
  // convert to physPermeability
//  output[0] = this->_converter.physPermeability(value);  // TODO converter MG
  // \todo Why do we need this???
  if (output[0] >= 42 && output[0] <= 42 && output[0] != 42) {
    output[0] = 999999;
  }
  if (std::isinf(output[0])) {
    output[0] = 1e6;
  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::BlockLatticePhysCroppedPermeability3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1)
{
  this->getName() = "cropped_permeability";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T value;
  // get porosity
  this->_blockLattice.get(input[0], input[1], input[2]).computeExternalField(
    DESCRIPTOR<T>::ExternalField::porosityIsAt, DESCRIPTOR<T>::ExternalField::sizeOfPorosity, &value);
  // convert to physPermeability
//  output[0] = _converter.physPermeability(value); // TODO converter MG
  // \todo Why do we need this???
  if (output[0] >= 42 && output[0] <= 42 && output[0] != 42) {
    output[0] = 1;
  }
  if (std::isinf(output[0])) {
    output[0] = 1;
  }
  if (output[0] > 1.) {
    output[0] = 1.;
  }
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticePhysDarcyForce3D<T, DESCRIPTOR>::BlockLatticePhysDarcyForce3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "alphaU";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticePhysDarcyForce3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  BlockLatticePhysPermeability3D<T, DESCRIPTOR> permeability(this->_blockLattice, this->_converter);
  BlockLatticeVelocity3D<T, DESCRIPTOR> velocity(this->_blockLattice);

  T nu = this->_converter.getPhysViscosity();
  permeability(output,input);
  T K = output[0];
  velocity(output,input);

  output[0] *= -nu / K;
  output[1] *= -nu / K;
  output[2] *= -nu / K;

  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeAverage3D<T, DESCRIPTOR>::BlockLatticeAverage3D(BlockLatticeF3D<T, DESCRIPTOR>& f,
    BlockGeometry3D<T>& blockGeometry, int material, T radius)
  : BlockLatticeF3D<T, DESCRIPTOR>(f.getBlockLattice(), f.getTargetDim()),
    _f(f),
    _blockGeometry(blockGeometry),
    _material(material),
    _radius(radius)
{
  this->getName() = "Average(" + f.getName() + ")";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeAverage3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  //  CuboidGeometry3D<T>& cGeometry = f.getBlockLattice().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice().get_load();

  //  //create boolean indicator functor isInSphere
  //  std::vector<T> center(3,0);
  //  center[0] = (int)cGeometry.get(load.glob(input[0])).get_globPosX() + input[1];
  //  center[1] = (int)cGeometry.get(load.glob(input[0])).get_globPosY() + input[2];
  //  center[2] = (int)cGeometry.get(load.glob(input[0])).get_globPosZ() + input[3];
  //  SphereAnalyticalF3D<bool,T> isInSphere(center,radius);

  // iterate over all cuboids & points and test for material && isInSphere
  output[0]=0;
  //  int numVoxels(0);
  //  if (this->blockGeometry.getMaterial(center[0],center[1],center[2]) == material) {
  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      for (int iX=0; iX<nX; ++iX) {
  //        for (int iY=0; iY<nY; ++iY) {
  //          for (int iZ=0; iZ<nZ; ++iZ) {
  //            std::vector<T> glob(3,0);
  //            glob[0] = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            glob[1] = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            glob[2] = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->blockGeometry.getMaterial(glob[0],glob[1],glob[2]) == material
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
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockEuklidNorm3D<T, DESCRIPTOR>::BlockEuklidNorm3D(BlockF3D<T>& f)
  : BlockF3D<T>(f.getBlockStructure(), 1),
    _f(f)
{
  this->getName() = "EuklidNorm(" + f.getName() + ")";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockEuklidNorm3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();  // flash output, since this methods adds values.
  T data[_f.getTargetDim()];
  _f(data,input);
  for (int i = 0; i < _f.getTargetDim(); ++i) {
    output[0] += data[i] * data[i];
  }
  output[0] = sqrt(output[0]);
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, UnitConverter<T,DESCRIPTOR> const& converter, Cuboid3D<T>* c, int overlap)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _cuboid(c),
    _overlap(overlap)
{
  this->getName() = "BlockLatticeInterpVelocity3D";
}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3D(
  const BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& rhs) :
  BlockLatticePhysF3D<T, DESCRIPTOR>(rhs._blockLattice, rhs._converter, 3),
  _cuboid(rhs._cuboid)
{
}

template<typename T, template<typename U> class DESCRIPTOR>
void BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::operator()(T output[3], const T input[3])
{
  T u[3], rho, volume;
  T d[3], e[3];
  int latIntPos[3] = {0};
  T latPhysPos[3] = {T()};
  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
  _cuboid->getPhysR(latPhysPos, latIntPos);

  T deltaRinv = 1. / _cuboid->getDeltaR();
  d[0] = (input[0] - latPhysPos[0]) * deltaRinv;
  d[1] = (input[1] - latPhysPos[1]) * deltaRinv;
  d[2] = (input[2] - latPhysPos[2]) * deltaRinv;

  e[0] = 1. - d[0];
  e[1] = 1. - d[1];
  e[2] = 1. - d[2];

  latIntPos[0]+=_overlap;
  latIntPos[1]+=_overlap;
  latIntPos[2]+=_overlap;

  this->_blockLattice.get(latIntPos[0], latIntPos[1],
                          latIntPos[2]).computeRhoU(rho, u);
  volume = e[0] * e[1] * e[2];
  output[0] = u[0] * volume;
  output[1] = u[1] * volume;
  output[2] = u[2] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1] + 1,
                          latIntPos[2]).computeRhoU(rho, u);
  volume = e[0] * d[1] * e[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1],
                          latIntPos[2]).computeRhoU(rho, u);
  volume = d[0] * e[1] * e[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1] + 1,
                          latIntPos[2]).computeRhoU(rho, u);
  volume = d[0] * d[1] * e[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1],
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = e[0] * e[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1] + 1,
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = e[0] * d[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1],
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = d[0] * e[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1] + 1,
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = d[0] * d[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  output[0] = this->_converter.getPhysVelocity(output[0]);
  output[1] = this->_converter.getPhysVelocity(output[1]);
  output[2] = this->_converter.getPhysVelocity(output[2]);
}

//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3Degree3D(
//  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, UnitConverter<T>& conv, Cuboid3D<T>* c, int overlap)
//  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
//    _conv(conv),
//    _cuboid(c),
//    _overlap(overlap)
//{
//  this->getName() = "BlockLatticeInterpVelocity3Degree3D";
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3Degree3D(
//  const BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>& rhs) :
//  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
//  _conv(rhs._conv),
//  _cuboid(rhs._cuboid)
//{
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//void BlockLatticeInterpPhysVelocity3Degree3D<T, DESCRIPTOR>::operator()(T output[3], const T input[3])
//{
//  T u[3], rho, volume;
//  int latIntPos[3] = {0};
//  T latPhysPos[3] = {T()};
//  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
//  _cuboid->getPhysR(latPhysPos, latIntPos);
//
//  latIntPos[0]+=_overlap;
//  latIntPos[1] += _overlap;
//  latIntPos[2] += _overlap;
//
//  volume=T(1);
//  for (int i = -1; i <= 2; ++i) {
//    for (int j = -1; j <= 2; ++j) {
//      for (int k = -1; k <= 2; ++k) {
//
//        this->_blockLattice.get(latIntPos[0]+i, latIntPos[1]+j, latIntPos[2]+k).computeRhoU(
//          rho, u);
//        for (int l = -1; l <= 2; ++l) {
//          if (l != i) {
//            volume *= (input[0] - (latPhysPos[0]+ l *_cuboid->getDeltaR()))
//                      / (latPhysPos[0] + i *_cuboid->getDeltaR()
//                         - (latPhysPos[0] + l *_cuboid->getDeltaR()));
//          }
//        }
//        for (int m = -1; m <= 2; ++m) {
//          if (m != j) {
//            volume *= (input[1]
//                       - (latPhysPos[1] + m *_cuboid->getDeltaR()))
//                      / (latPhysPos[1] + j * _cuboid->getDeltaR()
//                         - (latPhysPos[1] + m * _cuboid->getDeltaR()));
//          }
//        }
//        for (int n = -1; n <= 2; ++n) {
//          if (n != k) {
//            volume *= (input[2]
//                       - (latPhysPos[2] + n * _cuboid->getDeltaR()))
//                      / (latPhysPos[2] + k * _cuboid->getDeltaR()
//                         - (latPhysPos[2] + n * _cuboid->getDeltaR()));
//          }
//        }
//        output[0] += u[0] * volume;
//        output[1] += u[1] * volume;
//        output[2] += u[2] * volume;
//        volume=T(1);
//      }
//    }
//  }
//
//  output[0] = _conv.physVelocity(output[0]);
//  output[1] = _conv.physVelocity(output[1]);
//  output[2] = _conv.physVelocity(output[2]);
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpDensity3Degree3D(
//  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, UnitConverter<T>& conv, Cuboid3D<T>* c, int overlap)
//  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
//    _conv(conv),
//    _cuboid(c),
//    _overlap(overlap)
//{
//  this->getName() = "BlockLatticeInterpDensity3Degree3D";
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::BlockLatticeInterpDensity3Degree3D(
//  const BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>& rhs) :
//  BlockLatticeF3D<T, DESCRIPTOR>(rhs._blockLattice, 3),
//  _conv(rhs._conv),
//  _cuboid(rhs._cuboid)
//{
//}
//
//template<typename T, template<typename U> class DESCRIPTOR>
//void BlockLatticeInterpDensity3Degree3D<T, DESCRIPTOR>::operator()(T output[DESCRIPTOR<T>::q], const T input[3])
//{
//  T u[3], rho, volume, density;
//  int latIntPos[3] = {0};
//  T latPhysPos[3] = {T()};
//  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
//  _cuboid->getPhysR(latPhysPos, latIntPos);
//
//  latIntPos[0] += _overlap;
//  latIntPos[1] += _overlap;
//  latIntPos[2] += _overlap;
//
//  for (unsigned iPop = 0; iPop < DESCRIPTOR<T>::q; ++iPop) {
//    volume=T(1);
//    for (int i = -1; i <= 2; ++i) {
//      for (int j = -1; j <= 2; ++j) {
//        for (int k = -1; k <= 2; ++k) {
//
//          density = this->_blockLattice.get(latIntPos[0]+i, latIntPos[1]+j, latIntPos[2]+k).operator[](iPop);
//          for (int l = -1; l <= 2; ++l) {
//            if (l != i) {
//              volume *= (input[0] - (latPhysPos[0]+ l *_cuboid->getDeltaR()))
//                        / (latPhysPos[0] + i *_cuboid->getDeltaR()
//                           - (latPhysPos[0] + l *_cuboid->getDeltaR()));
//            }
//          }
//          for (int m = -1; m <= 2; ++m) {
//            if (m != j) {
//              volume *= (input[1]
//                         - (latPhysPos[1] + m *_cuboid->getDeltaR()))
//                        / (latPhysPos[1] + j * _cuboid->getDeltaR()
//                           - (latPhysPos[1] + m * _cuboid->getDeltaR()));
//            }
//          }
//          for (int n = -1; n <= 2; ++n) {
//            if (n != k) {
//              volume *= (input[2]
//                         - (latPhysPos[2] + n * _cuboid->getDeltaR()))
//                        / (latPhysPos[2] + k * _cuboid->getDeltaR()
//                           - (latPhysPos[2] + n * _cuboid->getDeltaR()));
//            }
//          }
//          output[iPop] += density * volume;
//
//          volume=T(1);
//        }
//      }
//    }
//  }
//}

template<typename T, template<typename U> class DESCRIPTOR>
BlockLatticeMomentumExchangeForce3D<T, DESCRIPTOR>::BlockLatticeMomentumExchangeForce3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>& blockGeometry,
  std::vector<ParticleIndicatorF3D<T,T>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter, int overlap)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 7*indicator.size()), _blockGeometry(blockGeometry), _vectorOfIndicator(indicator), _overlap(overlap)
{
  this->getName() = "physMomentumExchangeForce";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool BlockLatticeMomentumExchangeForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{

  // iterate over all particles in _indicator
  for (typename std::vector<ParticleIndicatorF3D<T,T>* >::size_type iInd=0; iInd!=_vectorOfIndicator.size(); iInd++) {
    int numVoxels = 0;
    // check for intersection of cubiod and indicator
    if (_blockGeometry.getOrigin()[0] <= _vectorOfIndicator[iInd]->getMax()[0]+_vectorOfIndicator[iInd]->getPos()[0]
        && _blockGeometry.getOrigin()[1] <= _vectorOfIndicator[iInd]->getMax()[1]+_vectorOfIndicator[iInd]->getPos()[1]
        && _blockGeometry.getOrigin()[2] <= _vectorOfIndicator[iInd]->getMax()[2]+_vectorOfIndicator[iInd]->getPos()[2]
        && _vectorOfIndicator[iInd]->getMin()[0]+_vectorOfIndicator[iInd]->getPos()[0] <= _blockGeometry.getOrigin()[0] + _blockGeometry.getExtend()[0] * _blockGeometry.getDeltaR()
        && _vectorOfIndicator[iInd]->getMin()[1]+_vectorOfIndicator[iInd]->getPos()[1] <= _blockGeometry.getOrigin()[1] + _blockGeometry.getExtend()[1] * _blockGeometry.getDeltaR()
        && _vectorOfIndicator[iInd]->getMin()[2]+_vectorOfIndicator[iInd]->getPos()[2] <= _blockGeometry.getOrigin()[2] + _blockGeometry.getExtend()[2] * _blockGeometry.getDeltaR() ) {

      // compute size of intersection for iteration
      int start[3] = {0}, span[3] = {0};
      T invDeltaR = 1./_blockGeometry.getDeltaR();
      for (int k=0; k<3; k++) {
        start[k] = int( (_vectorOfIndicator[iInd]->getPos()[k]+_vectorOfIndicator[iInd]->getMin()[k] - _blockGeometry.getOrigin()[k]) * invDeltaR );
        if (start[k] < 0) {
          start[k] = 0;
        }
        span[k] = int( (_vectorOfIndicator[iInd]->getMax()[k] - _vectorOfIndicator[iInd]->getMin()[k])*invDeltaR + 3 );
        if (span[k] + start[k] > _blockGeometry.getExtend()[k]) {
          span[k] = _blockGeometry.getExtend()[k] - start[k];
        }
      }

      // iterate over cells in the constructed intersection box
      for (int iX = start[0]; iX < start[0]+span[0]; iX++) {
        for (int iY = start[1]; iY < start[1]+span[1]; iY++) {
          for (int iZ = start[2]; iZ < start[2]+span[2]; iZ++) {

            // check if cell belongs to particle
            T inside[1] = {0.};
            T posIn[3] = {0.};
            _blockGeometry.getPhysR(posIn, iX, iY, iZ);
            (*(_vectorOfIndicator[iInd]))( inside, posIn);
            if ( !util::nearZero(inside[0]) && this->_blockGeometry.get(iX,iY,iZ)==1) {
              // compute momentum exchange force on particle
              T tmpForce[3] = {0.,0.,0.};
              for (int iPop = 1; iPop < DESCRIPTOR<T>::q ; ++iPop) {
                // Get direction
                const int* c = DESCRIPTOR<T>::c[iPop];
//                T posOut[3] = {0.};
//                T inside2[1] = {0.};
//                _blockGeometry.getPhysR(posOut, iX+c[0], iY+c[1], iZ+c[2]);
//                (*(_vectorOfIndicator[iInd]))( inside2, posOut);
                // if (util::nearZero(inside2[0])) {
                if (this->_blockGeometry.get(iX+c[0], iY+c[1], iZ+c[2])==1) {
                  // Get f_q of next fluid cell where l = opposite(q)
                  T f = this->_blockLattice.get(iX+_overlap + c[0], iY+_overlap + c[1], iZ+_overlap + c[2])[iPop];
                  // Get f_l of the boundary cell
                  // Add f_q and f_opp
                  f += this->_blockLattice.get(iX+_overlap, iY+_overlap, iZ+_overlap)[util::opposite<DESCRIPTOR<T> >(iPop)];
                  // Update force
                  tmpForce[0] -= c[0]*f;
                  tmpForce[1] -= c[1]*f;
                  tmpForce[2] -= c[2]*f;
                }
              }
              // convert force to SI units and compute torque
              numVoxels++;
              tmpForce[0] = this->_converter.getPhysForce(tmpForce[0]);
              tmpForce[1] = this->_converter.getPhysForce(tmpForce[1]);
              tmpForce[2] = this->_converter.getPhysForce(tmpForce[2]);
              output[0+iInd*7] += tmpForce[0];
              output[1+iInd*7] += tmpForce[1];
              output[2+iInd*7] += tmpForce[2];
              output[3+iInd*7] += (posIn[1]-_vectorOfIndicator[iInd]->getPos()[1])*tmpForce[2] - (posIn[2]-_vectorOfIndicator[iInd]->getPos()[2])*tmpForce[1];
              output[4+iInd*7] += (posIn[2]-_vectorOfIndicator[iInd]->getPos()[2])*tmpForce[0] - (posIn[0]-_vectorOfIndicator[iInd]->getPos()[0])*tmpForce[2];
              output[5+iInd*7] += (posIn[0]-_vectorOfIndicator[iInd]->getPos()[0])*tmpForce[1] - (posIn[1]-_vectorOfIndicator[iInd]->getPos()[1])*tmpForce[0];
            }

          }
        }
      }
    }

    output[6+iInd*7] = numVoxels;
  }
  return true;

}
template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
BlockLatticePhysTemperature3D<T,DESCRIPTOR,ThermalDESCRIPTOR>::BlockLatticePhysTemperature3D
(BlockLatticeStructure3D<T,ThermalDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF3D<T,DESCRIPTOR,ThermalDESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physTemperature";
}


template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
bool BlockLatticePhysTemperature3D<T,DESCRIPTOR,ThermalDESCRIPTOR>::operator() (T output[], const int input[])
{
  T latticeTemperature = this->_blockLattice.get( input[0] , input[1] , input[2]).computeRho();
  output[0] = this->_converter.getPhysTemperature(latticeTemperature);

  return true;
}

}  // end namespace olb

#endif
