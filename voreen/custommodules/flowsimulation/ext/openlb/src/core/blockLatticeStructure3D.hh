/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
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
 * Dynamics for a generic 3D block lattice view -- generic implementation.
 */
#ifndef BLOCK_LATTICE_STRUCTURE_3D_HH
#define BLOCK_LATTICE_STRUCTURE_3D_HH


#include "blockLatticeStructure3D.h"
#include "functors/lattice/indicator/blockIndicatorBaseF3D.hh"

namespace olb {

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure3D<T,Lattice>::defineRho(
  BlockGeometryStructure3D<T>& blockGeometry, int material, AnalyticalF3D<T,T>& rho)
{
  T rhoTmp;
  T physR[3] = {T(),T(),T()};
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (blockGeometry.getMaterial(iX, iY, iZ) == material) {
          blockGeometry.getPhysR(physR, iX,iY,iZ);
          rho(&rhoTmp,physR);
          get(iX,iY,iZ).defineRho(rhoTmp);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure3D<T,Lattice>::defineU(
  BlockGeometryStructure3D<T>& blockGeometry, int material, AnalyticalF3D<T,T>& u)
{
  T uTmp[3];
  T physR[3] = {T(),T(),T()};
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (blockGeometry.getMaterial(iX, iY, iZ) == material) {
          blockGeometry.getPhysR(physR, iX,iY,iZ);
          u(uTmp,physR);
          get(iX,iY,iZ).defineU(uTmp);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure3D<T,Lattice>::defineRhoU(
  BlockGeometryStructure3D<T>& blockGeometry, int material,
  AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  T rhoTmp;
  T uTmp[3];
  T physR[3] = {T(),T(),T()};
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (blockGeometry.getMaterial(iX, iY, iZ) == material) {
          blockGeometry.getPhysR(physR, iX,iY,iZ);
          rho(&rhoTmp,physR);
          u(uTmp,physR);
          get(iX,iY,iZ).defineRhoU(rhoTmp,uTmp);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure3D<T,Lattice>::definePopulations(
  BlockGeometryStructure3D<T>& blockGeometry, int material, AnalyticalF3D<T,T>& Pop)
{
  T physR[3] = {T(),T(),T()};
  T PopTmp[Lattice<T>::q];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (blockGeometry.getMaterial(iX, iY, iZ) == material) {
          blockGeometry.getPhysR(physR, iX,iY,iZ);
          Pop(PopTmp,physR);
          get(iX,iY,iZ).definePopulations(PopTmp);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure3D<T,Lattice>::defineExternalField(
  BlockGeometryStructure3D<T>& blockGeometry, int material, int fieldBeginsAt,
  int sizeOfField, AnalyticalF3D<T,T>& field)
{
  T physR[3] = {T(),T(),T()};
  T* fieldTmp = new T [sizeOfField];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (blockGeometry.getMaterial(iX, iY, iZ) == material) {
          blockGeometry.getPhysR(physR, iX,iY,iZ);
          field(fieldTmp,physR);
          get(iX,iY,iZ).defineExternalField(fieldBeginsAt,sizeOfField, fieldTmp);
        }
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure3D<T,Lattice>::defineExternalField(
  BlockIndicatorF3D<T>& indicator, int overlap,
  int fieldBeginsAt, int sizeOfField, AnalyticalF3D<T,T>& field)
{
  T physR[3] = {T(),T(),T()};
  T* fieldTmp = new T [sizeOfField];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        const int blockLocation[3] = { iX-overlap, iY-overlap, iZ-overlap };
        if (indicator(blockLocation)) {
          indicator.getBlockGeometry().getPhysR(physR,iX,iY,iZ);
          field(fieldTmp,physR);
          get(iX,iY,iZ).defineExternalField(fieldBeginsAt, sizeOfField, fieldTmp);
        }
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure3D<T,Lattice>::defineExternalField(
  BlockGeometryStructure3D<T>& blockGeometry, IndicatorF3D<T>& indicator, int fieldBeginsAt,
  int sizeOfField, AnalyticalF3D<T,T>& field)
{
  bool inside;
  T physR[3] = {T(),T(),T()};
  T* fieldTmp = new T [sizeOfField];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        blockGeometry.getPhysR(physR, iX,iY,iZ);
        indicator(&inside, physR);
        if (inside) {
          field(fieldTmp,physR);
          get(iX,iY,iZ).defineExternalField(fieldBeginsAt,sizeOfField, fieldTmp);
        }
      }
    }
  }
  delete[] fieldTmp;
}


template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure3D<T, Lattice>::setExternalParticleField(
  BlockGeometryStructure3D<T>& blockGeometry,
  AnalyticalF3D<T, T>& velocity, ParticleIndicatorF3D<T, T>& sIndicator)
{
  T foo[4] = { T(), T(), T(), T() }; /// Contains foo[0]=vel0; foo[1]=vel1; foo[2]=vel2; foo[3]=porosity
  T physR[3] = { T(), T(), T() };
  T porosity[1] = { T() };
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      for (int iZ = 0; iZ < this->_nz; iZ++) {
        blockGeometry.getPhysR(physR, iX, iY, iZ);
//        if (physR[0] > sIndicator.getMin()[0]
//            && physR[0] < sIndicator.getMax()[0]
//            && physR[1] > sIndicator.getMin()[1]
//            && physR[1] < sIndicator.getMax()[1]) {
        // TODO quick and dirty
        if (physR[0] < sIndicator.getPos()[0]+sIndicator.getMax()[0] &&
            physR[1] < sIndicator.getPos()[1]+sIndicator.getMax()[1] &&
            physR[2] < sIndicator.getPos()[2]+sIndicator.getMax()[2] &&
            physR[0] > sIndicator.getPos()[0]+sIndicator.getMin()[0] &&
            physR[1] > sIndicator.getPos()[1]+sIndicator.getMin()[1] &&
            physR[2] > sIndicator.getPos()[2]+sIndicator.getMin()[2]
           ) {
          sIndicator(porosity, physR);
//          std::cout << porosity[0] << std::endl;
          if (porosity[0] > 0.) {
            velocity(foo, physR);
            foo[0] *= porosity[0];
            foo[1] *= porosity[0];
            foo[2] *= porosity[0];
            foo[3] = porosity[0];
            get(iX, iY, iZ).addExternalField(1, 4, foo);
            porosity[0] = 1. - porosity[0];
            get(iX, iY, iZ).multiplyExternalField(0, 1, porosity);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure3D<T,Lattice>::iniEquilibrium(
  BlockGeometryStructure3D<T>& blockGeometry, int material,
  AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  T physR[3] = {T(),T(),T()};
  T uTmp[3] = {T(),T(),T()};
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (blockGeometry.getMaterial(iX, iY, iZ) == material) {
          blockGeometry.getPhysR(physR,iX,iY,iZ);
          u(uTmp,physR);
          rho(&rhoTmp,physR);
          get(iX,iY,iZ).iniEquilibrium(rhoTmp,uTmp);
        }
      }
    }
  }
}

}  // namespace olb

#endif
