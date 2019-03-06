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
 * Dynamics for a generic 2D block structure -- header file.
 */
#ifndef BLOCK_LATTICE_STRUCTURE_2D_HH
#define BLOCK_LATTICE_STRUCTURE_2D_HH

#include <vector>
#include "blockLatticeStructure2D.h"


namespace olb {


template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineDynamics(
  BlockGeometryStructure2D<T>& blockGeometry, int material,
  Dynamics<T,Lattice>* dynamics)
{
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        get(iX,iY).defineDynamics(dynamics);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineRho(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& rho)
{
  T rhoTmp;
  T physR[2]= {T(),T()};
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        blockGeometry.getPhysR(physR, iX,iY);
        rho(&rhoTmp,physR);
        get(iX,iY).defineRho(rhoTmp);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineU(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& u)
{
  T uTmp[2];
  T physR[2]= {T(),T()};
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        blockGeometry.getPhysR(physR, iX,iY);
        u(uTmp,physR);
        get(iX,iY).defineU(uTmp);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineRhoU(BlockGeometryStructure2D<T>& blockGeometry,
    int material, AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  T rhoTmp;
  T uTmp[2];
  T physR[2]= {T(),T()};
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        blockGeometry.getPhysR(physR, iX,iY);
        rho(&rhoTmp,physR);
        u(uTmp,physR);
        get(iX,iY).defineRhoU(rhoTmp,uTmp);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::definePopulations(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& Pop)
{
  T physR[2]= {T(),T()};
  T PopTmp[Lattice<T>::q];
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        blockGeometry.getPhysR(physR, iX,iY);
        Pop(PopTmp,physR);
        get(iX,iY).definePopulations(PopTmp);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineExternalField(
  BlockGeometryStructure2D<T>& blockGeometry, int material, int fieldBeginsAt,
  int sizeOfField, AnalyticalF2D<T,T>& field)
{
  T physR[2]= {T(),T()};
  T* fieldTmp = new T [sizeOfField];
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        blockGeometry.getPhysR(physR, iX,iY);
        field(fieldTmp,physR);
        get(iX,iY).defineExternalField(fieldBeginsAt,sizeOfField, fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::defineExternalField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator, int fieldBeginsAt,
  int sizeOfField, AnalyticalF2D<T,T>& field)
{
  T* fieldTmp = new T [sizeOfField];
  bool inside;
  T physR[2]= {T(),T()};
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      blockGeometry.getPhysR(physR, iX,iY);
      indicator(&inside, &physR[0]);
      if (inside) {
        field(fieldTmp,physR);
        get(iX,iY).defineExternalField(fieldBeginsAt,sizeOfField, fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::addExternalField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator, int fieldBeginsAt,
  int sizeOfField, AnalyticalF2D<T,T>& field)
{
  T* fieldTmp = new T [sizeOfField];
  bool inside;
  T physR[2]= {T(),T()};
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      blockGeometry.getPhysR(physR, iX,iY);
      indicator(&inside, &(physR[0]));
      if (inside) {
        field(fieldTmp,physR);
        get(iX,iY).addExternalField(fieldBeginsAt,sizeOfField, fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::addExternalField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator, int fieldBeginsAt,
  int sizeOfField, AnalyticalF2D<T,T>& field, AnalyticalF2D<T,T>& porous)
{
  T* fieldTmp = new T [sizeOfField];
  bool inside;
  T physR[2]= {T(),T()};
  T porousA[1] = {T()};
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      blockGeometry.getPhysR(physR, iX,iY);
      indicator(&inside, physR);
      if (inside) {
        porous(porousA, physR);
        field(fieldTmp,physR);
        for (int i = 0; i<sizeOfField; ++i) {
          fieldTmp[i] *= porousA[0];
        }
        get(iX,iY).addExternalField(fieldBeginsAt,sizeOfField, fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::resetExternalParticleField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator)
{
  T physR[2]= {T(),T()};
  T fieldTmp[4] = {T(1), T(), T(), T()};
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      blockGeometry.getPhysR(physR, iX,iY);
      if (indicator.getMin()[0] < physR[0] &&
          indicator.getMin()[1] < physR[1] &&
          indicator.getMax()[0] > physR[0] &&
          indicator.getMax()[1] > physR[1]) {
        get(iX,iY).defineExternalField(0,4, fieldTmp);
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::setExternalParticleField(BlockGeometryStructure2D<T>& blockGeometry, AnalyticalF2D<T,T>& velocity, ParticleIndicatorF2D<T,T>& sIndicator)
{

  // check for intersection of cuboid and indicator
  if (blockGeometry.getOrigin()[0] <= sIndicator.getCircumRadius() + sIndicator.getPos()[0]
      && blockGeometry.getOrigin()[1] <= sIndicator.getCircumRadius() + sIndicator.getPos()[1]
      && sIndicator.getPos()[0] - sIndicator.getCircumRadius() <= blockGeometry.getOrigin()[0] + blockGeometry.getExtend()[0] * blockGeometry.getDeltaR()
      && sIndicator.getPos()[1] - sIndicator.getCircumRadius() <= blockGeometry.getOrigin()[1] + blockGeometry.getExtend()[1] * blockGeometry.getDeltaR()) {

    int start[2] = {0}, span[2] = {0};
    T invDeltaR = 1./blockGeometry.getDeltaR();
    for (int k=0; k<2; k++) {
      start[k] = int( (sIndicator.getPos()[k]-sIndicator.getCircumRadius() - blockGeometry.getOrigin()[k]) * invDeltaR );
      if (start[k] < 0) {
        start[k] = 0;
      }
      span[k] = int( (2.*sIndicator.getCircumRadius())*invDeltaR + 3 );
      if (span[k] + start[k] > blockGeometry.getExtend()[k]) {
        span[k] = blockGeometry.getExtend()[k] - start[k];
      } else if ( (sIndicator.getPos()[k] - sIndicator.getCircumRadius() - blockGeometry.getOrigin()[k]) < 0. ) {
        span[k] = (sIndicator.getPos()[k] + sIndicator.getCircumRadius() - blockGeometry.getOrigin()[k]) * invDeltaR + 3;
      }
    }

    T foo[3] = {T(), T(), T()}; /// Contains foo[0]=vel0; foo[1]=vel1; foo[2]=porosity
    T physR[2]= {T(),T()};
    T porosity[1] = {T()};
    for (int iX = start[0]; iX < start[0]+span[0]; iX++) {
      for (int iY = start[1]; iY < start[1]+span[1]; iY++) {
        blockGeometry.getPhysR(physR, iX,iY);
        sIndicator(porosity, physR);
        if (!util::nearZero(porosity[0])) {
          velocity(foo,physR);
          foo[0] *= porosity[0];
          foo[1] *= porosity[0];
          foo[2] = porosity[0];
          get(iX,iY).addExternalField(1,3, foo);
          porosity[0] = 1.-porosity[0];
          get(iX,iY).multiplyExternalField(0,1, porosity);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::multiplyExternalField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator, int fieldBeginsAt,
  int sizeOfField, AnalyticalF2D<T,T>& field)
{
  T* fieldTmp = new T [sizeOfField];
  bool inside;
  T physR[3]= {T(),T(),T()};
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      blockGeometry.getPhysR(physR, iX,iY);
      indicator(&inside, &(physR[0]));
      if (inside) {
        field(fieldTmp,physR);
        get(iX,iY).multiplyExternalField(fieldBeginsAt,sizeOfField, fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeStructure2D<T,Lattice>::iniEquilibrium(
  BlockGeometryStructure2D<T>& blockGeometry, int material,
  AnalyticalF2D<T,T>& rho , AnalyticalF2D<T,T>& u)
{
  T physR[2]= {T(),T()};
  for (int iX = 0; iX < this->_nx; iX++) {
    for (int iY = 0; iY < this->_ny; iY++) {
      if (blockGeometry.getMaterial(iX, iY) == material) {
        blockGeometry.getPhysR(physR,iX,iY);
        T uTmp[] = {T(),T()};
        u(uTmp,physR);
        T rhoTmp = T();
        rho(&rhoTmp,physR);
        get(iX,iY).iniEquilibrium(rhoTmp,uTmp);
      }
    }
  }
}

}  // namespace olb

#endif
