/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2018 Jonas Latt, Adrian Kummerlaender
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

#include "blockLatticeStructure2D.h"
#include "functors/lattice/indicator/blockIndicatorBaseF2D.hh"
#include "functors/lattice/indicator/blockIndicatorF2D.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI2
#define M_PI2 1.57079632679489661923
#endif

namespace olb {


template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineRho(
  BlockIndicatorF2D<T>& indicator, AnalyticalF<2,T,T>& rho)
{
  T physR[2] = { };
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        rho(&rhoTmp, physR);
        get(iX, iY).defineRho(rhoTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineRho(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF<2,T,T>& rho)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  defineRho(indicator, rho);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineU(
  BlockIndicatorF2D<T>& indicator, AnalyticalF<2,T,T>& u)
{
  T physR[2] = { };
  T uTmp[2] = { };
  const auto& geometry = indicator.getBlockGeometryStructure();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        geometry.getPhysR(physR, iX, iY);
        u(uTmp, physR);
        get(iX, iY).defineU(uTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineU(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF<2,T,T>& u)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  defineU(indicator, u);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineRhoU(
  BlockIndicatorF2D<T>& indicator,
  AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u)
{
  T physR[2] = { };
  T uTmp[2] = { };
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        rho(&rhoTmp, physR);
        u(uTmp, physR);
        get(iX, iY).defineRhoU(rhoTmp, uTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineRhoU(
  BlockGeometryStructure2D<T>& blockGeometry, int material,
  AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  defineRhoU(indicator, rho, u);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::definePopulations(
  BlockIndicatorF2D<T>& indicator, AnalyticalF<2,T,T>& Pop)
{
  T physR[2] = { };
  T PopTmp[DESCRIPTOR::q];
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        Pop(PopTmp, physR);
        get(iX, iY).definePopulations(PopTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::definePopulations(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF<2,T,T>& Pop)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  definePopulations(indicator, Pop);
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineField(
  BlockIndicatorF2D<T>& indicator, AnalyticalF<2,T,T>& field)
{
  T* fieldTmp = new T[DESCRIPTOR::template size<FIELD>()];
  T physR[2] = { };
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        field(fieldTmp, physR);
        get(iX, iY).template defineField<FIELD>(fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineField(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF<2,T,T>& field)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  defineField<FIELD>(indicator, field);
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicatorF,
  AnalyticalF<2,T,T>& field)
{
  BlockIndicatorFfromIndicatorF2D<T> indicator(indicatorF, blockGeometry);
  defineField<FIELD>(indicator, field);
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::addField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator,
  AnalyticalF<2,T,T>& field)
{
  T* fieldTmp = new T[DESCRIPTOR::template size<FIELD>()];
  T physR[2] = { };
  bool inside;
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      blockGeometry.getPhysR(physR, iX, iY);
      indicator(&inside, &(physR[0]));
      if (inside) {
        field(fieldTmp,physR);
        get(iX, iY).template addField<FIELD>(fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::addField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator,
  AnalyticalF<2,T,T>& field, AnalyticalF<2,T,T>& porous)
{
  T* fieldTmp = new T[DESCRIPTOR::template size<FIELD>()];
  bool inside;
  T physR[2] = { };
  T porousA[1] = { };
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      blockGeometry.getPhysR(physR, iX, iY);
      indicator(&inside, physR);
      if (inside) {
        porous(porousA, physR);
        field(fieldTmp,physR);
        for (int i = 0; i < DESCRIPTOR::template size<FIELD>(); ++i) {
          fieldTmp[i] *= porousA[0];
        }
        get(iX, iY).template addField<FIELD>(fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}



template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::multiplyField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator,
  AnalyticalF<2,T,T>& field)
{
  T* fieldTmp = new T [DESCRIPTOR::template size<FIELD>()];
  bool inside;
  T physR[3] = { };
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      blockGeometry.getPhysR(physR, iX, iY);
      indicator(&inside, &(physR[0]));
      if (inside) {
        field(fieldTmp, physR);
        get(iX, iY).template multiplyField<FIELD>(fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::iniEquilibrium(
  BlockIndicatorF2D<T>& indicator,
  AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u)
{
  T physR[2] = { };
  T uTmp[2] = { };
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        u(uTmp, physR);
        rho(&rhoTmp, physR);
        get(iX, iY).iniEquilibrium(rhoTmp, uTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::iniEquilibrium(
  BlockGeometryStructure2D<T>& blockGeometry, int material,
  AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  iniEquilibrium(indicator, rho, u);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::iniRegularized(
  BlockIndicatorF2D<T>& indicator,
  AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u, AnalyticalF<2,T,T>& pi)
{
  T physR[2] = { };
  T uTmp[2] = { };
  T rhoTmp = T();
  T piTmp[util::TensorVal<DESCRIPTOR>::n] = { };
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        u(uTmp, physR);
        rho(&rhoTmp, physR);
        pi(piTmp, physR);
        get(iX, iY).iniRegularized(rhoTmp, uTmp, piTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::iniRegularized(
  BlockGeometryStructure2D<T>& blockGeometry, int material,
  AnalyticalF<2,T,T>& rho, AnalyticalF<2,T,T>& u, AnalyticalF<2,T,T>& pi)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  iniRegularized(indicator, rho, u, pi);
}



////////// FREE FUNCTIONS //////////

template <typename T>
bool getRangeBlockGeometrySmoothIndicatorIntersection2D(BlockGeometryStructure2D<T>& blockGeometry, SmoothIndicatorF2D<T,T,true>& sIndicator, T invDeltaX, std::vector<int>& start, std::vector<int>& end)
{
  T posXmin = sIndicator.getPos()[0] - sIndicator.getCircumRadius();
  T posXmax = sIndicator.getPos()[0] + sIndicator.getCircumRadius();
  T posYmin = sIndicator.getPos()[1] - sIndicator.getCircumRadius();
  T posYmax = sIndicator.getPos()[1] + sIndicator.getCircumRadius();

  start[0] = invDeltaX*(posXmin-blockGeometry.getOrigin()[0])+1;
  end[0] = std::ceil(invDeltaX*(posXmax-blockGeometry.getOrigin()[0]))-1;
  start[1] = invDeltaX*(posYmin-blockGeometry.getOrigin()[1])+1;
  end[1] = std::ceil(invDeltaX*(posYmax-blockGeometry.getOrigin()[1]))-1;

  start[0] = std::max(start[0],0);
  end[0] = std::min(end[0],blockGeometry.getNx()-1);
  start[1] = std::max(start[1],0);
  end[1] = std::min(end[1],blockGeometry.getNy()-1);

  if (!(end[0]>=start[0] && end[1]>=start[1])) {
    // for no intersection set size of range to -1
    start[0] = 1;
    end[0] = 0;
    start[1] = 1;
    end[1] = 0;
    return false;
  } else {
    // increase range of intersection by layer of one within
    // block boundaries to ensure complete object is inside
    for (int k=0; k<2; k++) {
      start[k] -= 1;
      if (start[k] < 0) {
        start[k] = 0;
      }
      end[k] += 2;
      if (end[k] > blockGeometry.getExtend()[k]) {
        end[k] = blockGeometry.getExtend()[k];
      }
    }
    return true;
  }
}

template<typename T>
void checkSmoothIndicatorOutOfGeometry( bool& outOfGeometry, Vector<T,2>& ghostPos,
                                        SmoothIndicatorF2D<T,T,true>& sIndicator, Vector<T,2> cellMin,
                                        Vector<T,2> cellMax, Vector<bool,2> periodic)
{
  Vector<bool,2> dir = Vector<bool,2> (false,false);
  for (int i=0; i<2; i++) {
    T posMin = sIndicator.getPos()[i] - sIndicator.getCircumRadius();
    T posMax = sIndicator.getPos()[i] + sIndicator.getCircumRadius();
    if (posMin < cellMin[i] && periodic[i]) {
      outOfGeometry = true;
      dir[i] = true;
      ghostPos[i] = cellMax[i] - (cellMin[i] - sIndicator.getPos()[i]);
    } else if (posMax > cellMax[i] && periodic[i]) {
      outOfGeometry = true;
      dir[i] = true;
      ghostPos[i] = cellMin[i] + (sIndicator.getPos()[i] - cellMax[i]);
    }
    if (!dir[i]) {
      ghostPos[i] = sIndicator.getPos()[i];
    }
  }
}

template<typename T, typename DESCRIPTOR>
void setBlockExternalParticleField(BlockGeometryStructure2D<T>& blockGeometry, AnalyticalF<2,T,T>& velocity, SmoothIndicatorF2D<T,T,true>& sIndicator, BlockLattice2D<T,DESCRIPTOR>& extendedBlockLattice)
{

  std::vector<int> start{0,0};
  std::vector<int> end{0,0};
  // check for intersection of cuboid and indicator
  if (getRangeBlockGeometrySmoothIndicatorIntersection2D(blockGeometry, sIndicator, 1./blockGeometry.getDeltaR(), start, end)) {
    T foo[3] = { }; /// Contains foo[0]=vel0; foo[1]=vel1; foo[2]=porosity
    T physR[2]= { };
    T porosity[1] = { };
    for (int iX = start[0]; iX < end[0]; ++iX) {
      for (int iY = start[1]; iY < end[1]; ++iY) {
        blockGeometry.getPhysR(physR, iX, iY);
        sIndicator(porosity, physR);
        if (!util::nearZero(porosity[0])) {
          // TODO: Check / adapt to use descriptor fields
          velocity(foo,physR);
          foo[0] *= porosity[0];
          foo[1] *= porosity[0];
          foo[2] = porosity[0];
          extendedBlockLattice.get(iX, iY).template addField<descriptors::VELOCITY_NUMERATOR>(foo);
          extendedBlockLattice.get(iX, iY).template addField<descriptors::VELOCITY_DENOMINATOR>(&foo[2]);
          porosity[0] = 1.-porosity[0];
          extendedBlockLattice.get(iX, iY).template getFieldPointer<descriptors::POROSITY>()[0] *= porosity[0];
        }
      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
void setBlockExternalParticleField( BlockGeometryStructure2D<T>& blockGeometry, AnalyticalF<2,T,T>& velocity,
                                    SmoothIndicatorF2D<T,T,true>& sIndicator, BlockLattice2D<T,DESCRIPTOR>& extendedBlockLattice,
                                    Vector<T,2> cellMin, Vector<T,2> cellMax, Vector<bool,2> periodic)
{
  // Checking if the Particle leaves the Domain
  bool outOfGeometry = false;
  Vector<T,2> ghostPos = Vector<T,2> (0.,0.);
  checkSmoothIndicatorOutOfGeometry(outOfGeometry, ghostPos, sIndicator, cellMin, cellMax, periodic);

  //Do the normal routine if the particle is in the geometry
  if (!outOfGeometry) {
    setBlockExternalParticleField( blockGeometry, velocity, sIndicator, extendedBlockLattice);
  } else {
    //sets the Particle to ghost position on the other side of the domain and sets the field
    Vector<T,2> particlePositionTmp = sIndicator.getPos();
    sIndicator.setPos(ghostPos);
    setBlockExternalParticleField( blockGeometry, velocity, sIndicator, extendedBlockLattice);
    //Reverting Particle to its Previous position and setting the field
    sIndicator.setPos(particlePositionTmp);
    setBlockExternalParticleField( blockGeometry, velocity, sIndicator, extendedBlockLattice);
  }
}

//Geng2019
template<typename T, typename DESCRIPTOR>
void setBlockZetaParticleField( BlockGeometryStructure2D<T>& blockGeometry, AnalyticalF<2,T,T>& velocity,
                                SmoothIndicatorF2D<T,T,true>& sIndicator,
                                BlockLattice2D<T,DESCRIPTOR>& extendedBlockLattice )
{

  int start[2] = {0};
  int end[2] = {0};
  // check for intersection of cuboid and indicator
  Cuboid2D<T> tmpCuboid(blockGeometry.getOrigin()[0], blockGeometry.getOrigin()[1], blockGeometry.getDeltaR(), blockGeometry.getNx(), blockGeometry.getNy());
  T posXmin = sIndicator.getPos()[0] - (sIndicator.getCircumRadius()+30);
  T posXmax = sIndicator.getPos()[0] + (sIndicator.getCircumRadius()+30);
  T posYmin = sIndicator.getPos()[1] - (sIndicator.getCircumRadius()+30);
  T posYmax = sIndicator.getPos()[1] + (sIndicator.getCircumRadius()+30);

  if (tmpCuboid.checkInters(posXmin, posXmax, posYmin, posYmax, start[0], end[0], start[1], end[1])) {

    for (int k=0; k<2; k++) {
      start[k] -= 1;
      if (start[k] < 0) {
        start[k] = 0;
      }
      end[k] += 2;
      if (end[k] > blockGeometry.getExtend()[k]) {
        end[k] = blockGeometry.getExtend()[k];
      }
    }

    T foo[3] = { }; /// Contains foo[0]=vel0; foo[1]=vel1; foo[2]=porosity
    T physR[2]= { };
    T porosity[1] = { };
    for (int iX = start[0]; iX < end[0]; ++iX) {
      for (int iY = start[1]; iY < end[1]; ++iY) {
        blockGeometry.getPhysR(physR, iX, iY);
        sIndicator(porosity, physR);
        if (!util::nearZero(porosity[0])) {
          // TODO: Check / adapt to use descriptor fields
          velocity(foo,physR);
          if (porosity[0]<0.5) {
            foo[0] = 0.;
            foo[1] = 0.;
          }

          foo[2] = porosity[0];
          extendedBlockLattice.get(iX, iY).template addField<descriptors::VELOCITY_NUMERATOR>(foo);
          extendedBlockLattice.get(iX, iY).template addField<descriptors::VELOCITY_DENOMINATOR>(&foo[2]);
          porosity[0] = 1.-porosity[0];


          double psi= min(porosity[0],extendedBlockLattice.get(iX, iY).template getField<descriptors::POROSITY>());
          if (psi < extendedBlockLattice.get(iX, iY).template getField<descriptors::POROSITY>()) {
            //note that xdis and ydis could be 0;
            double xdis=physR[0]-sIndicator.getPos()[0];
            double ydis=physR[1]-sIndicator.getPos()[1];
            double distToCenter=std::max(std::sqrt(std::pow(xdis,2)+std::pow(ydis,2)),1e-18);
            xdis/=distToCenter;
            ydis/=distToCenter;
            for (int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
              extendedBlockLattice.get(iX, iY).template getFieldPointer<descriptors::ZETA>()[tmp_iPop] = T(2.*(1.-psi)*abs(xdis*descriptors::c<DESCRIPTOR>(tmp_iPop,0)+ydis*descriptors::c<DESCRIPTOR>(tmp_iPop,1))/(sIndicator.getEpsilon()/blockGeometry.getDeltaR()));
            }
          }
          extendedBlockLattice.get(iX, iY).template setField<descriptors::POROSITY>(psi);
        }
      }
    }
  }
}

}  // namespace olb

#endif
