/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause, Adrian Kummerländer
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

#ifndef SUPER_LATTICE_INTEGRAL_F_2D_HH
#define SUPER_LATTICE_INTEGRAL_F_2D_HH

#include<vector>
#include<cmath>

#include "superLatticeIntegralF2D.h"
#include "blockLatticeIntegralF2D.h"
#include "indicator/superIndicatorF2D.h"
#include "utilities/vectorHelpers.h"
#include "core/vector.h"

using namespace olb::util;

namespace olb {


template <typename T>
SuperMax2D<T>::SuperMax2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Max("+_f.getName()+")";
}

template <typename T>
bool SuperMax2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i]=T();
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          if (_superGeometry.get(load.glob(iC), iX, iY) == _material) {
            T outputTmp[_f.getTargetDim()];
            _f(outputTmp,load.glob(iC),iX,iY);
            if (fabs(outputTmp[i]) > output[i]) {
              output[i] = fabs(outputTmp[i]);
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(output[i], MPI_MAX);
#endif
  }
  return true;
}

template <typename T>
SuperMin2D<T>::SuperMin2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Min("+_f.getName()+")";
}

template <typename T>
bool SuperMin2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i]=T();
    for (int iC = 0; iC < load.size(); ++iC) {
      int nX = cGeometry.get(load.glob(iC)).getNx();
      int nY = cGeometry.get(load.glob(iC)).getNy();
      for (int iX = 0; iX < nX; iX++) {
        for (int iY = 0; iY < nY; iY++) {
          if (_superGeometry.get(load.glob(iC), iX, iY) == _material) {
            T outputTmp[_f.getTargetDim()];
            _f(outputTmp,load.glob(iC),iX,iY);
            if (fabs(outputTmp[i]) < output[i]) {
              output[i] = fabs(outputTmp[i]);
            }
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(output[i], MPI_MAX);
#endif
  }
  return true;
}


template <typename T>
SuperSum2D<T>::SuperSum2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()+1),
    _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Sum("+_f.getName()+")";
}

template <typename T>
bool SuperSum2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  CuboidGeometry2D<T>& cGeometry = _f.getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  int numVoxels(0);
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i]=T();
  }
  for (int iC=0; iC<load.size(); ++iC) {
    int nX = cGeometry.get(load.glob(iC)).getNx();
    int nY = cGeometry.get(load.glob(iC)).getNy();
    for (int iX = 0; iX < nX; iX++) {
      for (int iY = 0; iY < nY; iY++) {
        if (this->_superGeometry.get(load.glob(iC), iX, iY) == _material) {
          T outputTmp[_f.getTargetDim()];
          _f(outputTmp,load.glob(iC),iX,iY);
          for (int i = 0; i < this->getTargetDim()-1 /*f.getTargetDim()*/; ++i) {
            output[i] += outputTmp[i];
          }
          numVoxels++;
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim()-1; ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  output[this->getTargetDim()-1] = numVoxels;
  return true;
}


template <typename T>
SuperSumIndicator2D<T>::SuperSumIndicator2D(SuperF2D<T>& f,
    SuperGeometry2D<T>& superGeometry, ParticleIndicatorF2D<T,T>& indicator)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()+1),
    _f(f), _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "Sum("+_f.getName()+")";
}

template <typename T>
bool SuperSumIndicator2D<T>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = 0.;
  }

  T physR[2];
  T inside[1];
  int numVoxels(0);
  T outputTmp[_f.getTargetDim()];
  Cuboid2D<T>* cub = nullptr;
  int start[2] = {0}, span[2] = {0};
  for (int iC = 0; iC < load.size(); ++iC) {
    int globiC = load.glob(iC);
    cub = &_superGeometry.getCuboidGeometry().get(globiC);
    if (! (cub->get_globPosX() > _indicator.getPos()[0]+_indicator.getCircumRadius() ||
           cub->get_globPosY() > _indicator.getPos()[1]+_indicator.getCircumRadius() ||
           _indicator.getPos()[0]-_indicator.getCircumRadius() > cub->get_globPosX() + cub->getExtend()[0] * cub->getDeltaR() ||
           _indicator.getPos()[1]-_indicator.getCircumRadius() > cub->get_globPosY() + cub->getExtend()[1] * cub->getDeltaR())) {
      for (int k=0; k<2; k++) {
        start[k] = int( (_indicator.getPos()[k]-_indicator.getCircumRadius() - cub->getOrigin()[k]) /
                        cub->getDeltaR() );
        if (start[k] < 0) {
          start[k] = 0;
        }
        span[k] = int( (2.*_indicator.getCircumRadius())/cub->getDeltaR() + 3 );
        if (span[k] + start[k] > cub->getExtend()[k]) {
          span[k] = int( cub->getExtend()[k] - start[k] );
        }
      }
      for (int iX = start[0]; iX < start[0]+span[0]; iX++) {
        for (int iY = start[1]; iY < start[1]+span[1]; iY++) {
          if (_superGeometry.get(globiC, iX, iY) == 1) {
            cub->getPhysR(physR,iX,iY);
            _indicator(inside, physR);
            if ( !util::nearZero(inside[0]) ) {
              _f(outputTmp,globiC,iX,iY);
              for (int iDim = 0; iDim < this->getTargetDim()-1; ++iDim) {
                output[iDim] += outputTmp[iDim];
              }
              numVoxels++;
            }
          }
        }
      }
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim()-1; ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
#endif
  output[this->getTargetDim()-1] = numVoxels;
  return true;
}


template <typename T>
SuperIntegral2D<T>::SuperIntegral2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry, const int material)
  : SuperF2D<T>(f.getSuperStructure(),f.getTargetDim()), _f(f), _superGeometry(superGeometry), _material(material)
{
  this->getName() = "Integral("+_f.getName()+")";
}

template <typename T>
bool SuperIntegral2D<T>::operator() (T output[], const int input[])
{
  //  f.getSuperStructure().communicate();
  //  CuboidGeometry2D<T>& cGeometry = f.getSuperStructure().getCuboidGeometry();
  //  LoadBalancer<T>& load = f.getSuperStructure().getLoadBalancer();

  //  std::vector<T> tmp(this->_n, T() );
  //  for (int i=0; i<this->_n; ++i) {
  //    for (int iC=0; iC<load.size(); ++iC) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  ////      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      T weight = pow(this->superGeometry.getDeltaR(),3);
  //      for (int iX=0; iX<nX; iX++) {
  //        for (int iY=0; iY<nY; iY++) {
  ////          for (int iZ=0; iZ<nZ; iZ++) {
  //            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  ////            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  ////            if (this->superGeometry.getMaterial(globX, globY) == material) {
  ////              tmp[i]+=f(load.glob(iC),iX,iY)[i]*weight;
  ////            }
  ////            if (this->superGeometry.getMaterial(globX, globY, globZ) == material) {
  ////              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*weight;
  ////            }

  ////          }
  //        }
  //      }
  //    }
  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
  //#endif
  //  }
  //  return tmp;
  return false;
}


template <typename T>
template<template<typename U> class DESCRIPTOR>
SuperGeometryFaces2D<T>::SuperGeometryFaces2D(SuperGeometry2D<T>& superGeometry,
    const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : GenericF<T,int>(7,3), _superGeometry(superGeometry), _material(material), _latticeL(converter.getConversionFactorLength())
{
  this->getName() = "superGeometryFaces";
}

template <typename T>
SuperGeometryFaces2D<T>::SuperGeometryFaces2D(SuperGeometry2D<T>& superGeometry,
    const int material, T latticeL)
  : GenericF<T,int>(7,3), _superGeometry(superGeometry), _material(material), _latticeL(latticeL)
{
  this->getName() = "superGeometryFaces";
}



template <typename T>
bool SuperGeometryFaces2D<T>::operator() (T output[], const int input[])
{
  for (int iDim = 0; iDim < 7; ++iDim) {
    output[iDim]=T();
  }
  _superGeometry.communicate();
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    BlockGeometryFaces2D<T> f(_superGeometry.getBlockGeometry(iC), _material, _latticeL);
    T outputTmp[f.getTargetDim()];
    f(outputTmp,input);
    for (int iDim = 0; iDim < 7; ++iDim) {
      output[iDim] += outputTmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 7; ++iDim) {
    singleton::mpi().reduceAndBcast( output[iDim], MPI_SUM);
  }
#endif
  return true;
}


template <typename T>
template<template<typename U> class DESCRIPTOR>
SuperGeometryFacesIndicator2D<T>::SuperGeometryFacesIndicator2D(SuperGeometry2D<T>& superGeometry,
    SmoothIndicatorCircle2D<T,T>& indicator, const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : GenericF<T,int>(7,0), _superGeometry(superGeometry), _indicator(indicator), _material(material),
    _latticeL(converter.getConversionFactorLength())
{
  this->getName() = "superGeometryFacesInd";
}
template <typename T>
SuperGeometryFacesIndicator2D<T>::SuperGeometryFacesIndicator2D(SuperGeometry2D<T>& superGeometry,
    SmoothIndicatorCircle2D<T,T>& indicator, const int material, T latticeL)
  : GenericF<T,int>(7,0), _superGeometry(superGeometry), _indicator(indicator), _material(material), _latticeL(latticeL)
{
  this->getName() = "superGeometryFacesInd";
}

template <typename T>
bool SuperGeometryFacesIndicator2D<T>::operator() (T output[], const int input[])
{
  _superGeometry.communicate();
  for (int iDim = 0; iDim < 7; ++iDim) {
    output[iDim]=T();
  }
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    BlockGeometryFacesIndicator2D<T> f(_superGeometry.getBlockGeometry(iC), _indicator, _material, _latticeL);
    T outputTmp[f.getTargetDim()];
    f(outputTmp,input);
    for (int iDim = 0; iDim < 7; ++iDim) {
      output[iDim] += outputTmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 7; ++iDim) {
    singleton::mpi().reduceAndBcast( output[iDim], MPI_SUM);
  }
#endif
  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDrag2D<T,DESCRIPTOR>::SuperLatticePhysDrag2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physDrag";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  SuperGeometryFaces2D<T> faces(_superGeometry, _material, this->_converter.getConversionFactorLength());
  SuperLatticePhysBoundaryForce2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _material, this->_converter);
  SuperSum2D<T> sumF(f, _superGeometry, _material);

  T factor = 2. / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity() * this->_converter.getCharPhysVelocity());

  T facesTmp[faces.getTargetDim()], sumFTmp[sumF.getTargetDim()];
  sumF(sumFTmp, input);
  faces(facesTmp, input);
  output[0] = factor * sumFTmp[0] / facesTmp[0];
  output[1] = factor * sumFTmp[1] / facesTmp[1];

  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDragIndicator2D<T,DESCRIPTOR>::SuperLatticePhysDragIndicator2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 ParticleIndicatorF2D<T,T>& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "physDragIndicator";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysDragIndicator2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  SuperLatticePhysBoundaryForceIndicator2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _indicator, this->_converter);
  SuperSumIndicator2D<T> sumF(f, _superGeometry, _indicator);

  T factor = 2. / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity() * this->_converter.getCharPhysVelocity());

  T sumFTmp[sumF.getTargetDim()];
  sumF(sumFTmp, input);
  //output[0] = factor * sumFTmp[0] / _indicator.getDiam();
  //output[1] = factor * sumFTmp[1] / _indicator.getDiam();

  return true;
}

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysDragIndicator2D_2<T,DESCRIPTOR>::SuperLatticePhysDragIndicator2D_2
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 SmoothIndicatorF2D<T,T>& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "physDragIndicator";
}

/*
template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysDragIndicator2D_2<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  SuperLatticePhysBoundaryForceIndicator2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _indicator, this->_converter);
//  SuperLatticePhysVolumeForceIndicator2D<T,DESCRIPTOR> f(this->_sLattice, _superGeometry, _indicator, this->_converter);
  SuperSumIndicator2D<T> sumF(f, _superGeometry, _indicator);

  T sumFTmp[sumF.getTargetDim()];
  sumF(output, input);

  return true;
}
*/

template <typename T, template <typename U> class DESCRIPTOR>
SuperLatticePhysCorrDrag2D<T,DESCRIPTOR>::SuperLatticePhysCorrDrag2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physCorrDrag";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool SuperLatticePhysCorrDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  SuperGeometryFaces2D<T> faces(superGeometry, material, this->converter);

  //  SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR> f(this->sLattice, superGeometry, material, this->converter);
  //  SuperSum2D<T> sumF(f, superGeometry, material);

  //  T factor = 2. / (this->converter.getCharRho() * this->converter.getCharU() * this->converter.getCharU());

  //  std::vector<T> drag(2,T());
  //  drag[0] = factor * sumF(input)[0] / faces(input)[0];
  //  drag[1] = factor * sumF(input)[1] / faces(input)[1];
  ////  drag[2] = factor * sumF(input)[2] / faces(input)[2];
  //  return drag;
  return false;
}


} //end namespace olb

#endif
