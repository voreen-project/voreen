/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Fabian Klemens, Benjamin Förster, Adrian Kummerländer
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

#ifndef SUPER_LATTICE_INTEGRAL_F_3D_HH
#define SUPER_LATTICE_INTEGRAL_F_3D_HH

#include <cmath>
#include <algorithm>

#include "superLatticeIntegralF3D.h"
#include "indicator/indicatorBaseF3D.hh"
#include "utilities/vectorHelpers.h"
#include "io/ostreamManager.h"
#include "utilities/functorPtr.hh"

using namespace olb::util;

namespace olb {


template <typename T, typename W>
SuperMin3D<T,W>::SuperMin3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                            FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperF3D<T,W>(f->getSuperStructure(), f->getTargetDim()),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "Min("+_f->getName()+")";

  LoadBalancer<T>&     load   = _f->getSuperStructure().getLoadBalancer();
  CuboidGeometry3D<T>& cuboid = _f->getSuperStructure().getCuboidGeometry();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockMin3D<T,W>(_f->getBlockF(iC),
                            _indicatorF->getBlockIndicatorF(iC),
                            cuboid.get(load.glob(iC)))
      );
    }
  }
}

template <typename T, typename W>
SuperMin3D<T,W>::SuperMin3D(FunctorPtr<SuperF3D<T,W>>&& f,
                            SuperGeometry3D<T>& superGeometry,
                            const int material)
  : SuperMin3D(
      std::forward<decltype(f)>(f),
      std::unique_ptr<SuperIndicatorF3D<T>>(
        new SuperIndicatorMaterial3D<T>(superGeometry,
                                        std::vector<int>(1, material))
      ))
{ }

template <typename T, typename W>
bool SuperMin3D<T,W>::operator() (W output[], const int input[])
{
  _f->getSuperStructure().communicate();
  CuboidGeometry3D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&     load     = _f->getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = std::numeric_limits<W>::max();
  }

  if (this->_blockF.empty()) {
    W outputTmp[_f->getTargetDim()];
    int inputTmp[_f->getSourceDim()];

    for (int iC = 0; iC < load.size(); ++iC) {
      const Cuboid3D<T> cuboid = geometry.get(load.glob(iC));
      inputTmp[0] = load.glob(iC);
      for (inputTmp[1] = 0; inputTmp[1] < cuboid.getNx(); ++inputTmp[1]) {
        for (inputTmp[2] = 0; inputTmp[2] < cuboid.getNy(); ++inputTmp[2]) {
          for (inputTmp[3] = 0; inputTmp[3] < cuboid.getNz(); ++inputTmp[3]) {
            if (_indicatorF(inputTmp)) {
              _f(outputTmp,inputTmp);
              for (int i = 0; i < this->getTargetDim(); ++i) {
                if (outputTmp[i] < output[i]) {
                  output[i] = outputTmp[i];
                }
              }
            }
          }
        }
      }
    }
  }
  else {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->getBlockF(iC)(output, input);
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim(); ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_MIN);
  }
#endif
  return true;
}


template <typename T, typename W>
SuperMax3D<T,W>::SuperMax3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                            FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperF3D<T,W>(f->getSuperStructure(), f->getTargetDim()),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "Max("+_f->getName()+")";

  LoadBalancer<T>&     load   = _f->getSuperStructure().getLoadBalancer();
  CuboidGeometry3D<T>& cuboid = _f->getSuperStructure().getCuboidGeometry();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockMax3D<T,W>(_f->getBlockF(iC),
                            _indicatorF->getBlockIndicatorF(iC),
                            cuboid.get(load.glob(iC)))
      );
    }
  }
}

template <typename T, typename W>
SuperMax3D<T,W>::SuperMax3D(FunctorPtr<SuperF3D<T,W>>&& f,
                            SuperGeometry3D<T>& superGeometry,
                            const int material)
  : SuperMax3D(
      std::forward<decltype(f)>(f),
      std::unique_ptr<SuperIndicatorF3D<T>>(
        new SuperIndicatorMaterial3D<T>(superGeometry,
                                        std::vector<int>(1, material))
      ))
{ }

template <typename T, typename W>
bool SuperMax3D<T,W>::operator() (W output[], const int input[])
{
  _f->getSuperStructure().communicate();
  CuboidGeometry3D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&     load     = _f->getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = std::numeric_limits<W>::min();
  }

  if (this->_blockF.empty()) {
    W outputTmp[_f->getTargetDim()];
    int inputTmp[_f->getSourceDim()];

    for (int iC = 0; iC < load.size(); ++iC) {
      const Cuboid3D<T> cuboid = geometry.get(load.glob(iC));
      inputTmp[0] = load.glob(iC);
      for (inputTmp[1] = 0; inputTmp[1] < cuboid.getNx(); ++inputTmp[1]) {
        for (inputTmp[2] = 0; inputTmp[2] < cuboid.getNy(); ++inputTmp[2]) {
          for (inputTmp[3] = 0; inputTmp[3] < cuboid.getNz(); ++inputTmp[3]) {
            if (_indicatorF(inputTmp)) {
              _f(outputTmp,inputTmp);
              for (int i = 0; i < this->getTargetDim(); ++i) {
                if (outputTmp[i] > output[i]) {
                  output[i] = outputTmp[i];
                }
              }
            }
          }
        }
      }
    }
  }
  else {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->getBlockF(iC)(output, input);
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim(); ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_MAX);
  }
#endif
  return true;
}


template <typename T, typename W>
SuperAverage3D<T,W>::SuperAverage3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                                    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperF3D<T,W>(f->getSuperStructure(), f->getTargetDim()+1),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "Average("+_f->getName()+")";

  LoadBalancer<T>&     load   = _f->getSuperStructure().getLoadBalancer();
  CuboidGeometry3D<T>& cuboid = _f->getSuperStructure().getCuboidGeometry();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockAverage3D<T,W>(_f->getBlockF(iC),
                                _indicatorF->getBlockIndicatorF(iC),
                                cuboid.get(load.glob(iC)))
      );
    }
  }
}

template <typename T, typename W>
SuperAverage3D<T,W>::SuperAverage3D(FunctorPtr<SuperF3D<T,W>>&& f,
                                    SuperGeometry3D<T>& superGeometry,
                                    const int material)
  : SuperAverage3D(
      std::forward<decltype(f)>(f),
      std::unique_ptr<SuperIndicatorF3D<T>>(
        new SuperIndicatorMaterial3D<T>(superGeometry,
                                        std::vector<int>(1, material))
      ))
{ }

template <typename T, typename W>
bool SuperAverage3D<T,W>::operator() (W output[], const int input[])
{
  _f->getSuperStructure().communicate();
  CuboidGeometry3D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&     load     = _f->getSuperStructure().getLoadBalancer();

  for (int i = 0; i < _f->getTargetDim(); ++i) {
    output[i] = W(0);
  }

  W outputTmp[_f->getTargetDim()];
  int inputTmp[_f->getSourceDim()];
  std::size_t voxels(0);

  for (int iC = 0; iC < load.size(); ++iC) {
    const Cuboid3D<T> cuboid = geometry.get(load.glob(iC));
    inputTmp[0] = load.glob(iC);
    for (inputTmp[1] = 0; inputTmp[1] < cuboid.getNx(); ++inputTmp[1]) {
      for (inputTmp[2] = 0; inputTmp[2] < cuboid.getNy(); ++inputTmp[2]) {
        for (inputTmp[3] = 0; inputTmp[3] < cuboid.getNz(); ++inputTmp[3]) {
          if (_indicatorF(inputTmp)) {
            _f(outputTmp,inputTmp);
            for (int i = 0; i < _f->getTargetDim(); ++i) {
              output[i] += outputTmp[i];
            }
            voxels += 1;
          }
        }
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < _f->getTargetDim(); ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(voxels, MPI_SUM);
#endif

  output[_f->getTargetDim()] = voxels;
  for (int i = 0; i < _f->getTargetDim(); ++i) {
    output[i] /= output[_f->getTargetDim()];
  }

  return true;
}


template <typename T, typename W>
SuperSumIndicator3D<T,W>::SuperSumIndicator3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, ParticleIndicatorF3D<T,T>& indicator)
  : SuperF3D<T,W>(f.getSuperStructure(),f.getTargetDim()+1),
    _f(f), _superGeometry(superGeometry), _indicator(indicator)
{
  this->getName() = "Sum("+_f.getName()+")";
}


template <typename T, typename W>
bool SuperSumIndicator3D<T,W>::operator() (T output[], const int input[])
{
  _f.getSuperStructure().communicate();
  LoadBalancer<T>& load = _f.getSuperStructure().getLoadBalancer();

  int numVoxels(0);
  T outputTmp[_f.getTargetDim()];
  T inside[1];
  Cuboid3D<T>* cub = nullptr;
  int start[3] = {0}, span[3] = {0};

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T(0);
  }

  for (int iC = 0; iC < load.size(); ++iC) {
    int globiC = load.glob(iC);
    cub = &_superGeometry.getCuboidGeometry().get(globiC);

    // check for intersection of cubiod and indicator
    if (cub->getOrigin()[0] <= _indicator.getMax()[0]+_indicator.getPos()[0]
        && cub->getOrigin()[1] <= _indicator.getMax()[1]+_indicator.getPos()[1]
        && cub->getOrigin()[2] <= _indicator.getMax()[2]+_indicator.getPos()[2]
        && _indicator.getMin()[0]+_indicator.getPos()[0] <= cub->getOrigin()[0] + cub->getExtend()[0] * cub->getDeltaR()
        && _indicator.getMin()[1]+_indicator.getPos()[1] <= cub->getOrigin()[1] + cub->getExtend()[1] * cub->getDeltaR()
        && _indicator.getMin()[2]+_indicator.getPos()[2] <= cub->getOrigin()[2] + cub->getExtend()[2] * cub->getDeltaR() ) {

      // compute size of intersection for iteration
      T invDeltaR = 1./cub->getDeltaR();
      for (int k=0; k<3; k++) {
        start[k] = int( (_indicator.getPos()[k]+_indicator.getMin()[k] - cub->getOrigin()[k]) * invDeltaR );
        if (start[k] < 0) {
          start[k] = 0;
        }
        span[k] = int( (_indicator.getMax()[k] - _indicator.getMin()[k])*invDeltaR + 3 );
        if (span[k] + start[k] > cub->getExtend()[k]) {
          span[k] = cub->getExtend()[k] - start[k];
        }
      }

      for (int iX = start[0]; iX < start[0]+span[0]; iX++) {
        for (int iY = start[1]; iY < start[1]+span[1]; iY++) {
          for (int iZ = start[2]; iZ < start[2]+span[2]; iZ++) {

            _indicator( inside, &(_superGeometry.getPhysR(load.glob(iC), iX, iY, iZ)[0]) );
            if ( !util::nearZero(inside[0]) ) {
              _f(outputTmp,load.glob(iC),iX,iY,iZ);
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
  output[this->getTargetDim() - 1] = numVoxels;
  return true;
}


template<typename T>
template<template<typename U> class DESCRIPTOR>
SuperGeometryFaces3D<T>::SuperGeometryFaces3D(SuperGeometry3D<T>& superGeometry,
    const int material,
    const UnitConverter<T,DESCRIPTOR>& converter)
  : GenericF<T, int>(7, 0),
    _superGeometry(superGeometry),
    _material(material),
    _latticeL(converter.getConversionFactorLength())
{
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    _blockGeometryFaces.push_back(
      new BlockGeometryFaces3D<T>(_superGeometry.getBlockGeometry(iC),
                                  _material, _latticeL));
  }
  this->getName() = "superGeometryFaces";
}

template<typename T>
SuperGeometryFaces3D<T>::SuperGeometryFaces3D(SuperGeometry3D<T>& superGeometry,
    const int material,
    T latticeL)
  : GenericF<T, int>(7, 0),
    _superGeometry(superGeometry),
    _material(material),
    _latticeL(latticeL)
{
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    _blockGeometryFaces.push_back(
      new BlockGeometryFaces3D<T>(_superGeometry.getBlockGeometry(iC),
                                  _material, _latticeL));
  }
  this->getName() = "superGeometryFaces";
}

template<typename T>
bool SuperGeometryFaces3D<T>::operator()(T output[], const int input[])
{
  _superGeometry.communicate();
  T tmp[7] = { T() };
  for (int iDim = 0; iDim < 7; ++iDim) {
    output[iDim] = T();
  }
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    (*(_blockGeometryFaces[iC]))(tmp, input);
    for (int iDim = 0; iDim < 7; ++iDim) {
      output[iDim] += tmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 7; ++iDim) {
    singleton::mpi().reduceAndBcast(output[iDim], MPI_SUM);
  }
#endif
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysDrag3D<T, DESCRIPTOR>::SuperLatticePhysDrag3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry),
    _material(material),
    _faces(_superGeometry, _material, this->_converter),
    _pBoundForce(this->_sLattice, _superGeometry, _material,
                 this->_converter),
    _sumF(_pBoundForce, _superGeometry, _material),
    _factor(
      2.
      / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity()
         * this->_converter.getCharPhysVelocity()))
{
  this->getName() = "physDrag";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysDrag3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  T faces[7] = { 0 };
  T sumF[4] = { 0 };
  _sumF(sumF, input);
  _faces(faces, input);
  output[0] = _factor * sumF[0] / faces[0];
  output[1] = _factor * sumF[1] / faces[1];
  output[2] = _factor * sumF[2] / faces[2];
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysDragIndicator3D<T, DESCRIPTOR>::SuperLatticePhysDragIndicator3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  ParticleIndicatorSphere3D<T, T>& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry),
    _indicator(indicator)
{
  this->getName() = "physDrag";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysDragIndicator3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  SuperLatticePhysBoundaryForceIndicator3D<T,DESCRIPTOR> pBoundForce(this->_sLattice, _superGeometry, _indicator, this->_converter);
  SuperSumIndicator3D<T,T> sumF(pBoundForce, _superGeometry, _indicator);

  T factor = 2. / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity() * this->_converter.getCharPhysVelocity());
  T tmp[4] = {};
  sumF(tmp,input);
  output[0] = factor * tmp[0] / _indicator.getDiam();//faces(input)[0];
  output[1] = factor * tmp[1] / _indicator.getDiam();//faces(input)[1];
  output[2] = factor * tmp[2] / _indicator.getDiam();//faces(input)[2];
  return true;
}

template<typename T, template<typename U> class DESCRIPTOR>
SuperLatticePhysCorrDrag3D<T, DESCRIPTOR>::SuperLatticePhysCorrDrag3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "physCorrDrag";
}

template<typename T, template<typename U> class DESCRIPTOR>
bool SuperLatticePhysCorrDrag3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  SuperGeometryFaces3D<T> faces(_superGeometry, _material, this->_converter);
  SuperLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR>  pBoundForce(this->_sLattice, _superGeometry,
      _material, this->_converter);
  SuperSum3D<T,T> sumF(pBoundForce, _superGeometry, _material);

  T factor = 2. / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity() * this->_converter.getCharPhysVelocity());

  T sum_tmp[4] = {};
  T face_tmp[7] = {};
  sumF(sum_tmp,input);
  faces(face_tmp,input);
  output[0] = factor * sum_tmp[0] / face_tmp[0];
  output[1] = factor * sum_tmp[1] / face_tmp[1];
  output[2] = factor * sum_tmp[2] / face_tmp[2];
  return true;
}


}  // end namespace olb

#endif
