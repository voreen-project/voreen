/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef SUPER_INTEGRAL_F_2D_HH
#define SUPER_INTEGRAL_F_2D_HH

#include <cmath>
#include <vector>

#include "superIntegralF2D.h"
#include "blockIntegralF2D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.hh"
#include "utilities/functorPtr.hh"

namespace olb {


template <typename T, typename W>
SuperSum2D<T,W>::SuperSum2D(FunctorPtr<SuperF2D<T,W>>&&        f,
                            FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF)
  : SuperF2D<T,W>(f->getSuperStructure(), f->getTargetDim()+1),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "Sum("+_f->getName()+")";

  LoadBalancer<T>& load = _f->getSuperStructure().getLoadBalancer();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockSum2D<T,W>(_f->getBlockF(iC),
                            _indicatorF->getBlockIndicatorF(iC))
      );
    }
  }
}

template <typename T, typename W>
SuperSum2D<T,W>::SuperSum2D(FunctorPtr<SuperF2D<T,W>>&& f,
                            SuperGeometry<T,2>& superGeometry,
                            const int material)
  : SuperSum2D(
      std::forward<decltype(f)>(f),
      superGeometry.getMaterialIndicator(material))
{ }

template <typename T, typename W>
bool SuperSum2D<T,W>::operator() (W output[], const int input[])
{
  _f->getSuperStructure().communicate();
  CuboidGeometry2D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&     load     = _f->getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = W(0);
  }

  if (this->_blockF.empty()) {
    W outputTmp[_f->getTargetDim()];
    int inputTmp[_f->getSourceDim()];
    std::size_t voxels(0);

    for (int iC = 0; iC < load.size(); ++iC) {
      const Cuboid2D<T> cuboid = geometry.get(load.glob(iC));
      inputTmp[0] = load.glob(iC);
      for (inputTmp[1] = 0; inputTmp[1] < cuboid.getNx(); ++inputTmp[1]) {
        for (inputTmp[2] = 0; inputTmp[2] < cuboid.getNy(); ++inputTmp[2]) {
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
    output[_f->getTargetDim()] = voxels;
  }
  else {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->getBlockF(iC)(output, input);
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim(); ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
#endif
  return true;
}


template <typename T, typename W>
SuperIntegral2D<T,W>::SuperIntegral2D(FunctorPtr<SuperF2D<T,W>>&&        f,
                                      FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF)
  : SuperF2D<T,W>(f->getSuperStructure(), f->getTargetDim()),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "Integral("+_f->getName()+")";

  LoadBalancer<T>& load = _f->getSuperStructure().getLoadBalancer();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockIntegral2D<T,W>(_f->getBlockF(iC),
                                 _indicatorF->getBlockIndicatorF(iC))
      );
    }
  }
}

template <typename T, typename W>
SuperIntegral2D<T,W>::SuperIntegral2D(FunctorPtr<SuperF2D<T,W>>&& f,
                                      SuperGeometry<T,2>& superGeometry,
                                      const int material)
  : SuperIntegral2D(
      std::forward<decltype(f)>(f),
      superGeometry.getMaterialIndicator(material))
{ }

template <typename T, typename W>
bool SuperIntegral2D<T,W>::operator() (W output[], const int input[])
{
  _f->getSuperStructure().communicate();
  CuboidGeometry2D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&     load     = _f->getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = W(0);
  }

  if (this->_blockF.empty()) {
    W outputTmp[_f->getTargetDim()];
    int inputTmp[_f->getSourceDim()];

    for (int iC = 0; iC < load.size(); ++iC) {
      const Cuboid2D<T> cuboid = geometry.get(load.glob(iC));
      const W weight = pow(cuboid.getDeltaR(), 2);
      inputTmp[0] = load.glob(iC);
      for (inputTmp[1] = 0; inputTmp[1] < cuboid.getNx(); ++inputTmp[1]) {
        for (inputTmp[2] = 0; inputTmp[2] < cuboid.getNy(); ++inputTmp[2]) {
            if (_indicatorF(inputTmp)) {
              _f(outputTmp,inputTmp);
              for (int i = 0; i < this->getTargetDim(); ++i) {
                output[i] += outputTmp[i] * weight;
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
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
#endif
  return true;
}


}

#endif
