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

#ifndef SUPER_GEOMETRY_FACES_2D_HH
#define SUPER_GEOMETRY_FACES_2D_HH

#include "superGeometryFaces2D.h"
#include "geometry/superGeometry.h"

namespace olb {


template<typename T>
SuperGeometryFaces2D<T>::SuperGeometryFaces2D(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF, T latticeL)
  : SuperF2D<T>(indicatorF->getSuperStructure(), 5),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "superGeometryFaces";
  for (int iC = 0; iC < this->getSuperStructure().getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockGeometryFaces2D<T>(indicatorF->getBlockIndicatorF(iC), latticeL));
  }
}

template<typename T>
SuperGeometryFaces2D<T>::SuperGeometryFaces2D(
  SuperGeometry<T,2>& superGeometry, const int material, T latticeL)
  : SuperGeometryFaces2D(superGeometry.getMaterialIndicator(material), latticeL)
{ }

template<typename T>
bool SuperGeometryFaces2D<T>::operator()(T output[], const int input[])
{
  this->getSuperStructure().communicate();

  T blockOutput[5] = { };
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = T();
  }

  for (int iC = 0; iC < this->getSuperStructure().getLoadBalancer().size(); ++iC) {
    this->getBlockF(iC)(blockOutput, input);
    for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
      output[iDim] += blockOutput[iDim];
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    singleton::mpi().reduceAndBcast(output[iDim], MPI_SUM);
  }
#endif
  return true;
}

template <typename T, bool HLBM>
SuperGeometryFacesIndicator2D<T,HLBM>::SuperGeometryFacesIndicator2D(
  SuperGeometry<T,2>& superGeometry,
  SmoothIndicatorF2D<T,T,HLBM>& indicator,
  const int material, T deltaX)
  : GenericF<T,int>(5,0), _superGeometry(superGeometry), _indicator(indicator),
    _material(material), _latticeL(deltaX)
{
  this->getName() = "superGeometryFacesInd";
}

template <typename T, bool HLBM>
bool SuperGeometryFacesIndicator2D<T,HLBM>::operator() (T output[], const int input[])
{
  _superGeometry.communicate();
  for (int iDim = 0; iDim < 5; ++iDim) {
    output[iDim]=T();
  }
  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    BlockGeometryFacesIndicator2D<T,HLBM> f(_superGeometry.getBlockGeometry(iC),
                                            _indicator, _material, _latticeL);
    T outputTmp[f.getTargetDim()];
    f(outputTmp,input);
    for (int iDim = 0; iDim < 5; ++iDim) {
      output[iDim] += outputTmp[iDim];
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < 5; ++iDim) {
    singleton::mpi().reduceAndBcast( output[iDim], MPI_SUM);
  }
#endif
  return true;
}

}

#endif
