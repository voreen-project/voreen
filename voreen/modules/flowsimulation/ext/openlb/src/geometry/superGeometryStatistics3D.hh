/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Mathias J. Krause
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
 * Representation of a statistic for a parallel 3D geometry -- generic implementation.
 */

#include <iostream>
#include "utilities/omath.h"
#include <math.h>
#include <sstream>

#include "geometry/superGeometry.h"
#include "geometry/superGeometryStatistics3D.h"
#include "utilities/vectorHelpers.h"
#include "core/vector.h"


namespace olb {

template<typename T>
SuperGeometryStatistics3D<T>::SuperGeometryStatistics3D(SuperGeometry<T,3>* superGeometry)
  : _superGeometry(superGeometry), _statisticsUpdateNeeded(true),
    _nMaterials(0), _minOverMaterial(0), _maxOverMaterial(0),
    clout(std::cout,"SuperGeometryStatistics3D")
{
}

template<typename T>
SuperGeometryStatistics3D<T>::SuperGeometryStatistics3D(SuperGeometryStatistics3D const& rhs)
  : _superGeometry(rhs._superGeometry), _statisticsUpdateNeeded(true),
    _nMaterials(rhs._nMaterials),
    _minOverMaterial(rhs._minOverMaterial), _maxOverMaterial(rhs._maxOverMaterial),
    clout(std::cout,"SuperGeometryStatistics3D")
{
}

template<typename T>
SuperGeometryStatistics3D<T>& SuperGeometryStatistics3D<T>::operator=(SuperGeometryStatistics3D const& rhs)
{
  _superGeometry = rhs._superGeometry;
  _statisticsUpdateNeeded = true;
  _nMaterials = rhs._nMaterials;
  _minOverMaterial = rhs._minOverMaterial;
  _maxOverMaterial = rhs._maxOverMaterial;
  return *this;
}


template<typename T>
bool& SuperGeometryStatistics3D<T>::getStatisticsStatus()
{
  return _statisticsUpdateNeeded;
}

template<typename T>
bool const & SuperGeometryStatistics3D<T>::getStatisticsStatus() const
{
  return _statisticsUpdateNeeded;
}


template<typename T>
void SuperGeometryStatistics3D<T>::update(bool verbose)
{
  const_this = const_cast<const SuperGeometryStatistics3D<T>*>(this);

#ifdef PARALLEL_MODE_MPI
  // This needs to be done with integers due to undesired behaviour with bool and LOR implicated by the MPI standard
  int updateReallyNeededGlobal = 0;
  if (_statisticsUpdateNeeded) {
    updateReallyNeededGlobal = 1;
  }
  singleton::mpi().reduceAndBcast(updateReallyNeededGlobal, MPI_SUM);

  if (updateReallyNeededGlobal>0) {
    _statisticsUpdateNeeded = true;
  }
  //singleton::mpi().reduceAndBcast(_statisticsUpdateNeeded, MPI_LOR);
#endif

  // check if update is really needed
  if (_statisticsUpdateNeeded ) {
    int updateReallyNeeded = 0;
    for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
      if (_superGeometry->getBlockGeometry(iCloc).getStatistics().getStatisticsStatus() ) {
        auto& blockGeometry = const_cast<BlockGeometry<T,3>&>(_superGeometry->getBlockGeometry(iCloc));
        blockGeometry.getStatistics(false).update(false);
        updateReallyNeeded++;
      }
    }


#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(updateReallyNeeded, MPI_SUM);
#endif

    if (updateReallyNeeded==0) {
      _statisticsUpdateNeeded = false;
      //      clout << "almost updated" << std::endl;
      return;
    }


    // get total number of different materials right
    _material2n = std::map<int, std::size_t>();
    _nMaterials = int();
    for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
      if (_superGeometry->getBlockGeometry(iCloc).getStatistics(false).getNmaterials() > 0) {
        _nMaterials = _superGeometry->getBlockGeometry(iCloc).getStatistics(false).getNmaterials();
      }
    }


#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(_nMaterials, MPI_SUM);
#endif

    // store the number and min., max. possition for each rank
    for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
      std::map<int, std::size_t> material2n = _superGeometry->getBlockGeometry(iCloc).getStatistics(false).getMaterial2n();
      std::map<int, std::size_t>::iterator iter;

      for (iter = material2n.begin(); iter != material2n.end(); iter++) {
        if (iter->second!=0) {
          std::vector<T> minPhysR = _superGeometry->getBlockGeometry(iCloc).getStatistics(false).getMinPhysR(iter->first);
          std::vector<T> maxPhysR = _superGeometry->getBlockGeometry(iCloc).getStatistics(false).getMaxPhysR(iter->first);
          if (_material2n.count(iter->first) == 0) {
            _material2n[iter->first] = iter->second;
            _material2min[iter->first] = minPhysR;
            _material2max[iter->first] = maxPhysR;
          }
          else {
            _material2n[iter->first] += iter->second;
            for (int iDim=0; iDim<3; iDim++) {
              if (_material2min[iter->first][iDim] > minPhysR[iDim]) {
                _material2min[iter->first][iDim] = minPhysR[iDim];
              }
              if (_material2max[iter->first][iDim] < maxPhysR[iDim]) {
                _material2max[iter->first][iDim] = maxPhysR[iDim];
              }
            }
          }
        }
      }
    }

    // store the number and min., max. possition for all ranks
#ifdef PARALLEL_MODE_MPI
    std::ptrdiff_t materials[_nMaterials];
    std::ptrdiff_t materialsInBuf[_nMaterials];
    std::ptrdiff_t materialCount[_nMaterials];
    std::ptrdiff_t materialCountInBuf[_nMaterials];
    T materialMinR[3*_nMaterials];
    T materialMaxR[3*_nMaterials];
    T materialMinRinBuf[3*_nMaterials];
    T materialMaxRinBuf[3*_nMaterials];

    for (int iM=0; iM<_nMaterials; iM++) {
      materials[iM]=-1;
      materialCount[iM]=0;
      for (int iDim=0; iDim<3; iDim++) {
        materialMinRinBuf[3*iM + iDim] = T();
        materialMaxRinBuf[3*iM + iDim] = T();
      }
    }
    std::size_t counter = 0;
    std::map<int, std::size_t>::iterator iMaterial;
    for (iMaterial = _material2n.begin(); iMaterial != _material2n.end(); iMaterial++) {
      materials[counter] = iMaterial->first;
      materialCount[counter] = iMaterial->second;
      for (int iDim=0; iDim<3; iDim++) {
        materialMinR[3*counter + iDim] = _material2min[iMaterial->first][iDim];
        materialMaxR[3*counter + iDim] = _material2max[iMaterial->first][iDim];
      }
      counter++;
    }

    for (int iRank=1; iRank<singleton::mpi().getSize(); iRank++) {
      int myRank = singleton::mpi().getRank();
      singleton::mpi().sendRecv(materials, materialsInBuf, _nMaterials,
                                (myRank+iRank)%singleton::mpi().getSize(),
                                (myRank-iRank+singleton::mpi().getSize())%singleton::mpi().getSize(), 0);
      singleton::mpi().sendRecv(materialCount, materialCountInBuf, _nMaterials,
                                (myRank+iRank)%singleton::mpi().getSize(),
                                (myRank-iRank+singleton::mpi().getSize())%singleton::mpi().getSize(), 1);
      singleton::mpi().sendRecv(materialMinR, materialMinRinBuf, 3*_nMaterials,
                                (myRank+iRank)%singleton::mpi().getSize(),
                                (myRank-iRank+singleton::mpi().getSize())%singleton::mpi().getSize(), 2);
      singleton::mpi().sendRecv(materialMaxR, materialMaxRinBuf, 3*_nMaterials,
                                (myRank+iRank)%singleton::mpi().getSize(),
                                (myRank-iRank+singleton::mpi().getSize())%singleton::mpi().getSize(), 3);
      for (int iM=0; iM<_nMaterials; iM++) {
        if (materialsInBuf[iM]!=-1) {
          std::vector<T> minPhysR(3,T());
          std::vector<T> maxPhysR(3,T());
          for (int iDim=0; iDim<3; iDim++) {
            minPhysR[iDim] = materialMinRinBuf[3*iM + iDim];
            maxPhysR[iDim] = materialMaxRinBuf[3*iM + iDim];
          }
          if (_material2n.count(materialsInBuf[iM]) == 0) {
            _material2n[materialsInBuf[iM]] = materialCountInBuf[iM];
            _material2min[materialsInBuf[iM]] = minPhysR;
            _material2max[materialsInBuf[iM]] = maxPhysR;
          }
          else {
            _material2n[materialsInBuf[iM]] += materialCountInBuf[iM];
            for (int iDim=0; iDim<3; iDim++) {
              if (_material2min[materialsInBuf[iM]][iDim] > minPhysR[iDim]) {
                _material2min[materialsInBuf[iM]][iDim] = minPhysR[iDim];
              }
              if (_material2max[materialsInBuf[iM]][iDim] < maxPhysR[iDim]) {
                _material2max[materialsInBuf[iM]][iDim] = maxPhysR[iDim];
              }
            }
          }
        }
      }
    }
#endif

    // update _minOverMaterial and _maxOverMaterial
    typename std::map< int, olb::Vector<T, 3> >::const_iterator iter;
    // find componentwise minimal extension over all material numbers
    for ( int iDim = 0; iDim < 3; iDim++ ) {
      // minimum
      for ( iter = _material2min.begin(); iter != _material2min.end(); iter++ ) {
        if ( iter->first != 0 ) { // only relevant material number are considered
          if ( _minOverMaterial[iDim] > iter->second[iDim] ) {
            _minOverMaterial[iDim] = iter->second[iDim];
          }
        }
      }
      // maximum
      for ( iter = _material2max.begin(); iter != _material2max.end(); iter++ ) {
        if ( iter->first != 0 ) { // only relevant material number are considered
          if ( _maxOverMaterial[iDim] < iter->second[iDim] ) {
            _maxOverMaterial[iDim] = iter->second[iDim];
          }
        }
      }
    }

    //clout.setMultiOutput(true);
//    print();
    //clout.setMultiOutput(false);

    if (verbose) {
      clout << "updated" << std::endl;
    }
    _statisticsUpdateNeeded = false;
  }
}

template<typename T>
int SuperGeometryStatistics3D<T>::getNmaterials()
{
  update();
  return const_this->getNmaterials();
}

template<typename T>
int SuperGeometryStatistics3D<T>::getNmaterials() const
{
  return _nMaterials;
}

template<typename T>
std::size_t SuperGeometryStatistics3D<T>::getNvoxel(int material)
{
  update(true);
  return const_this->getNvoxel(material);
}

template<typename T>
std::size_t SuperGeometryStatistics3D<T>::getNvoxel(int material) const
{
  try {
    return _material2n.at(material);
  }
  catch (std::out_of_range& ex) {
    return 0;
  }
}

template<typename T>
std::size_t SuperGeometryStatistics3D<T>::getNvoxel()
{
  update();
  return const_this->getNvoxel();
}

template<typename T>
std::size_t SuperGeometryStatistics3D<T>::getNvoxel() const
{
  std::size_t total = 0;
  for (const auto& material : _material2n) {
    if (material.first!=0) {
      total += material.second;
    }
  }
  return total;
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getMinPhysR(int material)
{
  update();
  return const_this->getMinPhysR(material);
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getMinPhysR(int material) const
{
  try {
    return _material2min.at(material);
  }
  catch (std::out_of_range& ex) {
    return {0, 0, 0};
  }
}


template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getMinPhysR()
{
  update();
  return const_this->getMinPhysR();
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getMinPhysR() const
{
  return _minOverMaterial;
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getMaxPhysR(int material)
{
  update();
  return const_this->getMaxPhysR(material);
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getMaxPhysR(int material) const
{
  try {
    return _material2max.at(material);
  }
  catch (std::out_of_range& ex) {
    return {0, 0, 0};
  }
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getMaxPhysR()
{
  update();
  return const_this->getMaxPhysR();
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getMaxPhysR() const
{
  return _maxOverMaterial;
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getPhysExtend(int material)
{
  update();
  return const_this->getPhysExtend(material);
}

template<typename T>
std::vector<T> SuperGeometryStatistics3D<T>::getPhysExtend(int material) const
{
  try {
    std::vector<T> extend;
    for (int iDim = 0; iDim < 3; iDim++) {
      extend.push_back(_material2max.at(material)[iDim] - _material2min.at(material)[iDim]);
    }
    return extend;
  }
  catch (std::out_of_range& ex) {
    return { 0, 0, 0 };
  }
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getPhysRadius(int material)
{
  update();
  return const_this->getPhysRadius(material);
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getPhysRadius(int material) const
{
  olb::Vector<T, 3> radius;
  for (int iDim=0; iDim<3; iDim++) {
    //radius.push_back((getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim])/2.);
    radius[iDim] = (getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim])/2.;
  }
  return radius;
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getCenterPhysR(int material)
{
  update();
  return const_this->getCenterPhysR(material);
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::getCenterPhysR(int material) const
{
  olb::Vector<T, 3> center;
  for (int iDim=0; iDim<3; iDim++) {
    center[iDim] = getMinPhysR(material)[iDim] + getPhysRadius(material)[iDim];
  }
  return center;
}

template<typename T>
olb::Vector<int, 3> SuperGeometryStatistics3D<T>::getType(int iC, int iX, int iY, int iZ)
{
  update();
  return const_this->getType(iC, iX, iY, iZ);
}

template<typename T>
olb::Vector<int, 3> SuperGeometryStatistics3D<T>::getType(int iC, int iX, int iY, int iZ) const
{
  int iCloc=_superGeometry->getLoadBalancer().loc(iC);
  olb::Vector<int, 3> discreteNormal = _superGeometry->getBlockGeometry(iCloc).getStatistics(false).getType(iX, iY, iZ);
  return discreteNormal;
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::computeNormal(int material)
{
  update();
  return const_this->computeNormal(material);
}

template<typename T>
olb::Vector<T, 3> SuperGeometryStatistics3D<T>::computeNormal(int material) const
{
  std::vector<T> normal (3,int());
  for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
    for (int iDim=0; iDim<3; iDim++) {
      if (_superGeometry->getBlockGeometry(iCloc).getStatistics(false).getNvoxel(material)!=0) {
        normal[iDim] += _superGeometry->getBlockGeometry(iCloc).getStatistics(false).computeNormal(material)[iDim]*_superGeometry->getBlockGeometry(iCloc).getStatistics(false).getNvoxel(material);
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int iDim=0; iDim<3; iDim++) {
    singleton::mpi().reduceAndBcast((normal[iDim]), MPI_SUM);
  }
#endif

  if (getNvoxel(material) == 0) {
    std::cerr << "Unkown material number: " << material << std::endl;
    std::exit(-1);
  }

  for (int iDim=0; iDim<3; iDim++) {
    normal[iDim] /= getNvoxel(material);
  }

  T norm = util::sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  if (norm>0.) {
    normal[0]/=norm;
    normal[1]/=norm;
    normal[2]/=norm;
  }
  return normal;
}

template<typename T>
olb::Vector<int, 3> SuperGeometryStatistics3D<T>::computeDiscreteNormal(int material, T maxNorm)
{
  update();
  return const_this->computeDiscreteNormal(material, maxNorm);
}

template<typename T>
olb::Vector<int, 3> SuperGeometryStatistics3D<T>::computeDiscreteNormal(int material, T maxNorm) const
{
  olb::Vector<T,3> normal = computeNormal(material);
  olb::Vector<int, 3> discreteNormal(0);

  T smallestAngle = T(0);
  for (int iX = -1; iX<=1; iX++) {
    for (int iY = -1; iY<=1; iY++) {
      for (int iZ = -1; iZ<=1; iZ++) {
        T norm = util::sqrt(iX*iX+iY*iY+iZ*iZ);
        if (norm>0.&& norm<maxNorm) {
          T angle = (iX*normal[0] + iY*normal[1] + iZ*normal[2])/norm;
          if (angle>=smallestAngle) {
            smallestAngle=angle;
            discreteNormal[0] = iX;
            discreteNormal[1] = iY;
            discreteNormal[2] = iZ;
          }
        }
      }
    }
  }
  return discreteNormal;
}

template<typename T>
T SuperGeometryStatistics3D<T>::computeMaxPhysDistance( int material ) const
{
  Vector<T,3> vec(getMaxPhysR(material)[0] -getMinPhysR(material)[0], getMaxPhysR(material)[1] -getMinPhysR(material)[1], getMaxPhysR(material)[2] -getMinPhysR(material)[2]);
  return norm(vec);
}

template<typename T>
T SuperGeometryStatistics3D<T>::computeMaxPhysDistance() const
{
  Vector<T,3> vec(getMaxPhysR()[0] -getMinPhysR()[0], getMaxPhysR()[1] -getMinPhysR()[1], getMaxPhysR()[2] -getMinPhysR()[2]);
  return norm(vec);
}

template<typename T>
void SuperGeometryStatistics3D<T>::print()
{
  update();
  return const_this->print();
}

template<typename T>
void SuperGeometryStatistics3D<T>::print() const
{
  try {
    for (const auto& material : _material2n) {
      clout << "materialNumber=" << material.first
            << "; count=" << material.second
            << "; minPhysR=(" << _material2min.at(material.first)[0] <<","<< _material2min.at(material.first)[1] <<","<< _material2min.at(material.first)[2] <<")"
            << "; maxPhysR=(" << _material2max.at(material.first)[0] <<","<< _material2max.at(material.first)[1] <<","<< _material2max.at(material.first)[2] <<")"
            << std::endl;
    }
  }
  catch (std::out_of_range& ex) {
  }
}



} // namespace olb
