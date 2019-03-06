/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Mathias J. Krause,
 *  Albert Mink, Adrian Kummeränder
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

#ifndef BLOCK_LATTICE_INTEGRAL_F_3D_HH
#define BLOCK_LATTICE_INTEGRAL_F_3D_HH

#include <vector>
#include <cmath>

#include "blockLatticeIntegralF3D.h"
#include "blockLatticeLocalF3D.h"
#include "blockCalcF3D.h" // for IdentityF
#include "core/blockLattice3D.h"
#include "core/olbDebug.h"
#include "indicator/blockIndicatorBaseF3D.h"

namespace olb {


template <typename T, typename W>
BlockMin3D<T,W>::BlockMin3D(BlockF3D<W>&          f,
                            BlockIndicatorF3D<T>& indicatorF,
                            Cuboid3D<T>&          cuboid)
  : BlockF3D<W>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF),
    _cuboid(cuboid)
{
  this->getName() = "BlockMin("+_f.getName()+")";
}

template <typename T, typename W>
bool BlockMin3D<T,W>::operator() (W output[], const int input[])
{
  OLB_ASSERT(_f.getSourceDim() == _indicatorF.getSourceDim(),
             "functor source dimension equals indicator source dimension");

  W outputTmp[_f.getTargetDim()];
  int inputTmp[_f.getSourceDim()];

  for (inputTmp[0] = 0; inputTmp[0] < _cuboid.getNx(); ++inputTmp[0]) {
    for (inputTmp[1] = 0; inputTmp[1] < _cuboid.getNy(); ++inputTmp[1]) {
      for (inputTmp[2] = 0; inputTmp[2] < _cuboid.getNz(); ++inputTmp[2]) {
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

  return true;
}


template <typename T, typename W>
BlockMax3D<T,W>::BlockMax3D(BlockF3D<W>&          f,
                            BlockIndicatorF3D<T>& indicatorF,
                            Cuboid3D<T>&          cuboid)
  : BlockF3D<W>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF),
    _cuboid(cuboid)
{
  this->getName() = "BlockMax("+_f.getName()+")";
}

template <typename T, typename W>
bool BlockMax3D<T,W>::operator() (W output[], const int input[])
{
  OLB_ASSERT(_f.getSourceDim() == _indicatorF.getSourceDim(),
             "functor source dimension equals indicator source dimension");

  W outputTmp[_f.getTargetDim()];
  int inputTmp[_f.getSourceDim()];

  for (inputTmp[0] = 0; inputTmp[0] < _cuboid.getNx(); ++inputTmp[0]) {
    for (inputTmp[1] = 0; inputTmp[1] < _cuboid.getNy(); ++inputTmp[1]) {
      for (inputTmp[2] = 0; inputTmp[2] < _cuboid.getNz(); ++inputTmp[2]) {
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

  return true;
}


template <typename T, typename W>
BlockAverage3D<T,W>::BlockAverage3D(BlockF3D<W>&          f,
                                    BlockIndicatorF3D<T>& indicatorF,
                                    Cuboid3D<T>&          cuboid)
  : BlockSum3D<T,W>(f, indicatorF, cuboid)
{
  this->getName() = "BlockAverage("+f.getName()+")";
}

template <typename T, typename W>
bool BlockAverage3D<T,W>::operator() (W output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = W(0);
  }

  const bool result = BlockSum3D<T,W>::operator()(output, input);

  for (int i = 0; i < this->getTargetDim()-1; ++i) {
    output[i] /= output[this->getTargetDim()-1];
  }

  return result;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockL1Norm3D<T,DESCRIPTOR>::BlockL1Norm3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
    BlockGeometry3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "L1("+_f.getName()+")";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockL1Norm3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  BlockIdentity3D<T> ff(_f); // exists only to prevent f from being deleted
  T outputTmp[this->getTargetDim()];
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T(0);
    for (int iX = 0; iX < _f.getBlockLattice().getNx(); ++iX) {
      for (int iY = 0; iY < _f.getBlockLattice().getNy(); ++iY) {
        for (int iZ = 0; iZ < _f.getBlockLattice().getNz(); ++iZ) {
          if (this->_blockGeometry.getMaterial(iX, iY, iZ) == _material) {
            _f(outputTmp,iX, iY, iZ);
            T tmp = fabs(outputTmp[i]);
            if (tmp > output[i]) {
              output[i] = tmp;
            }
          }
        }
      }
    }
  }
  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockL223D<T,DESCRIPTOR>::BlockL223D(BlockLatticeF3D<T,DESCRIPTOR>& f,
                                     BlockGeometry3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T,DESCRIPTOR>(f.getBlockLattice(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "L22("+f.getName()+")";
}


template <typename T, template <typename U> class DESCRIPTOR>
bool BlockL223D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  f.getBlockLattice().communicate();
  //  CuboidGeometry3D<T>& cGeometry = f.getBlockLattice().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice().get_load();

  output[0]=0;
  //  for (int i=0; i<this->n; i++) {

  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      T weight = pow(this->blockGeometry.getDeltaR(),3);
  //      for (int iX=0; iX<nX; ++iX) {
  //        for (int iY=0; iY<nY; ++iY) {
  //          for (int iZ=0; iZ<nZ; ++iZ) {
  //            int globX = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            int globY = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            int globZ = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->blockGeometry.get_material(globX, globY, globZ) == material) {
  //              tmp[i]+=f(load.glob(iC),iX,iY,iZ)[i]*f(load.glob(iC),iX,iY,iZ)[i]*weight;
  //            }
  //          }
  //        }
  //      }
  //    }
  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(tmp[i], MPI_SUM);
  //#endif
  //  }
  return true;
}


template <typename T>
BlockGeometryFaces3D<T>::BlockGeometryFaces3D(BlockGeometryStructure3D<T>& blockGeometry, int material, T latticeL)
  : GenericF<T,int>(7,0), _blockGeometry(blockGeometry), _material(material), _latticeL(latticeL)
{
  this->getName() = "blockGeometryFaces";
}

template <typename T>
bool BlockGeometryFaces3D<T>::operator() (T output[], const int input[])
{
  int counter[7] = {0};
  if (_blockGeometry.getStatistics().getNvoxel(_material)!=0) {
    const int x0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[0];
    const int y0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[1];
    const int z0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[2];
    const int x1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[0];
    const int y1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[1];
    const int z1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[2];

    // Iterate over all cells and count the cells of the face
    for (int iX = x0; iX <= x1; ++iX) {
      for (int iY = y0; iY <= y1; ++iY) {
        for (int iZ = z0; iZ <= z1; ++iZ) {
          // Lock at solid nodes only
          if (_blockGeometry.getMaterial(iX, iY, iZ) == _material) {
            if (_blockGeometry.getMaterial(iX-1, iY, iZ) == 1) {
              counter[0]++;
            }
            if (_blockGeometry.getMaterial(iX, iY-1, iZ) == 1) {
              counter[1]++;
            }
            if (_blockGeometry.getMaterial(iX, iY, iZ-1) == 1) {
              counter[2]++;
            }
            if (_blockGeometry.getMaterial(iX+1, iY, iZ) == 1) {
              counter[3]++;
            }
            if (_blockGeometry.getMaterial(iX, iY+1, iZ) == 1) {
              counter[4]++;
            }
            if (_blockGeometry.getMaterial(iX, iY, iZ+1) == 1) {
              counter[5]++;
            }
          }
        }
      }
    }

    T dx2 = _latticeL*_latticeL;
    output[6] = T();
    for (int i=0; i<6; ++i) {
      output[i] = (T)counter[i]*dx2;
      output[6] += (T) counter[i] * dx2;
    }
  } else {
    for (int i=0; i<7; ++i) {
      output[i] = T();
    }
  }
  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysDrag3D<T,DESCRIPTOR>::BlockLatticePhysDrag3D
(BlockLattice3D<T,DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physDrag";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysDrag3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  /*BlockGeometryFaces3D<T> faces(_blockGeometry, _material, this->_converter.getConversionFactorLength());
  BlockLatticePhysBoundaryForce3D<T,DESCRIPTOR> fTemp(this->_blockLattice, _blockGeometry, _material, this->_converter);
  BlockSum3D<T> sumF(fTemp, _blockGeometry, _material);

  T factor = 2. / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity() * this->_converter.getCharPhysVelocity());

  T outputSumF[3] = { T() };
  sumF(outputSumF,input);
  T outputFaces[3] = { T() };
  faces(outputFaces,input);

  output[0] = factor * outputSumF[0] / outputFaces[0];
  output[1] = factor * outputSumF[1] / outputFaces[1];
  output[2] = factor * outputSumF[2] / outputFaces[2];*/

  return true;
}


template <typename T, template <typename U> class DESCRIPTOR>
BlockLatticePhysCorrDrag3D<T,DESCRIPTOR>::BlockLatticePhysCorrDrag3D
(BlockLattice3D<T,DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,3),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physCorrDrag";
}

template <typename T, template <typename U> class DESCRIPTOR>
bool BlockLatticePhysCorrDrag3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  /*BlockGeometryFaces3D<T> faces(_blockGeometry, _material, this->_converter.getConversionFactorLength());

  BlockLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR> fTemp(this->_blockLattice, _blockGeometry, _material, this->_converter);
  BlockSum3D<T> sumF(fTemp, _blockGeometry, _material);

  T factor = 2. / (this->_converter.getPhysDensity() * this->_converter.getCharPhysVelocity() * this->_converter.getCharPhysVelocity());

  T outputSumF[3] = { T() };
  sumF(outputSumF,input);
  T outputFaces[3] = { T() };
  faces(outputFaces,input);

  output[0] = factor * outputSumF[0] / outputFaces[0];
  output[1] = factor * outputSumF[1] / outputFaces[1];
  output[2] = factor * outputSumF[2] / outputFaces[2];*/

  return true;
}


} // end namespace olb

#endif
