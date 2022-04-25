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

#ifndef BLOCK_LATTICE_INTEGRAL_F_2D_HH
#define BLOCK_LATTICE_INTEGRAL_F_2D_HH

#include <vector>
#include <cmath>

#include "blockLatticeIntegralF2D.h"
#include "blockGeometryFaces2D.h"
#include "blockCalcF2D.h" // for IdentityF
#include "core/olbDebug.h"

namespace olb {


template <typename T, typename DESCRIPTOR>
BlockL1Norm2D<T,DESCRIPTOR>::BlockL1Norm2D(BlockLatticeF2D<T,DESCRIPTOR>& f,
    BlockGeometry<T,2>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlock(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "L1("+_f.getName()+")";
}

template <typename T, typename DESCRIPTOR>
bool BlockL1Norm2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  BlockIdentity2D<T> ff(_f); // exists only to prevent f from being deleted
  T outputTmp[this->getTargetDim()];
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T(0);
    for (int iX = 0; iX < _f.getBlock().getNx(); ++iX) {
      for (int iY = 0; iY < _f.getBlock().getNy(); ++iY) {
          if (this->_blockGeometry.getMaterial(iX, iY) == _material) {
            _f(outputTmp,iX, iY);
            T tmp = fabs(outputTmp[i]);
            if (tmp > output[i]) {
              output[i] = tmp;
            }
          }
      }
    }
  }
  return true;
}


template <typename T, typename DESCRIPTOR>
BlockL222D<T,DESCRIPTOR>::BlockL222D(BlockLatticeF2D<T,DESCRIPTOR>& f,
                                     BlockGeometry<T,2>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlock(),f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "L22("+f.getName()+")";
}


template <typename T, typename DESCRIPTOR>
bool BlockL222D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  f.getBlock().communicate();
  //  CuboidGeometry2D<T>& cGeometry = f.getBlock().get_cGeometry();
  //  loadBalancer& load = f.getBlock().get_load();

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
  //            if (this->blockGeometry.get_material({globX, globY, globZ}) == material) {
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


template <typename T, typename DESCRIPTOR>
BlockLatticePhysDrag2D<T,DESCRIPTOR>::BlockLatticePhysDrag2D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockIndicatorF2D<T>&                  indicatorF,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2),
    _indicatorF(indicatorF),
    _facesF(indicatorF, converter.getConversionFactorLength()),
    _pBoundForceF(blockLattice, indicatorF, converter),
    _sumF(_pBoundForceF, indicatorF),
    _factor(2./( converter.getPhysDensity()*converter.getCharPhysVelocity()*converter.getCharPhysVelocity() ))
{
  this->getName() = "physDrag";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T faces[5] = { };
  T sum[3]   = { };
  _sumF(sum, input);
  _facesF(faces, input);

  output[0] = _factor * sum[0] / faces[0];
  output[1] = _factor * sum[1] / faces[1];

  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysCorrDrag2D<T,DESCRIPTOR>::BlockLatticePhysCorrDrag2D(
  BlockLattice<T,DESCRIPTOR>& blockLattice,
  BlockIndicatorF2D<T>&                  indicatorF,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2),
    _indicatorF(indicatorF),
    _facesF(indicatorF, converter.getConversionFactorLength()),
    _pBoundForceF(blockLattice, indicatorF, converter),
    _sumF(_pBoundForceF, indicatorF),
    _factor(2./( converter.getPhysDensity()*converter.getCharPhysVelocity()*converter.getCharPhysVelocity() ))
{
  this->getName() = "physCorrDrag";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysCorrDrag2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T faces[5] = { };
  T sum[3]   = { };
  _facesF(faces, input);
  _sumF(sum, input);

  output[0] = _factor * sum[0] / faces[0];
  output[1] = _factor * sum[1] / faces[1];

  return true;
}


} // end namespace olb

#endif
