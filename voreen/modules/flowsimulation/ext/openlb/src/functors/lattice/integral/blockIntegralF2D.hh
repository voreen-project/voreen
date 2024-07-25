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

#ifndef BLOCK_INTEGRAL_F_2D_HH
#define BLOCK_INTEGRAL_F_2D_HH

#include "blockIntegralF2D.h"
#include "core/olbDebug.h"
#include "functors/lattice/indicator/blockIndicatorBaseF2D.h"

namespace olb {


template <typename T, typename W>
BlockSum2D<T,W>::BlockSum2D(BlockF2D<W>&          f,
                            BlockIndicatorF2D<T>& indicatorF)
  : BlockF2D<W>(f.getBlockStructure(), f.getTargetDim()+1),
    _f(f),
    _indicatorF(indicatorF)
{
  this->getName() = "BlockSum("+_f.getName()+")";
}

template <typename T, typename W>
bool BlockSum2D<T,W>::operator() (W output[], const int input[])
{
  OLB_ASSERT(_f.getSourceDim() == _indicatorF.getSourceDim(),
             "functor source dimension equals indicator source dimension");

  W outputTmp[_f.getTargetDim()];
  int inputTmp[_f.getSourceDim()];
  std::size_t voxels(0);

  const auto& blockStructure = this->getBlockStructure();

  for (inputTmp[0] = 0; inputTmp[0] < blockStructure.getNx(); ++inputTmp[0]) {
    for (inputTmp[1] = 0; inputTmp[1] < blockStructure.getNy(); ++inputTmp[1]) {
        if (_indicatorF(inputTmp)) {
          _f(outputTmp,inputTmp);
          for (int i = 0; i < _f.getTargetDim(); ++i) {
            output[i] += outputTmp[i];
          }
          voxels += 1;
        }
    }
  }
  output[_f.getTargetDim()] += voxels;

  return true;
}


template <typename T, typename W>
BlockIntegral2D<T,W>::BlockIntegral2D(BlockF2D<W>&          f,
                                      BlockIndicatorF2D<T>& indicatorF)
  : BlockF2D<W>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF)
{
  this->getName() = "BlockIntegral("+_f.getName()+")";
}

template <typename T, typename W>
bool BlockIntegral2D<T,W>::operator() (W output[], const int input[])
{
  OLB_ASSERT(_f.getSourceDim() == _indicatorF.getSourceDim(),
             "functor source dimension equals indicator source dimension");

  const W weight = pow(_indicatorF.getBlockGeometry().getDeltaR(), 2);// was 3

  W outputTmp[_f.getTargetDim()];
  int inputTmp[_f.getSourceDim()];

  const auto& blockStructure = this->getBlockStructure();

  for (inputTmp[0] = 0; inputTmp[0] < blockStructure.getNx(); ++inputTmp[0]) {
    for (inputTmp[1] = 0; inputTmp[1] < blockStructure.getNy(); ++inputTmp[1]) {
        if (_indicatorF(inputTmp)) {
          _f(outputTmp,inputTmp);
          for (int i = 0; i < this->getTargetDim(); ++i) {
            output[i] += outputTmp[i] * weight;
          }
        }
    }
  }

  return true;
}


}

#endif
