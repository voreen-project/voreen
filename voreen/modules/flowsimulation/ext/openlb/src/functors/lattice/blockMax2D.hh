/*  This file is part of the OpenLB library
 *
 *  Copyright (C) Adrian Kummerlaender
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

#ifndef BLOCK_MAX_2D_HH
#define BLOCK_MAX_2D_HH

#include "blockMax2D.h"

namespace olb {


template <typename T, typename W>
BlockMax2D<T,W>::BlockMax2D(BlockF2D<W>&          f,
                            BlockIndicatorF2D<T>& indicatorF,
                            Cuboid2D<T>&          cuboid)
  : BlockF2D<W>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF),
    _cuboid(cuboid)
{
  this->getName() = "BlockMax("+_f.getName()+")";
}

template <typename T, typename W>
bool BlockMax2D<T,W>::operator() (W output[], const int input[])
{
  OLB_ASSERT(_f.getSourceDim() == _indicatorF.getSourceDim(),
             "functor source dimension equals indicator source dimension");

  W outputTmp[_f.getTargetDim()];
  int inputTmp[_f.getSourceDim()];

  for (inputTmp[0] = 0; inputTmp[0] < _cuboid.getNx(); ++inputTmp[0]) {
    for (inputTmp[1] = 0; inputTmp[1] < _cuboid.getNy(); ++inputTmp[1]) {
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

  return true;
}


}

#endif
