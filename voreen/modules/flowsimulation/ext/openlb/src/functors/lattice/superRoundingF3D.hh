/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Tom Braun
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

#ifndef SUPER_ROUNDING_F_3D_HH
#define SUPER_ROUNDING_F_3D_HH

#include "superRoundingF3D.h"
#include "blockRoundingF3D.h"

namespace olb {


template <typename T, typename W>
SuperRoundingF3D<T, W>::SuperRoundingF3D(FunctorPtr<SuperF3D<T, W>>&& f,
    RoundingMode rounding)
  :SuperF3D<T, W>(f->getSuperStructure(),f->getTargetDim()), _f(f), _rounding(rounding)
{
  this->getName()="Rounding(" + _f->getName() + ")";

  LoadBalancer<T>& load = _f->getSuperStructure().getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockRoundingF3D<W>(_f->getBlockF(iC),_rounding));
  }
}

}
#endif
