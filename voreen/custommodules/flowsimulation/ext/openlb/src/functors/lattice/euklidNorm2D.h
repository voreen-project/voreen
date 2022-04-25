/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef EUKLID_NORM_2D_H
#define EUKLID_NORM_2D_H

#include <vector>

#include "superBaseF2D.h"
#include "core/superLattice2D.h"
#include "indicator/superIndicatorBaseF2D.h"
#include "utilities/functorPtr.h"
#include "blockBaseF2D.h"
#include "geometry/blockGeometry.h"
#include "indicator/blockIndicatorF2D.h"
#include "dynamics/porousBGKdynamics.h"

namespace olb {

/// functor that returns pointwise the l2-norm, e.g. of a velocity
template <typename T, typename DESCRIPTOR>
class SuperEuklidNorm2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
private:
  SuperLatticeF2D<T,DESCRIPTOR>& _f;
public:
  SuperEuklidNorm2D(SuperLatticeF2D<T,DESCRIPTOR>& f);
};

///  BlockL2Norm2D returns pointwise the l2-norm, e.g. of a velocity.
template <typename T, typename DESCRIPTOR>
class BlockEuklidNorm2D final : public BlockF2D<T> {
private:
  BlockF2D<T>& _f;
public:
  BlockEuklidNorm2D(BlockF2D<T>& f);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
