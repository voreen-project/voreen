/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Clara Schragmann
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

#ifndef LATTICE_DISCRETE_NORMAL_3D_H
#define LATTICE_DISCRETE_NORMAL_3D_H

#include<vector>

#include "functors/analytical/indicator/indicatorBaseF3D.h"

#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/blockIndicatorBaseF3D.h"


/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

/// functor to get pointwise the discrete normal vector of local lattice boundary cells
template <typename T, typename DESCRIPTOR>
class SuperLatticeDiscreteNormal3D final :public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperGeometry<T,3>& _superGeometry;
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  SuperLatticeDiscreteNormal3D(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry,
                               FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
};

/// BlockLatticeDiscreteNormal3D returns pointwise the discrete normal vector of the local lattice boundary cells
template <typename T, typename DESCRIPTOR>
class BlockLatticeDiscreteNormal3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
  BlockGeometry<T,3>& _blockGeometry;
  BlockIndicatorF3D<T>& _indicatorF;
public:
  BlockLatticeDiscreteNormal3D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                               BlockGeometry<T,3>& blockGeometry,
                               BlockIndicatorF3D<T>& indicatorF);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise the type of a discrete normal vector
template <typename T, typename DESCRIPTOR>
class SuperLatticeDiscreteNormalType3D final :public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperGeometry<T,3>& _superGeometry;
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  SuperLatticeDiscreteNormalType3D(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry,
                                   FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
};

/// BlockLatticeDiscreteNormalType3D returns pointwise the type of a discrete normal vector
template <typename T, typename DESCRIPTOR>
class BlockLatticeDiscreteNormalType3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
  BlockGeometry<T,3>& _blockGeometry;
  BlockIndicatorF3D<T>& _indicatorF;
public:
  BlockLatticeDiscreteNormalType3D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                   BlockGeometry<T,3>& blockGeometry,
                                   BlockIndicatorF3D<T>& indicatorF);
  bool operator() (T output[1], const int input[]) override;
};

}
#endif
