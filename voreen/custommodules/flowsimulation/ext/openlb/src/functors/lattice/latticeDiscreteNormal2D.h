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

#ifndef LATTICE_DISCRETE_NORMAL_2D_H
#define LATTICE_DISCRETE_NORMAL_2D_H

#include <vector>

#include "core/superLattice2D.h"
#include "utilities/functorPtr.h"
#include "indicator/blockIndicatorF2D.h"


namespace olb {

template<typename T, unsigned D> class SuperGeometry;

/// functor to get pointwise the discrete normal vector of local lattice boundary cells
template <typename T, typename DESCRIPTOR>
class SuperLatticeDiscreteNormal2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
private:
  SuperGeometry<T,2>& _superGeometry;
  FunctorPtr<SuperIndicatorF2D<T>> _indicatorF;
public:
  SuperLatticeDiscreteNormal2D(SuperLattice<T,DESCRIPTOR>& sLattice,SuperGeometry<T,2>& superGeometry,
                               FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);
};

template <typename T, typename DESCRIPTOR>
SuperLatticeDiscreteNormal2D(SuperLattice<T,DESCRIPTOR>&,
                             SuperGeometry<typename SuperLattice<T,DESCRIPTOR>::value_t,2>&,
                             FunctorPtr<SuperIndicatorF2D<typename SuperLattice<T,DESCRIPTOR>::value_t>>&&)
  -> SuperLatticeDiscreteNormal2D<T,DESCRIPTOR>;

/// BlockLatticeDiscreteNormal2D returns pointwise the discrete normal vector of the local lattice boundary cells
template <typename T, typename DESCRIPTOR>
class BlockLatticeDiscreteNormal2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockGeometry<T,2>& _blockGeometry;
  BlockIndicatorF2D<T>& _indicatorF;
public:
  BlockLatticeDiscreteNormal2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                               BlockGeometry<T,2>& blockGeometry,
                               BlockIndicatorF2D<T>& indicatorF);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise the type of a discrete normal vector
template <typename T, typename DESCRIPTOR>
class SuperLatticeDiscreteNormalType2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
private:
  SuperGeometry<T,2>& _superGeometry;
  FunctorPtr<SuperIndicatorF2D<T>> _indicatorF;
public:
  SuperLatticeDiscreteNormalType2D(SuperLattice<T,DESCRIPTOR>& sLattice,SuperGeometry<T,2>& superGeometry,
                                   FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);
};

template <typename T, typename DESCRIPTOR>
SuperLatticeDiscreteNormalType2D(SuperLattice<T,DESCRIPTOR>&,
                                 SuperGeometry<typename SuperLattice<T,DESCRIPTOR>::value_t,2>&,
                                 FunctorPtr<SuperIndicatorF2D<typename SuperLattice<T,DESCRIPTOR>::value_t>>&&)
  -> SuperLatticeDiscreteNormalType2D<T,DESCRIPTOR>;


/// BlockLatticeDiscreteNormalType2D returns pointwise the type of a discrete normal vector
template <typename T, typename DESCRIPTOR>
class BlockLatticeDiscreteNormalType2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockGeometry<T,2>& _blockGeometry;
  BlockIndicatorF2D<T>& _indicatorF;
public:
  BlockLatticeDiscreteNormalType2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                   BlockGeometry<T,2>& blockGeometry,
                                   BlockIndicatorF2D<T>& indicatorF);
  bool operator() (T output[1], const int input[]) override;
};

}
#endif
