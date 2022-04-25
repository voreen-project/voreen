/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Mathias J. Krause,
 *  Albert Mink, Adrian Kummerlaender
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

#ifndef BLOCK_LATTICE_INTEGRAL_F_2D_H
#define BLOCK_LATTICE_INTEGRAL_F_2D_H

#include "functors/genericF.h"
#include "blockBaseF2D.h"
#include "geometry/blockGeometry.h"
#include "latticePhysBoundaryForce2D.h"
#include "latticePhysCorrBoundaryForce2D.h"
#include "blockGeometryFaces2D.h"
#include "integral/blockIntegralF2D.h"
#include "indicator/blockIndicatorBaseF2D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF2D.h"

/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


template<typename T> class BlockIndicatorF2D;

template <typename T, typename DESCRIPTOR>
class BlockL1Norm2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
protected:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
  BlockGeometry<T,2>& _blockGeometry;
  int _material;
public:
  BlockL1Norm2D(BlockLatticeF2D<T,DESCRIPTOR>& f, BlockGeometry<T,2>& blockGeometry, int material);
  bool operator() (T output[], const int input[]) override;
};


/// BlockL223D returns componentwise the squared l2-norm
template <typename T, typename DESCRIPTOR>
class BlockL222D final : public BlockLatticeF2D<T,DESCRIPTOR> {
protected:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
  BlockGeometry<T,2>& _blockGeometry;
  int _material;
public:
  BlockL222D(BlockLatticeF2D<T,DESCRIPTOR>& f,
             BlockGeometry<T,2>& blockGeometry,
             int material);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys force acting on a indicated boundary on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysDrag2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockIndicatorF2D<T>&                         _indicatorF;
  BlockGeometryFaces2D<T>                       _facesF;
  BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR> _pBoundForceF;
  BlockSum2D<T>                                 _sumF;

  const T _factor;
public:
  BlockLatticePhysDrag2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                         BlockIndicatorF2D<T>&                  indicatorF,
                         const UnitConverter<T,DESCRIPTOR>&     converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys force acting on a indicated boundary on local lattice
/**
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysCorrDrag2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockIndicatorF2D<T>&                             _indicatorF;
  BlockGeometryFaces2D<T>                           _facesF;
  BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR> _pBoundForceF;
  BlockSum2D<T>                                     _sumF;

  const T _factor;
public:
  BlockLatticePhysCorrDrag2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                             BlockIndicatorF2D<T>&                  indicatorF,
                             const UnitConverter<T,DESCRIPTOR>&     converter);
  bool operator() (T output[], const int input[]) override;
};


} // end namespace olb

#endif
