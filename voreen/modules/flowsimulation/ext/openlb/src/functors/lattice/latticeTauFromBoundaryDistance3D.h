/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef LATTICE_TAU_FROM_BOUNDARY_DISTANCE_3D_H
#define LATTICE_TAU_FROM_BOUNDARY_DISTANCE_3D_H

#include<vector>

#include "superBaseF3D.h"
#include "superCalcF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"

#include "blockBaseF3D.h"
#include "geometry/blockGeometry.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/blockIndicatorBaseF3D.h"
#include "dynamics/smagorinskyBGKdynamics.h"
#include "dynamics/porousBGKdynamics.h"
#include "latticePhysBoundaryDistance3D.h"


/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

/// functor returns pointwise pore radius (in m) for packings of spheres given by an xmlReader
/// returns NAN for non-pore voxels
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticePhysTauFromBoundaryDistance3D final : public SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
public:
  SuperLatticePhysTauFromBoundaryDistance3D(SuperLattice<T,DESCRIPTOR>& sLattice,
      SuperGeometry<T,3>& sGeometry,
      XMLreader const& xmlReader,
      ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter,
      const T p, const T T_avg, const T c_p, const T beta, const T lambda_0, const T sigma, const T p_0, const T n_0);
};

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class BlockLatticePhysTauFromBoundaryDistance3D final  : public BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
private:
  BlockGeometry<T,3>& _blockGeometry;
  BlockLatticePhysBoundaryDistance3D<T,DESCRIPTOR> _distanceFunctor;
  const T _tmp1, _tmp2;
public:
  BlockLatticePhysTauFromBoundaryDistance3D(BlockLattice<T,DESCRIPTOR>& blockLattice,
      BlockGeometry<T,3>& blockGeometry,
      XMLreader const& xmlReader,
      ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter,
      const T p, const T T_avg, const T c_p, const T beta, const T lambda_0, const T sigma, const T p_0, const T n_0);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
