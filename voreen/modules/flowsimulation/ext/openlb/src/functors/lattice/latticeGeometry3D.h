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

#ifndef LATTICE_GEOMETRY_3D_H
#define LATTICE_GEOMETRY_3D_H

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


/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

/// functor to get pointwise the material no. presenting the geometry on local lattice
template <typename T, typename DESCRIPTOR=void>
class SuperLatticeGeometry3D final : public SuperF3D<T> {
private:
  SuperGeometry<T,3>& _superGeometry;
  const int _material;
public:
  SuperLatticeGeometry3D(SuperLattice<T,DESCRIPTOR>& sLattice,
                         SuperGeometry<T,3>& superGeometry, const int material = -1);
  SuperLatticeGeometry3D(SuperGeometry<T,3>& superGeometry, const int material = -1);
};

/// functor returns pointwise the material no. presenting the geometry on local lattice
template <typename T, typename DESCRIPTOR=void>
class BlockLatticeGeometry3D final : public BlockF3D<T> {
  BlockGeometry<T,3>& _blockGeometry;
  const int _material;
public:
  BlockLatticeGeometry3D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                         BlockGeometry<T,3>& blockGeometry, int material = -1);
  BlockLatticeGeometry3D(BlockGeometry<T,3>& blockGeometry, int material = -1);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
