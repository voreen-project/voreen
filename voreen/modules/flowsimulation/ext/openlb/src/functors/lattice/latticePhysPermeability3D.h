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

#ifndef LATTICE_PHYS_PERMEABILITY_3D_H
#define LATTICE_PHYS_PERMEABILITY_3D_H

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

/// functor to get pointwise mesh-independent permeability values in (0,inf) in combination with (Extended)PorousBGKdynamics
/// note: result is cropped to 999999
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysPermeability3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysPermeability3D(SuperLattice<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter);
};

/**
 *  functor to get pointwise mesh-independent permeability values in (0,inf)
 *  in combination with (Extended)PorousBGKdynamics
 *  note: result is cropped to 999999
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysPermeability3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysPermeability3D(BlockLattice<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

//TODO: consistency with 2D (181219)
/*
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysPermeability3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometry<T,3>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysPermeability3D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                 BlockGeometry<T,3>& blockGeometry,
                                 int material, const UnitConverter<T>& converter);
  bool operator() (T output[], const int input[]);
};*/

}
#endif
