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

#ifndef LATTICE_INTERP_PHYS_VELOCITY_3D_H
#define LATTICE_INTERP_PHYS_VELOCITY_3D_H

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

template <typename T, typename DESCRIPTOR>
class SuperLatticeInterpPhysVelocity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticeInterpPhysVelocity3D(SuperLattice<T,DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR> const& converter);
  bool operator()(T output[], const int input[]) override;
  void operator()(T output[], const T input[], const int iC);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeInterpPhysVelocity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
protected:
  Cuboid3D<T>* _cuboid;
public:
  BlockLatticeInterpPhysVelocity3D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                                   const UnitConverter<T,DESCRIPTOR>& conv, Cuboid3D<T>* c);
  BlockLatticeInterpPhysVelocity3D(const BlockLatticeInterpPhysVelocity3D<T,DESCRIPTOR>& rhs);
  bool operator() (T output[3], const int input[3]) override
  {
    return false;
  }
  void operator() (T output[3], const T input[3]);
};

}
#endif
