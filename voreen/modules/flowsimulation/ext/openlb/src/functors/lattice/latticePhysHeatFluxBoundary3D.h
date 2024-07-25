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

#ifndef LATTICE_PHYS_HEAT_FLUX_BOUNDARY_3D_H
#define LATTICE_PHYS_HEAT_FLUX_BOUNDARY_3D_H

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

/// functor to get pointwise phys heat flux on a boundary with a given material on local lattice
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticePhysHeatFluxBoundary3D final : public SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
private:
  SuperGeometry<T,3>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysHeatFluxBoundary3D(SuperLattice<T,TDESCRIPTOR>& sLattice,
                                     SuperGeometry<T,3>& superGeometry, const int material,
                                     const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter,
                                     IndicatorF3D<T>& indicator);
};

/// functor returns pointwise phys heat flux on a boundary with a given material on local lattice
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class BlockLatticePhysHeatFluxBoundary3D final : public BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
private:
  BlockGeometry<T,3>& _blockGeometry;
  const int _material;
  std::vector<std::vector<std::vector<std::vector<int>>>> _discreteNormal;
  std::vector<std::vector<std::vector<std::vector<T>>>> _normal;
public:
  BlockLatticePhysHeatFluxBoundary3D(BlockLattice<T,TDESCRIPTOR>& blockLattice,
                                     BlockGeometry<T,3>& blockGeometry,
                                     int material,
                                     const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter,
                                     IndicatorF3D<T>& indicator);
  bool operator() (T output[], const int input[]) override;
};

}
#endif
